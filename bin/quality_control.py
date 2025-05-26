import math
import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew
from typing import List, Tuple
from itertools import groupby
from operator import itemgetter
from config import AVERAGING_PERIOD_MINUTES, FREQUENCY_HZ, ANEMOMETER_FILTER, IRGA_FILTER


def instrument(df: pd.DataFrame) -> Tuple[float, float]:
    """Instrument diagnostics

    If available records for current time period is less than a threshold, or is more than 18000, a flag is raised
    """
    SECONDS_IN_MINUTE = 60
    total_records = AVERAGING_PERIOD_MINUTES * SECONDS_IN_MINUTE * FREQUENCY_HZ

    return len(df[ANEMOMETER_FILTER]) / total_records, len(df[IRGA_FILTER]) / total_records


def spikes(df: pd.DataFrame) -> Tuple[pd.Series, float]:
    """Detect and remove spikes

    The method computes the mean and standard deviation for a series of moving windows of length L1. The window moves
    one point at a time through the series. Any point in the window that is more than 3.5 standard deviations from the
    window mean is considered a spike. The point is replaced using linear interpolation between data points. When four
    or more consecutive points are detected, they are not considered spikes and are not replaced. The entire process is
    repeated until no more spikes are detected. During the second pass, when the standard deviations may be smaller if
    spikes were replaced on the previous pass, the threshold for spike detection increases to 3.6 standard deviations
    and a like amount for each subsequent pass. The record is hard flagged when the total number of spikes replaced
    exceeds 1% of the total number of data points.
    """
    MAX_WIDTH_ELECTRONIC_SPIKES: int = 3
    L1_SECONDS: int = 300
    window_position: int = 0
    spike_detection_threshold: float = 4.5
    threshold_increment: float = 0.1
    pass_spikes: list = []
    all_spikes: list = []

    L1: int = L1_SECONDS * FREQUENCY_HZ

    df['spike'] = False
    col: str = df.columns[0]
    spikes_found: list = []

    while True:
        # Loop through each window is very computationally expensive. Instead, the evaluations are performed using a
        # vectorized method.
        # First calculate the rolling mean and standard deviation of each window with a width of L1. Then for each
        # window, calculate the acceptable range of values,
        #   left: window_mean - spike_detection_threshold * window_std, and
        #   right: window_mean + spike_detection_thresold * window_std].
        # Because the rolling method labels each window at the right edge of the window index, value at index k should
        # be evaluated in those windows labeled k - L1 to k. Calculate the rolling min and max of those window's left
        # and right, and align with the values to be evaluated.

        # Set min_periods = L1 // 2 to force calculating means and standard deviations when missing data are present.
        df['rolling_mean'] = df[col].rolling(L1, min_periods=L1 // 2).mean()
        df['rolling_std'] = df[col].rolling(L1, min_periods=L1 // 2).std()
        # Set means and standard deviations to NaN inside the first moving window to ensure the first window size to be
        # L1
        df.loc[:L1 - 2, ['rolling_mean', 'rolling_std']] = np.nan, np.nan

        # Flag spikes
        df['max'] = df['rolling_mean'] + spike_detection_threshold * df['rolling_std']
        df['min'] = df['rolling_mean'] - spike_detection_threshold * df['rolling_std']
        df['rolling_max'] = df['max'].iloc[::-1].rolling(L1, min_periods=1).min().iloc[::-1]
        df['rolling_min'] = df['min'].iloc[::-1].rolling(L1, min_periods=1).max().iloc[::-1]

        df['spike'] = (df[col] < df['rolling_min']) | (df[col] > df['rolling_max'])

        # Return if no spikes can be found any more
        if not df['spike'].any():
            return df[col], len(spikes_found) / len(df)

        # Find consecutive spikes and remove them
        for _, g in groupby(enumerate(df[df['spike']].index), lambda x: x[0] - x[1]):
            indices = list(map(itemgetter(1), g))
            if len(indices) > MAX_WIDTH_ELECTRONIC_SPIKES:
                df.loc[indices, 'spike'] = False

        # Record all spikes that are found
        if not spikes_found:
            spikes_found = df[df['spike']].index.to_list()

        # Note that the code below may fail to interpolate at some spikes, e.g., when the data records before and/or
        # after a spike are NaNs. However, this should be the expected behavior since NaNs indicate bad data
        df.loc[df['spike'], col] = np.interp(
            df[df['spike']].index.values.astype(float),
            df[~df['spike']].index.values.astype(float),
            df[~df['spike']][col].values,
        )

        # Reset spike flags for next pass
        df['spike'] = False
        spike_detection_threshold += threshold_increment
        threshold_increment += 0.1


def amplitude_resolution_dropouts(df: pd.DataFrame) -> Tuple[float, float, float]:
    """Detect resolution problems and dropouts

    An amplitude resolution problem is detected by computing a series of discrete frequency distributions for half-
    overlapping windows of length 1000 data points. These windows move one-half the window width at a time through the
    series. For each window position, the number of bins is set to 100 and the interval for the distribution is taken as
    the smaller of seven standard deviations and the range. When the number of empty bins in the discrete frequency
    distribution exceeds a critical threshold value, the record is hard flagged as a resolution problem.

    Dropouts are identified using the same window and frequency distributions used for the resolution problem.
    Consecutive points that fall into the same bin of the frequency distribution are tentatively identified as dropouts.
    When the total number of dropouts in the record exceeds a threshold value, the record is flagged for dropouts.
    """
    L1: int = 1000
    N_BINS: int = 100
    empty_bin_ratio: float = 0.0
    dropout_ratio: float = 0.0
    extreme_dropout_ratio: float = 0.0

    for window_position in range(0, len(df) - L1 + 1, L1 //2):
        # If missing half of data in the window, skip
        if df[window_position:window_position + L1].count()[0] < L1 // 2:
            continue

        array: np.array = df[window_position:window_position + L1].values

        window_mean: float = np.nanmean(array)
        distribution_interval: float = min(7.0 * np.nanstd(array), np.nanmax(array) - np.nanmin(array))
        bins: np.array = np.arange(
            window_mean - distribution_interval / 2.0,
            window_mean + distribution_interval / 2.0 + distribution_interval / N_BINS,
            distribution_interval / N_BINS,
        )

        # Compute the histogram of variables, and count empty bins
        hist, _ = np.histogram(array, bins=bins)
        empty_bin_ratio = max(empty_bin_ratio, np.count_nonzero(hist==0) / N_BINS)

        # Find consecutive points that fall into the same bin of the frequency distribution, following
        # https://stackoverflow.com/questions/6352425/whats-the-most-pythonic-way-to-identify-consecutive-duplicates-in-a-list
        indices = np.digitize(array, bins)
        grouped_indices = [(k, sum(1 for _ in g)) for k, g in groupby(indices)]

        max_dropouts = max(grouped_indices, key=lambda x: x[1])
        dropout_ratio = max(dropout_ratio, max_dropouts[1] / len(df))

        # Find extreme dropouts which are less than the 10th or greater than the 90th percentile values of the
        # distribution
        if  (max_dropouts[0] < 0.1 * N_BINS or max_dropouts[0] > 0.9 * N_BINS):
            extreme_dropout_ratio = max(extreme_dropout_ratio, max_dropouts[1] / len(df))

    return empty_bin_ratio, dropout_ratio, extreme_dropout_ratio


def detrend(df: pd.DataFrame) -> pd.Series:
    """Linear de-trend
    """
    col = df.columns[0]
    mask = ~df[col].isna()

    s = np.polyfit(df[mask].index, df[mask][col].values, 1)

    linear_fit = s[0] * df[mask].index + s[1]
    df.loc[mask, f'{col}_'] = df[mask][col].values - linear_fit

    return df[f'{col}_']


def higher_moment_statistics(array: np.array) -> Tuple[float, float]:
    """Detect possible instrument or recording problems and physical but unusual behavior using higher-moment statistics

    The skewness and kurtosis of the fields are computed for the entire record. The record is hard flagged when the
    skewness is outside the range (-2, 2) or the kurtosis is outside the range (1, 8).
    """
    return skew(array, nan_policy='omit'), kurtosis(array, fisher=False, nan_policy='omit')


def discontinuities(df: pd.DataFrame) -> Tuple[float, float]:
    """Detect discontinuities in the data using the Haar transform

    The transform is computed for a series of moving windows of width L1 and then normalized by the smaller of the
    standard deviation for the entire record and one-fourth the range for the entire record. The record is hard flagged
    if the absolute value of any single normalized transform exceeds 3 and soft flagged at 2. To identify coherent
    changes over the window width L1 in the intensity of the fluctuations, we compute the variance for each half-window
    and then compute the difference normalized by the variance over the entire record. The record is hard flagged if the
    absolute value of any single normalized transform exceeds 3 and soft flagged at 2.
    """
    L1_SECONDS: int = 300
    L1: int = L1_SECONDS * FREQUENCY_HZ

    col: str = df.columns[0]

    record_std: float = np.nanstd(df[col].values)
    record_range: float = np.nanmax(df[col].values) - np.nanmin(df[col].values)
    norm: float = min(record_std, record_range / 4.0)

    haar_mean_max: float = 0.0
    haar_variance_max: float = 0.0

    for window_position in range(0, len(df) - L1 + 1, L1 // 4):
        half_point = L1 // 2
        half1 = slice(0, half_point)
        half2 = slice(half_point, L1)

        if df[half1][col].count() < L1 // 4 or df[half2][col].count() < L1 // 4:
            continue

        half1_mean = np.nanmean(df[half1][col].values)
        half2_mean = np.nanmean(df[half2][col].values)

        half1_variance = np.nanvar(df[half1][col].values)
        half2_variance = np.nanvar(df[half2][col].values)

        haar_mean_max = max(haar_mean_max, abs((half2_mean - half1_mean) / norm))
        haar_variance_max = max(haar_variance_max, abs((half2_variance - half1_variance) / (record_std * record_std)))

    return haar_mean_max, haar_variance_max


def nonstationary(df: pd.DataFrame) -> Tuple[float, float, float, float]:
    """Identify nonstaionarity of horizontal wind

    The wind speed reduction is defined as the ratio of the speed of the vector averaged wind to the averaged
    instantaneous speed. When this ratio falls below 0.9, there is some cancellation in the vector average of the wind
    components and a soft flag is raised.
    The alongwind relative nonstationarity is calculated using linear regression to estimate the difference in the
    alongwind component between the beginning and end of the record. This difference normalized by the record mean of
    the alongwind component is used to compute the relative nonstationarity.
    The crosswind relative nonstationarity is computed from the difference based on the regression of the crosswind
    component.
    The flow is classified as nonstationary if RNu, RNv, or RNS > 0.50
    """
    # Wind speed reduction
    df.dropna(inplace=True)

    u = df['u'].values
    v = df['v'].values

    vector_average = math.sqrt(u.mean() * u.mean() + v.mean() * v.mean())
    instant_average = np.sqrt(u * u + v * v).mean()
    wind_speed_reduction = vector_average / instant_average

    # Relative nonstationarity
    time = np.array(df.index)
    su = np.polyfit(time, u, 1)
    du = su[0] * (time[-1] - time[0])
    sv = np.polyfit(time, v, 1)
    dv = sv[0] * (time[-1] - time[0])

    rnu = du / u.mean()
    rnv = dv / u.mean()
    rns = math.sqrt(du * du + dv * dv) / u.mean()

    return wind_speed_reduction, rnu, rnv, rns
