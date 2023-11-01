#!/usr/bin/env python3
import math
import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew
from itertools import groupby
from operator import itemgetter
from site_parameters import *

SECONDS_IN_MINUTE = 60.0
VARIABLES = ['u', 'v', 'w', 'tsonic', 'co2', 'h2o']

filter = {
    'u': ANEMOMETER_FILTER,
    'v': ANEMOMETER_FILTER,
    'w': ANEMOMETER_FILTER,
    'tsonic': ANEMOMETER_FILTER,
    'co2': IRGA_FILTER,
    'h2o': IRGA_FILTER,
}

def quality_control(unit_i, unit_j, unit_k, df):
    """Quality control for eddy covariance flux data
    Vickers, D., and L. Mahrt, 1997: Quality control and flux sampling problems for tower and aircraft data.
    J. Atmos. Oceanic tech., 14, 512-526
    """
    diagnostics = {}

    diagnostics['3dua'], diagnostics['irga'] = instrument(df)

    for var in VARIABLES:
        if len(df[filter[var]]) == 0:
            diagnostics[f'{var}_spikes'] = None
        else:
            diagnostics[f'{var}_spikes'] = spikes(df[filter[var]], var)

    # Rotate to the natural wind coordinate system after de-spiking
    if len(df[ANEMOMETER_FILTER]) > 0:
        _u = df[ANEMOMETER_FILTER]['u'].values
        _v = df[ANEMOMETER_FILTER]['v'].values
        _w = df[ANEMOMETER_FILTER]['w'].values
        u = _u * unit_i[0] + _v * unit_i[1] + _w * unit_i[2]
        v = _u * unit_j[0] + _v * unit_j[1] + _w * unit_j[2]
        w = _u * unit_k[0] + _v * unit_k[1] + _w * unit_k[2]
        [df.loc[ANEMOMETER_FILTER, 'u'], df.loc[ANEMOMETER_FILTER, 'v'], df.loc[ANEMOMETER_FILTER, 'w']] = [u, v, w]

    for var in VARIABLES:
        if len(df[filter[var]]) == 0:
            diagnostics[f'{var}_resolution'] = diagnostics[f'{var}_dropouts'] = diagnostics[f'{var}_extreme_dropouts'] = None
        else:
            diagnostics[f'{var}_resolution'], diagnostics[f'{var}_dropouts'], diagnostics[f'{var}_extreme_dropouts'] = amplitude_resolution_dropouts(df[filter[var]], var)

    if len(df[ANEMOMETER_FILTER]) == 0:
        diagnostics['u_max'] = diagnostics['v_max'] = diagnostics['w_max'] = diagnostics['tsonic_min'] = diagnostics['tsonic_max'] = None
    else:
        diagnostics['u_max'] = abs(df[ANEMOMETER_FILTER]['u'].values).max()
        diagnostics['v_max'] = abs(df[ANEMOMETER_FILTER]['v'].values).max()
        diagnostics['w_max'] = abs(df[ANEMOMETER_FILTER]['w'].values).max()
        diagnostics['tsonic_min'] = df[ANEMOMETER_FILTER]['tsonic'].values.min()
        diagnostics['tsonic_max'] = df[ANEMOMETER_FILTER]['tsonic'].values.max()
    if len(df[IRGA_FILTER]) == 0:
        diagnostics['co2_min'] = diagnostics['co2_max'] = diagnostics['h2o_min'] = diagnostics['h2o_max'] = None
    else:
        diagnostics['co2_min'] = df[IRGA_FILTER]['co2'].values.min()
        diagnostics['co2_max'] = df[IRGA_FILTER]['co2'].values.max()
        diagnostics['h2o_min'] = df[IRGA_FILTER]['h2o'].values.min()
        diagnostics['h2o_max'] = df[IRGA_FILTER]['h2o'].values.max()

    for var in VARIABLES:
        if len(df[filter[var]]) == 0:
            diagnostics[f'{var}_skewness'] = diagnostics[f'{var}_kurtosis'] = None
        else:
            detrend(df, var)
            diagnostics[f'{var}_skewness'], diagnostics[f'{var}_kurtosis'] = higher_moment_statistics(df[filter[var]], var)

    for var in VARIABLES:
        if len(df[filter[var]]) == 0:
            diagnostics[f'{var}_haar_mean'] = diagnostics[f'{var}_haar_variance'] = None
        else:
            diagnostics[f'{var}_haar_mean'], diagnostics[f'{var}_haar_variance'] = discontinuities(df[filter[var]], var)

    if len(df[ANEMOMETER_FILTER]) == 0:
        diagnostics['wind_speed_reduction'] = diagnostics['rnu'] = diagnostics['rnv'] = diagnostics['rns'] = None
    else:
        diagnostics['wind_speed_reduction'], diagnostics['rnu'], diagnostics['rnv'],  diagnostics['rns'] = nonstationary(df[ANEMOMETER_FILTER])

    return diagnostics


def instrument(df):
    """Instrument diagnostics
    If available records for current time period is less than a threshold, or is more than 18000, a flag is raised
    """
    return float(len(df[ANEMOMETER_FILTER])) / (AVERAGING_PERIOD_MINUTES * SECONDS_IN_MINUTE * FREQUENCY_HZ), \
        float(len(df[IRGA_FILTER])) / (AVERAGING_PERIOD_MINUTES * SECONDS_IN_MINUTE * FREQUENCY_HZ)


def spikes(df, var):
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
    CONSECUTIVE_SPIKES = 3
    L1_SECONDS = 300
    window_position = 0
    std_threshold = 4.5
    threshold_increment = 0.1
    pass_spikes = []
    all_spikes = []

    L1 = min(int(L1_SECONDS * FREQUENCY_HZ), len(df))

    while True:
        while window_position + L1 <= len(df):
            sub_df = df.loc[window_position:window_position + L1, :]
            array = sub_df[var].values

            # Get the index of spikes
            pass_spikes += sub_df.index[abs(array - array.mean()) > std_threshold * array.std()].tolist()

            window_position += 1

        # Remove duplicate spikes from this pass
        pass_spikes = np.sort(np.unique(pass_spikes))

        if len(pass_spikes) > 0:
            # Filter out consecutive spikes
            # df[TIME].values * 1E9 converts to seconds, and (df[TIME].values - df[TIME0].value) / 1E9 * FREQUENCY_HZ)
            # should be equivalent to df.index, but the latter does not take into account the time gaps between records.
            consecutive_spikes = []
            for c in [list(map(itemgetter(0), g)) for _, g in groupby(enumerate((df.loc[pass_spikes, 'time'].values.astype(int) - df.loc[0, 'time'].value) / 1E9 * FREQUENCY_HZ), lambda x: x[0] - x[1])]:
                if len(c) > CONSECUTIVE_SPIKES: consecutive_spikes += c
            if len(consecutive_spikes) > 0:
                pass_spikes = np.delete(pass_spikes, consecutive_spikes)

        if len(pass_spikes) > 0:
            all_spikes += pass_spikes.tolist()

            # Replace spikes with linear interpolation
            spike_filter = np.full(len(df), False)
            spike_filter[pass_spikes] = True
            df.loc[spike_filter, var] = np.interp(df.loc[spike_filter, 'time'].values.astype(float), df.loc[~spike_filter, 'time'].values.astype(float), df.loc[~spike_filter, var].values)

            window_position = 0
            std_threshold += threshold_increment
            threshold_increment += 0.1
            pass_spikes = []
        else:
            return float(len(np.unique(all_spikes))) / float(len(df))


def amplitude_resolution_dropouts(df, var):
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
    L1 = 1000
    N_BINS = 100.0

    L1 = min(L1, len(df))
    empty_bins = 0
    window_position = 0
    dropout_ratio = 0
    extreme_dropout_ratio = 0

    while (window_position + L1 <= len(df)):
        array = df.loc[window_position:window_position + L1, var].values
        window_mean = array.mean()
        window_std = array.std()
        distribution = min(7.0 * window_std, array.ptp())
        bin_edges = np.arange(
            window_mean - distribution / 2.0,
            window_mean + distribution / 2.0 + distribution / N_BINS,
            distribution / N_BINS,
        )

        hist, _ = np.histogram(array, bins=bin_edges)
        empty_bins = max(empty_bins, float(np.count_nonzero(hist==0)) / N_BINS)

        #https://stackoverflow.com/questions/6352425/whats-the-most-pythonic-way-to-identify-consecutive-duplicates-in-a-list
        bins = np.digitize(array, bin_edges)
        grouped_bins = [(k, sum(1 for _ in g)) for k, g in groupby(bins)]

        max_dropouts = max(grouped_bins, key=lambda x: x[1])
        dropout_ratio = max(dropout_ratio, float(max_dropouts[1]) / float(len(df)))
        if  (max_dropouts[0] < 10 or max_dropouts[0] > 90):
            extreme_dropout_ratio = max(extreme_dropout_ratio, float(max_dropouts[1]) / float(len(df)))

        window_position += L1 / 2

    return empty_bins, dropout_ratio, extreme_dropout_ratio


def detrend(df, var):
    """Linear de-trend
    """
    time_array = df[filter[var]]['time'].values.astype(float) / 1.0E9 # Convert to seconds
    time_array -= time_array[0]

    s = np.polyfit(time_array, df[filter[var]][var].values, 1)

    trend = s[0] * time_array + s[1]
    df.loc[filter[var], f'{var}_fluct'] = df.loc[filter[var], var].values - trend


def higher_moment_statistics(df, var):
    """Detect possible instrument or recording problems and physical but unusual behavior using higher-moment statistics
    The skewness and kurtosis of the fields are computed for the entire record. The record is hard flagged when the
    skewness is outside the range (-2, 2) or the kurtosis is outside the range (1, 8).
    """
    return skew(df[f'{var}_fluct'].values), kurtosis(df[f'{var}_fluct'].values, fisher=False)


def discontinuities(df, var):
    """Detect discontinuities in the data using the Haar transform
    The transform is computed for a series of moving windows of width L1 and then normalized by the smaller of the
    standard deviation for the entire record and one-fourth the range for the entire record. The record is hard flagged
    if the absolute value of any single normalized transform exceeds 3 and soft flagged at 2. To identify coherent
    changes over the window width L1 in the intensity of the fluctuations, we compute the variance for each half-window
    and then compute the difference normalized by the variance over the entire record. The record is hard flagged if the
    absolute value of any single normalized transform exceeds 3 and soft flagged at 2.
    """
    L1_SECONDS = 300
    L1 = min(int(L1_SECONDS * FREQUENCY_HZ), len(df))

    record_std = df[var].values.std()
    record_range = df[var].values.ptp()
    normal = min(record_std, record_range / 4.0)

    window_position = 0
    haar_mean_max = 0
    haar_variance_max = 0

    while (window_position + L1 <= len(df)):
        half_point = int(L1 / 2)

        half1_mean = df.loc[0:half_point, var].values.mean()
        half2_mean = df.loc[half_point:L1, var].values.mean()

        half1_variance = df.loc[0:half_point, var].values.var()
        half2_variance = df.loc[half_point:L1, var].values.var()

        haar_mean_max = max(haar_mean_max, abs((half2_mean - half1_mean) / normal))
        haar_variance_max = max(haar_variance_max, abs((half2_variance - half1_variance) / (record_std * record_std)))

        window_position += L1 / 4

    return haar_mean_max, haar_variance_max


def nonstationary(df):
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
    u = df['u'].values
    v = df['v'].values
    vector_average = math.sqrt(u.mean() * u.mean() + v.mean() * v.mean())
    instant_average = np.sqrt(u * u + v * v).mean()
    wind_speed_reduction = vector_average / instant_average

    # Relative nonstationarity
    time_array = df['time'].values.astype(float) / 1.0E9 # Convert to seconds
    time_array -= time_array[0]
    su = np.polyfit(time_array, u, 1)
    du = su[0] * (time_array[-1] - time_array[0])
    sv = np.polyfit(time_array, v, 1)
    dv = sv[0] * (time_array[-1] - time_array[0])

    rnu = du / u.mean()
    rnv = dv / u.mean()
    rns = math.sqrt(du * du + dv * dv) / u.mean()

    return wind_speed_reduction, rnu, rnv, rns
