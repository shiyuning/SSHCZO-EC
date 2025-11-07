import math
import numpy as np
import pandas as pd
from dataclasses import dataclass
from datetime import datetime
from dateutil.relativedelta import relativedelta
from itertools import groupby
from operator import itemgetter
from scipy.stats import kurtosis, skew
from typing import List
from config import AVERAGING_PERIOD_MINUTES, FREQUENCY_HZ
from config import CSAT3_AZIMUTH, ANEMOMETER_FLAGS, ANEMOMETER_FILTER, IRGA_FLAGS, IRGA_FILTER
from config import BADVAL
from quality_control import instrument, spikes, amplitude_resolution_dropouts, detrend, higher_moment_statistics
from quality_control import discontinuities, nonstationary
from unit_vectors import unit_vector_ij

@dataclass
class InstantaneousVariable:
    name: str
    unit: str

    def filter(self, x):
        if self.name in ['u', 'v', 'w', 'tsonic']:
            return ANEMOMETER_FILTER(x)
        else:
            return IRGA_FILTER(x)

INSTANTANEOUS_VARIABLES = (
    InstantaneousVariable('u', 'm/s'),
    InstantaneousVariable('v', 'm/s'),
    InstantaneousVariable('w', 'm/s'),
    InstantaneousVariable('tsonic', 'deg C'),
    InstantaneousVariable('co2', 'g/m3'),
    InstantaneousVariable('h2o', 'mg/m3'),
)

@dataclass
class OutputVariable:
    name: str
    unit: str
    format: str
    flag: List[str]

OUTPUT_VARIABLES = (
    OutputVariable('USTAR', 'm/s', '%.2f', ['u_flag', 'v_flag', 'w_flag']),
    OutputVariable('WD', 'degree', '%.1f', ['u_flag', 'v_flag']),
    OutputVariable('WS', 'm/s', '%.2f', ['u_flag', 'v_flag']),
    OutputVariable('FC', 'umolCO2/m2/s', '%.2f', ['w_flag', 'co2_flag']),
    OutputVariable('H', 'W/m2', '%.2f', ['w_flag', 'tsonic_flag']),
    OutputVariable('LE', 'W/m2', '%.2f', ['w_flag', 'h2o_flag']),
    OutputVariable('CO2', 'Mg CO2/m3', '%.2f', ['co2_flag']),
    OutputVariable('H2O', 'g H2O/m3', '%.2f', ['h2o_flag']),
    OutputVariable('PA', 'kPa', '%.2f', []),
    OutputVariable('T_SONIC', 'deg C', '%.2f', ['tsonic_flag']),
    OutputVariable('TA', 'deg C', '%.2f', []),
)

class EddyCovariance:
    def __init__(self, *, time: datetime, unit_k: np.array, data: pd.DataFrame, pressure_data: pd.DataFrame):
        self.start: datetime = time
        self.end: datetime = time + relativedelta(minutes=AVERAGING_PERIOD_MINUTES)
        self.unit_k: np.array = unit_k

        self.diagnostics: dict = {
            '3dua': np.nan,
            'irga': np.nan,
            **{f'{v.name}_spikes': np.nan for v in INSTANTANEOUS_VARIABLES},
            **{f'{v.name}_resolution': np.nan for v in INSTANTANEOUS_VARIABLES},
            **{f'{v.name}_dropouts': np.nan for v in INSTANTANEOUS_VARIABLES},
            **{f'{v.name}_extreme_dropouts': np.nan for v in INSTANTANEOUS_VARIABLES},
            'u_max': np.nan,
            'v_max': np.nan,
            'w_max': np.nan,
            'tsonic_min': np.nan,
            'tsonic_max': np.nan,
            'co2_min': np.nan,
            'co2_max': np.nan,
            'h2o_min': np.nan,
            'h2o_max': np.nan,
            **{f'{v.name}_skewness': np.nan for v in INSTANTANEOUS_VARIABLES},
            **{f'{v.name}_kurtosis': np.nan for v in INSTANTANEOUS_VARIABLES},
            **{f'{v.name}_haar_mean': np.nan for v in INSTANTANEOUS_VARIABLES},
            **{f'{v.name}_haar_variance': np.nan for v in INSTANTANEOUS_VARIABLES},
            'wind_speed_reduction': np.nan,
            'rnu': np.nan,
            'rnv': np.nan,
            'rns': np.nan,
        }

        self.output: dict = {
            'TIMESTAMP_START': self.start.strftime("%Y%m%d%H%M"),
            'TIMESTAMP_END': self.end.strftime("%Y%m%d%H%M"),
            **{v.name: np.nan for v in OUTPUT_VARIABLES},
        }

        self.data = data[(data['time'] >= self.start) & (data['time'] < self.end)].copy()

        if not self.data.empty:
            self.data['time'] = self.data['time'].values.astype(int) // 1.0E8 - (time - datetime(1970,1,1)).total_seconds() * 10
            self.data['time'] = self.data['time'].astype(int)
            self.data.set_index('time', inplace=True)

            # Re-index the dataframe to make sure it has the right number of rows. Filled rows are labeled with flags
            try:
                self.data = self.data.reindex(list(range(AVERAGING_PERIOD_MINUTES * 60 * FREQUENCY_HZ)))
            except:
                print(f'Warning: raw data from {self.start} to {self.end} may have invalid timestamps.')
                self.data = pd.DataFrame()

        if not self.data.empty:
            self.data['diag'] = self.data['diag'].fillna(ANEMOMETER_FLAGS + IRGA_FLAGS).astype(int)

            if self.data[ANEMOMETER_FILTER].empty:
                self.unit_i, self.unit_j = None, None
            else:
                self.unit_i, self.unit_j = unit_vector_ij(self.data['u'].values, self.data['v'].values, self.data['w'].values, unit_k)

            print(f'\n{self.end.strftime("%Y-%m-%d %H:%M")}')

            # Get air pressure and temperature
            if pressure_data.empty:
                self.pressure = np.nan
                self.tair = np.nan
            else:
                self.pressure = pressure_data[(pressure_data['time'] >= self.start) & (pressure_data['time'] < self.end)].mean()['pressure']
                self.tair = pressure_data[(pressure_data['time'] >= self.start) & (pressure_data['time'] < self.end)].mean()['tair']


    def quality_control(self):
        """Quality control for eddy covariance flux data

        Vickers, D., and L. Mahrt, 1997: Quality control and flux sampling problems for tower and aircraft data.
        J. Atmos. Oceanic tech., 14, 512-526
        """
        # Instrument
        self.diagnostics['3dua'], self.diagnostics['irga'] = instrument(self.data)

        # Spikes
        for v in INSTANTANEOUS_VARIABLES:
            if self.data[v.name].isna().all(): continue

            self.data[v.name], self.diagnostics[f'{v.name}_spikes'] = spikes(self.data[[v.name]].copy())

        # Rotate to the natural wind coordinate system after de-spiking
        if not self.data[ANEMOMETER_FILTER].empty:
            _u = self.data['u'].values
            _v = self.data['v'].values
            _w = self.data['w'].values
            self.data['u'] = _u * self.unit_i[0] + _v * self.unit_i[1] + _w * self.unit_i[2]
            self.data['v'] = _u * self.unit_j[0] + _v * self.unit_j[1] + _w * self.unit_j[2]
            self.data['w'] = _u * self.unit_k[0] + _v * self.unit_k[1] + _w * self.unit_k[2]

        # Resolution problems and dropouts
        for v in INSTANTANEOUS_VARIABLES:
            if self.data[v.name].isna().all(): continue

            self.diagnostics[f'{v.name}_resolution'], self.diagnostics[f'{v.name}_dropouts'], self.diagnostics[f'{v.name}_extreme_dropouts'] = amplitude_resolution_dropouts(self.data[[v.name]])

        # Absolute limits
        if not self.data['u'].isna().all():
            self.diagnostics['u_max'] = np.nanmax(abs(self.data['u'].values))
        if not self.data['v'].isna().all():
            self.diagnostics['v_max'] = np.nanmax(abs(self.data['v'].values))
        if not self.data['w'].isna().all():
            self.diagnostics['w_max'] = np.nanmax(abs(self.data['w'].values))
        if not self.data['tsonic'].isna().all():
            self.diagnostics['tsonic_min'] = np.nanmin(self.data['tsonic'].values)
            self.diagnostics['tsonic_max'] = np.nanmax(self.data['tsonic'].values)
        if not self.data['co2'].isna().all():
            self.diagnostics['co2_min'] = np.nanmin(self.data['co2'].values)
            self.diagnostics['co2_max'] = np.nanmax(self.data['co2'].values)
        if not self.data['h2o'].isna().all():
            self.diagnostics['h2o_min'] = np.nanmin(self.data['h2o'].values)
            self.diagnostics['h2o_max'] = np.nanmax(self.data['h2o'].values)

        ## Higher moment statistics
        for v in INSTANTANEOUS_VARIABLES:
            if self.data[v.name].isna().all():
                self.data[f'{v.name}_'] = np.nan
                continue

            # Linear detrend
            self.data[f'{v.name}_'] = detrend(self.data[[v.name]].copy())
            self.diagnostics[f'{v.name}_skewness'], self.diagnostics[f'{v.name}_kurtosis'] = higher_moment_statistics(self.data[f'{v.name}_'].values)

        # Discontinuities
        for v in INSTANTANEOUS_VARIABLES:
            if self.data[v.name].isna().all(): continue

            self.diagnostics[f'{v.name}_haar_mean'], self.diagnostics[f'{v.name}_haar_variance'] = discontinuities(self.data[[v.name]])

        if not (self.data['u'].isna() & self.data['v'].isna()).all():
            self.diagnostics['wind_speed_reduction'], self.diagnostics['rnu'], self.diagnostics['rnv'],  self.diagnostics['rns'] = nonstationary(self.data[['u', 'v']].copy())


    def calculate_fluxes(self) -> None:
        RD = 287.05
        LV = 2503000.0
        C_AIR = 1004.0
        RHO_AIR = 1.20

        u_ = self.data['u_'].values
        v_ = self.data['v_'].values
        w_ = self.data['w_'].values
        co2_ = self.data['co2_'].values
        h2o_ = self.data['h2o_'].values
        tsonic_ = self.data['tsonic_'].values

        eta = _wind_direction(self.unit_i)

        if not (np.isnan(u_) & np.isnan(w_) & np.isnan(v_)).all():
            ust = math.sqrt(math.sqrt(np.nanmean(u_ * w_) * np.nanmean(u_ * w_) + np.nanmean(v_ * w_) * np.nanmean(v_ * w_)))
            print(f'  Friction velocity: {ust:.3f} m/s')
        else:
            ust = np.nan

        sh = RHO_AIR * C_AIR * np.nanmean(w_ * tsonic_) if not (np.isnan(w_) & np.isnan(tsonic_)).all() else np.nan
        fc = np.nanmean(w_ * co2_) if not (np.isnan(w_) & np.isnan(co2_)).all() else np.nan
        e = np.nanmean(w_ * h2o_) / 1000.0 if not (np.isnan(w_) & np.isnan(h2o_)).all() else np.nan

        if not np.isnan(sh): print(f'  Sensible heat flux: {sh:.2f} umol/m2/s')
        if not np.isnan(e): print(f'  Latent heat flux: {e * LV:.2f} W/m2')
        if not np.isnan(fc): print(f'  CO2 flux: {fc * 1000.0 / 44.0:.2f} W/m2')

        if (not np.isnan(self.pressure) and not np.isnan(self.tair) and
            not self.data['h2o'].isna().all() and not self.data['co2'].isna().all()):
            wpl = True

            ta = self.tair + 273.15
            rho_v = np.nanmean(self.data['h2o'].values) / 1000.0
            rho_c = np.nanmean(self.data['co2'].values)

            vp = rho_v * RD * ta / 0.622
            q = 0.622 * vp / (self.pressure - 0.378 * vp)
            rho = self.pressure / (RD * (1.0 + 0.608 * q) * ta)
            rho_d = rho - rho_v

            print(f'  Air pressure = {self.pressure:.2f} Pa')
            print(f'  Air density = {rho:.2f} kg/m3')

            sh *= rho / RHO_AIR

            e += sh / rho / C_AIR * rho_v / ta
            e *= 1.0 + 1.6077 * rho_v / rho_d

            fc += 1.6077 * e / rho_d * rho_c / (1.0 + 1.6077 * (rho_v / rho_d)) + sh / rho / C_AIR * rho_c / ta

            print(f'  After WPL correction:')
            print(f'  Sensible heat flux after WPL corerction = {sh:.2f} W/m2')
            print(f'  Latent flux after WPL correction = {e * LV:.2f} W/m2')
            print(f'  CO2 flux after WPL correction = {fc * 1000.0 / 44.0:.2f} umol/m2/s')
        else:
            wpl = False

        self.output['USTAR'] = ust
        self.output['WD'] = eta
        self.output['WS'] = np.nan if self.data['u'].isna().all() or self.data['v'].isna().all() else math.sqrt(np.nanmean(self.data['u'].values * self.data['u'].values + self.data['v'].values * self.data['v'].values))
        self.output['FC'] = fc * 1000.0 / 44.0
        self.output['H'] = sh
        self.output['LE'] = e * LV
        self.output['CO2'] = np.nan if self.data['co2'].isna().all() else np.nanmean(self.data['co2'].values)
        self.output['H2O'] = np.nan if self.data['h2o'].isna().all() else np.nanmean(self.data['h2o'].values)
        self.output['PA'] = self.pressure
        self.output['T_SONIC'] = np.nan if self.data['tsonic'].isna().all() else np.nanmean(self.data['tsonic'].values)
        self.output['TA'] = self.tair


    def write_diag_file(self, *, first: bool, fn: str) -> None:
        df = pd.DataFrame.from_dict([{**self.output, **self.diagnostics}])

        for v in OUTPUT_VARIABLES:
            df[v.name] = df[v.name].map(lambda x: v.format % x)

        for col in self.diagnostics:
            df[col] = df[col].map(lambda x: '%.3f' % x)

        df.replace('nan', str(BADVAL), inplace=True)

        with open(fn, 'w' if first else 'a') as f:
            df.to_csv(
                f,
                index=False,
                header=True if first else False,
            )


def _wind_direction(unit_i):
    if unit_i is None: return np.nan

    ce = unit_i[0] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])
    se = unit_i[1] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])

    eta = math.acos(ce) / math.pi * 180.0
    eta = 360.0 - eta if se < 0 else eta
    eta = (360.0 - eta + CSAT3_AZIMUTH) % 360.0
    print(f'  Wind direction: {eta:.3f} degree')

    return eta
