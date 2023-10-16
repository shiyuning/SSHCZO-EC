#!/usr/bin/env python3
import numpy as np
import math
from site_parameters import *

RD = 287.05
LV = 2503000.0
C_AIR = 1004.0

def calculate_fluxes(unit_i, pressure_kpa, ta_c, df):
    # Calculate wind direction and speed
    ce = unit_i[0] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])
    se = unit_i[1] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])

    eta = math.acos(ce) / math.pi * 180.0
    eta = 360.0 - eta if se < 0 else eta
    eta = (360.0 - eta + CSAT3_AZIMUTH) % 360.0
    print(f'  Wind direction: {eta:.3f} degree')

    u = df[f'{U}_fluct']
    v = df[f'{V}_fluct']
    w = df[f'{W}_fluct']

    ust = math.sqrt(math.sqrt(np.mean(u * w) * np.mean(u * w) + np.mean(v * w) * np.mean(v * w)))
    print(f'  Friction velocity: {ust:.3f} m/s')

    rho = 1.20

    f0 = np.mean(w * df[f'{CO2}_fluct'])
    e0 = np.mean(w * df[f'{H2O}_fluct']) / 1000.0
    sh0 = rho * C_AIR * np.mean(w * df[f'{TS}_fluct'])

    if pressure_kpa is not None:
        p = pressure_kpa * 1000.0
        ta = ta_c + 273.15 if ta_c is not None else np.mean(df[TS]) + 273.15
        rho_v = np.mean(df[H2O]) / 1000.0
        rho_c = np.mean(df[CO2])

        vp = rho_v * RD * ta / 0.622
        q = 0.622 * vp / (p - 0.378 * vp)
        rho = p / (RD * (1.0 + 0.608 * q) * ta)
        rho_d = rho - rho_v

        sh = rho * C_AIR * np.mean(w * df[f'{TS}_fluct'])

        # Sonic correction below is commented out because it is already done in LI-7500
        # H = H + rho * C_AIR * (-0.51 * ta * np.mean(w * df['h2o_fluct'] / 1000.0) / rho

        e = (1.0 + 1.6077 * rho_v / rho_d) * (e0 + sh / rho / C_AIR * rho_v / ta)
        fc = f0 + 1.6077 * e / rho_d * rho_c / (1.0 + 1.6077 * (rho_v / rho_d)) + sh / rho / C_AIR * rho_c / ta

        print(f'  Air pressure = {p:.2f} Pa')
        print(f'  Air density = {rho:.2f} kg/m3')
        print(f'  Sensible heat flux = {sh:.2f} W/m2')
        print(f'  Latent flux before WPL correction = {e0 * LV:.2f} W/m2')
        print(f'  Latent flux after WPL correction = {e * LV:.2f} W/m2')
        print(f'  CO2 flux before WPL correction = {f0 * 1000.0 / 44.0:.2f} umol/m2/s')
        print(f'  CO2 flux after WPL correction = {fc * 1000.0 / 44.0:.2f} umol/m2/s')
    else:
        sh = sh0
        e = e0
        fc = f0
        print(f'  Sensible heat flux: {sh:.2f} umol/m2/s')
        print(f'  Latent heat flux: {e * LV:.2f} W/m2')
        print(f'  CO2 flux: {fc * 1000.0 / 44.0:.2f} W/m2')

    fluxes = {
        'USTAR': ust,
        'WD': eta,
        'WS': math.sqrt(np.mean(df[U]) * np.mean(df[U]) + np.mean(df[V]) * np.mean(df[V])),
        'FC': fc * 1000.0 / 44.0,
        'H': sh,
        'LE': e * LV,
        'CO2': np.mean(df[CO2]),
        'H2O': np.mean(df[H2O]),
        'PA': pressure_kpa * 1000.0 if pressure_kpa is not None else np.nan,
        'T_SONIC': np.mean(df[TS]),
        'TA': ta_c if ta_c is not None else np.nan,
    }

    return fluxes
