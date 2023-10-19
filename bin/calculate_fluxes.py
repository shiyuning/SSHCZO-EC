#!/usr/bin/env python3
import numpy as np
import math
from site_parameters import *

RD = 287.05
LV = 2503000.0
C_AIR = 1004.0

def calculate_fluxes(unit_i, pressure_pa, tair_c, df):
    # Calculate wind direction and speed
    ce = unit_i[0] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])
    se = unit_i[1] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])

    eta = math.acos(ce) / math.pi * 180.0
    eta = 360.0 - eta if se < 0 else eta
    eta = (360.0 - eta + CSAT3_AZIMUTH) % 360.0
    print(f'  Wind direction: {eta:.3f} degree')

    u = df['u_fluct'].values
    v = df['v_fluct'].values
    w = df['w_fluct'].values

    ust = math.sqrt(math.sqrt(np.mean(u * w) * np.mean(u * w) + np.mean(v * w) * np.mean(v * w)))
    print(f'  Friction velocity: {ust:.3f} m/s')

    rho = 1.20

    f0 = (w * df['co2_fluct'].values).mean()
    e0 = (w * df['h2o_fluct'].values).mean() / 1000.0
    sh0 = rho * C_AIR * (w * df['tsonic_fluct'].values).mean()

    if pressure_pa is not None:
        ta = tair_c + 273.15 if tair_c is not None else df['tsonic'].values.mean() + 273.15
        rho_v = df['h2o'].values.mean() / 1000.0
        rho_c = df['co2'].values.mean()

        vp = rho_v * RD * ta / 0.622
        q = 0.622 * vp / (pressure_pa - 0.378 * vp)
        rho = pressure_pa / (RD * (1.0 + 0.608 * q) * ta)
        rho_d = rho - rho_v

        sh = rho * C_AIR * (w * df['tsonic_fluct'].values).mean()

        # Sonic correction below is commented out because it is already done in LI-7500
        # H = H + rho * C_AIR * (-0.51 * ta * np.mean(w * df['h2o_fluct'] / 1000.0) / rho

        e = (1.0 + 1.6077 * rho_v / rho_d) * (e0 + sh / rho / C_AIR * rho_v / ta)
        fc = f0 + 1.6077 * e / rho_d * rho_c / (1.0 + 1.6077 * (rho_v / rho_d)) + sh / rho / C_AIR * rho_c / ta

        print(f'  Air pressure = {pressure_pa:.2f} Pa')
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
        'WS': math.sqrt(df['u'].values.mean() * df['u'].values.mean() + df['v'].values.mean() * df['v'].values.mean()),
        'FC': fc * 1000.0 / 44.0,
        'H': sh,
        'LE': e * LV,
        'CO2': df['co2'].values.mean(),
        'H2O': df['h2o'].values.mean(),
        'PA': pressure_pa if pressure_pa is not None else np.nan,
        'T_SONIC': df['tsonic'].values.mean(),
        'TA': tair_c if tair_c is not None else np.nan,
    }

    return fluxes
