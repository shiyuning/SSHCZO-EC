#!/usr/bin/env python3
import pandas as pd
import numpy as np
import math
from site_parameters import *

RD = 287.05
LV = 2503000.0
C_AIR = 1004.0

def calculate_fluxes(unit_k, pressure_kpa, df):
    # Transformation to the natural wind coordinate system
    # (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 4)
    # Calculate unit vector i and j
    u_vector = np.array([np.mean(df[U]), np.mean(df[V]), np.mean(df[W])])
    unit_j = np.cross(unit_k, u_vector)
    unit_j /= np.linalg.norm(unit_j)
    unit_i = np.cross(unit_j, unit_k)

    # Calculate wind direction and speed
    ce = unit_i[0] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])
    se = unit_i[1] / math.sqrt(unit_i[0] * unit_i[0] + unit_i[1] * unit_i[1])

    eta = math.acos(ce) / math.pi * 180.0
    eta = 360.0 - eta if se < 0 else eta
    eta = (360.0 - eta + CSAT3_AZIMUTH) % 360.0
    print(f'  Wind direction: {eta:.3f} degree')

    u_bar = np.dot(unit_i, u_vector)
    v_bar = np.dot(unit_j, u_vector)

    u = df[f'{U}_fluct'] * unit_i[0] + df[f'{V}_fluct'] * unit_i[1] + df[f'{W}_fluct'] * unit_i[2]
    v = df[f'{U}_fluct'] * unit_j[0] + df[f'{V}_fluct'] * unit_j[1] + df[f'{W}_fluct'] * unit_j[2]
    w = df[f'{U}_fluct'] * unit_k[0] + df[f'{V}_fluct'] * unit_k[1] + df[f'{W}_fluct'] * unit_k[2]

    ust = math.sqrt(math.sqrt(np.mean(u * w) * np.mean(u * w) + np.mean(v * w) * np.mean(v * w)))
    print(f'  Friction velocity: {ust:.3f} m/s')

    rho = 1.20

    f0 = np.mean(w * df[f'{CO2}_fluct'])
    e0 = np.mean(w * df[f'{H2O}_fluct']) / 1000.0
    sh0 = rho * C_AIR * np.mean(w * df[f'{TS}_fluct'])

    if pressure_kpa is not None:
        p = pressure_kpa * 1000.0
        ta = np.mean(df[TS]) + 273.15
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
        print(f'  H2O flux before WPL correction = {e0 * LV:.2f} kg/m2')
        print(f'  H2O flux after WPL correction = {e * LV:.2f} kg/m2')
        print(f'  CO2 flux before WPL correction = {f0 * 1000.0 / 44.0:.2f} umol/m2/s')
        print(f'  CO2 flux after WPL correction = {fc * 1000.0 / 44.0:.2f} umol/m2/s')
    else:
        sh = sh0
        e = e0
        fc = f0
        print(f'  Sensible heat flux: {sh:.2f} umol/m2/s')
        print(f'  Latent heat flux: {e * LV:.2f} W/m2')
        print(f'  CO2 flux: {fc * 1000.0 / 44.0:.2f} W/m2')

    fluxes = pd.DataFrame(
        {
            'UST': ust,
            'TA': np.mean(df[TS]),
            'WD': eta,
            'WS': math.sqrt(u_bar * u_bar + v_bar * v_bar),
            'VWS': np.dot(unit_k, u_vector),
            'FC': fc * 1000.0 / 44.0,
            'H': sh,
            'LE': e * LV,
            'CO2': np.mean(df[CO2]),
            'H2O': np.mean(df[H2O]),
        },
        index=[df[TIME].iloc[0]])


    return fluxes
