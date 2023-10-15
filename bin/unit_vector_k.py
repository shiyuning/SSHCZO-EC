#!/usr/bin/env python3
import numpy as np

def unit_vector_k(u, v, w):
    """Determines unit vector k (parallel to new z-axis)
    k: unit vector parallel to new coordinate z axis
    b[0]: instrument offset in w1
    Adapted from Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 3, Section 9
    """
    su = np.mean(u)
    sv = np.mean(v)
    sw = np.mean(w)
    suv = np.mean(u * v)
    suw = np.mean(u * w)
    svw = np.mean(v * w)
    su2 = np.mean(u * u)
    sv2 = np.mean(v * v)

    h = np.array([
        [1, su, sv],
        [su, su2, suv],
        [sv, suv, sv2]
    ])
    g = np.array([sw, suw, svw])

    b = np.linalg.lstsq(h, g, rcond=None)[0]

    # Determine unit vector k
    unit_k = np.zeros(b.size)
    unit_k[2] = 1.0 / (1.0 + b[1] * b[1] + b[2] * b[2])
    unit_k[0] = -b[1] * unit_k[2]
    unit_k[1] = -b[2] * unit_k[2]
    unit_k /= np.linalg.norm(unit_k)

    print('\nUnit vector k of planar fit coordinate:')
    print(f'  k_vector = [{unit_k[0]:.3f}, {unit_k[1]:.3f}, {unit_k[2]:.3f}]')
    print(f'  b0 = {b[0]:.3f}')

    return unit_k, b[0]
