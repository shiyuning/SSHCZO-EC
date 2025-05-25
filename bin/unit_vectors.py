import numpy as np

def unit_vector_k(u: np.array, v: np.array, w: np.array) -> np.array:
    """Determines unit vectors
    Adapted from Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 3, Section 9
    """
    u_bar = np.nanmean(u)
    v_bar = np.nanmean(v)
    w_bar = np.nanmean(w)
    uv_bar = np.nanmean(u * v)
    uw_bar = np.nanmean(u * w)
    vw_bar = np.nanmean(v * w)
    u2_bar = np.nanmean(u * u)
    v2_bar = np.nanmean(v * v)

    h = np.array([
        [1, u_bar, v_bar],
        [u_bar, u2_bar, uv_bar],
        [v_bar, uv_bar, v2_bar]
    ])
    g = np.array([w_bar, uw_bar, vw_bar])

    b = np.linalg.lstsq(h, g, rcond=None)[0]

    # Determine unit vector k
    unit_k = np.zeros(b.size)
    unit_k[2] = 1.0 / (1.0 + b[1] * b[1] + b[2] * b[2])
    unit_k[0] = -b[1] * unit_k[2]
    unit_k[1] = -b[2] * unit_k[2]
    unit_k /= np.linalg.norm(unit_k)

    print(f'\nUnit vector k of planar fit coordinate: [{unit_k[0]:.3f}, {unit_k[1]:.3f}, {unit_k[2]:.3f}]')

    return unit_k


def unit_vector_ij(u: np.array, v: np.array, w: np.array, unit_k: np.array) -> tuple[np.array, np.array]:
    unit_j = np.cross(unit_k, [np.nanmean(u), np.nanmean(v), np.nanmean(w)])
    unit_j /= np.linalg.norm(unit_j)
    unit_i = np.cross(unit_j, unit_k)

    return unit_i, unit_j
