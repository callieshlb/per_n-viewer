"""geod_dist of E94 overlay_export.py (copied verbatim: importing that script would run it)."""
import numpy as np
def geod_dist(z, xi1, xi2):
    al, be = np.angle(xi1), np.angle(xi2); half = (be - al)/2
    if abs(np.cos(half)) < 1e-9: return np.arcsinh(2*np.abs((z*np.exp(-1j*al)).imag)/(1 - np.abs(z)**2))
    C = np.exp(1j*(al + be)/2)/np.cos(half); rho = abs(np.tan(half))
    return np.arcsinh(np.abs(np.abs(z - C)**2 - rho**2)/(rho*(1 - np.abs(z)**2)))
