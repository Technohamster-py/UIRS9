'''
Created on Nov 10 2024
@author: Khomchenko A.A.
'''

import re
import numpy as np

# Config
LATITUDE = '55N'
LONGITUDE = '83E'

c = 2.99792458e8
deg2semi =  1./180.
semi2rad =  np.pi;
deg2rad  =  np.pi/180.


def find_cords(filename: str, pos_lat: float, pos_long: float) -> list:
    with open(filename, 'r') as file:
        for line in file:
            if 'LAT/LON1/LON2/DLON/H' in line:
                lat, lon1, lon2, dlon, h = tuple(float(match) for match in re.findall(r'-?\d+\.\d+', line.strip()))

def io_delay(lat, long, delays):
    phi_pp = lat
    lambda_pp = long
    tau_v = delays

    W = list
    W[0] = x_pp * y_pp
    W[1] = (1 - x_pp) * y_pp
    W[2] = (1 - x_pp) * (1 - y_pp)
    W[3] = x_pp * (1 - y_pp)

    for k in range(3):
        tau_vpp = W[k] * tau_v[k]

    return tau_vpp

def klobuchar(fi, lamb, elev, azim, tow, alpha, beta):
    a = azim * deg2rad
    e = elev * deg2semi

    psi = 0.0137 / (e + 0.11) - 0.022

    lat_i = fi * deg2semi + psi * np.cos(a)
    if lat_i > 0.416: lat_i = 0.416
    elif lat_i < -0.416: lat_i = -0.416

    long_i = lamb * deg2semi + (psi * np.sin(a) / np.cos(lat_i * semi2rad))

    lat_m = lat_i + 0.064 * np.cos((long_i - 1.617) * semi2rad)

    t = 4.32e4 * long_i + tow
    t = t % 86400.
    if t > 86400.: t = t - 86400.

    sF = 1. + 16. * (0.53-e)**3

    PER = beta[0] + beta[1] * lat_m + beta[2] * lat_m ** 2 + beta[3] * lat_m ** 3
    if PER < 72000: PER = 72000.

    x = 2. * np.pi * (t - 50400.) / PER

    AMP = alpha[0] + alpha[1] * lat_m + alpha[2] * lat_m ** 2 + alpha[3] * lat_m ** 3
    if AMP < 0.: AMP = 0.

    if np.fabs(x) > 1.57: dIon = sF * (5.e-9)
    else: dIon = sF * (5.e-9 + AMP * (1. - x*x/2. + x*x*x*x/24.))

    return c * dIon


if __name__ == '__main__':
    LAT = LATITUDE[:-1] if LATITUDE[-1] == 'N' else '-' + LATITUDE[:-1]
    LONG = LONGITUDE[:-1] if LONGITUDE[-1] == 'E' else '-' + LONGITUDE[:-1]
    find_cords('data/igsg0010.18i')

