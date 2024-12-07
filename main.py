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
L1_freq = 1575.42e6
k = 40.308193
deg2semi =  1./180.
semi2rad =  np.pi;
deg2rad  =  np.pi/180.


def find_cords(filename: str, pos_lat: float, pos_long: float) -> tuple:
    lat_north = 90.
    lon_west = 180.
    lat_south = -lat_north
    lon_east = -lon_west

    with open(filename, 'r') as file:
        for line in file:
            if 'LAT/LON1/LON2/DLON/H' in line:
                LAT, LON1, LON2, DLON, H = tuple(float(match) for match in re.findall(r'-?\d+\.\d+', line.strip()))
                if LAT - pos_lat > 0 and lat_north > LAT: lat_north = LAT
                elif LAT - pos_lat < 0 and lat_south < LAT: lat_south = LAT

        for long in range(int(LON1), int(LON2) + 1, int(DLON)):
            if long - pos_long < 0 and np.fabs(pos_long - lon_west) > np.fabs(pos_long - long): lon_west = long
            elif long - pos_long > 0 and np.fabs(lon_east - pos_long) > np.fabs(pos_long - long): lon_east = long

        return lat_north, lat_south, float(lon_west), float(lon_east)


def find_TEC_delays(filename: str, latitude: float, longitude: float) -> dict:
    TEC_delays = dict()

    epoch = 0
    line_number = 0
    with open(filename, 'r') as file:
        lines = file.readlines()


    for line in lines:
        if 'EPOCH OF CURRENT MAP' in line:
            epoch = tuple(int(match) for match in re.findall(r'\d+', line.strip()))
            epoch_str = f"{epoch[2]}.{epoch[1]}.{epoch[0]} - {epoch[3]}"
        if 'START OF RMS MAP' in line: break

        if 'LAT/LON1/LON2/DLON/H' in line:
            LAT, LON1, LON2, DLON, H = tuple(float(match) for match in re.findall(r'-?\d+\.\d+', line.strip()))
            if LAT == latitude:
                delays_on_lat = [match for match in re.findall(r'\d+',''.join(lines[line_number + 1 : line_number + 6]))]
                if epoch_str not in TEC_delays.keys(): TEC_delays[epoch_str] = int(delays_on_lat[int(np.fabs(longitude - LON1) / DLON)])
        line_number += 1

    return TEC_delays


def io_delay(lat: float, long: float, delays: list, points: tuple) -> float:
    phi_pp = lat
    lambda_pp = long
    tau_v = delays

    lambda_1, lambda_2, phi_1, phi_2 = points

    x_pp = (lambda_pp - lambda_1) / (lambda_2 - lambda_1)
    y_pp = (phi_pp - phi_1) / (phi_2 - phi_1)

    W = [x_pp * y_pp, (1 - x_pp) * y_pp, (1 - x_pp) * (1 - y_pp), x_pp * (1 - y_pp)]

    for k in range(3):
        tau_vpp = W[k] * tau_v[k]

    return TECU_to_meters(tau_vpp)


def klobuchar(latitude, longitude, elev, azim, tow, alpha, beta):
    fi = latitude
    lamb = longitude

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


def io_delays_by_epoch(filename, LAT, LONG):
    lat_north, lat_south, lon_west, lon_east = find_cords(filename, float(LAT), float(LONG))

    delays_on_NW = find_TEC_delays(filename, lat_north, lon_west)
    delays_on_NE = find_TEC_delays(filename, lat_north, lon_east)
    delays_on_SW = find_TEC_delays(filename, lat_south, lon_west)
    delays_on_SE = find_TEC_delays(filename, lat_south, lon_east)

    io_delays = dict()
    for epoch in delays_on_NW.keys():
        delays = []
        delays.append(delays_on_NE[epoch])
        delays.append(delays_on_NW[epoch])
        delays.append(delays_on_SW[epoch])
        delays.append(delays_on_SE[epoch])
        io_delays[epoch] = io_delay(float(LAT), float(LONG), delays, (lon_west, lon_east, lat_south, lat_north))
    return io_delays


def TECU_to_meters(TECU):
    return -k * (TECU / L1_freq ** 2)

if __name__ == '__main__':
    LAT = LATITUDE[:-1] if LATITUDE[-1] == 'N' else '-' + LATITUDE[:-1]
    LONG = LONGITUDE[:-1] if LONGITUDE[-1] == 'E' else '-' + LONGITUDE[:-1]

    print(io_delays_by_epoch('data/igsg0010.18i', LAT, LONG))
    print(io_delays_by_epoch('data/igrg0010.18i', LAT, LONG))
