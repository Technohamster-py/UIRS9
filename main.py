'''
Created on Nov 10 2024
@author: Khomchenko A.A.
'''

import re

# Config
LATITUDE = '55N'
LONGITUDE = '83E'


def find_cords(filename: str, pos_lat: float, pos_long: float) -> list:
    with open(filename, 'r') as file:
        for line in file:
            if 'LAT/LON1/LON2/DLON/H' in line:
                lat, lon1, lon2, dlon, h = tuple(float(match) for match in re.findall(r'-?\d+\.\d+', line.strip()))

if __name__ == '__main__':
    LAT = LATITUDE[:-1] if LATITUDE[-1] == 'N' else '-' + LATITUDE[:-1]
    LONG = LONGITUDE[:-1] if LONGITUDE[-1] == 'E' else '-' + LONGITUDE[:-1]
    find_cords('data/igsg0010.18i')
