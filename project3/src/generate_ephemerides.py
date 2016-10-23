from jplephem.spk import SPK #https://pypi.python.org/pypi/jplephem
from astropy.time import Time #http://www.astropy.org/
from datetime import datetime
from collections import OrderedDict
import argparse
import numpy as np
import sys

AU = 149597871 #km per AU
AU_per_year = AU / 365.25 #km_per_AU/days_per_year

M = dict()
M['Sun'] = 2e30
M['Earth'] = 6e24
M['Jupiter'] = 1.9e27
M['Mars'] = 6.6e23
M['Venus'] = 4.9e24
M['Saturn'] = 5.5e26
M['Mercury'] = 3.3e23
M['Uranus'] = 8.8e25
M['Neptune'] = 1.03e26


def setup_parser():
    parser = argparse.ArgumentParser(description='Generate planet positions and velocities')
    parser.add_argument('year', metavar='yyyy', type=int, help='4 digit year for ephemerides')
    parser.add_argument('month', metavar='m', type=int, help='month[1-12] for ephemerides')
    parser.add_argument('day', metavar='d', type=int, help='day[1 - 31] for ephemerides')
    parser.add_argument('--sun_center', action="store_true", help='converts from solar system barycenter to sun center')
    return parser;
    

def parse_input(parser):
    args = parser.parse_args()
    datearg = Time(datetime(args.year, args.month, args.day))
    return datearg, args.sun_center

	
def print_ephemerides(solar_system, sun_center):

    print('#Ephemerides for %s' % solar_system.pop('date'))
    print('%8s %17s %17s %17s %17s %17s %17s %17s' % ('#Body'.ljust(8), 'x (AU)', 'y (AU)', 'z (AU)',
					'vel_x (AU/year)', 'vel_y (AU/year)', 'vel_z (AU/year)', 'mass (solar masses)'))

    if sun_center:
        for k, v in solar_system.items():
    	    solar_system[k] = v - solar_system['Sun']

    for k, v in solar_system.items():
        print(k.ljust(8), end=" ")
        print('%+12.10e %+12.10e %+12.10e ' % (v[0]/AU, v[1]/AU, v[2]/AU), end='')
        print('%+12.10e %+12.10e %+12.10e ' % (v[3]/AU_per_year, v[4]/AU_per_year, v[5]/AU_per_year), end='')
        print('%+12.10e' % (M[k]/M['Sun']))


def get_ephemerides(eph_date):
	#Can be downloaded from ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp
	eph = SPK.open('../input/de430.bsp')

	solar_system = OrderedDict()
	solar_system['date'] = eph_date.datetime.strftime('%Y-%m-%d')

	solar_system['Sun'] = np.hstack(eph[0,10].compute_and_differentiate(eph_date.jd))
	solar_system['Mercury'] = np.hstack(eph[0,1].compute_and_differentiate(eph_date.jd))
	solar_system['Venus'] = np.hstack(eph[0,2].compute_and_differentiate(eph_date.jd))
	solar_system['Earth'] = np.hstack(eph[0,3].compute_and_differentiate(eph_date.jd))
	solar_system['Mars'] = np.hstack(eph[0,4].compute_and_differentiate(eph_date.jd))
	solar_system['Jupiter'] = np.hstack(eph[0,5].compute_and_differentiate(eph_date.jd))	
	solar_system['Saturn'] = np.hstack(eph[0,6].compute_and_differentiate(eph_date.jd))
	solar_system['Uranus'] = np.hstack(eph[0,7].compute_and_differentiate(eph_date.jd))
	solar_system['Neptune'] = np.hstack(eph[0,8].compute_and_differentiate(eph_date.jd))

	return solar_system


def main():
    parser = setup_parser()
    eph_date, sun_center = parse_input(parser)
    solar_system = get_ephemerides(eph_date)		
    print_ephemerides(solar_system, sun_center)
	
if __name__ == '__main__':
	main()

