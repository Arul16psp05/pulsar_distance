#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 16:31:48 2023

@author: arul

Ref: https://stackoverflow.com/questions/68564393/astropy-ra-dec-to-altaz-conversion-in-degrees
"""

import numpy as np
from astropy import units as u
#from astropy.coordinates import SkyCoord
#from datetime import datetime
#j1713 = SkyCoord('17h13m49.5s','+07d47m37.4', frame='icrs')


from astropy.coordinates import EarthLocation,SkyCoord
from astropy.time import Time
#from astropy import units as u
from astropy.coordinates import AltAz

#from zoneinfo import ZoneInfo

#tz = ZoneInfo('Asia/Kolkata')

#from dateutil import tz
#from dateutil.tz import gettz

#dt = Time(datetime(2021, 4, 12, 20, 0, 0, tzinfo=tz))

#dt2 = dt.datetime(2015, 12, 21, 12, 0, tzinfo = NYC) 

#obs_time = Time(datetime.now(tz=ZoneInfo('Asia/Kolkata')))

#datetime.now(tz=gettz('Asia/Kolkata'))

"""
GMRT


latitude = 19.0912

longitude = 74.0432

altitude = 560

"""

def sperical_to_cartesian( psr_ra, psr_dec, dist):

    psr_ra_dec = psr_ra +' '+ psr_dec
    coord = SkyCoord(psr_ra_dec, unit=(u.hourangle, u.deg), distance=dist*u.kpc)
    cart = coord.cartesian
    x, y, z = cart.xyz      # kPC
    

    x = x.value  # lyr            # (x.si).value  ------> m
    y = y.value
    z = z.value
    """
    x = x.to(u.lyr).value  # lyr            # (x.si).value  ------> m
    y = y.to(u.lyr).value
    z = z.to(u.lyr).value
    """
    return x, y, z


def az_alt(latitude, longitude, altitude, psr_ra, psr_dec):
    
    observing_location = EarthLocation(lat = latitude, lon = longitude, height = altitude*u.m)  
    observing_time = Time.now() #Time(datetime(2000, 1, 2, 12, 0, 0), scale='utc')  #Time('2023-07-15 20:12:18')  
    aa = AltAz(location=observing_location, obstime=observing_time)
    
    psr_ra_dec = psr_ra +' '+ psr_dec
    coord = SkyCoord(psr_ra_dec, unit=(u.hourangle, u.deg))  # SkyCoord('4h42m', '-38d6m50.8s')
    az_alt = coord.transform_to(aa)
    return az_alt


#x.to(u.lyr).value

#az_alt(latitude, longitude, altitude, psr_ra[4], psr_dec[4])

"""
psr_ra_dec = psr_ra[0]+' '+ psr_dec[0]

c = SkyCoord('00:42.5 +41:12', unit=(u.hourangle, u.deg))

SkyCoord(psr_ra_dec, unit=(u.hourangle, u.deg))
"""
