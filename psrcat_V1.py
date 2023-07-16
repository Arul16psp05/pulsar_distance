#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 15:13:24 2023

@author: arul
"""

# Import pandas
#import pandas as pd
import numpy as np
import csv
from ra_dec_to_Alt_Az import az_alt
from ra_dec_to_Alt_Az import sperical_to_cartesian

data_file = "/home/arul/Documents/py_projects/pulsar_catalog/psrcat_db.csv"


data = open(data_file, "r")

psrcat = np.array(list(csv.reader(data, delimiter=",")))

"""
psr_file = "/home/arul/Documents/py_projects/pulsar_catalog/pulsars.csv"
psr_data = open(psr_file, "r")
psr = np.array(list(csv.reader(psr_data, delimiter=",")))
"""


file_name = "/home/arul/Documents/py_projects/pulsar_catalog/proposed_33_pulsars_for_cy45.txt"
psr = np.loadtxt(file_name,dtype='str')

psr_index = np.zeros(len(psr), dtype='int')

def find_idx(psr_name):
    index = int(np.where(psrcat[1:,0]== psr_name)[0][0])
    return index

def ra_dec_dist_idx():
    ra_idx = np.where(psrcat[0] == "RAJ")[0][0]
    dec_idx = np.where(psrcat[0] == "DECJ")[0][0]
    dist_idx = np.where(psrcat[0] == "DIST_DM")[0][0]
    return ra_idx, dec_idx, dist_idx

ra_idx, dec_idx, dist_idx =  ra_dec_dist_idx()

for idx in range(len(psr)):
    index = find_idx(psr[idx])
    psr_index[idx] = int(index)
    #print(psr_index)
    
    
def psr_ra_dec_dist(psr_index, ra_idx, dec_idx, dist_idx, idy):
    ra = str(psrcat[1:,ra_idx][psr_index[idy]])
    dec = str(psrcat[1:,dec_idx][psr_index[idy]])
    dist = str(psrcat[1:,dist_idx][psr_index[idy]])
    return ra, dec, dist
    
    
psr_ra  = np.zeros(len(psr), dtype='object')
psr_dec = np.zeros(len(psr), dtype='object')
psr_dist = np.zeros(len(psr), dtype='float')
    
for idy in range(len(psr)):
    psr_ra[idy], psr_dec[idy], psr_dist[idy] = psr_ra_dec_dist(psr_index, ra_idx, dec_idx, dist_idx, idy)


##------------------------------

#Convert RA DEC to AZ ALT

"""
GMRT
"""

latitude = 19.0912

longitude = 74.0432

altitude = 560


#psr_az  = np.zeros(len(psr), dtype='object')
psr_az_alt = np.zeros(len(psr), dtype='object')

for idz in range (len(psr)): 
    
    psr_az_alt[idz] = az_alt(latitude, longitude, altitude, psr_ra[idz], psr_dec[idz])

#psr_az[idz], 
#idz
# Converted into light years

psr_x_cord  = np.zeros(len(psr), dtype='float')
psr_y_cord  = np.zeros(len(psr), dtype='float')
psr_z_cord  = np.zeros(len(psr), dtype='float')

for ida in range (len(psr)):
    # Converted into light years
    psr_x_cord[ida], psr_y_cord[ida], psr_z_cord[ida] = sperical_to_cartesian( psr_ra[ida], psr_dec[ida], psr_dist[ida])


from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)
for idb in range(len(psr)):
    xline = [0,psr_x_cord[idb]]
    yline = [0,psr_y_cord[idb]]
    zline = [0,psr_z_cord[idb]]
    ax.plot3D(xline, yline, zline, 'blue')
    
"""
ax = plt.axes(projection='3d')
idb = 15
xline = [0,psr_x_cord[idb]]
yline = [0,psr_y_cord[idb]]
zline = [0,psr_z_cord[idb]]
ax.plot3D(xline, yline, zline, 'red')
"""

ax = plt.axes(projection='3d')
for idb in range(len(psr)):
    xline = [0,psr_x_cord[idb]]
    yline = [0,psr_y_cord[idb]]
    zline = [0,psr_z_cord[idb]]
    ax.plot3D(xline, yline, zline, 'blue')
    ax.plot3D(psr_x_cord[idb], psr_y_cord[idb], psr_z_cord[idb], 'red')

"""
ax.set_xlabel('Distance  (lYr)')
ax.set_ylabel('Distance  (lYr)')
ax.set_zlabel('Distance  (lYr)')
"""
#ax = plt.axes(projection='3d')
for idb in range(len(psr)):
    ax.scatter(psr_x_cord[idb], psr_y_cord[idb], psr_z_cord[idb], c='r')
    


ax.set_xlabel('Distance  (kpc)')
ax.set_ylabel('Distance  (kpc)')
ax.set_zlabel('Distance  (kpc)')


import scipy.io as sio

mdict = {"x":psr_x_cord, "y":psr_y_cord, "z":psr_z_cord, "psr":psr}
sio.savemat("psr_dist_XYZ.mat", mdict)

# plot ref: https://stackoverflow.com/questions/54722526/how-to-get-vertical-lines-in-a-3d-scatter-plot-in-matlab
"""
ax.set_xlabel('Distance  (lYr)')
ax.set_ylabel('Distance  (lYr)')
ax.set_zlabel('Distance  (lYr)')
"""

"""
for ii in range(0,360,1):
    ax.view_init(elev=10., azim=ii)
    plt.savefig("movie%d.png" % ii)

"""
"""
find_idx("J0007+7303")

ra = 

# reading csv file 
psrcat = pd.read_csv(data_file)

psr = psrcat.to_numpy()

ps = psrcat["PSRJ"].values

ps[np.where(ps == "J2355+0051")[0][0]]


#arr = np.genfromtxt(data_file, delimiter="," , dtype=str)
"""
