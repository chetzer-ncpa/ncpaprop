#!/usr/bin/env python
#
# G2S_GenGrid.py
#
# Usage: G2S_GenGrid.py <TimeStamp> <LatStart> <LatEnd> <LatCnt> <LongStart> <LongEnd>  <LongCnt>
# e.g., G2S_GenGrid.py 2013020512 33.0 34.0 2 -105.0 -109.0 10
#
# Philip Blom (pblom@lanl.gov)

import sys, os, pdb, math
import numpy as np

def OneLat2km(lat):
    rad = math.radians(lat)
    return 111.132 - 0.559*math.cos(2.0*rad) + 0.001175*math.cos(4.0*rad)

def OneLon2km(lat):
    rad = math.radians(lat)
    return math.pi/180.0 * 6378.137 * math.cos(rad) / math.sqrt(1.0 - 0.006694*math.sin(rad)**2)

def main():
    # Read in script inputs correcting order of upper and lower bounds
    DateTime =  sys.argv[1]

    LatStart =  min(float(sys.argv[2]), float(sys.argv[3]))
    LatEnd =    max(float(sys.argv[2]), float(sys.argv[3]))
    LatCnt =    int(sys.argv[4])

    LonStart = min(float(sys.argv[5]), float(sys.argv[6]))
    LonEnd =   max(float(sys.argv[5]), float(sys.argv[6]))
    LonCnt =   int(sys.argv[7])

    # Set up latitude and longitude values and generate profiles from G2S .bin file
    LatArray = np.linspace(LatStart, LatEnd, num=LatCnt)
    LonArray = np.linspace(LonStart, LonEnd, num=LonCnt)

    for nlat in range(LatCnt):
        for nlon in range(LonCnt):
            lat = LatArray[nlat]
            lon = LonArray[nlon]
        
            Command = './extractpro -v -t %s -k G2SGCSE -l %s %s > Profile%s.met' % (DateTime, lat, lon, nlon + nlat*LonCnt)

            print 'Extracting profile at %s degrees latitude, %s degrees longitude into Profile%s.met' % (lat, lon, nlon + nlat*LonCnt)
            os.system(Command)

    print ''

    # Set region origin as mid point of latitudes and longitudes
    # Then output the locx and locy file information
    print 'Writing locx and locy files...'
    midlat = LatArray.mean()
    midlon = LonArray.mean()

    locx_out = open('Profile.locx', 'w')
    locy_out = open('Profile.locy', 'w')
           
    for nlat in range(LatCnt):
        print >> locy_out, OneLat2km(midlat) * (LatArray[nlat] - midlat)

    for nlon in range(LonCnt):
        print >> locx_out, OneLon2km(midlat) * (LonArray[nlon] - midlon)

    locx_out.close()
    locy_out.close()

    print 'Writing locx and locy files...'

    loclat_out = open('Profile.loclat', 'w')
    loclon_out = open('Profile.loclon', 'w')

    for nlat in range(LatCnt):
        print >> loclat_out, LatArray[nlat]

    for nlon in range(LonCnt):
        print >> loclon_out, LonArray[nlon]

    loclat_out.close()
    loclon_out.close()

    # Clean the various Profile{Number}.met files
    print 'Cleaning profiles...'
    for i in range(LatCnt*LonCnt):
        f = open('Profile' + str(i) + '.met','r')
        f_out = open('Profile' + str(i) + '.out','w')
    
        inData = 0
        for line in f:
            if (line.split()[0] == 'Z'):
                inData = 1
                continue
        
            if (inData == 1):
                # Extract the relevant fields:
                z = line[0:12]
                t = line[12:24]
                u = line[24:36]
                v = line[36:48]
                d = line[48:60]
                p = line[60:72]
            
                print >> f_out, z, t, u, v, d, p

        f.close()
        f_out.close()

        os.system('rm Profile' + str(i) + '.met')
        os.system('mv Profile' + str(i) + '.out Profile' + str(i) + '.met')

main()
