#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 16:12:13 2022

@author: Zeephod
"""

# Python 3 program for the
# haversine formula
import math

def spherical_distance(lat1, long1, lat2, long2):
    phi1 = 0.5*math.pi - lat1
    phi2 = 0.5*math.pi - lat2
    r = 0.5*(6378137 + 6356752) # mean radius in meters
    t = math.sin(phi1)*math.sin(phi2)*math.cos(long1-long2) + math.cos(phi1)*math.cos(phi2)
    return r * math.acos(t) 

def ellipsoidal_distance(lat1, long1, lat2, long2):

    a = 6378137.0 # equatorial radius in meters 
    f = 1/298.257223563 # ellipsoid flattening 
    b = (1 - f)*a 
    tolerance = 1e-11 # to stop iteration

    phi1, phi2 = lat1, lat2
    U1 = math.atan((1-f)*math.tan(phi1))
    U2 = math.atan((1-f)*math.tan(phi2))
    L1, L2 = long1, long2
    L = L2 - L1

    lambda_old = L + 0

    while True:
    
        t = (math.cos(U2)*math.sin(lambda_old))**2
        t += (math.cos(U1)*math.sin(U2) - math.sin(U1)*math.cos(U2)*math.cos(lambda_old))**2
        sin_sigma = t**0.5
        cos_sigma = math.sin(U1)*math.sin(U2) + math.cos(U1)*math.cos(U2)*math.cos(lambda_old)
        sigma = math.atan2(sin_sigma, cos_sigma) 
    
        sin_alpha = math.cos(U1)*math.cos(U2)*math.sin(lambda_old) / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2
        cos_2sigma_m = cos_sigma - 2*math.sin(U1)*math.sin(U2)/cos_sq_alpha
        C = f*cos_sq_alpha*(4 + f*(4-3*cos_sq_alpha))/16
    
        t = sigma + C*sin_sigma*(cos_2sigma_m + C*cos_sigma*(-1 + 2*cos_2sigma_m**2))
        lambda_new = L + (1 - C)*f*sin_alpha*t
        if abs(lambda_new - lambda_old) <= tolerance:
            break
        else:
            lambda_old = lambda_new

    u2 = cos_sq_alpha*((a**2 - b**2)/b**2)
    A = 1 + (u2/16384)*(4096 + u2*(-768+u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
    t = cos_2sigma_m + 0.25*B*(cos_sigma*(-1 + 2*cos_2sigma_m**2))
    t -= (B/6)*cos_2sigma_m*(-3 + 4*sin_sigma**2)*(-3 + 4*cos_2sigma_m**2)
    delta_sigma = B * sin_sigma * t
    s = b*A*(sigma - delta_sigma)

    return s
 

# Python 3 program for the
# haversine formula
def haversine(lat1, lon1, lat2, lon2):
     
    # distance between latitudes
    # and longitudes
    dLat = (lat2 - lat1) * math.pi / 180.0
    dLon = (lon2 - lon1) * math.pi / 180.0
 
    # convert to radians
    lat1 = (lat1) * math.pi / 180.0
    lat2 = (lat2) * math.pi / 180.0
 
    # apply formulae
    a = (pow(math.sin(dLat / 2), 2) +
         pow(math.sin(dLon / 2), 2) *
             math.cos(lat1) * math.cos(lat2));
    rad = 6371
    c = 2 * math.asin(math.sqrt(a))
    return rad * c
 
# Driver code
if __name__ == "__main__":
    lat1 = 51.5081
    lon1 = 0.0759
    lat2 = 48.8738
    lon2 = -2.2950
     
    print("Haversine distance in miles ",haversine(lat1, lon1,lat2, lon2)/1.6)
    
#print out the distance using Vincenty's ellipsoid formula
print("distance in miles vincenty ", ellipsoidal_distance(lat1, lon1, lat2, lon2 )/160000 )  
print("Spherical distance ",spherical_distance(lat1, lon1, lat2, lon2 )/160000)  

 
# This code is contributed
# by ChitraNayal

LocationOneV = lon1
LocationTwoV = lon2
LocationOneH = lat1
LocationTwoH = lat2
#My own version of distance created using the flat earth model
    #Distance calculation
DistanceV = (LocationOneV - LocationTwoV) *54.6
DistanceH = (LocationOneH - LocationTwoH) *68.7
DistanceCalc = (DistanceV*DistanceV)+(DistanceH*DistanceH)
ActualDistance = math.sqrt(DistanceCalc)
print("Distance in Miles using flat earth calculation ", ActualDistance)
