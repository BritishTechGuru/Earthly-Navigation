#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 12:04:26 2022

@author: zephod
"""
import math
DiameterOfEarth = 7917.5

#Assumung two fixed point - just for the sake of making life easy...
#North = positive values. South = negative values
#West = Positive values, East = negative values
LocationOneV = 51.5081
LocationOneH = 0.0759
#Tower of London

LocationTwoV = 48.8738
LocationTwoH = -2.2950
#Arc de Triomph

X = math.cos(LocationTwoV)*math.sin(LocationOneH-LocationTwoH)
Y = (math.cos(LocationOneV)*math.sin(LocationTwoV)-math.sin(LocationOneV)*
    math.cos(LocationTwoV)*math.cos(LocationOneH-LocationTwoH))

Z = math.atan2(X,Y)  #radians
print("radians ", Z)

Degrees = math.degrees(Z)
print ("Degrees", Degrees)

#Distance calculation
DistanceV = (LocationOneV - LocationTwoV) *54.6
DistanceH = (LocationOneH - LocationTwoH) *68.7
DistanceCalc = (DistanceV*DistanceV)+(DistanceH*DistanceH)
ActualDistance = math.sqrt(DistanceCalc)
print("Distance in Miles ", ActualDistance)
