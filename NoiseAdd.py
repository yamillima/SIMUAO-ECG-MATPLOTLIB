# -*- coding: utf-8 -*-
"""
Created on Thu May 23 11:09:32 2019

@author: Juan Camilo

"""
import numpy as np
import math

def NoiseAdd(ECG,Noise,t):
    ECGN=np.zeros((15,len(t)))
    pi=3.1416
    for i in range (0,len(t)):
            ECGN[0,i]=ECG[0,i]+Noise[0]*2*math.sin(2*pi*Noise[1]*t[i])
            ECGN[1,i]=ECG[1,i]+Noise[0]*0.03*math.sin(2*pi*Noise[1]*t[i])
            ECGN[8,i]=ECG[8,i]+Noise[0]*0.02*math.sin(2*pi*Noise[1]*t[i])
            ECGN[9,i]=ECG[9,i]+Noise[0]*0.025*math.sin(2*pi*Noise[1]*t[i])
            ECGN[10,i]=ECG[10,i]+Noise[0]*0.022*math.sin(2*pi*Noise[1]*t[i])
            ECGN[11,i]=ECG[11,i]+Noise[0]*0.021*math.sin(2*pi*Noise[1]*t[i])
            ECGN[2,i]=ECGN[8,i]+(1/3)*(ECGN[0,i]+ECGN[1,i])
            ECGN[3,i]=ECGN[9,i]+(1/3)*(ECGN[0,i]+ECGN[1,i])
            ECGN[4,i]=ECGN[10,i]+(1/3)*(ECGN[0,i]+ECGN[1,i])
            ECGN[5,i]=ECGN[11,i]+(1/3)*(ECGN[0,i]+ECGN[1,i])
            ECGN[6,i]=ECGN[12,i]+(1/3)*(ECGN[0,i]+ECGN[1,i])
            ECGN[7,i]=ECGN[13,i]+(1/3)*(ECGN[0,i]+ECGN[1,i])           
    
    return ECGN