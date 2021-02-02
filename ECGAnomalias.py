# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 13:13:13 2018

@author: Juan Camilo
ANOMALIAS DE ECG

"""
import numpy as np
import math
import ECGGenAn as ECGAn

############## VALORES DETERMINADOS POR EL USUARIO ######
def AnomFunc(N,T,t,Anomalia):
    print(Anomalia)

    if (Anomalia==11):
        BPM=45
        F=BPM/60
        T=1/F
        t=np.linspace(0,3*T,3*T*100)
        
    elif (Anomalia==12):
        BPM=160
        F=BPM/60
        T=1/F
        t=np.linspace(0,3*T,3*T*100)
    
    elif (Anomalia==13):
        BPM=150
        F=BPM/60
        T=1/F
        t=np.linspace(0,3*T,3*T*100)
    
    elif (Anomalia==14):
        t=np.linspace(0,9*T,9*T*100)
        BPM=100
        F=BPM/60
        T=1/F
        t=np.linspace(0,3*T,3*T*100)
    
    elif (Anomalia==16):
        BPM=200
        F=BPM/60
        T=1/F
        t=np.linspace(0,3*T,3*T*100)
        
    pi=3.1416
    
    ############### DURACIÓN DE ONDAS ########################
    
    TPW=0.37*(math.sqrt(T))-0.22*T-0.06
    TTW=1.06*(math.sqrt(T))-0.51*T-0.33
    TQRS=0.25*(math.sqrt(T))-0.16*T-0.01 #0.02
    TPQ=0.33*(math.sqrt(T))-0.18*T-0.08
    TST=-0.09*(math.sqrt(T))+0.13*T+0.04
    
    ############## DESFASES ONDA P Y T #######################
    
    DesfaseP=-0.5*TPW-TPQ-0.5*TQRS
    DesfaseT=0.5*TPW+TST+0.5*TQRS
    
    ############# INTERVALOS DE TIEMPO #######################
    
    A=-TPQ-0.5*TQRS-TPW
    B=-TPQ-0.5*TQRS
    C=TST+0.5*TQRS
    D=TTW+TST+0.5*TQRS
    
    ############# CREACIÓN VECTORES DERIVADAS ##################
    
    D1=np.zeros([len(t),1],dtype=np.float64)
    D2=np.zeros([len(t),1],dtype=np.float64)
    D3=np.zeros([len(t),1],dtype=np.float64)
    avR=np.zeros([len(t),1],dtype=np.float64)
    avL=np.zeros([len(t),1],dtype=np.float64)
    avF=np.zeros([len(t),1],dtype=np.float64)
    V1=np.zeros([len(t),1],dtype=np.float64)
    V2=np.zeros([len(t),1],dtype=np.float64)
    V3=np.zeros([len(t),1],dtype=np.float64)
    V4=np.zeros([len(t),1],dtype=np.float64)
    V5=np.zeros([len(t),1],dtype=np.float64)
    V6=np.zeros([len(t),1],dtype=np.float64)
    LV1=np.zeros([len(t),1],dtype=np.float64)
    LV2=np.zeros([len(t),1],dtype=np.float64)
    LV3=np.zeros([len(t),1],dtype=np.float64)
    LV4=np.zeros([len(t),1],dtype=np.float64)
    LV5=np.zeros([len(t),1],dtype=np.float64)
    LV6=np.zeros([len(t),1],dtype=np.float64)
    
    
    if (Anomalia==2):
        print('Fibrilación auricular')
        Aqrs=0.8
        D1=ECGAn.ECG_NSFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,N)
        Aqrs=1.0
        D2=ECGAn.ECG_NSFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,N)
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar1=0.2
        As1=0.6
        V1=ECGAn.ECG_SFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ar1,As1,N)
        Ar1=0.4
        As1=0.9
        V2=ECGAn.ECG_SFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ar1,As1,N)
        Ar1=1.0
        As1=0.6
        V3=ECGAn.ECG_SFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ar1,As1,N)
        Ar1=1.1
        As1=0.3
        V4=ECGAn.ECG_SFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ar1,As1,N)
        Aqrs=1.1
        V5=ECGAn.ECG_NSFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,N)
        Aqrs=0.8
        V6=ECGAn.ECG_NSFA(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,N)
    
    elif (Anomalia==3):
        print('Fibrilación ventricular')
    
    elif (Anomalia==4):
        print('Contracción ventricular prematura')
    
    elif (Anomalia==5):
        print('Asistole')
        ####### Vectores en cero #####################
        
    elif (Anomalia==6):
        print('Bigeminismo')
    
    elif (Anomalia==7):
        print('Trigeminismo')
    
    elif (Anomalia==8):
        print('Perdida de latido')
        Ap=0.1
        At=0.15
        Aqrs=0.6
        D1=ECGAn.ECG_NSPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.15
        At=0.25
        Aqrs=0.9
        D2=ECGAn.ECG_NSPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar=0.2
        As=0.5
        At=0.1
        Ap=0.05
        V1=ECGAn.ECG_SPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=0.4
        As=0.7
        At=0.2
        Ap=0.1
        V2=ECGAn.ECG_SPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=1.0
        As=0.5
        At=0.3
        Ap=0.2
        V3=ECGAn.ECG_SPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ap=0.1
        At=0.25
        Aqrs=1.1
        V4=ECGAn.ECG_NSPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.2
        Aqrs=0.9
        V5=ECGAn.ECG_NSPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.15
        Aqrs=0.7
        V6=ECGAn.ECG_NSPL(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        ################
    
    elif (Anomalia==9):
        print('Contracción ventricular prematura multifocal')
        
    elif (Anomalia==10):
        print('Bloqueo cardiaco de tercer grado')
        
    elif (Anomalia==11):
        print('Bradicarida extrema')
        Ap=0.1
        At=0.15
        Aqrs=0.6
        D1=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.15
        At=0.25
        Aqrs=0.9
        D2=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar=0.2
        As=0.5
        At=0.1
        Ap=0.05
        V1=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=0.4
        As=0.7
        At=0.2
        Ap=0.1
        V2=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=1.0
        As=0.5
        At=0.3
        Ap=0.2
        V3=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ap=0.1
        At=0.25
        Aqrs=1.1
        V4=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.2
        Aqrs=0.9
        V5=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.15
        Aqrs=0.7
        V6=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
    
    elif (Anomalia==12):
        print('Taquicardia extrema')
        Ap=0.1
        At=0.15
        Aqrs=0.6
        D1=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.15
        At=0.25
        Aqrs=0.9
        D2=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar=0.2
        As=0.5
        At=0.1
        Ap=0.05
        V1=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=0.4
        As=0.7
        At=0.2
        Ap=0.1
        V2=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=1.0
        As=0.5
        At=0.3
        Ap=0.2
        V3=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ap=0.1
        At=0.25
        Aqrs=1.1
        V4=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.2
        Aqrs=0.9
        V5=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.15
        Aqrs=0.7
        V6=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
    
    elif (Anomalia==13):
        print('Flúter auricular')
        Ap=0.1
        At=0.1
        Aqrs=0.7
        D1=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.15
        At=0.2
        Aqrs=1
        D2=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar=0.2
        As=0.6
        At=0.05
        Ap=0.02
        V1=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=0.4
        As=0.9
        At=0.15
        Ap=0.1
        V2=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=1.0
        As=0.7
        At=0.3
        Ap=0.2
        V3=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ap=0.2
        At=0.15
        Aqrs=1.2
        V4=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.15
        At=0.1
        Aqrs=1
        V5=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.15
        Aqrs=0.6
        V6=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        
    elif (Anomalia==14):
        print('Taquicardia paroxistica')
        Ap=0.1
        At=0.15
        Aqrs=0.6
        D1=ECGAn.ECG_NSTP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.15
        At=0.25
        Aqrs=0.9
        D2=ECGAn.ECG_NSTP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar=0.2
        As=0.5
        At=0.1
        Ap=0.05
        V1=ECGAn.ECG_STP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=0.4
        As=0.7
        At=0.2
        Ap=0.1
        V2=ECGAn.ECG_STP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=1.0
        As=0.5
        At=0.3
        Ap=0.2
        V3=ECGAn.ECG_STP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ap=0.1
        At=0.25
        Aqrs=1.1
        V4=ECGAn.ECG_NSTP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.2
        Aqrs=0.9
        V5=ECGAn.ECG_NSTP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.15
        Aqrs=0.7
        V6=ECGAn.ECG_NSTP(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
    
    elif (Anomalia==15):
        print('Ritmo nodal')
        Ap=0.2
        At=0.2
        Aqrs=0.6
        D1=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=-0.1
        At=0.2
        Aqrs=0.95
        D2=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        D2=D2+0.05
        
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar=0.2
        As=0.5
        At=0.1
        Ap=0.05
        V1=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=0.4
        As=0.7
        At=0.2
        Ap=0.1
        V2=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=1.0
        As=0.5
        At=0.3
        Ap=0.2
        V3=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ap=0.1
        At=0.25
        Aqrs=1.1
        V4=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.2
        Aqrs=0.9
        V5=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.15
        Aqrs=0.7
        V6=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        
    elif (Anomalia==16):
        print('Taquicardia supraventricular')
        Ap=0.1
        At=0.15
        Aqrs=0.6
        D1=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.15
        At=0.25
        Aqrs=0.9
        D2=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        D3=D2-D1
        ##########################################################################
        avR=(-1/2)*(D1+D2)
        avL=D1-(1/2)*D2
        avF=D2-(1/2)*D1
        ##########################################################################
        Ar=0.2
        As=0.5
        At=0.1
        Ap=0.05
        V1=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=0.4
        As=0.7
        At=0.2
        Ap=0.1
        V2=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ar=1.0
        As=0.5
        At=0.3
        Ap=0.2
        V3=ECGAn.ECG_S(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Ap,Ar,As,At,N)
        Ap=0.1
        At=0.25
        Aqrs=1.1
        V4=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.2
        Aqrs=0.9
        V5=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
        Ap=0.1
        At=0.15
        Aqrs=0.7
        V6=ECGAn.ECG_NS(t,DesfaseP,DesfaseT,A,B,C,D,pi,TPW,TTW,TQRS,TPQ,TST,T,Aqrs,Ap,At,N)
    
    elif (Anomalia==17):
        print('Taquicardia ventricular')
    
    Signal=np.zeros((16,len(t)))
    Signal[0,0:len(t)]=D1[:,0]
    Signal[1,0:len(t)]=D2[:,0]
    Signal[2,0:len(t)]=LV1[:,0]
    Signal[3,0:len(t)]=LV2[:,0]
    Signal[4,0:len(t)]=LV3[:,0]
    Signal[5,0:len(t)]=LV4[:,0]
    Signal[6,0:len(t)]=LV5[:,0]
    Signal[7,0:len(t)]=LV6[:,0]
    Signal[8,0:len(t)]=V1[:,0]
    Signal[9,0:len(t)]=V2[:,0]
    Signal[10,0:len(t)]=V3[:,0]
    Signal[11,0:len(t)]=V4[:,0]
    Signal[12,0:len(t)]=V5[:,0]
    Signal[13,0:len(t)]=V6[:,0]
    Signal[14,0:len(t)]=t[:]
    
    return Signal

######################### PLOTEO DE FUNCIONES ###############################

#fig, aD1 = plt.subplots()
#aD1.plot(t, D1,'b')
#aD1.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD DI')
#aD1.grid()
#plt.plot(t,D1)
#plt.xlim(TTW+TST+TQRS+0.2,1+TTW+TST+TQRS)
#plt.show()
#
#fig, aD2 = plt.subplots()
#aD2.plot(t, D2,'b')
#aD2.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD DII')
#aD2.grid()
#plt.plot(t,D2[:])
#plt.xlim(TTW+TST+TQRS+0.2,1+TTW+TST+TQRS)
#plt.show()
#
#fig, aD3 = plt.subplots()
#aD3.plot(t, D3)
#aD3.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD DIII')
#aD3.grid()
#plt.plot(t,D3)
#plt.show()
#
#fig, aavR = plt.subplots()
#aavR.plot(t, avR)
#aavR.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD avR')
#aavR.grid()
#plt.plot(t,avR)
#plt.show()
#
#fig, aavL = plt.subplots()
#aavL.plot(t, avL)
#aavL.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD avL')
#aavL.grid()
#plt.plot(t,avL)
#plt.show()
#
#fig, aavF = plt.subplots()
#aavF.plot(t, avF)
#aavF.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD avF')
#aavF.grid()
#plt.plot(t,avF)
#plt.show()
#
#fig, aV1 = plt.subplots()
#aV1.plot(t, V1)
#aV1.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD V1')
#aV1.grid()
#plt.plot(t,V1)
#plt.show()
#
#fig, aV2 = plt.subplots()
#aV2.plot(t, V2)
#aV2.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD V2')
#aV2.grid()
#plt.plot(t,V2)
#plt.show()
#
#fig, aV3 = plt.subplots()
#aV3.plot(t, V3)
#aV3.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD V3')
#aV3.grid()
#plt.plot(t,V3)
#plt.show()
#
#fig, aV4 = plt.subplots()
#aV4.plot(t, V4)
#aV4.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD V4')
#aV4.grid()
#plt.plot(t,V4)
#plt.show()
#
#fig, aV5 = plt.subplots()
#aV5.plot(t, V5)
#aV5.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD V5')
#aV5.grid()
#plt.plot(t,V5)
#plt.show()
#
#fig, aV6 = plt.subplots()
#aV6.plot(t, V6)
#aV6.set(xlabel='Tiempo (s)', ylabel='Voltaje (mV)', title='LEAD V6')
#aV6.grid()
#plt.plot(t,V6)
#plt.show()