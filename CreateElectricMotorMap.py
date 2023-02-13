# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 16:07:24 2023

@author: Leo Laine,

Edited: Sachin Janardhanan, 
"""

import numpy as np
import scipy.constants as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d
from scipy import interpolate

def CreateElectricMotorMap(pltnr,MotorName,MotorType,PMax,TMax,OmegaMax,FieldWR):
    """ Create Electric Motor Maps
   Input=Motortype=AM (induction) or PM (permanent magnet)
         PowerContinuousMax electric Motor PMAX     [kW] 
         TorqueContinuousMax electric Motor         [Nm]
         OmegaMax electric Motor                    [rpm]
         Field weakening ratio,                     [-]
         i.e. ratio between max speed and base
   Output=TorqueContinuousMax Electric Motor TMAX (or given as input =/0)
          Efficiency map (Omega, Torque)
          MaxTorque (Omega)
          MinTorque (Omega)
    """
    #
    # General Properties
    #
    Epsilon = np.finfo(float).eps # (small value to avoid /0..)

        
    #
    # Zero Variables
    #
    i=0
    j=0
    Psi=0
    Isx=0 
    Isy=0 
    IsyMAX=0
    Curr=0 
    Tmax=0 
    T=0 
    Pout=0 
    Pcopper=0 
    Piron=0 
    Pwindage=0 
    Ploss=0 
    Eta=0


    def MaxContinuousTorque(MaxContinuousPower,FWR,OmegaMax):
        """Calculate Peak Continuous Torque"""
        return MaxContinuousPower/(OmegaMax/FWR)
#    
    def MaxContinuousPower(MaxContinuousTorque,FWR,OmegaMax):
        """Calculate Peak Continuous Power"""
        return MaxContinuousTorque*(OmegaMax*FWR)

    #
    # Motor Properties	                           Unit
    #
    FWR = FieldWR 
    OmegaMax = OmegaMax*np.pi/30                   # [rpm] -> [rad/s]
    if PMax == 0:
        PowerMax = MaxContinuousPower(TMax, FWR, OmegaMax)
    else:
        PowerMax = PMax*1000                       # [kW] ->  [W]                            
    
    print('PowerMax [kW] =', PowerMax/1000)                                               # [-]
    if TMax == 0:
        TorqueMax = MaxContinuousTorque(PowerMax, FWR, OmegaMax)
    else:
        TorqueMax = TMax

    print('TorqueMax [Nm] =', TorqueMax)

    ElectricMotorType=MotorType


    #def MotorMaps(PowerMax,TorqueMax,ElectricMotorType):
    #    """"Electric Motor Maps"""
    StepsTorque=100 #i loop [25]
    StepsOmega=100 #j loop [25]
    OmegaVector = np.linspace(0,FWR,StepsOmega)
    TrqRefVector= np.linspace(-1,1,StepsTorque)
    
    
    # Pre-allocate Matrices
    Psi = np.zeros((StepsTorque,StepsOmega))      # Pre-allocate matrix
    Isx = np.zeros((StepsTorque,StepsOmega))      # Pre-allocate matrix
    Isy = np.zeros((StepsTorque,StepsOmega))      # Pre-allocate matrix
    IsyMAX = np.zeros((StepsTorque,StepsOmega))   # Pre-allocate matrix
    Curr = np.zeros((StepsTorque,StepsOmega))     # Pre-allocate matrix
    Tmax = np.zeros((StepsTorque,StepsOmega))     # Pre-allocate matrix
    Torque = np.zeros((StepsTorque,StepsOmega))   # Pre-allocate matrix
    PowerOut = np.zeros((StepsTorque,StepsOmega)) # Pre-allocate matrix
    Pcopper = np.zeros((StepsTorque,StepsOmega))  # Pre-allocate matrix
    Piron = np.zeros((StepsTorque,StepsOmega))    # Pre-allocate matrix
    Pwindage = np.zeros((StepsTorque,StepsOmega)) # Pre-allocate matrix
    Ploss =  np.zeros((StepsTorque,StepsOmega))   # Pre-allocate matrix
    Eta = np.zeros((StepsTorque,StepsOmega))      # Pre-allocate matrix

    for i in range(0,np.size(TrqRefVector),1):
        #print('i loop', i)
        for j in range(0,np.size(OmegaVector),1):
             #print('j loop', j)
             Psi[i,j]=min(1,1/(Epsilon+OmegaVector[j]))
             if ElectricMotorType == 'PermanentMagnetMotor':
                Isx[i,j]= min(1,max((OmegaVector[j]-1),0)/(max(OmegaVector)-1+Epsilon))
                Isy[i,j]= min(np.absolute(TrqRefVector[i])/Psi[i,j],np.sqrt(1-(Isx[i,j])**2))
                IsyMAX[i,j]= min(max(np.absolute(TrqRefVector))/Psi[i,j],np.sqrt(1-(Isx[i,j])**2))
             elif ElectricMotorType=='AsynchronousMotor': # induction motor
                Isx[i,j]= 0.15*Psi[i,j]
                Isy[i,j]= min(abs(TrqRefVector[i])/Psi[i,j],np.sqrt(1-(Isx[i,j])**2))
                IsyMAX[i,j]= min(max(np.absolute(TrqRefVector))/Psi[i,j],np.sqrt(1-(Isx[i,j])**2))
    # #       end
             Curr[i,j]= np.sqrt((Isx[i,j])**2+(Isy[i,j])**2)
             Tmax[i,j]= min(TorqueMax,PowerMax/(Epsilon+OmegaVector[j])) # % Pem_maxPsi(i,j)*IsyMAX(i,j);
             Torque[i,j] = min(TrqRefVector[i],Tmax[i,j])
             PowerOut[i,j]= OmegaVector[j]*Torque[i,j]
             if ElectricMotorType == 'PermanentMagnetMotor':
                   Pcopper[i,j]=(Curr[i,j])**2*0.025 
                   #print('Pcopper[i,j]=',Pcopper[i,j],Curr[i,j], j)
             elif ElectricMotorType == 'AsynchronousMotor':
                   Pcopper[i,j]=(Curr[i,j])**2*0.03 + (Isy[i,j])**2*0.035
    # #         end
             Piron[i,j] = OmegaVector[j]*0.01 + OmegaVector[j]*(Psi[i,j])**2*0.003 + \
             (OmegaVector[j]*Psi[i,j])**2*0.003
             Pwindage[i,j]=(OmegaVector[j]/max(OmegaVector))**3*0.02
             Ploss[i,j]= 0*0.005+Pcopper[i,j]+Piron[i,j]+Pwindage[i,j]
             if np.absolute(TrqRefVector[i]*OmegaVector[j])<=1:
                Eta[i,j] = min(1,(max(0.05,(abs(PowerOut[i,j])/(abs(PowerOut[i,j])+Ploss[i,j])))))
             else:
                Eta[i,j] = 0#np.NAN
    # #         end
             if Torque[i,j]==0:
                Eta[i,j] = 0.05
    #         end
    #     end
    # end
    
    # Plot 101
#    pltstringEMMAP101="plt "+str(pltnr+1)+" - Electric Motor Map scaled data",ElectricMotorType
#    print(pltstringEMMAP101)
    

      
    fig = plt.figure(pltnr+1,figsize=(12,6))
    # ax1 = fig.add_subplot(251, projection='3d')
    # ax2 = fig.add_subplot(252, projection='3d')
    # ax3 = fig.add_subplot(253, projection='3d')
    # ax4 = fig.add_subplot(254, projection='3d')
    # ax5 = fig.add_subplot(255, projection='3d')
    # ax6 = fig.add_subplot(256, projection='3d')
    # ax7 = fig.add_subplot(257, projection='3d')
    # ax8 = fig.add_subplot(258, projection='3d')
    # ax9 = fig.add_subplot(259, projection='3d')
    # ax10 = fig.add_subplot(2,5,10, projection='3d')
#    X,Y=np.meshgrid(OmegaVector,TrqRefVector)
#    plt.xlabel('OmegaVector [-]')
#    plt.ylabel('TorqueReferenceVector [-]')
#    fig.suptitle(pltstringEMMAP101) 
    #Z = X*np.exp(-X - Y)
    
    
    # Plot a basic wireframe
    # ax1.plot_surface(X, Y, Psi, rstride=1, cstride=1, cmap='viridis')
    # ax1.set_title('Psi')
    # ax2.plot_surface(X, Y, Isx,rstride=1, cstride=1, cmap='viridis')
    # ax2.set_title('Isx')
    # ax3.plot_surface(X, Y, Isy,rstride=1, cstride=1, cmap='viridis')
    # ax3.set_title('Isy')
    # ax4.plot_surface(X, Y, Torque,rstride=1, cstride=1, cmap='viridis')
    # ax4.set_title('Torque')
    # ax5.plot_surface(X, Y, PowerOut,rstride=1, cstride=1, cmap='viridis')
    # ax5.set_title('PowerOut')
    # ax6.plot_surface(X, Y, Pcopper,rstride=1, cstride=1, cmap='viridis')
    # ax6.set_title('Pcopper')
    # ax7.plot_surface(X, Y, Piron,rstride=1, cstride=1, cmap='viridis')
    # ax7.set_title('Piron')
    # ax8.plot_surface(X, Y, Pwindage,rstride=1, cstride=1, cmap='viridis')
    # ax8.set_title('Pwindage')
    # ax9.plot_surface(X, Y, Ploss,rstride=1, cstride=1, cmap='viridis')
    # ax9.set_title('Ploss')
    # ax10.plot_surface(X, Y, Eta,rstride=1, cstride=1, cmap='viridis')
    # ax10.set_title('Eta')
    
    PowerSteps=50
    TorqueEM = TorqueMax*TrqRefVector/max(TrqRefVector);
    OmegaEM = OmegaMax*OmegaVector/max(OmegaVector);
    PowerEM = np.linspace(-PowerMax,PowerMax,PowerSteps);
    # Pre Allocate
    TorqueLimResult = np.zeros(np.size(OmegaEM))
    TorqueLim = np.zeros(np.size(OmegaEM));
    EnergyConsEM = np.zeros((np.size(TrqRefVector),np.size(OmegaVector)))
    
    TorqueLim=np.minimum(TorqueMax*np.ones(np.size(OmegaEM)),\
              np.ones(np.size(OmegaEM))*PowerMax/(Epsilon+OmegaEM), TorqueLimResult)
    EtaEM=Eta;
    #EtaTemp=np.zeros((np.size(TrqRefVector),np.size(OmegaVector)))
    i=0
    j=0
    for i in range(0,np.size(TrqRefVector),1):
        #print('i loop', i)
        for j in range(0,np.size(OmegaVector),1):
            if OmegaEM[j]*TorqueEM[i] >= 0:
               EnergyConsEM[i,j] = min(PowerMax,OmegaEM[j]*TorqueEM[i])+Ploss[i,j]*PowerMax
            else:
               EnergyConsEM[i,j] = max(-PowerMax,OmegaEM[j]*TorqueEM[i])-Ploss[i,j]*PowerMax # Adding power-losses to Kinetic power instead of subtracting (2022-02-11)
               # EnergyConsEM[i,j] = max(-PowerMax,OmegaEM[j]*TorqueEM[i])+Ploss[i,j]*PowerMax

    
    PowerlossEM=Ploss*PowerMax
    
    
    X,Y=np.meshgrid(OmegaEM,TorqueEM)
    Z=EtaEM
    # old code removed 2021-02-07
    #tck = interpolate.bisplrep(X, Y, Z, s=0, kx=2,ky=2)
    #def givemeZ(x,y):
    #    return interpolate.bisplev(x,y,tck)
    givemeZ = interpolate.interp2d(X,Y,Z)
    
    
    PowerToTorqueEM=np.zeros((np.size(PowerEM),4))
    #EtaTemp=np.zeros(np.size(PowerEM))
    # Fix powertotorque columns remember 0 to 3 2020-11-05
    PowerToTorqueEM[:,0]=PowerEM
    for i in range(0,np.size(PowerEM),1):
        #print('i loop', i)
        if PowerEM[i]==0:
            PowerToTorqueEM[i,1]= 0
            PowerToTorqueEM[i,2]= Epsilon
            PowerToTorqueEM[i,3]= 0
        elif PowerEM[i]<0: #elif
            for j in range(0,np.int(np.floor(np.size(TorqueEM)/2)),1):
                #print('j loop PowerEM<0', j)
                TorqueTemp = TorqueEM[j]
                OmegaTemp  = PowerEM[i]/TorqueTemp
                # Get nearest id value of OmegaEM for OmegaTemp
                idx = (np.abs(OmegaEM-OmegaTemp)).argmin()  
                #EtaTemp = givemeZ(OmegaTemp, TorqueTemp)
                EtaTemp = givemeZ(OmegaEM[idx],TorqueEM[j])
                #EtaTemp=interpolate.griddata(X,Y,Z,(OmegaTemp,TorqueTemp),method='nearest')
                #print('EtaTemp =',idx, EtaTemp, OmegaTemp,TorqueTemp)
                if EtaTemp > PowerToTorqueEM[i,2]:
                         #print('banana')
                         PowerToTorqueEM[i,1]=TorqueEM[j]
                         PowerToTorqueEM[i,2]=EtaTemp
                         PowerToTorqueEM[i,3]=OmegaTemp
        elif PowerEM[i]>0: #elif
            for j in range(np.int(np.ceil(np.size(TorqueEM)/2)),np.size(TorqueEM),1):
                #print('PowerEm>0') 
                TorqueTemp = TorqueEM[j]
                OmegaTemp  = PowerEM[i]/TorqueTemp
                # Get nearest id value of OmegaEM for OmegaTemp
                idx = (np.abs(OmegaEM-OmegaTemp)).argmin()
                #EtaTemp = givemeZ(OmegaTemp, TorqueTemp)
                EtaTemp = givemeZ(OmegaEM[idx],TorqueEM[j])
                #print('d EtaTemp =',EtaTemp, OmegaTemp,TorqueTemp)
                if EtaTemp > PowerToTorqueEM[i,2]:
                        #print('lemon')
                        PowerToTorqueEM[i,1]=TorqueEM[j]
                        PowerToTorqueEM[i,2]=EtaTemp
                        PowerToTorqueEM[i,3]=OmegaTemp        
                         
    for i in range(0,np.size(PowerEM)-1,1):
        PowerToTorqueEM[i,1]=(PowerToTorqueEM[i-1,1]+PowerToTorqueEM[i,1]+PowerToTorqueEM[i+1,1])/3
        PowerToTorqueEM[i,3]=(PowerToTorqueEM[i-1,3]+PowerToTorqueEM[i,3]+PowerToTorqueEM[i+1,3])/3
    
    
    

    
    #Plot 102
    # pltstringEMMAP102="plt "+str(pltnr+2)+"- Electric Motor Map generated data Motortype="+str(ElectricMotorType)+", MaxPower [kW]= "+str(PowerMax/1000) 
    # print(pltstringEMMAP102)
    # fig = plt.figure(pltnr+2,figsize=plt.figaspect(0.5))
    # fig.suptitle(pltstringEMMAP101) 
    
    # ax1 = fig.add_subplot(121, projection='3d')
    # ax2 = fig.add_subplot(122, projection='3d')
    
    # X,Y=np.meshgrid(OmegaEM,TorqueEM)
    # ax1.plot_surface(X*30/np.pi, Y, EtaEM, rstride=1, cstride=1, cmap='viridis')
    # ax1.set_title('EtaEM')
    # ax2.plot_surface(X, Y, EnergyConsEM,rstride=1, cstride=1, cmap='viridis')
    # ax2.set_title('EnergyConsEM')
    
    # Plot 103
    # pltstringEMMAP103="plt "+str(pltnr+3)+"- Electric Motor Map generated data Motortype="+str(ElectricMotorType)+", MaxPower [kW]= "+str(PowerMax/1000) 
    # print(pltstringEMMAP103)
    
    # fig = plt.figure(pltnr+3,figsize=plt.figaspect(0.5))
    # fig.suptitle(pltstringEMMAP103) 
    # ax1 = fig.add_subplot(121)
    # ax2 = fig.add_subplot(122)
    
    # CS=ax1.contour(X*30/np.pi, Y, EtaEM, levels=[0.7,0.8,0.85,0.9,0.91,0.92,0.93,0.95,0.96], cmap="RdBu_r",linewidths=2)
    # ax1.clabel(CS, inline=1, fontsize=10)
    # ax1.set_title('Efficiency, Optimal points, and Torque limits')
    # ax1.set_ylim(-800,800)
    # ax1.plot(PowerToTorqueEM[:,3]*30/np.pi,PowerToTorqueEM[:,1],color='green', marker='o', linestyle='dashed')
    # ax1.plot(OmegaEM*30/np.pi,TorqueLim,color='black', linestyle='solid')
    # ax1.plot(OmegaEM*30/np.pi,-TorqueLim,color='black', linestyle='solid')
    # ax1.plot(PowerToTorqueEM[:,3]*30/np.pi,PowerToTorqueEM[:,1],color='green', marker='o', linestyle='dashed')
    # ax1.set_xlabel('OmegaEM [rpm]')
    # ax1.set_ylabel('TorqueEM [Nm]')
    # ax2.plot(PowerToTorqueEM[:,0]/1000,PowerToTorqueEM[:,2],color='green', marker='o', linestyle='dashed')
    # ax2.set_xlim(-300,300)
    # ax2.set_title('Efficiency vs Power')
    # ax2.set_xlabel('PowerEM [kW]')
    # ax2.set_ylabel('Efficiency[-]')
    
    class MoIDData:
        pass
    MotorData = MoIDData()
    
    
    vars()["{}Powerloss".format(MotorName)] =  PowerlossEM
    vars()["{}ZEfficiency".format(MotorName)] =  EtaEM 
    #Electric Motor Efficiency map eta (n*n) Omega (n*n)
    setattr(MotorData, "{}ZEfficiency".format(MotorName), EtaEM)
    setattr(MotorData, "{}XOmegaGrid".format(MotorName), X)   #OmegaEM X,Y=np.meshgrid(OmegaEM,TorqueEM)
    setattr(MotorData, "{}YTorqueGrid".format(MotorName), Y)  #TorqueEM
    setattr(MotorData, "{}TorqueLim".format(MotorName), TorqueLim)   #TorqueEM
    setattr(MotorData, "{}OmegaEM".format(MotorName), OmegaEM)      #OmegaEM
    setattr(MotorData, "{}PowerEM".format(MotorName), PowerEM)      #OmegaEM)
    setattr(MotorData, "{}EnergyConsEM".format(MotorName), EnergyConsEM)      #EnergyconsEM)
    setattr(MotorData, "{}Powerloss".format(MotorName), PowerlossEM)
    
    #MotorName='StrtAxleMotor'
    return MotorData
    ## Output Results
    ## Unit - MotorType - MotorSize Power
    #CreateElectricMotorMap(MotorName,MotorType,PMax,TMax,OmegaMax,FieldWR)  
    #class MoIDData:
    #    pass
    #global MotorData #=MoIDData()
#vars()["{}ZEfficiency".format(MotorName)] =  EtaEM    
    # Electric Motor Efficiency map eta (n*n) Omega (n*n)
    # setattr(MotorData, "{}ZEfficiency".format(MotorName), EtaEM)
    # setattr(MotorData, "{}XOmegaEfficiency".format(MotorName), X)   #OmegaEM X,Y=np.meshgrid(OmegaEM,TorqueEM)
    # setattr(MotorData, "{}YTorqueEfficiency".format(MotorName), Y)  #TorqueEM
    # setattr(MotorData, "{}TorqeLim".format(MotorName), TorqueLim)   #TorqueEM
    # setattr(MotorData, "{}OmegaEM".format(MotorName), OmegaEM)      #OmegaEM
    # setattr(MotorData, "{}PowerEM".format(MotorName), PowerEM)      #OmegaEM)

   

