# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:33:52 2022

Control allocation problem formulation solved using Quadratic programming. 
@ Editor: Sachin Janardhanan
Contributors: Leo Laine
------------------------------------------------------------------------------
This script specifically solves a powerloss minimization problem as a QP problem. 

In specific torque request and service brake request for vehicle motion with propulsion distributed over
two axles is studied here. The request is coordinated equally on both the axles.

Vehicle used for this analysis is a representative of a 4X4 tractor, which has 4 service brakes- one on each wheel
and 2 electric machines one each per axle.

The cruise/front axle is assumed to efficient using PMSM machines on left and right wheel
while asynchronous machines are configured on the startability/ rear axle. 

But this script can be generalized for other applications and processes.
"""

' Declaration of packages to be used'
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy import interpolate
import time
import warnings
warnings.filterwarnings("ignore")

from CreateElectricMotorMap import CreateElectricMotorMap 


" Vehicle and road class description "
def VehicleDataBEVSemi():
    """ Vehicle Data for BEV Tractor SEmitrailer 35 ton """
#
# General Properties
# Vehicle Properties	Value	           Unit
#
GCWmass=                  9000            #	[kg] original value = 40000
FrontalArea=              9	               #    [m^2]
AirDragCoeff=	          0.59             #   0.59 - American cabs average #0.36    #    -
RollingResistance=	      0.005	           #   0.0032   #    - How much is bearings/tyre/air/drive?
WheelRadius=	            0.506            #    [m]	
		
	
#Performance Criterias		
Vxmax=	                  125/3.6	       #    [km/h->m/s] 
VxCruise=	              85/3.6	       #    [km/h->m/s] 
MaxCruiseSlope=	          2.5	               #    [%]
MinCruiseSlope=          -2.5                #    [%]
MaxStartabilitySlope=	  10	           #    [%] 
 
#Enivornment Properties		
AirDensity=               1.2              #    [kg/m^3]
Gravity_g=	              9.81	           #    [m/s^2]

# StartAxle and cruise axle power
TotalPeakPowerCapability= 600              #    [kW]
BatteryEnergyCapability=  6*90  #4*90      #    [kWh]
UsedSOC_MinCapability=    15               #    [%]
UsedSOC_MaxCapability=    85               #    [%]
kWhToMJ=                  3.6              #    [%]

RangeSlope=              0                 #    [%]
RangeVxAverage=          85/3.6            #    [km/h->m/s] 
RangeTarget=             800               #    [km]

GearRatio=               26

wheelbase=               3.7               #    [m]

cf=                      2*230e3           #    [N/rad] 
cr=                      4*200e3           #    [N/rad]

lf=                      1.32    #    [m]
lr=                      wheelbase-lf      #    [m]

VehicleDataBEVSemi()

###############################################################################
# Road Load forces
#
def FxAirDrag(Vx):
    """Airdrag resistance force FxAirDrag"""
    return 0.5*AirDragCoeff*FrontalArea*AirDensity*np.square(Vx)
def FxRollResistance(Mass,Vx,Gravity_g):
    """Rolling resistance force FxRollRes"""
    return Mass*Gravity_g*RollingResistance*(1-1/np.exp(0.5*np.abs(Vx)))*np.sign(Vx)
def FxSlope(Mass,Slope,Gravity_g):
    """Slope resistance force FxSlope"""
    return Mass*Gravity_g*np.sin(np.arctan(Slope/100))


###############################################################################################################################
# Solving using pure QP and powerlosses
###############################################################################################################################
class MoIDData:
        pass
MotorData1 = MoIDData()
MotorData2 = MoIDData()


def CruiseEaxleMaxForceLimits(vx,WheelRadius,CruiseAxleGearRatio):       
    CruiseOmega=vx/WheelRadius*CruiseAxleGearRatio
    f = interpolate.interp1d(MotorData1.CrsAxlOmegaEM, 
                              MotorData1.CrsAxlTorqueLim)
    CruiseTorquelim=f(CruiseOmega)
    CruiseForce=CruiseTorquelim/WheelRadius*CruiseAxleGearRatio
    return CruiseForce

def StartabilityEaxleMaxForceLimits(vx,WheelRadius,StartabilityAxleGearRatio):       
    StartOmega=vx/WheelRadius*StartabilityAxleGearRatio
    f = interpolate.interp1d(MotorData3.StrtAxlOmegaEM, 
                             MotorData3.StrtAxlTorqueLim)
    StartTorquelim=f(StartOmega)
    StartForce=StartTorquelim/WheelRadius*StartabilityAxleGearRatio
    return StartForce
    
def BrakeMaxForceLimits(ServiceBrakeMaxTorque,WheelRadius,mass):
    BrakeTorquelim= min(ServiceBrakeMaxTorque,mass*Gravity_g*WheelRadius)
    BrakeForcelim=BrakeTorquelim/WheelRadius
    return BrakeForcelim

def CruiseAxleEMPowerloss(WheelForceCruiseAxle,vx,GearCruise):
    
    #MotorTorque
    CruiseTorque=WheelForceCruiseAxle*WheelRadius/GearCruise
    CruiseOmega=vx/WheelRadius*GearCruise
    
    p=interpolate.interp2d(MotorData1.CrsAxlXOmegaGrid, 
                           MotorData1.CrsAxlYTorqueGrid,
                           MotorData1.CrsAxlPowerloss)
    
    PowerlossCruiseAxle=np.abs(float(p(CruiseOmega,CruiseTorque)))
    
    return CruiseOmega, CruiseTorque, PowerlossCruiseAxle

def StartabilityAxleEMPowerloss(WheelForceStartabilityAxle,vx,GearStart):
   
    #MotorTorque
    StartTorque=WheelForceStartabilityAxle*WheelRadius/GearStart
    StartOmega=vx/WheelRadius*GearStart

    
    p=interpolate.interp2d(MotorData3.StrtAxlXOmegaGrid, 
                           MotorData3.StrtAxlYTorqueGrid,
                           MotorData3.StrtAxlPowerloss)
   
    PowerlossStartAxle=np.abs(float(p(StartOmega,StartTorque)))
    
    return StartOmega, StartTorque, PowerlossStartAxle


"-----------------------------------------------------------------------------"

"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
"Initialization of EM - Axle 11"
"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
GearCruise=12

# Cruise Axle PM motors
pltnr=100
PMax_Ax11=300 # *Max power per EM  
OmegaMax_Ax11=10000 # rpm
FieldWR=2.5
TMax=0              # Calculated when set to zero 
MotorType='PermanentMagnetMotor'
MotorName='CrsAxl'

MotorData1=CreateElectricMotorMap(pltnr,MotorName,MotorType,PMax_Ax11,TMax,OmegaMax_Ax11,FieldWR)
MotorData2=MotorData1

"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
"Initialization of EM - Axle 12"
"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"

GearStart=23

PMax_Ax12=300 #*TwoMotors # calculate the power...              # Selected Value
OmegaMax_Ax12=13000 # rpm
FieldWR=2.5
TMax=0    # Calculated when set to zero 
MotorType='AsynchronousMotor'
MotorName='StrtAxl'

MotorData3=CreateElectricMotorMap(pltnr,MotorName,MotorType,PMax_Ax12,TMax,OmegaMax_Ax12,FieldWR)
MotorData4=MotorData3


'-----------------------------------------------------------------------------'
' Inputs to the CA problem'
'------------------------------------------------------------------------------'
###############################################################################
'Inputs to define different driving states or manouevre'


ServiceBrakeMaxTorque=20e3 # [Nm]  # per wheel

vx=np.array([10,70,50,40])/3.6 
Slope=np.array([5.0,-3.0,0.0,-2.0])# 
ax=np.array([0.15,-0.25,-0.3,-0.15])*9.81# 
ay=np.array([0.0,0.0,-0.3,0.1])*9.81#
Mu=np.array([0.7,0.5,0.6,0.3])

Marker_pattern=['o','o','v','v','P','D','D','s','<','^']
colour=['red','red','orange','orange','green','blue','blue','magenta']


'------------------------------------------------------------------------------------'

' Initialize loop varaibles'
Totalpropulsion_req=np.zeros(len(vx))
Totalpropulsion_req_lim=np.zeros(len(vx))
OmegaCrsSts=np.zeros(len(vx))
OmegaStrtSts=np.zeros(len(vx))
v_input=np.zeros(np.array([len(vx),3]))
u=np.zeros(np.array([len(vx),4]))
ax_Actual=np.zeros(len(vx))
Total_Powerloss=np.zeros(len(vx))

CruiseAxle_powerloss=np.zeros(len(vx))
StartabilityAxle_powerloss=np.zeros(len(vx))

Marker_pattern=['o','v','P','D','s','<','^']
colour=['red','orange','green','blue','magenta','grey']


for i in range(0,len(vx)):
    
    
    '------------------------------------------------------------------------------'
    'Study efficiency at specific speed'
    '------------------------------------------------------------------------------'
    OmegaCrsSts[i]=vx[i]*GearCruise/WheelRadius
    
    "-----------------------------------------------------------------------------"
    # Ax11
    "-----------------------------------------------------------------------------"
    
    # Capabilities
   
    umin_EM_Ax11=-CruiseEaxleMaxForceLimits(vx[i],WheelRadius,GearCruise)
    umax_EM_Ax11=CruiseEaxleMaxForceLimits(vx[i],WheelRadius,GearCruise)

 
      
    
    '------------------------------------------------------------------------------'
    'Study efficiency at specific speed'
    '------------------------------------------------------------------------------'
    
    OmegaStrtSts[i]=vx[i]*GearStart/WheelRadius
    
    
    "-----------------------------------------------------------------------------"
    # Ax12
    "-----------------------------------------------------------------------------"
    
    # Capabilities

    umin_EM_Ax12=-StartabilityEaxleMaxForceLimits(vx[i],WheelRadius,GearStart)
    umax_EM_Ax12=StartabilityEaxleMaxForceLimits(vx[i],WheelRadius,GearStart)
    
       
    
    "------------------------------------------------------------------------------"
    "------------------------------------------------------------------------------"
    "Capabilities of Brakes" 
    "------------------------------------------------------------------------------"
    "------------------------------------------------------------------------------"
    
    umin_brk_Ax11=-2*BrakeMaxForceLimits(ServiceBrakeMaxTorque,WheelRadius,GCWmass)
    umax_brk_Ax11=0;


    umin_brk_Ax12=-2*BrakeMaxForceLimits(ServiceBrakeMaxTorque,WheelRadius,GCWmass)
    umax_brk_Ax12=0;
    
    "------------------------------------------------------------------------------"
    " Capabilities of the actuators"
    "------------------------------------------------------------------------------"
    
    lb=np.array([umin_EM_Ax11,umin_EM_Ax12,umin_brk_Ax11,umin_brk_Ax12])
    ub=np.array([umax_EM_Ax11,umax_EM_Ax12,umax_brk_Ax11,umax_brk_Ax12])
    
    "---------------------------------------------------------------------------------------------------------------------------------------"
    " Preparation of the QP formulation using the power loss infomation of the actuators"
    "---------------------------------------------------------------------------------------------------------------------------------------"
    w=2.2
    
    aBrk=1e-5
    
    delta=1e-3
    
    Fz=np.array([GCWmass*lr/(wheelbase),GCWmass*lf/(wheelbase)])*9.81
    
    WheelOmega=vx[i]/WheelRadius
    
    "------------------------------------------------------------------------------"
    "Lateral capability of axles"
    "------------------------------------------------------------------------------"

    Fz=np.array([GCWmass*lr/(wheelbase),GCWmass*lf/(wheelbase)])*9.81
    
    Fy_f=((GCWmass*lr)/wheelbase)*(ay[i])  
    Fy_r=((GCWmass*lf)/wheelbase)*(ay[i])
   
    Fy_whl=np.array([Fy_f,Fy_r])
    
    Fz_Ax11=Mu[i]*(Fz[0])
    Fz_Ax12=Mu[i]*(Fz[1])

    Fy_Ax11=(Fy_whl[0])
    Fy_Ax12=(Fy_whl[1])
    
    Fx_lim_Ax11=(np.sqrt((Fz_Ax11**2)-(Fy_Ax11**2)))
    Fx_lim_Ax12=(np.sqrt((Fz_Ax12**2)-(Fy_Ax12**2)))


    " -----------------------------------------------------------------------------"
    " Check if the ax and ay inputs produce requests that are within friction limit"
    "------------------------------------------------------------------------------"
    if (np.abs(ax[i])>(Fx_lim_Ax11+Fx_lim_Ax12)/GCWmass) or (np.abs(ax[i])>(2*(min(Fx_lim_Ax11,Fx_lim_Ax12))*1/GCWmass)):
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        print('--------------------------Warning-----------------------------------')
        print('The inputs of ax,ay and mu are not feasible.')
        print('Kindly tune the values so that they are within the axle force capability limit.')
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        print('Exited the operating , due to infeasible solution')
        sys.exit()    
    
    " Virtual control input - Force request"
    
    if ax[i]!=0:
            MaxAccelerationCapability=((min(2*min(Fx_lim_Ax11,Fx_lim_Ax12),min(np.abs(ub[0]),Fx_lim_Ax11)+min(np.abs(ub[1]),Fx_lim_Ax12))
                                            -(FxSlope(GCWmass,Slope[i],Gravity_g)+FxRollResistance(GCWmass,vx[i],Gravity_g)+FxAirDrag(vx[i])))
                                            /GCWmass)
            ax_Actual[i]=np.sign(ax[i])*min(np.abs(ax[i]),MaxAccelerationCapability)
    else:
        ax_Actual[i]=0
   
    Totalpropulsion_req[i]=FxSlope(GCWmass,Slope[i],Gravity_g)+FxRollResistance(GCWmass,vx[i],Gravity_g)+FxAirDrag(vx[i])+(GCWmass*ax_Actual[i])
       
   
    # Calculate needed udes service brake from v_input and v_min_el
    v_max_total=(umax_EM_Ax11+umax_EM_Ax12+umax_brk_Ax11+umax_brk_Ax12)
    v_min_total=(umin_EM_Ax11+umin_EM_Ax12+umin_brk_Ax11+umin_brk_Ax12)
    v_min_el=(umin_EM_Ax11+umin_EM_Ax12)

    # Never violate vxmax=B*(umaxi) and vxmin=B*(umini) especially when used as equality constraint        
    if Totalpropulsion_req[i]>=0:
        Totalpropulsion_req_lim[i]=np.array([(np.min((Totalpropulsion_req[i],v_max_total)))])  # Virtual Force inputs
    else:
        Totalpropulsion_req_lim[i]=np.array([(np.max((Totalpropulsion_req[i],v_min_total)))])  
    
    v_input[i]=np.array([Totalpropulsion_req_lim[i],0,0]) 
    
    
    " Constraints"
    Flim_Ax11=Mu[i]*(Fz[0])
    Flim_Ax12=Mu[i]*(Fz[1])
    
    Fy_Ax11=(Fy_whl[0])
    Fy_Ax12=(Fy_whl[1])
    
    # Find the minimum between two axles
    Flim_vehicle=min((np.sqrt((Flim_Ax11**2)-(Fy_Ax11**2))),(np.sqrt((Flim_Ax12**2)-(Fy_Ax12**2))))
    
    
    if (v_input[i][0]<0)and(v_input[i][0]>=(umin_EM_Ax11+umin_EM_Ax12)):
        Fx_EM_A11=max(0.5*v_input[i][0],(np.sign(v_input[i][0])*(Flim_vehicle)))
        Fx_EM_A12=max(0.5*v_input[i][0],(np.sign(v_input[i][0])*(Flim_vehicle)))
        Fx_brk_A11=0
        Fx_brk_A12=0
    elif (v_input[i][0]<(umin_EM_Ax11+umin_EM_Ax12)):
            Fx_EM_A11=max(0.5*v_input[i][0],(np.sign(v_input[i][0])*(Flim_vehicle)))
            Fx_EM_A12=max(0.5*v_input[i][0],(np.sign(v_input[i][0])*(Flim_vehicle)))
            Fx_brk_A11=max(0.5*(v_min_el-v_input[i][0]),np.sign(v_input[i][0])*(Fx_EM_A11+Flim_vehicle))
            Fx_brk_A12=Fx_brk_A11
    else:
        Fx_EM_A11=0.5*v_input[i][0]
        Fx_EM_A12=Fx_EM_A11
        Fx_brk_A11=0
        Fx_brk_A12=0
 
   
    
    "-----------------------------------------------------------------------------"
    ' QP optimization'
    "-----------------------------------------------------------------------------"
    t0=time.time()
    
    [Fx_EM_A11,Fx_EM_A12,Fx_brk_A11,Fx_brk_A12]
    
    t1=time.time() 
    
    u[i,]=np.array([Fx_EM_A11,Fx_EM_A12,Fx_brk_A11,Fx_brk_A12])
    
    [CruiseOmega, CruiseTorque, PowerlossCruiseAxle]=CruiseAxleEMPowerloss(Fx_EM_A11,vx[i],GearCruise)
    CruiseAxle_powerloss[i]=PowerlossCruiseAxle+np.abs(Fx_brk_A11*vx[i])
    
    
    [StartOmega, StartTorque, PowerlossStartAxle]=StartabilityAxleEMPowerloss(Fx_EM_A12, vx[i], GearStart)
    StartabilityAxle_powerloss[i]=PowerlossStartAxle+np.abs(Fx_brk_A12*vx[i])
    
    Total_Powerloss[i]=CruiseAxle_powerloss[i]+StartabilityAxle_powerloss[i]
    
    print("--------------------------------------------------------------------")
    print("QP solution - Minimizing power loss ")
   
    
    print("Verification of the solution and constraints")
    print("--------------------------------------------------------------------")
    print('Lower Bounds of actuators:             ', lb)
    print("--------------------------------------------------------------------")
    print('Upper Bounds of actuators:             ', ub)
    print("--------------------------------------------------------------------")
    print('Aboslute Longitudinal force capability on Rear axle ', (np.sqrt((Flim_Ax11**2)-(Fy_Ax11**2))))
    print("--------------------------------------------------------------------")
    print('Absolute Longitudinal force capability on Front axle ',(np.sqrt((Flim_Ax12**2)-(Fy_Ax12**2))))
    print("--------------------------------------------------------------------")
    print('Longitudinal force from electric machines ', u[i,0]+u[i,1])
    print("--------------------------------------------------------------------")
    print('Longitudinal force from Brakes ', u[i,2]+u[i,3])
    print("--------------------------------------------------------------------")
    print('Total powerloss ',Total_Powerloss[i])
    print("--------------------------------------------------------------------")
    
    
    
    # Friction circle plot
    figure, axs = plt.subplots(figsize=(4,12))
    axs.set_xlim(0,6)
    axs.set_ylim(0,20)
    scaling_factor=7e-2
    
    xpos=3
    ypos_front=16
    ypos_rear=5
   
    Frct_cricle_Ax11 = plt.Circle((xpos, ypos_front),Mu[i]*(Fz[0])*1e-3*scaling_factor,fill=0,linewidth=2) 
    Frct_cricle_Ax12 = plt.Circle((xpos, ypos_rear),Mu[i]*(Fz[1])*1e-3*scaling_factor,fill=0,linewidth=2) 
  
    axs.add_patch(Frct_cricle_Ax12)
    axs.add_patch(Frct_cricle_Ax11)
    
    axs.arrow(xpos,ypos_rear,0,Fx_EM_A12*scaling_factor*1e-3,head_width=0.3, head_length=0.3,length_includes_head='true', color='g',width=0.06)
    axs.arrow(xpos,ypos_front,0,Fx_EM_A11*scaling_factor*1e-3,head_width=0.3, head_length=0.3,length_includes_head='true', color='g',width=0.06)
    
    
    # axs.arrow(1,1,0,5000*scaling_factor*1e-3,head_width=0.2, head_length=0.2, color='g',length_includes_head='true')
    # axs.text(1.25,1,'=5$kN$',fontsize=12)
    
    if ax[i]<0:
        if Fx_brk_A11<-1e-5:
            axs.arrow(xpos,ypos_front,0,(Fx_brk_A11)*1e-3*scaling_factor,head_width=0.3, head_length=0.3,length_includes_head='true',color='red',width=0.06)
        if Fx_brk_A12<-1e-5:
            axs.arrow(xpos,ypos_rear,0,(Fx_brk_A12)*1e-3*scaling_factor,head_width=0.3, head_length=0.3,length_includes_head='true',color='red',width=0.06)
    
    if ay[i]!=0:
        axs.arrow(xpos,ypos_rear,Fy_Ax12*scaling_factor*1e-3,0,head_width=0.3, head_length=0.3,length_includes_head='true',color='grey',width=0.06)
        axs.arrow(xpos,ypos_front,Fy_Ax11*scaling_factor*1e-3,0,head_width=0.3, head_length=0.3,length_includes_head='true',color='grey',width=0.06)

    axs.set_aspect( 1 )
    
    axs.plot([xpos-1,xpos+1],[ypos_front+((ub[0])*1e-3*scaling_factor),ypos_front+((ub[0])*1e-3*scaling_factor)],linewidth=3,color='orange')
    axs.plot([xpos-1,xpos+1],[ypos_front+((lb[0])*1e-3*scaling_factor),ypos_front+((lb[0])*1e-3*scaling_factor)],linewidth=3,color='orange',label='Front axle EM based Limit')
    axs.plot([xpos-1,xpos+1],[ypos_rear+((ub[1])*1e-3*scaling_factor),ypos_rear+((ub[1])*1e-3*scaling_factor)],linewidth=3,color='magenta')
    axs.plot([xpos-1,xpos+1],[ypos_rear+((lb[1])*1e-3*scaling_factor),ypos_rear+((lb[1])*1e-3*scaling_factor)],linewidth=3,color='magenta',label='Rear axle EM based Limit')

    
    # axs.legend(loc='upper center',fontsize=12)
    
   
    # axs.arrow(1,(ypos_rear+ypos_front)/2,0,10e3*scaling_factor*1e-3,head_width=0.7, head_length=0.7, color='black')
    # axs.text(0,ypos_rear,'Rear axle',fontsize=16)
    # axs.text(0,ypos_front,'Front axle',fontsize=16)
    # axs.text(6,0.5,'$OP$'+ str(i+1),fontsize=16)
    axs.text(2.5,ypos_front+4,'MER',fontsize=24)
    
    plt.axis('off')
    plt.show()
    # figure.tight_layout()
    plotname='MultipleEAxle_Baseline_OP'+str(i+1)+'.pdf'
    figure.savefig(plotname)
    
    
"------------------------------------------------Additional EM plots ----------------------------------------------------------"      
    
# Powerloss_1EM_Ax11=np.zeros(len(vx))    

# fig2 = plt.figure()
# host2 = fig2.add_subplot()    
# CS=host2.contour(MotorData1.CrsAxlXOmegaGrid*30/np.pi, 
#                 MotorData1.CrsAxlYTorqueGrid, 
#                 MotorData1.CrsAxlZEfficiency, \
#             levels=[0.7,0.8,0.9,0.95,0.98], cmap="copper",linewidths=0.75)
# host2.clabel(CS, inline=1, fontsize=8)
# host2.plot(MotorData1.CrsAxlOmegaEM*30/np.pi,MotorData1.CrsAxlTorqueLim,color='black',linewidth=2)
# host2.plot(MotorData1.CrsAxlOmegaEM*30/np.pi,-MotorData1.CrsAxlTorqueLim,color='black',linewidth=2)
# for j in range(len(vx)):
#     host2.plot(OmegaCrsSts[j]*30/np.pi,u[j,0]*WheelRadius/GearCruise,color=colour[j], 
#           marker=Marker_pattern[j], linestyle='dotted', label='$OP$'+str(j+1))
#     host2_legend = plt.legend(loc='upper right')  


# host2.set_ylabel('Torque [Nm]',fontsize=12)    
# host2.set_xlabel('speed [rpm]',fontsize=12)
# host2.set_ylim(-800,800)
# host2.grid()
# # fig2.savefig('Volverine_UC1_CruiseAxl_map.pdf')


# fig3 = plt.figure()
# host3 = fig3.add_subplot()
# #pltstring600='plt 600 - Start Axle efficiency, Operating points during driving cycle'
# #print(pltstring600)
# #host3.figure(600) 
# CS=host3.contour(MotorData3.StrtAxlXOmegaGrid*30/np.pi, 
#                 MotorData3.StrtAxlYTorqueGrid, 
#                 MotorData3.StrtAxlZEfficiency, \
#             levels=[0.7,0.8,0.9,0.95,0.98], cmap="copper",linewidths=0.5)
# host3.clabel(CS, inline=1, fontsize=8)
# host3.plot(MotorData3.StrtAxlOmegaEM *30/np.pi,MotorData3.StrtAxlTorqueLim,color='black',linewidth=2)
# host3.plot(MotorData3.StrtAxlOmegaEM*30/np.pi,-MotorData3.StrtAxlTorqueLim,color='black',linewidth=2)
# for k in range(len(vx)):
#     host3.plot(OmegaStrtSts[k]*30/np.pi,u[k,1]*WheelRadius/GearStart,color=colour[k], 
#           marker=Marker_pattern[k], linestyle='dotted', label='$OP$'+str(k+1))
#     host3_legend = plt.legend(loc='upper right')    
# host3.set_ylabel('Torque [Nm]',fontsize=12)
# host3.set_xlabel('speed [rpm]',fontsize=12)
# host3.set_ylim(-685,685)
# host3.grid()
# # fig3.savefig('Volverine_UC1_StartAxl_map.pdf')



    
# # vx_ref=np.linspace(1,50,200)

# # CrsAxlforce_lim=(MotorData1.CrsAxlTorqueLim*GearCruise/WheelRadius)+(MotorData2.CrsAxlTorqueLim*GearCruise/WheelRadius)
# # Vx_Crs=(MotorData1.CrsAxlOmegaEM*WheelRadius/GearCruise) 

# # StrtAxlforce_lim=(MotorData3.StrtAxlTorqueLim*GearStart/WheelRadius)+(MotorData4.StrtAxlTorqueLim*GearStart/WheelRadius)
# # Vx_Strt=(MotorData3.StrtAxlOmegaEM*WheelRadius/GearStart) # kmph

# # StrtAxlforce_interp=np.interp(vx_ref,Vx_Strt,StrtAxlforce_lim,right=0)

# # CrsAxlforce_interp=np.interp(vx_ref,Vx_Crs,CrsAxlforce_lim,right=0)

# # Total_EM_force=(CrsAxlforce_interp+StrtAxlforce_interp)


# # fig4 = plt.figure()
# # host4 = fig4.add_subplot() 
# # #pltstring700='plt 700 - Traction diagram with operating points'
# # #print(pltstring700)

# # host4.plot(vx_ref[vx_ref<Vx_Crs[-1]]*3.6,CrsAxlforce_interp[CrsAxlforce_interp>0]*1e-3,color='red',linewidth=1.5,label='Cruise Axle')
# # host4.plot(vx_ref[vx_ref<Vx_Strt[-1]]*3.6,StrtAxlforce_interp[StrtAxlforce_interp>0]*1e-3,color='blue',linewidth=1.5,label='Startability Axle')
# # host4.plot(vx_ref[vx_ref<Vx_Crs[-1]]*3.6,Total_EM_force[Total_EM_force>0]*1e-3,color='black',linewidth=1.5,label='Total Axle')
# # host4.legend(loc='lower right')
# # host4.plot(vx_ref[vx_ref<Vx_Crs[-1]]*3.6,-CrsAxlforce_interp[CrsAxlforce_interp>0]*1e-3,color='red',linewidth=1.5)
# # host4.plot(vx_ref[vx_ref<Vx_Strt[-1]]*3.6,-StrtAxlforce_interp[StrtAxlforce_interp>0]*1e-3,color='blue',linewidth=1.5)
# host4.plot(vx_ref[vx_ref<Vx_Crs[-1]]*3.6,-Total_EM_force[Total_EM_force>0]*1e-3,color='black',linewidth=1.5)
# for m in range(len(vx)):
#     host4.plot(vx[m]*3.6,Totalpropulsion_req[m]*1e-3,color=colour[m], 
#           marker=Marker_pattern[m], linestyle='dotted', label='Pos'+str(m))
#     # host4_legend = host4.legend(loc='upper right') 

# host4.set_xlim((0,Vx_Strt[-2]*3.6))
# host4.set_ylim((-50,50))
# host4.set_xlabel('Vehicle speed [km/h]',fontsize=12)
# host4.set_ylabel('Total Wheel force [kN]',fontsize=12)
# host4.grid()
# # fig4.savefig('Volverine_UC1_Wheelforce.pdf')





