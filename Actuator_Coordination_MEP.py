# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:33:52 2022

Control allocation problem formulation solved using Quadratic programming. 
@ Editor: Sachin Janardhanan
Contributors: Leo Laine, Esteban Gelso
------------------------------------------------------------------------------
This script specifically solves a powerloss minimization problem as a QP problem. 

In specific torque request and service brake request for vehicle motion with propulsion distributed over
two axles is studied here. 

Vehicle used for this analysis is a representative of a 4X4 tractor, which has 4 service brakes- one on each wheel
and 2 electric machines one each per axle.

The cruise/front axle is assumed to efficient using PMSM machine,
while asynchronous machines are configured on the startability/ rear axle. 

But this script can be generalized for other applications and processes.
"""

' Declaration of packages to be used'
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy import interpolate
from qpsolvers import solve_qp
from scipy.optimize import curve_fit
import time
import warnings
warnings.filterwarnings("ignore")


from CreateElectricMotorMap import CreateElectricMotorMap 


" Vehicle and road class description "
def VehicleDataBEVSemi():
    """ Vehicle Data for BEV Tractor """
# Vehicle Properties	Value	           Unit
#
GCWmass=                  9000             #	[kg] 
FrontalArea=              9	              #    [m^2]
AirDragCoeff=	          0.59            #   0.59 
RollingResistance=	      0.005	           #   0.0032  
WheelRadius=	          0.506            #    [m]	
			

#Enivornment Properties		
AirDensity=               1.2              #    [kg/m^3]
Gravity_g=	              9.81	           #    [m/s^2]

# StartAxle and cruise axle power
TotalPeakPowerCapability= 600              #    [kW]


wheelbase=               3.7               #    [m]

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


def MinMaxTorque(ActualOmega,OmegaVector,TorqueLimVector):
    givemeTorque=interpolate.interp1d(OmegaVector,TorqueLimVector,kind='linear')
    TorqueMax=givemeTorque(ActualOmega)
    TorqueMin=-TorqueMax
    return [TorqueMin,TorqueMax]

def Powerlosscurve(OmegaEffGrid,TorqeEffGrid,Powerloss,ActualOmega,TrqMin,TrqMax):
    """Effieciency curve for electric machine at a specific speed
    OmegaEffGrid : gridmap size x*y rotational speed
    TorqueEffGrid : gridmap size x*y Torque
    Efficiency : efficiency on gridmap x*y
    """
    X=OmegaEffGrid
    Y=TorqeEffGrid
    Z=Powerloss
    TrqMin=TrqMin    #
    TrqMax=TrqMax    #
    TorqueRange=np.linspace(TrqMin,TrqMax,101)
    TorqueRange=np.concatenate((TorqueRange[TorqueRange<0], [0], TorqueRange[TorqueRange>0]))#np.linspace(any(min(Y)),any(max(Y)),20)
    giveZ = interpolate.interp2d(X,Y,Z,kind='linear')
    ActualPowerloss=giveZ(ActualOmega,TorqueRange)           
    return [TorqueRange,ActualPowerloss]

def CruiseEaxleMaxForceLimits(vx,WheelRadius,CruiseAxleGearRatio):       
    CruiseOmega=vx/WheelRadius*CruiseAxleGearRatio
    f = interpolate.interp1d(MotorData1.CrsAxlOmegaEM, 
                              MotorData1.CrsAxlTorqueLim)
    CruiseTorquelim=f(CruiseOmega)
    CruiseForce=CruiseTorquelim/WheelRadius*CruiseAxleGearRatio
    return CruiseForce

def StartabilityEaxleMaxForceLimits(vx,WheelRadius,StartabilityAxleGearRatio):       
    StartOmega=vx/WheelRadius*StartabilityAxleGearRatio
    f = interpolate.interp1d(MotorData2.StrtAxlOmegaEM, 
                             MotorData2.StrtAxlTorqueLim)
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
    
    PowerlossCruiseEM=np.abs(float(p(CruiseOmega,CruiseTorque)))
    
    return CruiseOmega, CruiseTorque, PowerlossCruiseEM

def StartabilityAxleEMPowerloss(WheelForceStartabilityAxle,vx,GearStart):
    
   
    #MotorTorque
    StartTorque=WheelForceStartabilityAxle*WheelRadius/GearStart
    StartOmega=vx/WheelRadius*GearStart

    
    p=interpolate.interp2d(MotorData2.StrtAxlXOmegaGrid, 
                           MotorData2.StrtAxlYTorqueGrid,
                           MotorData2.StrtAxlPowerloss)
   
    PowerlossStartEM=np.abs(float(p(StartOmega,StartTorque)))
    
    return StartOmega, StartTorque, PowerlossStartEM

"------------------------------------------------------------------------------------------------"
" Quadratic Curve fitting funnction for the Power loss approaximation and R_2 value of the error"
"------------------------------------------------------------------------------------------------"

def objective(x, a, b, c):
 	return a * x **2 + b * x + c

def R2_Value(xdata, ydata, fdata):
    '''
    Parameters
    ----------
    xdata : vector length n
        original x values
    ydata : vector length n
        original y values
    fdata : vector length n
        curve_fit model y values for xdata

    Returns
    -------
    R2_Value : scalar
        Determination coefficient or R2 value for curve fit

    '''
    residuals=ydata-fdata
    SumOfSquares_residuals=np.sum(residuals**2)
    SumOfTotalSquares_total=np.sum((ydata-np.mean(ydata))**2)
    R2_Value=1-SumOfSquares_residuals/SumOfTotalSquares_total
    return R2_Value

"-----------------------------------------------------------------------------"

"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
"Configuration of EM on Axle 11"
"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
GearCruise=12

# Cruise Axle PM motors
pltnr=100
PMax_Ax11=300 # Max power per EM  
OmegaMax_Ax11=10000 # rpm
FieldWR=2.5
TMax=0              # Calculated when set to zero 
MotorType='PermanentMagnetMotor'
MotorName='CrsAxl'

MotorData1=CreateElectricMotorMap(pltnr,MotorName,MotorType,PMax_Ax11,TMax,OmegaMax_Ax11,FieldWR)


"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
"Configuration of EM on - Axle 12"
"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"

GearStart=23

# Startability Axle IM motors
PMax_Ax12=300 #
OmegaMax_Ax12=13000 # rpm
FieldWR=2.5
TMax=0    # Calculated when set to zero 
MotorType='AsynchronousMotor'
MotorName='StrtAxl'

MotorData2=CreateElectricMotorMap(pltnr,MotorName,MotorType,PMax_Ax12,TMax,OmegaMax_Ax12,FieldWR)



'-----------------------------------------------------------------------------'
' Inputs to the CA problem'
'------------------------------------------------------------------------------'
###############################################################################
'Inputs to define different driving states or manouevre'

ServiceBrakeMaxTorque=20e3 # [Nm]  # per wheel

vx=np.array([10,70,50,40])/3.6 # [m/s] - vehicle speed
Slope=np.array([5.0,-3.0,0.0,-2.0])# [%] - Gradient
ax=np.array([0.15,-0.25,-0.3,-0.15])*9.81# [m/s^2] Longitudinal acceleration request
ay=np.array([0.0,0.0,-0.3,0.1])*9.81# [m/s^2] Lateral acceleration input
Mu=np.array([0.7,0.5,0.6,0.3]) # Friction co-efficient

Marker_pattern=['o','o','v','v','P','D','D','s','<','^']
colour=['red','red','orange','orange','green','blue','blue','magenta']

 
'------------------------------------------------------------------------------------'
' Initialize loop varaibles'
'------------------------------------------------------------------------------------'
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
    'Ax11'
    '------------------------------------------------------------------------------'
    
    # Machine speed 
    OmegaCrsSts[i]=vx[i]*GearCruise/WheelRadius  
    
    # EM Capabilities
    [CrsTorqueMin_Ax11, CrsTorqueMax_Ax11]=MinMaxTorque(OmegaCrsSts[i], MotorData1.CrsAxlOmegaEM,MotorData1.CrsAxlTorqueLim)
    print('Min max torque cruise axle : ',CrsTorqueMin_Ax11, CrsTorqueMax_Ax11)
    [CrsTorqueRange_Ax11,ActualPowerloss_Ax11]=Powerlosscurve(MotorData1.CrsAxlXOmegaGrid
                                                          ,MotorData1.CrsAxlYTorqueGrid
                                                          ,MotorData1.CrsAxlPowerloss
                                                          ,OmegaCrsSts[i], CrsTorqueMin_Ax11, CrsTorqueMax_Ax11)
    
    umin_EM_Ax11=-CruiseEaxleMaxForceLimits(vx[i],WheelRadius,GearCruise)
    umax_EM_Ax11=CruiseEaxleMaxForceLimits(vx[i],WheelRadius,GearCruise)
    
    PowerLoss_EM_Ax11=np.reshape(np.absolute(ActualPowerloss_Ax11),np.size(ActualPowerloss_Ax11))
  
    
    '------------------------------------------------------------------------------'
    'Ax12'
    '------------------------------------------------------------------------------'
    # EM speed
    OmegaStrtSts[i]=vx[i]*GearStart/WheelRadius
    
    # Capabilities
    [StrtTorqueMin_Ax12, StrtTorqueMax_Ax12]=MinMaxTorque(OmegaStrtSts[i], MotorData2.StrtAxlOmegaEM,MotorData2.StrtAxlTorqueLim)
    print('Min max torque Start axle : ',StrtTorqueMin_Ax12, StrtTorqueMax_Ax12)
    
    [StrtTorqueRange_Ax12,ActualPowerloss_Ax12]=Powerlosscurve(MotorData2.StrtAxlXOmegaGrid
                                                          ,MotorData2.StrtAxlYTorqueGrid
                                                          ,MotorData2.StrtAxlPowerloss
                                                          ,OmegaStrtSts[i], StrtTorqueMin_Ax12, StrtTorqueMax_Ax12)
    
    
    umin_EM_Ax12=-StartabilityEaxleMaxForceLimits(vx[i],WheelRadius,GearStart)
    umax_EM_Ax12=StartabilityEaxleMaxForceLimits(vx[i],WheelRadius,GearStart)
    
    PowerLoss_EM_Ax12=np.reshape(np.absolute(ActualPowerloss_Ax12),np.size(ActualPowerloss_Ax12))
      
    
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
    
    "-----------------------------------------------------------------------------"
    " Curve fitting of powerlosses"
    "------------------------------------------------------------------------------"
    "-----------------------------"
    # Curve fit Ax11
    "------------------------------"
    poptCrs, pcovCrs = curve_fit(objective, CrsTorqueRange_Ax11, PowerLoss_EM_Ax11)
    # summarize the parameter values
    a_Ax11, b_Ax11, c_Ax11 = poptCrs
    #print('y = %.5f * x^2 + %.5f * x + %.5f' % (a_Ax11_W11, b_Ax11_W11, c_Ax11_W11))
    R2Crs=R2_Value(CrsTorqueRange_Ax11, PowerLoss_EM_Ax11,objective(CrsTorqueRange_Ax11,a_Ax11,b_Ax11,c_Ax11))
    #print('R2 value for fit :', R2Crs)
    
    
    "-----------------------------"
    # Curve fit Ax12
    "------------------------------"
    poptStrt, pcovStrt = curve_fit(objective, StrtTorqueRange_Ax12, PowerLoss_EM_Ax12)
    # summarize the parameter values
    a_Ax12, b_Ax12, c_Ax12 = poptStrt
    #print('y = %.5f * x^2 + %.5f * x + %.5f' % (a_Ax12_W22, b_Ax12_W22, c_Ax12_W22))
    R2Strt=R2_Value(StrtTorqueRange_Ax12, PowerLoss_EM_Ax12,objective(StrtTorqueRange_Ax12,a_Ax12,b_Ax12,c_Ax12))
    #print('R2 value for fit :', R2Strt)
    
    
    "---------------------------------------------------------------------------------------------------------------------------------------"
    " Preparation of the QP formulation using the power loss infomation of the actuators"
    "---------------------------------------------------------------------------------------------------------------------------------------"
      
    aBrk=1e-5
    
    delta=0
    
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
    if np.abs(ax[i])>(Fx_lim_Ax11+Fx_lim_Ax12)/GCWmass:
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        print('--------------------------Warning-----------------------------------')
        print('The inputs of ax,ay and mu are not feasible.')
        print('Kindly tune the values so that they are within the axle force capability limit.')
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        sys.exit()
    

    " Virtual control input - Force request generation and control"
    
    if ax[i]>0:
            MaxAccelerationCapability=np.abs((min(ub[0],Fx_lim_Ax11)+min(ub[1],Fx_lim_Ax12)
                                            -(FxSlope(GCWmass,Slope[i],Gravity_g)+FxRollResistance(GCWmass,vx[i],Gravity_g)+FxAirDrag(vx[i])))
                                            /GCWmass)
            ax_Actual[i]=min(np.abs(ax[i]),MaxAccelerationCapability)
    elif ax[i]<0:
            MaxAccelerationCapability=np.abs(((min(np.abs(np.sum(lb[0]+lb[2])),Fx_lim_Ax11)+(min(np.abs(np.sum(lb[1]+lb[3])),Fx_lim_Ax12)))
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
    
    P=2*np.array([[(a_Ax11*WheelRadius**2/GearCruise**2), 0, 0, 0 ], 
                  [0, (a_Ax12*WheelRadius**2/GearStart**2) , 0 , 0],
                  [0, 0 , aBrk, 0],
                  [0, 0 , 0, aBrk]])
    
    q=np.array([(b_Ax11*WheelRadius/GearCruise), (b_Ax12*WheelRadius/GearStart),-WheelOmega*WheelRadius,-WheelOmega*WheelRadius])
    
    
    " Constraints"
    
    Aeq=np.array([1, 1, 1, 1])    # Fx
    beq=v_input[i,0]
    
    G_ineq=np.array([[-1,0,-1,0],  # Fram axle
                    [0,-1,0,-1],
                    [1,0,1,0],
                    [0,1,0,1]]) # Bak axlen

    H_ineq=np.array([((np.sqrt((Fz_Ax11**2)-(Fy_Ax11**2)))-delta),((np.sqrt((Fz_Ax12**2)-(Fy_Ax12**2)))-delta),((np.sqrt((Fz_Ax11**2)-(Fy_Ax11**2)))-delta),((np.sqrt((Fz_Ax12**2)-(Fy_Ax12**2)))-delta)]) 
       
    "-----------------------------------------------------------------------------"
    ' QP optimization'
    "-----------------------------------------------------------------------------"
    t0=time.time()
    
    [Fx_EM_A11,Fx_EM_A12,Fx_brk_A11,Fx_brk_A12]=solve_qp(P,q, G=G_ineq, h=H_ineq, A=Aeq, b=beq,lb=lb,ub=ub,solver='quadprog')
    
    t1=time.time() 
    
    u[i,]=np.array([Fx_EM_A11,Fx_EM_A12,Fx_brk_A11,Fx_brk_A12])
    
    [CruiseOmega, CruiseTorque, PowerlossCruiseAxle]=CruiseAxleEMPowerloss(Fx_EM_A11,vx[i],GearCruise)
    CruiseAxle_powerloss[i]=PowerlossCruiseAxle+np.abs(Fx_brk_A11*vx[i])
    
    
    [StartOmega, StartTorque, PowerlossStartAxle]=StartabilityAxleEMPowerloss(Fx_EM_A12, vx[i], GearStart)
    StartabilityAxle_powerloss[i]=PowerlossStartAxle+np.abs(Fx_brk_A12*vx[i])
    
    Total_Powerloss[i]=CruiseAxle_powerloss[i]+StartabilityAxle_powerloss[i]
    
    print("--------------------------------------------------------------------")
    print("QP solution - Minimizing power loss ")
    print('QP solution:             ', solve_qp(P,q, G=G_ineq, h=H_ineq, A=Aeq, b=beq,lb=lb,ub=ub,solver='quadprog'))
    
    
    print("Verification of the solution and constraints")
    print("--------------------------------------------------------------------")
    print('Lower Bounds of actuators:             ', lb)
    print("--------------------------------------------------------------------")
    print('Upper Bounds of actuators:             ', ub)
    print("--------------------------------------------------------------------")
    print('Aboslute Longitudinal force capability on Rear axle ', H_ineq[1])
    print("--------------------------------------------------------------------")
    print('Absolute Longitudinal force capability on Front axle ', H_ineq[0])
    print("--------------------------------------------------------------------")
    print('Longitudinal force from electric machines ', u[i,0]+u[i,1])
    print("--------------------------------------------------------------------")
    print('Longitudinal force from Brakes ', u[i,2]+u[i,3])
    print("--------------------------------------------------------------------")
    print('Total powerloss ',Total_Powerloss[i])
    print("--------------------------------------------------------------------")
    
    " Friction circle plot"
    figure, axs = plt.subplots(figsize=(4,12)) # 
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
    axs.plot([xpos-1,xpos+1],[ypos_front+((lb[0])*1e-3*scaling_factor),ypos_front+((lb[0])*1e-3*scaling_factor)],linewidth=3,color='orange',label='Front axle EM-limit')
    axs.plot([xpos-1,xpos+1],[ypos_rear+((ub[1])*1e-3*scaling_factor),ypos_rear+((ub[1])*1e-3*scaling_factor)],linewidth=3,color='magenta')
    axs.plot([xpos-1,xpos+1],[ypos_rear+((lb[1])*1e-3*scaling_factor),ypos_rear+((lb[1])*1e-3*scaling_factor)],linewidth=3,color='magenta',label='Rear axle EM based Limit')

    axs.text(2.5,ypos_front+4,'MEP',fontsize=24)
    
    plt.axis('off')
    plt.show()
    # figure.tight_layout()
    # plotname='MultipleEAxle_OP'+str(i+1)+'.pdf'
    # figure.savefig(plotname)
    
    
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
# CS=host3.contour(MotorData2.StrtAxlXOmegaGrid*30/np.pi, 
#                 MotorData2.StrtAxlYTorqueGrid, 
#                 MotorData2.StrtAxlZEfficiency, \
#             levels=[0.7,0.8,0.9,0.95,0.98], cmap="copper",linewidths=0.5)
# host3.clabel(CS, inline=1, fontsize=8)
# host3.plot(MotorData2.StrtAxlOmegaEM *30/np.pi,MotorData2.StrtAxlTorqueLim,color='black',linewidth=2)
# host3.plot(MotorData2.StrtAxlOmegaEM*30/np.pi,-MotorData2.StrtAxlTorqueLim,color='black',linewidth=2)
# for k in range(len(vx)):
#     host3.plot(OmegaStrtSts[k]*30/np.pi,u[k,1]*WheelRadius/GearStart,color=colour[k], 
#           marker=Marker_pattern[k], linestyle='dotted', label='$OP$'+str(k+1))
#     host3_legend = plt.legend(loc='upper right')    
# host3.set_ylabel('Torque [Nm]',fontsize=12)
# host3.set_xlabel('speed [rpm]',fontsize=12)
# host3.set_ylim(-685,685)
# host3.grid()
# fig3.savefig('Volverine_UC1_StartAxl_map.pdf')



    
# vx_ref=np.linspace(1,50,200)

# CrsAxlforce_lim=(MotorData1.CrsAxlTorqueLim*GearCruise/WheelRadius)+(MotorData2.CrsAxlTorqueLim*GearCruise/WheelRadius)
# Vx_Crs=(MotorData1.CrsAxlOmegaEM*WheelRadius/GearCruise) 

# StrtAxlforce_lim=(MotorData3.StrtAxlTorqueLim*GearStart/WheelRadius)+(MotorData4.StrtAxlTorqueLim*GearStart/WheelRadius)
# Vx_Strt=(MotorData3.StrtAxlOmegaEM*WheelRadius/GearStart) # kmph

# StrtAxlforce_interp=np.interp(vx_ref,Vx_Strt,StrtAxlforce_lim,right=0)

# CrsAxlforce_interp=np.interp(vx_ref,Vx_Crs,CrsAxlforce_lim,right=0)

# Total_EM_force=(CrsAxlforce_interp+StrtAxlforce_interp)


# fig4 = plt.figure()
# host4 = fig4.add_subplot() 
# #pltstring700='plt 700 - Traction diagram with operating points'
# #print(pltstring700)

# host4.plot(vx_ref[vx_ref<Vx_Crs[-1]]*3.6,CrsAxlforce_interp[CrsAxlforce_interp>0]*1e-3,color='red',linewidth=1.5,label='Cruise Axle')
# host4.plot(vx_ref[vx_ref<Vx_Strt[-1]]*3.6,StrtAxlforce_interp[StrtAxlforce_interp>0]*1e-3,color='blue',linewidth=1.5,label='Startability Axle')
# host4.plot(vx_ref[vx_ref<Vx_Crs[-1]]*3.6,Total_EM_force[Total_EM_force>0]*1e-3,color='black',linewidth=1.5,label='Total Axle')
# host4.legend(loc='lower right')
# host4.plot(vx_ref[vx_ref<Vx_Crs[-1]]*3.6,-CrsAxlforce_interp[CrsAxlforce_interp>0]*1e-3,color='red',linewidth=1.5)
# host4.plot(vx_ref[vx_ref<Vx_Strt[-1]]*3.6,-StrtAxlforce_interp[StrtAxlforce_interp>0]*1e-3,color='blue',linewidth=1.5)
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





