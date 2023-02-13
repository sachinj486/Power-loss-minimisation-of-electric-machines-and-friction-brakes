"""
Created on Mon Oct 17 10:33:52 2022

Control allocation problem formulation solved using Quadratic programming. 
@ Editor: Sachin Janardhanan
Contributors: Leo Laine, Esteban Gelso
------------------------------------------------------------------------------
This script specifically solves a powerloss minimization problem as a QP problem. 

Vehicle used for this analysis is a representative of a 4X2 tractor, which has 4 service brakes- one on each wheel
and 4 electric machines centrally mounted and driving only the rear axle. 

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

t0=time.time()


" Vehicle and road class description "

def VehicleDataBEVSemi():
    """ Vehicle Data for BEV Tractor SEmitrailer 35 ton """
# Vehicle Properties	Value	           Unit

GCWmass=                  9000            #	[kg] original value = 40000
FrontalArea=              9	               #    [m^2]
AirDragCoeff=	          0.59             #    0.59 - American cabs average #0.36    #    -
RollingResistance=	      0.005	           #    0.0032   #    - How much is bearings/tyre/air/drive?
WheelRadius=	          0.506            #    [m]	
		

#Enivornment Properties		
AirDensity=               1.2              #    [kg/m^3]
Gravity_g=	              9.81	           #    [m/s^2]

# StartAxle and cruise axle power
TotalPeakPowerCapability= 600              #    [kW]



wheelbase=               3.7               #    [m]
lf=                      1.32   #    [m]  #2.12
lr=                      wheelbase-lf      #    [m]

VehicleDataBEVSemi()


###############################################################################
# Road Load forces
###############################################################################
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
# Powerloss formulation functions
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


"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
"Configuration of EM - Axle 11"
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


"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"
"Configuration of EM - Axle 12"
"------------------------------------------------------------------------------"
"------------------------------------------------------------------------------"

GearStart=23

PMax_Ax12=300 #*TwoMotors # calculate the power...              # Selected Value
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
##############################################################################

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
'------------------------------------------------------------------------------------'
Totalpropulsion_req=np.zeros(len(vx))
Totalpropulsion_req_lim=np.zeros(len(vx))
OmegaCrsSts=np.zeros(len(vx))
OmegaStrtSts=np.zeros(len(vx))
v_input=np.zeros(np.array([len(vx),3]))
u=np.zeros(np.array([len(vx),4]))
ax_Actual=np.zeros(len(vx))
Total_Powerloss=np.zeros(len(vx))

CruiseEM_powerloss=np.zeros(len(vx))
StartabilityEM_powerloss=np.zeros(len(vx))

Marker_pattern=['o','v','P','D','s','<','^']
colour=['red','orange','green','blue','magenta','grey']

for i in range(0,len(vx)):
    
   
    '------------------------------------------------------------------------------------'
    'Study efficiency at specific speed'
    '------------------------------------------------------------------------------------'
    OmegaCrsSts[i]=vx[i]*GearCruise/WheelRadius

    "-----------------------------------------------------------------------------"
    # Ax12_EM1
    "-----------------------------------------------------------------------------"
    
    # Capabilities
    [EM1_TorqueMin, EM1_TorqueMax]=MinMaxTorque(OmegaCrsSts[i],  MotorData1.CrsAxlOmegaEM,MotorData1.CrsAxlTorqueLim)
    print('Min max torque cruise axle : ',EM1_TorqueMin, EM1_TorqueMax)
    [EM1TorqueRange_Ax12,EM1ActualPowerloss_Ax12]=Powerlosscurve(MotorData1.CrsAxlXOmegaGrid
                                                          ,MotorData1.CrsAxlYTorqueGrid
                                                          ,MotorData1.CrsAxlPowerloss
                                                          ,OmegaCrsSts[i], EM1_TorqueMin, EM1_TorqueMax)
    
    umin_EM1_Ax12=-CruiseEaxleMaxForceLimits(vx[i],WheelRadius,GearCruise)
    umax_EM1_Ax12= CruiseEaxleMaxForceLimits(vx[i],WheelRadius,GearCruise)
    
    PowerLoss_EM1_Ax12=np.reshape(np.absolute(EM1ActualPowerloss_Ax12),np.size(EM1ActualPowerloss_Ax12))
    

    "-----------------------------------------------------------------------------"
    # Ax12_EM2
    "-----------------------------------------------------------------------------"
    
    '------------------------------------------------------------------------------'
    'Study efficiency at specific speed'
    '------------------------------------------------------------------------------'
    
    OmegaStrtSts[i]=vx[i]*GearStart/WheelRadius
    
    
    
    # Capabilities
    [EM2_TorqueMin, EM2_TorqueMax]=MinMaxTorque(OmegaStrtSts[i], MotorData2.StrtAxlOmegaEM,MotorData2.StrtAxlTorqueLim)
    print('Min max torque cruise axle : ',EM2_TorqueMin, EM2_TorqueMax)
    [EM2TorqueRange_Ax12,EM2ActualPowerloss_Ax12]=Powerlosscurve(MotorData2.StrtAxlXOmegaGrid
                                                          ,MotorData2.StrtAxlYTorqueGrid
                                                          ,MotorData2.StrtAxlPowerloss
                                                          ,OmegaStrtSts[i], EM2_TorqueMin, EM2_TorqueMax)
    
    umin_EM2_Ax12=-StartabilityEaxleMaxForceLimits(vx[i],WheelRadius,GearStart)
    umax_EM2_Ax12= StartabilityEaxleMaxForceLimits(vx[i],WheelRadius,GearStart)
    
    PowerLoss_EM2_Ax12=np.reshape(np.absolute(EM2ActualPowerloss_Ax12),np.size(EM2ActualPowerloss_Ax12))


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
        
    lb=np.array([umin_EM1_Ax12,umin_EM2_Ax12,umin_brk_Ax11,umin_brk_Ax12])
    ub=np.array([umax_EM1_Ax12,umax_EM2_Ax12,umax_brk_Ax11,umax_brk_Ax12])
   

    "-----------------------------------------------------------------------------"
    " Curve fitting of powerlosses"
    "------------------------------------------------------------------------------"
    "-----------------------------"
    # Curve fit Ax12_EM1
    "------------------------------"
    poptCrs, pcovCrs = curve_fit(objective, EM1TorqueRange_Ax12, PowerLoss_EM1_Ax12)
    
    # summarize the parameter values
    a_Ax12_EM1, b_Ax12_EM1, c_Ax12_EM1 = poptCrs
    R2EM1=R2_Value(EM1TorqueRange_Ax12, PowerLoss_EM1_Ax12,objective(EM1TorqueRange_Ax12,a_Ax12_EM1,b_Ax12_EM1,c_Ax12_EM1))
    

    "-----------------------------"
    # Curve fit Ax12_EM2
    "------------------------------"
    poptCrs, pcovCrs = curve_fit(objective, EM2TorqueRange_Ax12, PowerLoss_EM2_Ax12)
    
    # summarize the parameter values
    a_Ax12_EM2, b_Ax12_EM2, c_Ax12_EM2 = poptCrs
    R2EM1=R2_Value(EM1TorqueRange_Ax12, PowerLoss_EM1_Ax12,objective(EM1TorqueRange_Ax12,a_Ax12_EM1,b_Ax12_EM1,c_Ax12_EM1))
    


    "---------------------------------------------------------------------------------------------------------------------------------------"
    " Preparation of the QP formulation using the power loss infomation of the actuators"
    "---------------------------------------------------------------------------------------------------------------------------------------"
    w=2.2
    
    aBrk=1e-5
    
    delta=1e-5
    
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
    if ax[i]>0 and ax[i]>(Fx_lim_Ax12)/GCWmass:
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        print('--------------------------Warning-----------------------------------')
        print('The inputs of ax,ay and mu are not feasible.')
        print('Kindly tune the values so that they are within the axle force capability limit.')
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
    elif np.abs(ax[i])>((Fx_lim_Ax11+Fx_lim_Ax12)/GCWmass):
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        print('--------------------------Warning-----------------------------------')
        print('The inputs of ax,ay and mu are not feasible.')
        print('Kindly tune the values so that they are within the axle force capability limit.')
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        sys.exit()
    


    " Virtual control input - Force request"
    
    if ax[i]>0:
            MaxAccelerationCapability=np.abs((min(np.sum(ub),Fx_lim_Ax12)
                                            -(FxSlope(GCWmass,Slope[i],Gravity_g)+FxRollResistance(GCWmass,vx[i],Gravity_g)+FxAirDrag(vx[i])))
                                            /GCWmass)
            ax_Actual[i]=min(np.abs(ax[i]),MaxAccelerationCapability)
    elif ax[i]<0:
            MaxAccelerationCapability=np.abs((min(np.abs(np.sum(lb)),(Fx_lim_Ax12+Fx_lim_Ax11))
                                            -(FxSlope(GCWmass,Slope[i],Gravity_g)+FxRollResistance(GCWmass,vx[i],Gravity_g)+FxAirDrag(vx[i])))
                                            /GCWmass)
            ax_Actual[i]=np.sign(ax[i])*min(np.abs(ax[i]),MaxAccelerationCapability)
    else:
          ax_Actual[i]=0
          
    Totalpropulsion_req[i]=FxSlope(GCWmass,Slope[i],Gravity_g)+FxRollResistance(GCWmass,vx[i],Gravity_g)+FxAirDrag(vx[i])+(GCWmass*ax_Actual[i])
    
    
    
    # Calculate needed udes service brake from v_input and v_min_el
    v_max_total=(umax_EM1_Ax12+umax_EM2_Ax12+umax_brk_Ax11+umax_brk_Ax12)
    v_min_total=(umin_EM1_Ax12+umin_EM2_Ax12+umin_brk_Ax11+umin_brk_Ax12)
    v_min_el=(umin_EM1_Ax12+umin_EM2_Ax12)

    # Never violate vxmax=B*(umaxi) and vxmin=B*(umini) especially when used as equality constraint        
    if Totalpropulsion_req[i]>=0:
        Totalpropulsion_req_lim[i]=np.array([(np.min((Totalpropulsion_req[i],v_max_total)))])  # Virtual Force inputs
    else:
        Totalpropulsion_req_lim[i]=np.array([(np.max((Totalpropulsion_req[i],v_min_total)))])  
    
    v_input[i]=np.array([Totalpropulsion_req_lim[i],0,0]) 
    
    P=2*np.array([[(a_Ax12_EM1*WheelRadius**2/GearCruise**2), 0, 0, 0], 
                  [0,(a_Ax12_EM2*WheelRadius**2/GearStart**2), 0, 0 ], 
                  [ 0 , 0 , aBrk, 0 ],
                  [ 0, 0 , 0, aBrk]])

    q=np.array([((b_Ax12_EM1*WheelRadius)/GearCruise),((b_Ax12_EM2*WheelRadius)/GearStart),-WheelOmega*WheelRadius,-WheelOmega*WheelRadius])
         
    " Constraints"
    Aeq=np.array([1, 1, 1, 1])    # Fx
    beq=v_input[i,0]
    
    G_ineq=np.array([[0,0,-1,0],  # Fram axle
                    [-1,-1,0,-1],   # Bak axle
                    [0,0,1,0],     # Fram axle
                    [1,1,0,1]])    # Bak axle

    Fz_Ax11=Mu[i]*(Fz[0])
    Fz_Ax12=Mu[i]*(Fz[1])

    Fy_Ax11=(Fy_whl[0])
    Fy_Ax12=(Fy_whl[1])

    H_ineq=np.array([((np.sqrt((Fz_Ax11**2)-(Fy_Ax11**2)))+delta),((np.sqrt((Fz_Ax12**2)-(Fy_Ax12**2)))+delta),((np.sqrt((Fz_Ax11**2)-(Fy_Ax11**2)))+delta),((np.sqrt((Fz_Ax12**2)-(Fy_Ax12**2)))+delta)]) 

    "-----------------------------------------------------------------------------"
    ' QP optimization'
    "-----------------------------------------------------------------------------"
    t0=time.time()
    
    [Fx_EM1_A12,Fx_EM2_A12,Fx_brk_A11,Fx_brk_A12]=solve_qp(P,q, G=G_ineq, h=H_ineq, A=Aeq, b=beq,lb=lb,ub=ub,solver='quadprog')
    
    t1=time.time() 
    
    u[i,]=np.array([Fx_EM1_A12,Fx_EM2_A12,Fx_brk_A11,Fx_brk_A12])

    [CruiseOmega, CruiseTorque, PowerlossCruiseEM]=CruiseAxleEMPowerloss(Fx_EM1_A12,vx[i],GearCruise)
    CruiseEM_powerloss[i]=PowerlossCruiseEM
    
    
    [StartOmega, StartTorque, PowerlossStartEM]=StartabilityAxleEMPowerloss(Fx_EM2_A12, vx[i], GearStart)
    StartabilityEM_powerloss[i]=PowerlossStartEM
   
    Total_Powerloss[i]=CruiseEM_powerloss[i]+StartabilityEM_powerloss[i]+np.abs(Fx_brk_A11*vx[i])+np.abs(Fx_brk_A12*vx[i])
    
    print("--------------------------------------------------------------------")
    print("QP solution - Minimizing power loss ")
    print('QP solution:             ', solve_qp(P,q, G=G_ineq, h=H_ineq, A=Aeq, b=beq,lb=lb,ub=ub,solver='quadprog'))
    #print(f"QP Solver time in Python GIL {(t1-t0)*1e3} milliseconds")
    #print("--------------------------------------------------------------------")
    
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
    axs.arrow(xpos,ypos_rear,0,(Fx_EM1_A12+Fx_EM2_A12)*scaling_factor*1e-3,head_width=0.3, head_length=0.3,length_includes_head='true', color='g',width=0.06)
    
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
    
    
    axs.plot([xpos-1,xpos+1],[ypos_rear+((ub[0]+ub[1])*1e-3*scaling_factor),ypos_rear+((ub[0]+ub[1])*1e-3*scaling_factor)],linewidth=3,color='magenta')
    axs.plot([xpos-1,xpos+1],[ypos_rear+((lb[0]+lb[1])*1e-3*scaling_factor),ypos_rear+((lb[0]+lb[1])*1e-3*scaling_factor)],linewidth=3,color='magenta',label='Rear axle EM based Limit')
   
    
    # axs.legend(loc='upper center',fontsize=16)
    
    # axs.arrow(1,(ypos_rear+ypos_front)/2,0,10e3*scaling_factor*1e-3,head_width=0.7, head_length=0.7, color='black')
    # axs.text(0,ypos_rear,'Rear axle',fontsize=16)
    # axs.text(0,ypos_front,'Front axle',fontsize=16)
    # axs.text(6,0.5,'$OP$'+ str(i+1),fontsize=16)
    axs.text(2.5,ypos_front+4,'SEP',fontsize=24)
    axs.grid()
    
    plt.axis('off')
    plt.show()
    # figure.tight_layout()
    plotname='SingleEAxle_OP'+str(i+1)+'.pdf'
    figure.savefig(plotname)
    





"------------------------------------------------Additional EM plots ----------------------------------------------------------"      
# fig2 = plt.figure()
# host2 = fig2.add_subplot()    

# CS=host2.contour(MotorData1.EACPAxlXOmegaGrid*30/np.pi, 
#                 MotorData1.EACPAxlYTorqueGrid, 
#                 MotorData1.EACPAxlZEfficiency, \
#             levels=[0.7,0.8,0.9,0.95,0.98], cmap="copper",linewidths=0.75)
# host2.clabel(CS, inline=1, fontsize=8)
# host2.plot(MotorData1.EACPAxlOmegaEM*30/np.pi,MotorData1.EACPAxlTorqueLim,color='black',linewidth=2)
# host2.plot(MotorData1.EACPAxlOmegaEM*30/np.pi,-MotorData1.EACPAxlTorqueLim,color='black',linewidth=2)
# for j in range(len(vx)):
#     host2.plot(OmegaEMSts[j]*30/np.pi,u[j,0]*WheelRadius/GearRatio,color=colour[j], 
#           marker=Marker_pattern[j], linestyle='dotted', label='OP'+str(j))
#     host2_legend = plt.legend(loc='upper right')
# host2.set_ylabel('Torque [Nm]',fontsize=12)
# host2.set_xlabel('speed [rpm]',fontsize=12)
# host2.set_ylim(-620,800)
# host2.grid()
# # fig2.savefig('Futuricum_UC1_EM_map.pdf')


    
# EM1Axlforce_lim=(MotorData1.EACPAxlTorqueLim*GearRatio/WheelRadius)
# Vx_EM=(MotorData1.EACPAxlOmegaEM*WheelRadius/GearRatio) 

# Total_EM_force=(EM1Axlforce_lim)*2


# fig4 = plt.figure()
# host4 = fig4.add_subplot() 

# host4.plot(Vx_EM*3.6,EM1Axlforce_lim*1e-3,color='red',linewidth=1.5,label='EM-1')
# host4.plot(Vx_EM*3.6,Total_EM_force*1e-3,color='black',linewidth=1.5,label='Total EM')
# host4.plot(Vx_EM*3.6,Fz_Ax12*np.ones(np.size(Vx_EM))*1e-3,color='orange',linewidth=1.5,linestyle='dashed',label='Friction limit')
# host4.legend(loc='lower right')
# host4.plot(Vx_EM*3.6,-EM1Axlforce_lim*1e-3,color='red',linewidth=1.5)
# host4.plot(Vx_EM*3.6,-Total_EM_force*1e-3,color='black',linewidth=1.5)
# host4.plot(Vx_EM*3.6,-Fz_Ax12*np.ones(np.size(Vx_EM))*1e-3,color='orange',linewidth=1.5,linestyle='dashed')
# for m in range(len(vx)):
#     host4.plot(vx[m]*3.6,Totalpropulsion_req[m]*1e-3,color=colour[m], 
#           marker=Marker_pattern[m], linestyle='dotted', label='Pos'+str(m))

# #host4.set_xlim((0,Vx_Strt[-2]*3.6))
# #host4.set_ylim((-50,50))
# host4.set_xlabel('Vehicle speed [km/h]',fontsize=12)
# host4.set_ylabel('Total Wheel force [kN]',fontsize=12)
# host4.grid()
# # fig4.savefig('Futuricum_UC1_Wheelforce.pdf')




