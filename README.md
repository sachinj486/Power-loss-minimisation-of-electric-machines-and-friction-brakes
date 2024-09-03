# Power-loss-minimisation-of-electric-machines-and-friction-brakes

This repository contains the python code that used in https://research.chalmers.se/en/publication/534516 for concept analysis of powertrain topologies. The repository has three files representing the topologies and configurations. In addition, a generic electric machine model function is available to configure the machine size. The code is well documented and vehicle engineers should be able to understand the code without much ambiguity. 
Inputs used to specify the operating points representing vehicle operation are fed in as vectors. The outputs are the coordinated actuator outputs by minimising power losses. Furthermore the plots of actuator coordination lumped per axle.

Prerequisites:
Python, qp solvers.
The implementation of qp solvers is done for version <= qp solvers 1.10.0. Hence installing the versions < qp solvers 1.10.0, guarantees the running of the script. Otherwise, the necessary changes are to be made to use latest versions. 

NOTE: The regulation of force request for the MER case has not been tested for all the range of inputs and may need improvements.	

Citation: Kindly refer the Licentiate thesis - On power loss minimisation for heavy vehicles with axle-wise and modular electrical propulsion and friction braking, in case of citation of this code.


