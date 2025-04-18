# CFD: 1D Euler Finite Volume Solver

This program can be used to solve both subsonic flow and supersonic flow in a variable area pipe. First and second order reconstruction is implemented in this code along with a minmod limiter for 2nd order reconstruction with discontinuities in the solution. 

## Inputs 

Xbounds - Domain Size
num_points - Number of solution points in the domain. These are evenly distributed. 
sigma - CFL Number
gamma - Specific Heat Ratios
order - Desired Order of Error
N - Number of time steps to run
islimiteron - Use minmod limiter? true/false
flow - Flow problem 1/2
isplot -  Plot during sim? true/false

There are two flow problems implemented into the code. 
1) Subsonic flow with an exit Mach number of 0.4. Assuming that the
density at the exit is 1, and pressure at the exit is 1/gamma and the speed of sound at the
exit is 1, and the exit velocity is 0.4.

2) Transonic flow with a shock wave down stream of the throat. Assuming that the inlet
Mach number is 0.2006533, inlet density 1 and pressure 1/gamma and exit pressure
0.6071752.
