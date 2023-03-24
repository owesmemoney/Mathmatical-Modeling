# RK2D.py: Plot out time series of integration steps of a 2D ODE
#      to illustrate the fourth-order Runge-Kutta method.
#
# For a 2D ODE
#     dx/dt = f(x,y)
#     dy/dt = g(x,y)
# See RKTwoD() below for how the fourth-order Rungle-Kutta method integrates.
#

# Import plotting routines
from pylab import *

# The van der Pol 2D ODE
#   x_dot = y
#   y_dot = u * (0.1 - x * x) * y - x
#    u in [0,10] is the control parameter
def VDPXDot(u,x,y):
	return x*(1-x-0.1*exp(0.1*y))

def VDPYDot(u,x,y):
	return y*(1-y-0.2*exp(0.3*x))

# 2D Fourth-Order Runge-Kutta Integrator
def RKTwoD(r,x,y,f,g,dt):
	k1x = dt * f(r,x,y)
	k1y = dt * g(r,x,y)
	k2x = dt * f(r,x + k1x / 2.0,y + k1y / 2.0)
	k2y = dt * g(r,x + k1x / 2.0,y + k1y / 2.0)
	k3x = dt * f(r,x + k2x / 2.0,y + k2y / 2.0)
	k3y = dt * g(r,x + k2x / 2.0,y + k2y / 2.0)
	k4x = dt * f(r,x + k3x,y + k3y)
	k4y = dt * g(r,x + k3x,y + k3y)
	x = x + ( k1x + 2.0 * k2x + 2.0 * k3x + k4x ) / 6.0
	y = y + ( k1y + 2.0 * k2y + 2.0 * k3y + k4y ) / 6.0
	return x,y

# Simulation parameters
# Integration time step
dt = 0.05
#
# Control parameter for the van der Pol 2D ODE:
u = 2.0
# Set up arrays of iterates for three different initial conditions
x1 = [ 0.1]
y1 = [ 0.1]
x2 = [ 0.8]
y2 = [ 0.8]
x3 = [ 1.5]
y3 = [ 1.5]
plot(x1,y1,'bo')
plot(x2,y2,'ro')
plot(x3,y3,'go')
# Time
t  = [ 0.0]
# The number of time steps to integrate over
N = 1000

# The main loop that generates the orbit, storing the states
for n in range(0,N):
  # at each time step calculate new x(t) and y(t)
  # and append to lists x1 and y1, x2 and y2
  x,y = RKTwoD(u,x1[n],y1[n],VDPXDot,VDPYDot,dt)
  x1.append(x)
  y1.append(y)
  x,y = RKTwoD(u,x2[n],y2[n],VDPXDot,VDPYDot,dt)
  x2.append(x)
  y2.append(y)
  x,y = RKTwoD(u,x3[n],y3[n],VDPXDot,VDPYDot,dt)
  x3.append(x)
  y3.append(y)
  t.append(t[n] + dt)

# Setup the parametric plot
xlabel('x(t)') # set x-axis label
ylabel('y(t)') # set y-axis label
title('4th order Runge-Kutta Method: van der Pol ODE at u = ' + str(u)) # set plot title
axis('equal')
axis([-0.5,2.0,-0.5,2.0])
# Plot the trajectory in the phase plane
plot(0.9,0,'b')
plot(0,0.8,'r')
plot(0.89,0.84,'g')

# Use command below to save figure
#savefig('RK2D', dpi=600)

# Display the plot in a window
show()
