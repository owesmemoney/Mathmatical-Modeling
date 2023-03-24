# MT2507 Computing Project 
# Part 2: Helper functions

from numpy import *
from pylab import *

# Assume we have the following coupled system:
# dx/dt = x*(1-y)
# dy/dt = y*(1-x)

# Define this as a vector valued function:
def f(x,y):
    return array([x*(1-y),y*(1-x)])

# By inspection we know the two steady states are
W1 = [0,0]
W2 = [1,1]
print('The steady states are:')
print(W1,W2)

# If the steady states were not known one could write a
# 2D Newton-Raphson scheme to find them (Chapter 2 of lecture notes)

print('---------------------------------------')


# Now take the same RHS to integrate the system
def f(x,y,t):
  return array([x*(1-y),y*(1-x)])


# Now define a second order Runge-Kutta scheme to
# integrate the function for one time step:
def RK2step(x,y,t,h,f):
  k1 = h*f(x,y,t)
  k2 = h*f(x+k1[0]/2.0,y+k1[1]/2.0,t+h/2.0)
  return [x+k2[0],y+k2[1]]

#Define the time step and number of steps to use
h=0.05
nsteps=25

# Take a certain initial condition
x = 1.5
y = 1.5
t = 0.0

# Now solve the coupled ODEs for this initial condition
# Do this by advancing the second order R-K scheme for several time steps

# Arrays in which to store solutions
xstore = []
ystore = []
xstore.append(x)
ystore.append(y)

# Use the second order Runge-Kutta method
for n in range(nsteps+1):
  [x,y] = RK2step(x,y,t,h,f)
  t = t+h
  xstore.append(x)
  ystore.append(y)

# Make a plot of this solution
plot(xstore, ystore, 'k-')
xlabel('x')
ylabel('y')
title('One solution')
xlim(0,2)
ylim(0,2)

plt.plot(0.0, 0.0, 'yo')

# This might be heading to the steady state at (1,1)?
# The steady states could be added to the plot (e.g. with coloured dots)









