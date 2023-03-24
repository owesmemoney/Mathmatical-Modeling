# MT2507 Computing Project 
# Part 1: Helper functions

# Import the packages we will use
from numpy import *
from pylab import *


# Let's imagine we want to find a root of the function
# f(x) = exp(2x+3sin(x)-1)-1

# First define the function
def f(x):
  return exp(2*x+3*sin(x)-1)-1

# Plot the function to find how many roots exist and roughly where they lie

# Range of x values (chosen by hand)
x = arange(-1.0,1.0,0.01)

# Make the plot
plot(x,f(x))
xlabel('x')
ylabel('f(x)')
title('f(x) on limited range')
# Display the figure.
# Useful to see the y=0 line too
plt.grid(True)
show()

# So we know there is one root for this function, somewhere in (0,0.25)


# Next define a simple Newton-Raphson scheme.Ωxc vbnm,
# This iterates towards a root, for a suitable starting value

# Firt define the derivative of the function, to be used in 
# the Newton-Raphson scheme
def dfdx(x):
  return (2+3*cos(x))*exp(2*x+3*sin(x)-1)
  

#Now make a first attempt at the scheme itself  
def NRroot(x,f,dfdx):
    for i in range(10):
        x = x - f(x)/dfdx(x)
        print(x)
    return x
    
# For initial guess x=0.7 the root is estimated as
r = NRroot(0.7,f,dfdx)
print('Root is at, roughly', r)
  
# The above Newton-Raphson scheme is flawed.
# For example, it does not include any way to determine if
# we have converged to a root, or other stopping criteria.


print('')
print('------------------------------------')
print('')

# Now assume we have a differential equation, 
# dy/dx= f(x,t) =  exp(2x+3sin(x)-1)-1
# i.e. the RHS is the above function (not depending explicitly on t)

# Define this function (could be adapted to include t if need be)
def f(x,t):
  return exp(2*x+3*sin(x)-1)-1


# Now write a second order Runge-Kutta scheme to 
# integrate our function for a single time step:
def RK2step(x,t,h,f):
  k1 = h*f(x,t)
  k2 = h*f(x+k1/2.0,t+h/2.0)
  return x+k2

# To integrate the function (solve the ODE) choose a suitable time step h,
# and number of steps to take. The maximum time t reached will be h*nsteps

# Choice of timestep, h, and total number of steps to take
h = 0.05
nsteps = 30

# Also take an initial condition, here x(0)=0.2
x = 0.2 
t = 0.0

# Create arrays in which to store solutions
xstore = []
tstore = []
xstore.append(x)
tstore.append(t)

# Now solve the ODE using the Runge-Kutta method
for n in range(nsteps+1):
  x = RK2step(x,t,h,f)
  t = t+h
  xstore.append(x)
  tstore.append(t)
  
# Plot the results  
plot(tstore,xstore,'k-')
xlabel('t')
ylabel('x')
title('x(t) for x(0)=0.2')
show()











