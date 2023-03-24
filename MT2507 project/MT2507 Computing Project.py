# -*- coding: utf-8 -*-
"""
ID: 190026112
MT2507: Mathematical Modelling Computing Project
"""
from numpy import *
from pylab import *
import numpy as np
import math

print('Part 1: Question 1')
# Let's imagine we want to find a root of the function
# First define the function
def f(x):
  return 0.01-0.4*x+x*x/(1+x*x)

# Plot the function to find how many roots exist and roughly where they lie

# Range of x values (chosen by hand)
x = arange(-3.0,3.0,0.01)

# Make the plot
plot(x,f(x))
xlabel('x')
ylabel('f(x)')
title('f(x) on limited range')
# Display the figure.
# Useful to see the y=0 line too
plt.grid(True)
show()
print('')
print('------------------------------------')
print('')

# So we know there were two roots for this function, somewhere between (0,0.5)
# and the third root is around to 2 in x-axis


# Next define a simple Newton-Raphson scheme.Ωxc vbnm,
# This iterates towards a root, for a suitable starting value

# Firt define the derivative of the function, to be used in 
# the Newton-Raphson scheme

def dfdx(x):
  return 2*x/(x*x+1)-2*x*x*x/((x*x+1)*(x*x+1))-0.4
  

#Now make a first attempt at the scheme itself  
def NRroot(x,f,dfdx):
    for i in range(10):
        x = x - f(x)/dfdx(x)
        print(x)
    return x
    
# For initial guess x=0.7 the root is estimated as
r = NRroot(0.7,f,dfdx)
print('Root is at, roughly', format(r, '.4f'))
  
# The above Newton-Raphson scheme is flawed.
# For example, it does not include any way to determine if
# we have converged to a root, or other stopping criteria.


print('')
print('------------------------------------')
print('')

# For initial guess x=1.7 the root is estimated as
r = NRroot(1.7,f,dfdx)
print('Root is at, roughly', format(r, '.4f'))
# try initial guesses 1.7 for root of f(x)

print('')
print('------------------------------------')
print('')

# For initial guess x=0.1 the root is estimated as
r = NRroot(0.1,f,dfdx)
print('Root is at, roughly', format(r, '.4f'))
  
print('')
print('------------------------------------')
print('')

print('Part 1: Question 2')
#Part 1: 2

# Now assume we have a differential equation, 
# Define this function (could be adapted to include t if need be)
def f(x,t):
  return 0.01-0.4*x+x*x/(1+x*x)


# Now write a fourth order Runge-Kutta scheme to 
# integrate our function for a single time step:
def RK4step(x,t,h,f):
  k1 = h*f(x,t)
  k2 = h*f(x+k1/2.0,t+h/2.0)
  k3 = h*f(x+k2/2.0,t+h/2.0)
  k4 = h*f(x+k3,t+h)
  return x+(k1+2*k2+2*k3+k4)/6.0

# To integrate the function (solve the ODE) choose a suitable time step h,
# and number of steps to take. The maximum time t reached will be h*nsteps

# Choice of timestep, h, and total number of steps to take
h = 0.05
nsteps = 50

# Also take an initial condition, here x(0)=0.45
x = 0.45
t = 0.0

# Create arrays in which to store solutions
xstore = []
tstore = []
xstore.append(x)
tstore.append(t)

# Now solve the ODE using the Runge-Kutta method
for n in range(nsteps+1):
  x = RK4step(x,t,h,f)
  t = t+h
  xstore.append(x)
  tstore.append(t)
  
# Plot the results  
plot(tstore,xstore,'k-')
xlabel('t')
ylabel('x')
title('x(t) for x(0)=0.45')
show()

# take another initial condition, here x(0)=0.5
x=0.5
t=0.0

xstore = []
tstore = []
xstore.append(x)
tstore.append(t)

# Now solve the ODE using the Runge-Kutta method
for n in range(nsteps+1):
  x = RK4step(x,t,h,f)
  t = t+h
  xstore.append(x)
  tstore.append(t)
  
# Plot the results  
plot(tstore,xstore,'k-')
xlabel('t')
ylabel('x')
title('x(t) for x(0)=0.5')
show()

print('')
print('------------------------------------')
print('')

print('Part 2: Question 1')
# Assume we have the following coupled system:
# dx/dt = x*(1-x-0.1*exp(0.1*y))
# dy/dt = y*(1-y-0.2*exp(0.3*x))
# define the domain of x, y
x = arange(0,100,0.01)
y = arange(0,100,0.01)

# Define this as a vector valued function:
def f(x,y):
    return array([x*(1-x-0.1*exp(0.1*y)),y*(1-y-0.2*exp(0.3*x))])
# To get J(x,y) take the derivative by hand 
# define fucntion g, h, and get the differentiating of g, h with respect to x, y
[x, y] = [1, 1]     #choose the initial condition of x, y
for i in range(10):
    g = 1-x-0.1*math.exp(0.1*y)     
    h = 1-y-0.2*math.exp(0.3*x)   
    gx = -1     
    gy = -(math.exp(y/10))/100            
    hx = -3*(math.exp(3*x/10))/50            
    hy = -1    
    # form J(x,y)
    J = np.array([[gx, gy], [hx, hy]]) 
    # get the inverse     
    invJ = np.linalg.inv(J) 
    # calculated using 2D Newton-Raphson Method:                  
    D = np.dot(invJ, np.array([[g, h]]).T)
    [x, y] = [x, y] - D.T[0]
    
    print([x, y])
    
[x, y] = np.around([x, y], decimals=3)
print('Root is at, roughly', [x,y])    

# By inspection we know the four steady states are
# the first three steady states could be calculated by hand 
W1 = [0,0]
W2 = [0,0.8]
W3 = [0.9,0]
W4 = [x,y]

print('The steady states are:')
print(W1,W2,W3,W4)


print('')
print('------------------------------------')
print('')

print('Part 2: Question 2')
# Now take the same RHS to integrate the system
def f(x,y,t):
  return array([x*(1-x-0.1*exp(0.1*y)),y*(1-y-0.2*exp(0.3*x))])

# Now define a fourth order Runge-Kutta scheme to
# integrate the function for one time step:
def RK4step(x,y,t,h,f):
  k1 = h*f(x,y,t)
  k2 = h*f(x+k1[0]/2.0,y+k1[1]/2.0,t+h/2.0)
  k3 = h*f(x+k2[0]/2.0,y+k2[1]/2.0,t+h/2.0)
  k4 = h*f(x+k3[0],y+k3[1],t+h)
  return [x+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6.0,y+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6.0]

#Define the time step and number of steps to use
h=0.05
nsteps=1000



# Now solve the coupled ODEs for this initial condition
# Do this by advancing the fourth order R-K scheme for several time steps

# define a solve function to simplify the process of taking various initial condition 
def solve(x0,y0):
    x=x0
    y=y0
    t=0.0
    # Arrays in which to store solutions
    xstore = []
    ystore = []
    xstore.append(x)
    ystore.append(y)
    # Use the fourth order Runge-Kutta method
    for n in range(nsteps+1):
        [x,y] = RK4step(x,y,t,h,f)
        t = t+h
        xstore.append(x)
        ystore.append(y)
        
    plot(xstore, ystore, 'k-')
        

# take x, y values close to the four boundries of the graph 
for x in arange(0.01,1.3,0.15):
    y=0.01
    solve(x,y)
    
for y in arange(0.01,1.3,0.15):
    x=0.01
    solve(x,y)
    
for y in arange(0.01,1.3,0.15):
    x=1.2
    solve(x,y)
    
for x in arange(0.01,1.3,0.15):
    y=1.2
    solve(x,y)

#put arrows on the lines
x_values,y_values=meshgrid(arange(0.0,1.2,0.01),arange(0.0,1.4,0.01))  
x=x_values*(1-x_values-0.1*exp(0.1*y_values))  
y=y_values*(1-y_values-0.2*exp(0.3*x_values)) 
streamplot(x_values,y_values,x,y)

# Make a plot of this solution
xlabel('x')
ylabel('y')
title('One solution')
xlim(0,1.2)
ylim(0,1.2)

# The steady states could be added to the plot  with coloured dots
plt.plot(0.89, 0.74, 'yo') 
plt.plot(0, 0.8, 'yo') 
plt.plot(0.9, 0, 'yo') 
plt.plot(0, 0, 'yo') 


