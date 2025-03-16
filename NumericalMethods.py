import matplotlib.pyplot as plt
import numpy as npy

#implicit euler -----------------------------------------------------
#xe = last time
#x0 = time start
#dx = interval
#y0 = velocity
def implicit_euler(x0, y0, xe, dx):
    
    #calculate steps
    stepCount = num_steps(x0, xe, dx)
    
    #create initial arrays
    time = npy.linspace(x0, xe, stepCount) 
    displacement = npy.zeros((len(time), len(y0))) 
    velocity = npy.zeros((len(time), len(y0))) 
    
    #set initial arrays
    velocity[0] = y0
    gravity = npy.array([0,9.81,0])
    
    #do implicit euler
    for i in range(1, len(displacement)):
        velocity[i] = velocity[i-1] - dx * gravity
        displacement[i] = displacement[i-1] + dx * velocity[i]
        
    return displacement, velocity, time
    
#rk2 -----------------------------------------------------
def rk2(f, x0, y0, xe, dx):
    stepCount = num_steps(x0, xe, dx)
    xx = npy.linspace(x0, xe, stepCount)
    yy = npy.zeros(len(xx))
    yy[0] = y0
    for i in range(1, len(xx)):
        k1 = dx * f(xx[i-1], yy[i-1]) 
        k2 = dx * f(xx[i-1] + dx/2, yy[i-1] + k1/2)
        yy[i] = yy[i -1] + k2
        
    return xx, yy
    
#rk4 -----------------------------------------------------
#xe = last time
#x0 = time start
#dx = interval
#y0 = value (anguler velocity)
#f = function
def rk4(f, x0, y0, xe, dx):
    stepCount = num_steps(x0, xe, dx)
    x_values = npy.linspace(x0, xe, stepCount) #create interval array
    y_values = npy.zeros((len(x_values), len(y0))) #2d array
    y_values[0] = y0 
    
    for i in range(1, len(x_values)):
        k1 = dx * f(x_values[i-1], y_values[i-1]) 
        k2 = dx * f(x_values[i-1] + dx/2, y_values[i-1] + k1/2)
        k3 = dx * f(x_values[i-1] + dx/2, y_values[i-1] + k2/2)
        k4 = dx * f(x_values[i-1] + dx, y_values[i-1] + k3)
        y_values[i] = y_values[i -1] + ((k1 + 2 * k2 + 2 * k3 + k4) / 6)
    return x_values, y_values
    
def num_steps(x0, xe, dx):
    return int((xe - x0) / dx) + 1
