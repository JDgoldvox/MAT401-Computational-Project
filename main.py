import NumericalMethods
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as npy
import math

mass = 5
radius = 2
height = 6
w0 = 1, 2, 3
timeStart, timeEnd = 0, 20
steps = 0.01;

velocity0 = npy.array([0,0,200])

def solve_intertia_tensor(mass, radius, height):
    I1 = (3/20) * mass * ((radius**2) + (1/4 * height**2))
    I2 = (3/20) * mass * ((radius**2) + (1/4 * height**2))
    I3 = (3/20) * mass * (2*radius**2)
    
    tensor = I1, I2, I3
    return tensor

tensor = solve_intertia_tensor(mass, radius, height)

#returns array of euler euqation results
def solve_euler_equation(t,  w0):
    I1, I2, I3 = tensor
    wx0, wy0, wz0 = w0
    
    gamma1 = (I3 - I2) / I1
    gamma2 = (I1 - I3) / I2
    gamma3 = (I2 - I1) / I3
    
    Wx = -gamma1 * wy0 * wz0
    Wy = -gamma2 * wx0 * wz0
    Wz = -gamma3 * wx0 * wy0
    
    return npy.array([Wx, Wy, Wz])

def calculate_angular_speed(wx, wy, wz):
    
    #Create an array of zeros with the same length as wx
    angular_speed = npy.zeros(len(wx)) 
    
    #Solve for magnitude
    for i in range(0, len(wx)):
        x = wx[i] ** 2
        y = wy[i] ** 2
        z = wz[i] ** 2
        
        angular_speed[i] = math.sqrt(x + y + z)
       
    return angular_speed

def main():
    
    #Q1, Q2 ---------------------------------------------------------------
    
    #f = solve_euler_equation
    #x0 = timeStart
    #y0 = w0
    #xe = timeEnd
    #dx = steps
    xx, yy = NumericalMethods.rk4(solve_euler_equation, timeStart, w0, timeEnd, steps);
    
    #Grab values from each axis
    wx_values = yy[:, 0]  #wx
    wy_values = yy[:, 1]  #wy
    wz_values = yy[:, 2]  #wz
    
    angular_speed = calculate_angular_speed(wx_values,wy_values,wz_values)
    
    plt.figure()
    plt.plot(xx, wx_values, label="wx")
    plt.plot(xx, wy_values, label="wy")
    plt.plot(xx, wz_values, label="wz")
    plt.plot(xx, angular_speed, label="angular velocity")
    plt.xlabel("Time")
    plt.ylabel("Angular Velocity")
    
    plt.legend()
    plt.title("Task 2 - Anguler velocity over t [0,20]")
    plt.show()
    
    #Q3, Q4 ---------------------------------------------------------------
    
    displacement, velocity, time = NumericalMethods.implicit_euler(timeStart, velocity0, timeEnd, steps)
    
    #xDisplacement = displacement[:, 0]  #wx
    yDisplacement = displacement[:, 1]  #wy
    #zDisplacement = displacement[:, 2]  #wz
    
    #xVelocity = velocity[:, 0]  #wx
    yVelocity = velocity[:, 1]  #wy
    #zVelocity = velocity[:, 2]  #wz
    
    plt.figure()
    #plt.plot(time, xDisplacement, label="x Displacement")
    plt.plot(time, yDisplacement, label="y Displacement")
    #plt.plot(time, zDisplacement, label="z Displacement")

    #plt.plot(time, xVelocity, label="x velocity")
    plt.plot(time, yVelocity, label="y velocity")
    #plt.plot(time, zVelocity, label="z velocity")
    
    plt.legend()
    plt.title("Vertical velocity and displacement with a = 9.8 and z velocity = 200ms")
    plt.show()
    
    #Q5 ---------------------------------------------------------------
    
    
    
if __name__ == "__main__":
    main()