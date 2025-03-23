import NumericalMethods
import matplotlib.pyplot as plt
import numpy as npy
import math

#Declare variables to question
mass = 5
radius = 2
height = 6
w0 = 1, 2, 3
timeStart, timeEnd = 0, 20
stepSize = 0.01;
velocity0 = npy.array([0,0,200])
position0 = npy.array([0, 3 * radius/4, 0])

def solve_intertia_tensor(mass, radius, height):
    
    #Solve Intertia given in Q
    I1 = (3/20) * mass * ((radius**2) + (1/4 * height**2))
    I2 = (3/20) * mass * ((radius**2) + (1/4 * height**2))
    I3 = (3/20) * mass * (2*radius**2)
    
    #unify tensor
    tensor = I1, I2, I3
    
    return tensor

tensor = solve_intertia_tensor(mass, radius, height)

#returns array of euler euqation results
def solve_euler_equation(t,  w0):
    I1, I2, I3 = tensor
    wx0, wy0, wz0 = w0
    
    #Solve euler equation
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
    xx, yy = NumericalMethods.rk4(solve_euler_equation, timeStart, w0, timeEnd, stepSize);
    
    #Grab values from each axis
    wx_values = yy[:, 0]  #wx
    wy_values = yy[:, 1]  #wy
    wz_values = yy[:, 2]  #wz
    
    #calculate angular speed
    angular_speed = calculate_angular_speed(wx_values,wy_values,wz_values)
    
    plt.figure()
    plt.plot(xx, wx_values, label="wx")
    plt.plot(xx, wy_values, label="wy")
    plt.plot(xx, wz_values, label="wz")
    plt.plot(xx, angular_speed, label="angular velocity")
    plt.xlabel("Time(s)")
    plt.ylabel("Angular Velocity(w)")
    
    plt.legend()
    plt.title("Task 2 - Anguler velocity over t [0,20]")
    plt.show()
    
    #Q3, Q4 ---------------------------------------------------------------
    
    displacement, velocity, time = NumericalMethods.implicit_euler(timeStart, velocity0, timeEnd, stepSize)
    
    zDisplacement = displacement[:, 2] #Z Axis
    zVelocity = velocity[:, 2] #Z Axis
    
    plt.figure()
    plt.plot(time, zDisplacement, label="z Displacement")
    
    plt.legend()
    plt.title("Vertical displacement of cone vs time")
    plt.xlabel("Time(s)")
    plt.ylabel("Displacement(m)")
    plt.show()
    
    plt.plot(time, zVelocity, label="z velocity")
    
    plt.legend()
    plt.title("Vertical velocity of cone vs time")
    plt.xlabel("Time(s)")
    plt.ylabel("Linear Velocity(ms^-1)")
    plt.show()
    
    #Q5 ---------------------------------------------------------------
    #The full motion of the cone will be a combination of the rotational
    #motion (Tasks 1 & 2) and the translational motion (Tasks 3 & 4).\
        
    #Use your results to determine the trajectory of the point P with initial
    #position (0, 3r/4, 0) as a function of time. Produce 2D projection with times
    #0 < t < 20 s. onto the:
    #(i) x-y-plane, 
    #(ii), x-z-plane
    #(iii) y-z-plane of the position of P for all
    
    #Add all values of rotation into a cumulative summation
    # and multiply by steps to get fraction of this
    # Give us total rotation at each point
    # theta = wt, increased by time
    theta_x = npy.cumsum(wx_values * stepSize) 
    theta_y = npy.cumsum(wy_values * stepSize)
    theta_z = npy.cumsum(wz_values * stepSize)
    
    #declare position arrays
    x_positions = npy.zeros(len(time))
    y_positions = npy.zeros(len(time))
    z_positions = npy.zeros(len(time))
    
    #Do rotational matrices
    #Figure out each Axis of rotation
    for i in range(len(time)):      
        
        #Rotate about X
        Rx = npy.array([[1, 0, 0], 
                        [0, npy.cos(theta_x[i]), -npy.sin(theta_x[i])], 
                        [0, npy.sin(theta_x[i]), npy.cos(theta_x[i])]])
        
        #Rotate about Y
        Ry = npy.array([[npy.cos(theta_y[i]), 0, npy.sin(theta_y[i])], 
                        [0, 1, 0], 
                        [-npy.sin(theta_y[i]), 0, npy.cos(theta_y[i])]])
        
        #Rotate about Z
        Rz = npy.array([[npy.cos(theta_z[i]), -npy.sin(theta_z[i]), 0], 
                        [npy.sin(theta_z[i]), npy.cos(theta_z[i]), 0], 
                        [0, 0, 1]])
    
        #Combine rotation matrices into 1
        R = Rz @ Ry @ Rx  # Dot product all matrices in euler standard order
        
        #Now we found rotation matrix, now dot product translation
        P0_rotated = R @ position0  #Apply the rotation in correct order, Rotation then translation
        
        #Set our positions to rotated t0 position
        x_positions[i] = P0_rotated[0]
        y_positions[i] = P0_rotated[1]
        z_positions[i] = P0_rotated[2]

    #Add displacement from past calculation to my new rotated t0 positions
    #to get the final displaced rotated point
    x_positions += displacement[:, 0]  #x displacement
    y_positions += displacement[:, 1]  #y displacement
    z_positions += displacement[:, 2]  #z displacement
    
    plt.figure()
    plt.plot(x_positions, y_positions)
    plt.xlabel("X position")
    plt.ylabel("Y position")
    plt.title("Trajectory in X-Y Plane")
    plt.show()
    
    plt.figure()
    plt.plot(x_positions, z_positions)
    plt.xlabel("X position")
    plt.ylabel("Z position")        
    plt.title("Trajectory in X-Z Plane")
    plt.show()
    
    plt.figure()
    plt.plot(y_positions, z_positions)
    plt.xlabel("Y position(m)")
    plt.ylabel("Z position")
    plt.title("Trajectory in Y-Z Plane")
    plt.show()
    
if __name__ == "__main__":
    main()