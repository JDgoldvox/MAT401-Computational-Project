import matplotlib.pyplot as plt
import numpy as npy

stepSize = 40.50;
total = 100.0;

def func(x,y):
    return 5.0 * x ** 4

def ans(x,y):
    return x**5

#euler -----------------------------------------------------
#Euler's method 
#dy = m * dx
#y(n+1) = y(n) + m * dx
def euler(f, x0, y0, xe, dx):
    stepCount = num_steps(x0, xe, dx)
    xx = npy.linspace(x0, xe, stepCount) #Add + dx to include last element
    yy = npy.zeros(len(xx))
    yy[0] = y0
    for i in range(1, len(xx)):
        yy[i] = yy[i-1] + dx * f(xx[i-1], yy[i-1])
    return xx, yy
    
#rk2 -----------------------------------------------------
# k1 = dx * f(x,y)
# k2 = dx * f(x(n) + dx/2, y(n) + k1/2)
# y(n) + 1 = y(n) + k2
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
# k1 = dx * f(x(n), y(n))
# k2 = dx * f(x(n) + dx/2, y(n) + k1/2)
# k3 = dx * f(x(n) + dx/2, y(n) + k2/2)
# k4 = dx * f(x(n) + dx, y(n) + k3)
# y(n+1) = y(n) + (k1 + 2k2, + 2k3 + k4 / 6)
#xe = last time, dx = interval
#x0 = time
def rk4(f, x0, y0, xe, dx):
    stepCount = num_steps(x0, xe, dx)
    xx = npy.linspace(x0, xe, stepCount)
    yy = npy.zeros(len(xx))
    yy[0] = y0
    for i in range(1, len(xx)):
        k1 = dx * f(xx[i-1], yy[i-1]) 
        k2 = dx * f(xx[i-1] + dx/2, yy[i-1] + k1/2)
        k3 = dx * f(xx[i-1] + dx/2, yy[i-1] + k2/2)
        k4 = dx * f(xx[i-1] + dx, yy[i-1] + k3)
        yy[i] = yy[i -1] + ((k1 + 2 * k2 + 2 * k3 + k4) / 6)
    return xx, yy
    
def num_steps(x0, xe, dx):
    return int((xe - x0) / dx) + 1
    
euler_result_x, euler_result_y = euler(func, 0, 0, total, stepSize)
rk2_result_x, rk2_result_y = rk2(func, 0, 0, total, stepSize)
rk4_result_x, rk4_result_y = rk4(func, 0, 0, total, stepSize)

stepCount =  num_steps(0, total, stepSize)
xa = npy.linspace(0, total, stepCount)
ya = xa**5.0
    
#Plot -----------------------------------------------------
#plt.plot(xa, ya, label = "real") #real answer
#plt.plot(euler_result_x, euler_result_y, label = "euler")
#plt.plot(rk2_result_x, rk2_result_y, label = "rk2")
#plt.plot(rk4_result_x, rk4_result_y, label = "rk4")

#print(".", ya[-1])
#print(".", 20**5)
#print(".", xa[-1])
#print("Euler error: ", abs(ya[-1] - euler_result_y[-1]))
#print("rk2 error: ", abs(ya[-1] - rk2_result_y[-1]))
#print("rk4 error: ", abs(ya[-1] - rk4_result_y[-1]))

#plt.legend()
#plt.title(label="Euler vs rk2 vs rk4 results")
#plt.show()

#plt.plot(ya, abs(ya-euler_result_y), label = "euler")
#plt.plot(ya, abs(ya-rk2_result_y), label = "rk2")
#plt.plot(ya, abs(ya-rk4_result_y), label = "rk4")

#plt.legend()
#plt.title(label="Euler vs rk2 vs rk4 error difference")
#plt.show()