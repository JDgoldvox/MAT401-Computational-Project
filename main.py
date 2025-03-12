import NumericalMethods

mass = 5
radius = 2
height = 6
wx0, wy0, wz0 = 1, 2, 3
timeStart, timeEnd = 0, 20
steps = 2000;

def main():
    Ix, Iy, Iz = solve_intertia_tensor(mass, radius, height)
    print(Ix)
    print(Iy)
    print(Iz)
    print("---------------------------")
    Wx, Wy, Wz = solve_euler_equation(Ix, Iy, Iz, wx0, wy0, wz0)
    print(Wx)
    print(Wy)
    print(Wz)
    
    #f = Wx
    #x0 = timeStart
    #y0 = wx0
    #xe = timeEnd
    #dx = steps
    xx, yy = NumericalMethods.rk4(Wx, timeStart, wx0, timeEnd, steps);
    
    
    
    xx, yy = NumericalMethods.rk4(Wx, timeStart, wx0, timeEnd, steps);
    xx, yy = NumericalMethods.rk4(Wx, timeStart, wx0, timeEnd, steps);
    

def solve_intertia_tensor(mass, radius, height):
    Ix = (3/20) * mass * ((radius**2) + (1/4 * height**2))
    Iy = (3/20) * mass * ((radius**2) + (1/4 * height**2))
    Iz = (3/20) * mass * (2*radius**2)
    return Ix, Iy, Iz
    

def solve_euler_equation(I1, I2, I3,  wx, wy, wz):
    gamma1 = (I3 - I2) / I1
    gamma2 = (I1 - I3) / I2
    gamma3 = (I2 - I1) / I3
    
    Wx = -gamma1 * wy * wz
    Wy = -gamma2 * wx * wz
    Wz = -gamma3 * wx * wy
    return Wx, Wy, Wz

if __name__ == "__main__":
    main()