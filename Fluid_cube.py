import numpy as np
import math as m
import matplotlib.cm as cm
import time
import matplotlib.pyplot as plt


class FluidCube:
    def __init__(self, size, dt, diff, visc):
        self.size = size                           #length in one dimension
        self.dt = dt                               #timestep
        self.diff = diff                           #diffusion constant
        self.visc = visc                           #viscosity
        self.s = np.zeros((size,size,size))        #unsure
        self.dens = np.zeros((size,size,size))     #density
        self.vx = np.random.rand(size,size,size)*-.1       #velocity in x
        self.vy = np.random.rand(size,size,size)*-.1      #velocity in y
        self.vz = np.random.rand(size,size,size)*-.1       #velocity in z
        self.vx0 = np.random.rand(size,size,size)      #previous velocity in x
        self.vy0 = np.random.rand(size,size,size)      #previous velocity in y
        self.vz0 = np.random.rand(size,size,size)      #previous velocity in z
    def set_velocities(self, vx, vy, vz, x, y, z):
        self.vx0[x,y,z] = self.vx[x,y,z]
        self.vy0[x,y,z] = self.vy[x,y,z]
        self.vz0[x,y,z] = self.vz[x,y,z]
        self.vx[x,y,z] += vx
        self.vy[x,y,z] += vy
        self.vz[x,y,z] += vz
    def set_density(self, x, y, z, dens):
        self.dens[x,y,z] += dens
    
def step(cube):
    size = cube.size
    dt = cube.dt
    diffusion = cube.diff
    visc = cube.visc
    s = cube.s
    dens = cube.dens
    vx = cube.vx
    vy = cube.vy
    vz = cube.vz
    vx0 = cube.vx0
    vy0 = cube.vy0
    vz0 = cube.vz0
    # print(dens[:,:,4])
    diffuse(1, vx0, vx, visc, dt, 4, size)
    diffuse(2, vy0, vy, visc, dt, 4, size)
    diffuse(3, vz0, vz, visc, dt, 4, size)

    project(vx0, vy0, vz0, vx, vy, 4, size)

    advect(1, vx, vx0, vx0, vy0, vz0, dt, size)
    advect(2, vy, vy0, vx0, vy0, vz0, dt, size)
    advect(3, vz, vz0, vx0, vy0, vz0, dt, size)

    project(vx0, vy0, vz0, vx, vy, 4, size)

    diffuse(0, s, dens, diffusion, dt, 4, size)
    advect(0, dens, s, vx, vy, vz, dt, size)
    # print(dens[:,:,4])

    cube.dens = dens
    cube.s = s
    cube.vx = vx
    cube.vy = vy
    cube.vz = vz
    cube.vx0 = vx0
    cube.vy0 = vy0
    cube.vz0 = vz0
    

def set_bnd(b, x, N):
    for j in range(1, N-1):
        for i in range(1, N-1):
            if b == 3:
                x[i, j, 0] = -x[i, j, 1]
                x[i, j, N-1] = -x[i, j, N-2]
            else:
                x[i, j, 0] = x[i, j, 1]
                x[i, j, N-1] = x[i, j, N-2]
    for k in range(1, N-1):
        for i in range(1, N-1):
            if b == 2:
                x[i, 0  , k] = -x[i, 1  , k]
                x[i, N-1, k] = -x[i, N-2, k]
            else:
                x[i, 0  , k] = x[i, 1  , k]
                x[i, N-1, k] = x[i, N-2, k]
    for k in range(1, N-1):
        for i in range(1, N-1):
            if b ==1:
                x[0  , j, k] = -x[1  , j, k]
                x[N-1, j, k] = -x[N-2, j, k]
            else:
                x[0, j, k] = x[1, j, k]
                x[N-1, j, k] = x[N-2, j, k]
    x[0,0,0] = (x[1,0,0]+ x[0,1,0]+ x[0,0,1])/3
    x[0,N-1,0] = (x[1,N-1,0]+ x[0,N-2,0]+ x[0,N-1,1])/3
    x[0, 0, N-1] = (x[1,0,N-1]+ x[0,1,N-1]+ x[0,0,N-2])/3
    x[0, N-1, N-1] = (x[1,N-1,N-1]+ x[0,N-2,N-1]+ x[0,N-1,N-2])/3
    x[N-1, 0, 0] = (x[N-2,0,0]+ x[N-1,1,0]+ x[N-1,0,1])/3
    x[N-1,N-1,0] = (x[N-2,N-1,0]+ x[N-1,N-2,0]+ x[N-1,N-1,1])/3
    x[N-1,N-1,N-1] = (x[N-2,N-1,N-1]+ x[N-1,N-2,N-1]+ x[N-1,N-1,N-2])/3
    x[N-1,0,N-1] = (x[N-2,0,N-1]+ x[N-1,1,N-1]+ x[N-1,0,N-2])/3

def lin_solve(b, x, x0, a, c, iter, N):
    for k in range(iter):
        for m in range(1,N-1):
            for j in range(1,N-1):
                for i in range(1, N-1):
                    x[i,j,m] = (x0[i, j, m]+ a*(x[i+1, j, m]+x[i-1, j, m]+x[i, j-1, m]+ x[i, j-1, m]+x[i, j, m+1]+x[i, j, m-1]))/c 
        set_bnd(b, x, N)
    
    

def diffuse(b, x, x0, diff, dt, iter, N):
    a = dt * diff* (N-2) *(N-2)
    lin_solve(b, x, x0, a, 1 + 6 *a, iter, N)

def project(velocX, velocY, velocZ, p, div, iter, N):
    for k in range(N-1):
        for j in range(N-1):
            for i in range(N-1):
                div[i, j, k] = -0.5*(velocX[i+1, j, k]-velocX[i-1, j, k] + velocY[i, j+1, k] - velocY[i, j-1, k] + velocZ[i, j, k+1] - velocZ[i, j, k-1])/N
                p[i,j,k] = 0
    set_bnd(0, div, N)
    set_bnd(0, p, N)
    lin_solve(0, p, div, 1, 6, iter, N) #Do you have to return this for numpy arrays or are the local operations on them saved? let's see # no return needed
    for k in range(N-1):
        for j in range(N-1):
            for i in range(N-1):
                velocX[i, j, k] = -0.5*(p[i+1, j, k]-p[i-1,j,k])
                velocY[i, j, k] = -0.5*(p[i, j+1, k]-p[i,j+1,k])
                velocZ[i, j, k] = -0.5*(p[i, j, k+1]-p[i,j,k-1])
    set_bnd(1, velocX, N)
    set_bnd(2, velocY, N)
    set_bnd(3, velocZ, N)

def advect(b, d, d0, velocX, velocY, velocZ, dt, N):
    dtx = dt *(N-2)
    dty = dt *(N-2)
    dtz = dt *(N-2)

    for k in range(N-1):
        for j in range(N-1):
            for i in range(N-1):
                tmp1 = dtx * velocX[i, j, k]
                tmp2 = dty * velocY[i, j, k]
                tmp3 = dtz * velocZ[i, j, k]
                x = i - tmp1
                y = j - tmp2
                z = k - tmp3
                if x <0.5:
                    x = 0.5
                if x > (N + 0.5):
                    x = N + 0.5
                i0 = m.floor(x)
                i1 = i0 + 1
                if y <0.5:
                    y = 0.5
                if y > (N + 0.5):
                    y = N + 0.5
                j0 = m.floor(y)
                j1 = j0 + 1
                if z <0.5:
                    z = 0.5
                if z > (N + 0.5):
                    z = N + 0.5
                k0 = m.floor(z)
                k1 = k0 + 1

                s1 = x - i0
                s0 = 1 - s1
                t1 = y- j0
                t0 = 1 - t1
                u1 = z -k0
                u0 = 1 - u1
                # double check this
                d[i,j,k] = s0*(t0*(u0*d0[i0, j0, k0] + u1 * d0[i1, j1, k1]) + (t1* (u0 * d0 [i0, j1, k0]+ u1* d0[i0, j1, k1]))) +s1 * (t0 *(u0* d0[i1, j0, k0] + u1* d0[i1, j0, k1])+ t1 * (u0*d0[i1,j1,k] + u1*d0[i1,j1,k1]))

    set_bnd(b, d, N)

def visualize(density, zaxis,index):
    fig, ax = plt.subplots()
    im = ax.imshow(density[:,:,zaxis], cmap=cm.RdYlGn)
    plt.savefig("tryout_%g.png"%(index))


def main():
    size = 64 #number of blocks in one direction
    viscosity = 10 #kg/m*s (viscosity of cream)
    diffusion = 10**-4 #diffusion constant, m^2/s
    dt = .1 #seconds
    middle = size //2
    fluid = FluidCube(size, dt, diffusion, viscosity)
    fluid.set_density(middle, middle, middle, 0.0000001)
    for k in range(100):
        # start = time.perf_counter()
        step(fluid)
        # end = time.perf_counter()
        # print(end-start)
        visualize(fluid.dens, middle, k)
    
if __name__=='__main__':
    main()