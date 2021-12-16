import numpy as np





class FluidCube:
    def __init__(self, size, dt, diff, visc):
        self._size = size                           #length in one dimension
        self._dt = dt                               #timestep
        self._diff = diff                           #diffusion constant
        self._visc = visc                           #viscosity
        self._s = np.zeros((size,size,size))        #unsure
        self._dens = np.zeros((size,size,size))     #density
        self._vx = np.zeros((size,size,size))       #velocity in x
        self._vy = np.zeros((size,size,size))       #velocity in y
        self._vz = np.zeros((size,size,size))       #velocity in z
        self._vx0 = np.zeros((size,size,size))      #previous velocity in x
        self._vy0 = np.zeros((size,size,size))      #previous velocity in y
        self._vz0 = np.zeros((size,size,size))      #previous velocity in z
    def set_velocities(self, vx, vy, vz, x, y, z):
        self._vx0[x,y,z] = self._vx[x,y,z]
        self._vy0[x,y,z] = self._vy[x,y,z]
        self._vz0[x,y,z] = self._vz[x,y,z]
        self._vx[x,y,z] += vx
        self._vy[x,y,z] += vy
        self._vz[x,y,z] += vz
    def set_density(self, x, y, z, dens):
        self._dens+= dens[x,y,z]
    
def setting_steps(cube):
    size = cube._size
    dt = cube._dt
    diffusion = cube._diff
    visc = cube._visc
    s = cube._s
    dens = cube._dens
    vx = cube._vx
    vy = cube._vy
    vz = cube._vz
    vx0 = cube._vx0
    vy0 = cube._vy0
    vz0 = cube._vz0

    vx = diffuse(1, vx0, vx, visc, dt, 4, size)
    vy = diffuse(2, vy0, vy, visc, dt, 4, size)
    vz = diffuse(3, vz0, vz, visc, dt, 4, size)

    vx, vy, vz = project(vx0, vy0, vz0, vx, vy, 4, size)

    vx = advect(1, vx, vx0, vx0, vy0, vz0, dt, size)
    vy = advect(1, vy, vy0, vx0, vy0, vz0, dt, size)
    vz = advect(1, vz, vz0, vx0, vy0, vz0, dt, size)

    vx, vy, vz = project(vx0, vy0, vz0, vx, vy, 4, size)

    diffuse(0, s, dens, diffusion, dt, 4, size)
    advect(0, dens, s, vx, vy, vz, dt, size)

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
    return x

def lin_solve(b, x, x0, a, c, iter, N):
    for k in range(iter):
        for m in range(1,N-1):
            for j in range(1,N-1):
                for i in range(1, N-1):
                    x[i,j,m] = (x0[i, j, m]+ a*(x[i+1, j, m]+x[i-1, j, m]+x[i, j-1, m]+ x[i, j-1, m]+x[i, j, m+1]+x[i, j, m-1]))/c 
        x = set_bnd(b, x, N)
    return x
    

def diffuse(b, x, x0, diff, dt, iter, N):
    a = dt * diff* (N-2) *(N-2)
    return lin_solve(b, x, x0, a, 1 + 6 *a, iter, N)


def main():
    size = 512 #number of blocks in one direction
    viscosity = 10 #kg/m*s (viscosity of cream)
    diffusion = 10**-12 #diffusion constant, m^2/s
    dt = 0.001 #seconds