import numpy

size = 0.01 #meter
viscosity = 10 #kg/m*s
dt = 0.001 #seconds





class FluidCube:
    def __init__(self, size, dt, diff, visc, s, density, vx, vy, vz, vx0,vy0, vz0):
        self._size = size               #length in one dimension
        self._dt = dt                   #timestep
        self._diff = diff               #diffusion constant
        self._visc = visc               #viscosity
        self._s = s                     #unsure
        self._density = density         #density
        self._vx = vx                   #velocity in x
        self._vy = vy                   #velocity in y
        self._vz = vz                   #velocity in z
        self._vx0 = vx0                 #previous velocity in x
        self._vy0 = vy0                 #previous velocity in y
        self._vz0 = vz0                 #previous velocity in z

FluidCube = FluidCube(s)