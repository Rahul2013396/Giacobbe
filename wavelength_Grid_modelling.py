import numpy as np

vellow = -225
velhigh = 225
velgap  = 2.7


radial_vel = np.arange(vellow,velhigh,velgap)


c= 3*10e5
def grid_maker(wavelength):
    grid = []
    for i in range(len(radial_vel)):
        shiftfactor = ((1+(radial_vel[i]/c))/(1-(radial_vel[i]/c)))**0.5
        shifted_wavelength = np.array(wavelength*shiftfactor)
        grid.append(shifted_wavelength)
    return grid
