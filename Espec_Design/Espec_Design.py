#Imports
from __future__ import print_function #In case you're using python 2
from math import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#Import custom modules
from physics import *

#%%
#==============================================================================
# Here we define the analytical expressions. These have been derived assuming a 
# uniform magnetic field and relativistic electrons. The radius of the electron
# trajectory in this case is given by

# $$ r = \frac{m_e c}{e B_0} \sqrt{ \left( \frac{E}{m_e c^2} + 1 \right)^2 - 1 }$$

# where $B_0$ is the magnetic field strength and $E$ is the electron energy. In
# the simplest case, the detector is placed at the edge of the magnet. The 
# final electron position along the detector as a function of energy is then 
# given by

# $$ z = \sqrt{ r^2 - (y_f-r)^2)} $$

# where $y_f$ is the distance between the electron injection position (the 
# slit) and the screen (the final electron y-position).
#==============================================================================

#Defined analytic functions
def KEcalc(Z,A,zf,yf):
    """Returns KE in J given final z-position in m"""
    Q = Z*q #Total charge
    M = A*mp #Total mass (technically this should use the atomic mass unit instead of the proton mass)
    return M*c**2*(np.sqrt((Q*B0/(M*c))**2*((zf**2+yf**2)/(2*yf))**2+1)-1)

def Radius(Z,A,KE):
    """Returns radius of electron orbit in m given KE in keV"""
    Q = Z*q #Total charge
    M = A*mp #Total mass (technically this should use the atomic mass unit instead of the proton mass)
    return M*c/(Q*B0)*np.sqrt((KE*1000*Q/(M*c**2)+1)**2-1)

def zfcalc(Z,A,KE,yf):
    """Returns the z-position at the screen in inches given the KE in keV and the distance between the 
    slit and the screen (yf) in m"""
    R = Radius(Z,A,KE)
    return np.sqrt(R**2 - (yf-R)**2)
    
    
#%%
#==============================================================================
# Input constants
#==============================================================================

#Espec parameters
B0 = 2136.0 #magnetic field strength in Gauss
yM = 0.5 #Distance between slit and magnet edge (screen location) in inches. Should be no greater than the magnet width.
maglength = 3.0 #magnet length in inches
magwidth = 2.0 #magnet width in inches

#Ion species parameters
Z = 1 #atomic number
A = me/mp #atomic mass

#Conversion to base units
B0 = B0*10**-4
yM = yM*0.0254
maglength = maglength*0.0254
magwidth = magwidth*0.0254


#%%
#==============================================================================
# Calculate final screen position as a function of energy
#==============================================================================

NumPoints = 1000
zposition = np.linspace(0,maglength,NumPoints)
Energy = KEcalc(Z,A,zposition,yM)

#Convert back to standard units
zposition = zposition*10**3 #Decided to use mm
Energy = Energy/(1000*q)

print('At z=',zposition[0],'mm, KE =',Energy[0],'keV')
print('At z=',zposition[-1],'mm, KE =',Energy[-1],'keV')

#%%
#==============================================================================
# Plot the results
#==============================================================================

mpl.rcParams.update({'font.size': 18, 'font.family': 'serif'})
plt.figure(figsize=(10,7))
plt.plot(zposition,Energy)
plt.xlabel('Position (mm)')
plt.ylabel('Energy (keV)')
plt.tight_layout()
plt.show()
