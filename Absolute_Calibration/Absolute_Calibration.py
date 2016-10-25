#Imports
from __future__ import print_function
from math import *
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

#Matplotlib settings
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

sns.set(font_scale=2.0)
sns.set_style("darkgrid")
sns.set_palette(palette='deep')
sns.set_color_codes(palette='deep')

#%%
#==============================================================================
# Functions for determining the background and true signal
#==============================================================================
def Background(data,interval):
    """Calculates the background of a data set given an interval over which there is no true signal"""
    return np.mean(np.mean(data[:,interval[0]:interval[1]],axis=0))

def Signal(data,interval):
    """Calculates the true signal by subtracting the background"""
    return data - Background(data,interval)

def AverageSignal(data,interval):
    """Calculates the average signal"""
    return np.mean(Signal(data,interval),axis=0)


#%%
#==============================================================================
# Load in the data
#==============================================================================
filepath1MeV = './ND_Data/0.541758/traces.csv'
data1MeV = np.genfromtxt(filepath1MeV, delimiter=',')

filepath2MeV = './ND_Data/0.039960/traces.csv'
data2MeV = np.genfromtxt(filepath2MeV, delimiter=',')

filepath2p5MeV = './ND_Data/0.597246/traces.csv'
data2p5MeV = np.genfromtxt(filepath2p5MeV, delimiter=',')


#%%
#==============================================================================
# Determine signal from data
#==============================================================================
interval = (0,500) #there appears to be no signal over this interval
Signal1MeV = Signal(data1MeV,interval)
AverageSignal1MeV = AverageSignal(data1MeV,interval)

Signal2MeV = Signal(data2MeV,interval)
AverageSignal2MeV = AverageSignal(data2MeV,interval)

Signal2p5MeV = Signal(data2p5MeV,interval)
AverageSignal2p5MeV = AverageSignal(data2p5MeV,interval)


#%%
#==============================================================================
# Mesh plot of signal
#==============================================================================
fig1 = plt.figure(figsize=(12,8))
ax1 = fig1.add_subplot(111)
ax1.set_xlim(0,3648)
ax1.set_xlabel('Pixel')
ax1.set_ylabel('Instance')


mesh1 = ax1.pcolormesh(Signal1MeV, cmap='viridis',vmin=0, vmax=np.max(Signal1MeV))
colorbar = fig1.colorbar(mesh1)
colorbar.set_label('Pixel Value')
plt.show()


#%%
#==============================================================================
# Plot of average signals at each energy
#==============================================================================
plt.figure()
plt.plot(AverageSignal1MeV,label='1.0 MeV')
plt.plot(AverageSignal2MeV,label='2.0 MeV')
plt.plot(AverageSignal2p5MeV,label='2.5 MeV')
plt.xlim(0,3648)
plt.ylim(0)
plt.xlabel('Pixel')
plt.ylabel('Pixel Value')
plt.legend()
plt.show()

#%%
#==============================================================================
# Integrate over the signals at each energy
#==============================================================================
Integral1MeV = np.sum(AverageSignal1MeV[1200:2000])
Integral2MeV = np.sum(AverageSignal2MeV[1000:2000])
Integral2p5MeV = np.sum(AverageSignal2p5MeV[1000:2000])
print('1 MeV integral =',Integral1MeV)
print('2 MeV integral =',Integral2MeV)
print('2.5 MeV integral =',Integral2p5MeV)
