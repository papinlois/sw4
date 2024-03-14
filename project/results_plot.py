#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 18:35:58 2024

@author: sydneydybing
"""
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('vci-results/ALBH4.txt')

plt.plot(data[:,0],data[:,1], label = 'E-W')
plt.plot(data[:,0],data[:,2], label = 'N-S')
plt.plot(data[:,0],data[:,3], label = 'Vertical')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.legend()