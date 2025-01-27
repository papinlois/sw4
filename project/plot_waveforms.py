#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:06:03 2024

@author: lpapin
"""

import sw4_tools
from obspy.core.utcdatetime import UTCDateTime

## Case you want to plot
# Family of LFEs
base_dir = 'fam099/'
# Stations
# stas = ['GLBC','LZB'] # Choosen stations
stas = sw4_tools.find_all_stations(base_dir) # All stations
print(stas)
# Moment tensor components
comps = ['mxx', 'mxy', 'mxz', 'myy', 'myz', 'mzz']
# Initial time (to plot the 15s window of simulated wveform)
starttime = UTCDateTime("2000-01-01")

# Read Green's functions as streams
gfs_streams, short = sw4_tools.read_gfs(base_dir, stas, comps, starttime)
# print(gfs_streams["GLBC"]["mxx"]["north"]) #test
if short:
    print("Warning: The following Green's functions are missing or incomplete:")
    print("\n".join(short))

# Calculate final displacement streams
displacement_streams = sw4_tools.simul_waveforms(gfs_streams, stas)

# Plot the waveforms
sw4_tools.plot_displacement_streams_indiv(base_dir, displacement_streams, stas)

# =============================================================================
# # Print and plot the results
# for sta in stas:
#     print(f"Displacement Stream for Station {sta}:")
#     print(displacement_streams[sta])
#     displacement_streams[sta].plot()
# =============================================================================

##

sw4_tools.plot_displacement_streams(base_dir, displacement_streams, stas)
