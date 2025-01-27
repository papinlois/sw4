#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:31:57 2024

@author: lpapin
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, Stream, Trace
from pyrocko import moment_tensor as pMT

# ========== Notebook Problem_Setup.ipynb ==========

def calculate_moment_tensor(strike, dip, rake, magnitude):
    """
    Calculate and format the normalized moment tensor components.

    Parameters:
    - strike (float): Strike angle (degrees)
    - dip (float): Dip angle (degrees)
    - rake (float): Rake angle (degrees)
    - magnitude (float): Magnitude of the earthquake (Mw)

    Returns:
    - moment (list): Formatted list of moment tensor components as [mxx, myy, mzz, mxy, mxz, myz]
    """
    # Convert magnitude to scalar moment
    M0 = pMT.magnitude_to_moment(magnitude)

    # Create moment tensor
    MT = pMT.MomentTensor(strike=strike, dip=dip, rake=rake, scalar_moment=M0)

    # Extract moment tensor components
    M6 = [MT.mnn, MT.mee, MT.mdd, MT.mne, MT.mnd, MT.med]

    # Normalize moment tensor components
    M6_normalized = np.array(M6) / MT.scalar_moment()

    # Set values below 1e-16 to zero
    M6_normalized = np.where(np.abs(M6_normalized) < 1e-16, 0, M6_normalized)

    # Function to format values with + for positive numbers
    def format_value(val):
        if val > 0:
            return f"+{val:.3f}"
        if val == 0:
            return f" {val:.2f}"
        return f"{val:.2f}"

    # Format the components
    moment_tensor = [format_value(val) for val in M6_normalized]

    # Map to the desired variable names
    moment_labels = ['mxx', 'myy', 'mzz', 'mxy', 'mxz', 'myz']
    moment_tensor = dict(zip(moment_labels, moment_tensor))

    return moment_tensor, M0

# ========== Script plot_waveforms.py ==========

def find_all_stations(base_dir):
    """
    Identify all unique station names in the directory structure.
    Station names are assumed to be 3-4 characters before `_m` in filenames.

    Args:
        base_dir (str): Path to the base directory containing parent_dirs.

    Returns:
        list: A sorted list of unique station names.
    """
    station_set = set()  # Use a set to avoid duplicates

    # Walk through all subdirectories and files
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            # Match station names (3-4 characters before "_m")
            match = re.match(r"([A-Z]{3,4})_m", file)
            if match:
                station_set.add(match.group(1))  # Add the station name

    return sorted(station_set)  # Return sorted list of unique station names

def read_gfs(base_dir, stas, comps, starttime):
    """
    Dynamically reads Green's functions from all subdirectories in base_dir.

    Parameters:
    - base_dir (str): Base directory containing multiple parent directories (e.g., "fam099").
    - stas (list): List of station names.
    - comps (list): List of moment tensor components ['mxx', 'mxy', 'mxz', 'myy', 'myz', 'mzz'].
    - starttime (UTCDateTime): Start time to assign to the traces.

    Returns:
    - gfs_streams (dict): Nested dictionary of Green's function streams:
      Format: {station: {component: {'north': Trace, 'east': Trace, 'vertical': Trace}}}.
    - short (list): List of stations/components missing **across all parent_dirs**.
    """
    
    # Final green functions streams
    gfs_streams = {}
    
    # Tracks pairs found and missing
    found_files = set()
    missing_files = set()

    # Dynamically find all parent directories in base_dir
    parent_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

    # Go through all the folders that contain results
    for parent in parent_dirs:
        parent_path = os.path.join(base_dir, parent)

        # Find all folders that contain "results_" and any station number
        result_dirs = [
            d for d in os.listdir(parent_path)
            if os.path.isdir(os.path.join(parent_path, d)) and "results_" in d
        ]

        # Extract the component from the folder name
        for result_dir in result_dirs:
            comp = result_dir.split("_")[-1]

            # Skip irrelevant directories
            if comp not in comps:
                continue

            result_path = os.path.join(parent_path, result_dir)

            for sta in stas:
                # Initialize station and component in the dictionary
                if sta not in gfs_streams:
                    gfs_streams[sta] = {}
                if comp not in gfs_streams[sta]:
                    gfs_streams[sta][comp] = {}

                # Attempt to read the files for the current station/component
                try:
                    north_path = os.path.join(result_path, f"{sta}_{comp}.n")
                    temp_n = read(north_path)
                    temp_n[0].stats.station = sta
                    temp_n[0].stats.starttime = starttime
                    gfs_streams[sta][comp]['north'] = temp_n[0]

                    east_path = os.path.join(result_path, f"{sta}_{comp}.e")
                    temp_e = read(east_path)
                    temp_e[0].stats.station = sta
                    temp_e[0].stats.starttime = starttime
                    gfs_streams[sta][comp]['east'] = temp_e[0]

                    vertical_path = os.path.join(result_path, f"{sta}_{comp}.u")
                    temp_z = read(vertical_path)
                    temp_z[0].stats.station = sta
                    temp_z[0].stats.starttime = starttime
                    gfs_streams[sta][comp]['vertical'] = temp_z[0]

                    # Record this station/component as found
                    found_files.add((sta, comp))

                except FileNotFoundError:
                    continue

    # Identify files truly missing across all folders
    for sta in stas:
        for comp in comps:
            if (sta, comp) not in found_files:
                missing_files.add((sta, comp))

    return gfs_streams, sorted(missing_files)

def simul_waveforms(gfs_streams, stas):
    """
    Calculate the final displacement as a Stream for each station and direction 
    by summing Green's functions.
    
    Parameters:
    - gfs_streams (dict): Green's function streams from `read_gfs`
    - stas (list): List of stations

    Returns:
    - displacement_streams (dict): Displacement Stream for each station
      Format: {station: Stream (with Traces for N, E, Z)}
      
    NB: If you have normalized moment tensor components as inputs and a moment 
    magnitude M0, use this function. Otherwise, check the displacement expression 
    and add the necessary values.
    """
    displacement_streams = {}

    for sta in stas:
        # Initialize empty Traces for each direction and station
        north_data = None
        east_data = None
        vertical_data = None

        # Sum the contributions of each component
        for comp in ['mxx', 'mxy', 'mxz', 'myy', 'myz', 'mzz']:
            if north_data is None:
                # Initialize with zeros the same shape as the Green's function data
                north_data = np.zeros_like(gfs_streams[sta][comp]['north'].data)
                east_data = np.zeros_like(gfs_streams[sta][comp]['east'].data)
                vertical_data = np.zeros_like(gfs_streams[sta][comp]['vertical'].data)

            # Add contributions to each direction
            north_data += gfs_streams[sta][comp]['north'].data
            east_data += gfs_streams[sta][comp]['east'].data
            vertical_data += gfs_streams[sta][comp]['vertical'].data

        # Create new traces for the final displacement
        north_tr = Trace(data=north_data, header=gfs_streams[sta]['mxx']['north'].stats)
        north_tr.stats.channel = 'N'
        east_tr = Trace(data=east_data, header=gfs_streams[sta]['mxx']['east'].stats)
        east_tr.stats.channel = 'E'
        vertical_tr = Trace(data=vertical_data, header=gfs_streams[sta]['mxx']['vertical'].stats)
        vertical_tr.stats.channel = 'Z'

        # Create a Stream for the station
        displacement_stream = Stream(traces=[north_tr, east_tr, vertical_tr])
        displacement_streams[sta] = displacement_stream

    return displacement_streams

def plot_displacement_streams_indiv(displacement_streams, stas):
    """
    Custom plot for displacement streams.

    Parameters:
    - displacement_streams (dict): Displacement Stream for each station
      Format: {station: Stream (with Traces for N, E, Z)}
    - stas (list): List of stations
    """
    for sta in stas:
        if sta not in displacement_streams:
            print(f"Station {sta} not found in displacement streams.")
            continue

        stream = displacement_streams[sta]

        # Create a figure for the station
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(9, 8), sharex=True)
        fig.suptitle(f"Displacements for Station {sta}", fontsize=14, y=0.96)

        # Plot each component
        components = ['N', 'E', 'Z']
        colors = ['b', 'g', 'r']
        for i, (comp, color) in enumerate(zip(components, colors)):
            trace = next((tr for tr in stream if tr.stats.channel == comp), None)
            if trace:
                axes[i].plot(trace.times(), trace.data, color=color, label=f"{comp} Component")
                axes[i].set_ylabel("Amplitude")
                axes[i].legend()
            else:
                axes[i].text(0.5, 0.5, f"{comp} Component Missing", ha='center', va='center')
            axes[i].grid(True)

        axes[-1].set_xlabel("Time (s)")
        plt.tight_layout()
        plt.show()

def plot_displacement_streams(displacement_streams, stas):
    """
    Custom plot for displacement streams with ?3 columns for N, E, Z and a row for each station.
    Divides the streams into multiple figures, each containing up to ? stations.

    Parameters:
    - displacement_streams (dict): Displacement Stream for each station
      Format: {station: Stream (with Traces for N, E, Z)}
    - stas (list): List of stations
    """
    components = ['N', 'E', 'Z']
    colors = ['b', 'g', 'r']
    
    # Split the stations into chunks
    chunk_size = 3
    stas_chunks = [stas[i:i + chunk_size] for i in range(0, len(stas), chunk_size)]

    # Loop over each chunk of stations and plot
    for stas_chunk in stas_chunks:
        n_stas = len(stas_chunk)

        # Create a new figure for each chunk
        fig, axes = plt.subplots(nrows=n_stas, ncols=3, figsize=(10, 2 * n_stas), sharex=True, sharey=True)
        fig.suptitle("Displacement Components by Station", fontsize=14)

        # Ensure axes is always a 2D array, even for a single station
        if n_stas == 1:
            axes = [axes]

        for row_idx, sta in enumerate(stas_chunk):
            print(f"Plotting {sta}")  # Debugging line
            if sta not in displacement_streams:
                print(f"Station {sta} not found in displacement streams.")
                continue

            stream = displacement_streams[sta]

            for col_idx, (comp, color) in enumerate(zip(components, colors)):
                ax = axes[row_idx][col_idx] if n_stas > 1 else axes[col_idx]
                trace = next((tr for tr in stream if tr.stats.channel == comp), None)

                if trace:
                    ax.plot(trace.times(), trace.data, color=color, label=f"{sta} - {comp}")
                else:
                    ax.text(0.5, 0.5, f"{comp} Component Missing", transform=ax.transAxes,
                            ha='center', va='center', fontsize=10, color='gray')

                ax.set_title(f"{comp} Component" if row_idx == 0 else "", fontsize=12)
                ax.set_ylabel(f"Station: {sta}" if col_idx == 0 else "")
                ax.grid(True)

        # Set x-axis label only for the bottom row
        if n_stas > 1:
            axes[-1][1].set_xlabel("Time (s)", fontsize=14)
        else:
            axes[0][1].set_xlabel("Time (s)", fontsize=14)

        # Adjust layout to make space for the title
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()

        # Save each figure with a filename that indicates the range of stations
        plt.savefig(f"{stas_chunk[0]}_to_{stas_chunk[-1]}.png", dpi=450)

