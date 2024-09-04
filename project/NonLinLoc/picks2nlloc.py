#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 3 19:55:57 2022 

@original author:Tim Lin
@email:jiunting@uoregon.edu

Code adapted from get_arrivals.py to use on the stacked families of the full
catalog of Bostock for P and S picks.

@author: lpapin
"""

import glob
import json
import numpy as np
import pandas as pd
from obspy import UTCDateTime
from datetime import datetime, timedelta

#---------------parameter setting-----------------
N_min = 3 # at least N station detect arrival in the same 15s time window
threshold = 0.1 # decision threshold

use_P = True # add P arrival?
N_min_P = 3
threshold_P = 0.5

#excluded_sta = ['PO.KLNB', ] # excluded stations
excluded_sta = []

fileout = "./arrivals_sta%d_y%.1f.csv"%(N_min,threshold)
#fileout = "./arrivals_sta%d_y%.1f_PS.csv"%(N_min,threshold)
base_dir ='/Users/lpapin/Downloads/picker_full-2' ###

with open("templates_dict_bostock.json", 'r') as json_file:
    templates_dict_bostock = json.load(json_file)
families=list(templates_dict_bostock.keys())

#---------------parameter setting END-----------------

def dect_time(file_path: str, thresh=0.1, is_catalog: bool=False, EQinfo: str=None) -> np.ndarray:
    """
    Get detection time for 1). ML detection in .csv file or 2). catalog file.
    If input is .csv file, return the `window` starttime so that is easier to match with other stations
    
    INPUTs
    ======
    file_path: str
        Path of the detection file
    is_catalog: bool
        Input is catalog or not
    thresh: float
        Threshold of selection. Only apply when is_catalog==False (i.e. using ML detection csv files)
    EQinfo: str
        Path of the EQinfo file i.e. "sav_family_phases.npy"
    
    OUTPUTs
    =======
    sav_OT: np.array[datetime] or None if empty detection
        Detection time in array. For ML detection, time is the arrival at the station. For catalog, time is the origin time
    
    EXAMPLEs
    ========
    #T1 = dect_time('./Data/total_mag_detect_0000_cull_NEW.txt',None,True,'./Data/sav_family_phases.npy')  #for catalog
    #T2 = dect_time('./results/Detections_S_small/cut_daily_PO.GOWB.csv')                             # for ML detection output
    """
    if is_catalog:
        EQinfo = np.load(EQinfo,allow_pickle=True) #Note! detection file has a shift for each family
        EQinfo = EQinfo.item()
        sav_OT = []
        with open(file_path,'r') as IN1:
            for line in IN1.readlines():
                line = line.strip()
                ID = line.split()[0] #family ID
                OT = UTCDateTime('20'+line.split()[1]) #YYMMDD
                HH = int(line.split()[2])*3600
                SS = float(line.split()[3])
                OT = OT + HH + SS + EQinfo[ID]['catShift'] #detected origin time
                sav_OT.append(OT.datetime)
        sav_OT.sort()
        return pd.Series(sav_OT)
    else:
        csv = pd.read_csv(file_path)
        if len(csv)==0:
            return None, None
        csv=modify_cut_format(csv)
        T = csv[csv['y']>=thresh].starttime.values
        T_arr = np.array([i.split('_')[1] for i in csv[csv['y']>=thresh].id.values])
        fam = np.sort(csv[csv['y']>=thresh].lfe_family.values)
        return T, T_arr, fam

def modify_cut_format(csv): ###
    """
    Modify the format of a seismic detection CSV file to standardize time and create unique IDs.

    This function adjusts the time format in the CSV file, standardizes the start time, and modifies
    the `id` column by appending a timestamp based on the time column in the original data. The function 
    also retains only the necessary columns for further analysis.
    """
    print('Now modifying :',csv)
    # Identify the time column
    time_col = None
    for col in csv.columns:
        if col.endswith('_times'):
            time_col = col
            break
    if time_col is None:
        time_col = csv.columns[1]
    fixed_starttime = '2000-01-01T00:00:00.000000Z'
    csv['starttime'] = fixed_starttime
    base_time = datetime.strptime('2000-01-01T00:00:00.000000Z', '%Y-%m-%dT%H:%M:%S.%fZ')
    csv['id'] = csv.apply(lambda row: f"{row['id']}_{(base_time + timedelta(seconds=row[time_col])).strftime('%Y-%m-%dT%H:%M:%S.%f')[:-4]}", axis=1)
    csv = csv[['starttime', 'lfe_family', time_col, 'nb_events', 'y', 'id']]
    # csv=csv.sort_values(by='lfe_family')
    return csv

#STEP 1. =====Load all the detections from ML=====

csvs = glob.glob(base_dir+"/Detections_S_stack/*.csv")
csvs.sort()
# csvs=csvs[0:3]

sav_T = {}
sav_T_arr = {} # link starttime to arrival
sav_fam = {}
for csv in csvs:
    print('Now loading :',csv)
    net_sta = csv.split('/')[-1].split('_')[-1].replace('.csv','')
    if net_sta in excluded_sta:
        print(' Skip station:',net_sta)
        continue
    T,T_arr, fam = dect_time(csv, thresh=threshold)
    if T is None:
        continue
    # print(T,'\n',T_arr)
    sav_T[net_sta] = T
    sav_T_arr[net_sta] = T_arr#{i:j for i,j in zip(T,T_arr)}
    sav_fam[net_sta]=fam
    print(sav_T[net_sta].size,sav_T_arr[net_sta].size,sav_fam[net_sta].size)

#STEP 2. =====Count the number of detections in the same starttime window. For simplicity, this doesnt consider nearby windows i.e. +-15 s.
print('-----Initialize all_fam, number of detections for each family, (Min=1, Max=%d)'%(len(sav_fam)))
all_fam = {}
for station, families in sav_fam.items():
    print('Processing station:', station)
    for family in families:
        if family not in all_fam:
            all_fam[family] = {'num': 1, 'sta': [station]}
        else:
            all_fam[family]['num'] += 1
            all_fam[family]['sta'].append(station)

#STEP 3. =====Apply min N stations filter=====
#N_min = 3 #define in the begining
sav_k = [] # keys that pass the filter
for k in all_fam.keys():
    if all_fam[k]['num']>=N_min:
        sav_k.append(k)

sav_k.sort()
#sav_k = np.array([UTCDateTime(i) for i in sav_k])
print('Total of candidate templates=%d after N_min=%d filter'%(len(sav_k), N_min))

#STEP 4. ======Add P arrival?, repeat STEP1-3 with P wave======

if use_P:
    csvs = glob.glob(base_dir + "/Detections_P_stack/*.csv")
    csvs.sort()
    sav_T_P = {}
    sav_T_arr_P = {}  # link starttime to arrival
    sav_fam_P = {}

    # Load detections from P wave files
    for csv in csvs:
        print('Now loading :', csv)
        net_sta = csv.split('/')[-1].split('_')[-1].replace('.csv', '')
        if net_sta in excluded_sta:
            print(' Skip station:', net_sta)
            continue
        T, T_arr, fam = dect_time(csv, thresh=threshold_P)
        if T is None:
            continue
        sav_T_P[net_sta] = T
        sav_T_arr_P[net_sta] = T_arr#{i: j for i, j in zip(T, T_arr)}
        sav_fam_P[net_sta] = fam

    # Initialize counting for P wave detections by lfe_family
    print('-----Initialize all_fam_P, number of detections for each family, (Min=1, Max=%d)' % (len(sav_fam_P)))
    all_fam_P = {}
    for station, families in sav_fam_P.items():
        print('Processing station:', station)
        for family in families:
            if family not in all_fam_P:
                all_fam_P[family] = {'num': 1, 'sta': [station]}
            else:
                all_fam_P[family]['num'] += 1
                all_fam_P[family]['sta'].append(station)
'''
    # Apply min N stations filter for P wave detections
    sav_k_P = []  # keys that pass the filter
    for k in all_fam_P.keys():
        if all_fam_P[k]['num'] >= N_min_P:
            sav_k_P.append(k)

    sav_k_P.sort()
    print('Total of candidate P wave templates=%d after N_min=%d filter' % (len(sav_k_P), N_min_P))
'''

#STEP final. =====Find the corresponding arrival time and write the output=====
OUT1 = open(fileout,'w')
NID = 1

for fam in sav_k:  # Loop through each family that passed the filter
    results = {}  # Collect results to write later
    
    # For each station that detected this family, retrieve the corresponding detection time
    for sta in all_fam[fam]['sta']:  # Stations that detected this family
        try:
            # Find the detection time that corresponds to the current family
            Sarr = next(time for time, family in zip(sav_T_arr[sta], sav_fam[sta]) if family == fam)
            # Sarr = sav_T_arr[sta][detection_time]
            results[sta] = {'S': Sarr}
            print(f'S arrival for station {sta} and family {fam}: {Sarr}')
        except StopIteration:
            print(f"Family {fam} not found in station {sta}'s data.")
            continue
    # If P arrivals are being used, retrieve those as well
    if use_P:
        if fam in all_fam_P:
            if all_fam_P[fam]['num'] >= N_min_P:
                for sta in all_fam_P[fam]['sta']:
                    try:
                        # Find the detection time that corresponds to the current family
                        Parr = next(time for time, family in zip(sav_T_arr_P[sta], sav_fam_P[sta]) if family == fam)
                        # Parr = sav_T_arr_P[sta][detection_time]
                        if sta in results:
                            results[sta]['P'] = Parr  # This station has both S and P
                        else:
                            results[sta] = {'P': Parr}
                        print(f'P arrival for station {sta} and family {fam}: {Parr}')
                    except StopIteration:
                        print(f"P wave family {fam} not found in station {sta}'s data.")
                        continue

    # Write to NLLoc format
    OUT1.write(f'#ID:{fam:03}\n')
    NID += 1
    for sta, phases in results.items():
        for phs, arr in phases.items():
            arr_time = UTCDateTime(arr)
            OUT1.write(f'{sta.split(".")[1]:6s} ?  ?    ? {phs}      ? {arr_time.strftime("%Y%m%d %H%M %S.%f")[:19]} GAU  1.00e-01 -1.00e+00 -1.00e+00 -1.00e+00\n')

OUT1.close()
