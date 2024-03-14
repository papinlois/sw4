import xarray as xr
import pandas as pd
import numpy as np
from scipy.interpolate import griddata

def write_topo(file2load, name2write, lonMin, lonMax, latMin, latMax):
    """ Takes file path for netcdf file which should end in '.nc', and writes to 'file2write' 
    which should end in '.topo' for the region specified by lonMin,lonMax, latMin, latMax"""
    suffix=".topo"
    file2write=name2write+suffix
    # read topo from netcdf
    grid_big = xr.open_dataset(file2load)
    grid=grid_big.loc[dict(lat=slice(latMin, latMax),lon=slice(lonMin, lonMax))]
    #convert to pandas DF
    # sort by lat and lon
    df_grid = grid.to_dataframe().reset_index().sort_values(['lat', 'lon'])
                                           
    # add nb lon nb lat as first line
    custom_first_line = f"{len(grid.lon.data)} {len(grid.lat.data)}"


    # Open the file in write mode and add the custom first line
    with open(file2write, 'w') as file:
        file.write(custom_first_line + '\n')

    # Append the DataFrame to the CSV file
    df_grid.to_csv(file2write, mode='a', index=False, header=False, float_format='%.5f', sep=' ')

def vp2rho(pvel):
    """ This function takes a p velocity in km/s and returns a density in
    g/cc using Brocher 2005 eqn 1 """
    # check if your velocities are in bounds of the parameterization
    if np.min(pvel)<1.5:
        print("Proceed with caution lowest Vp is outside of parameterization")
    if np.max(pvel)>8.5:
        print("Proceed with caution highest Vp is outside of parameterization")

    # use equation
    rho=1.6612*pvel-0.4721*pvel**2+0.0671*pvel**3-0.0043*pvel**4+0.000106*pvel**5
    # #convert from g/cc to kg/m**3
    # rho=rho*1000
    rho=np.round(rho,2) #round to 2 decimal places
    return rho

def generate_pfile(vmod, name, lonmin, lonmax, latmin, latmax):
    """Input includes a pandas dataframe of velocity model data and the boundaries that the files will be interpolated to, velocities should have units of km/s and depths should have units of km, name should be the prefix for the file without the '.ppmod' suffix"""
    suffix=".ppmod"
    filename=name+suffix
    delta=0.25

    lats=np.round(np.arange(latmin,latmax,0.025),decimals=3)
    lons=np.round(np.arange(lonmin,lonmax,0.025),decimals=3)
    depths= vmod['depth'].unique()
    grid=np.zeros((len(lats)*len(lons)*len(depths),3))
    ii = 0
    for lon in lons:
        for lat in lats:
            for dep in depths:
                grid[ii,0], grid[ii,1], grid[ii,2]=lon, lat, dep
                ii+=1
        
    grid_vp = griddata((vmod['lon'],vmod['lat'],vmod['depth']), vmod['vp'], 
                    (grid), method='linear')
    grid_vs = griddata((vmod['lon'],vmod['lat'],vmod['depth']), vmod['vs'], 
                    (grid), method='linear')
    vmod_interp = pd.DataFrame(np.hstack((grid,
                                      grid_vp.reshape(len(lats)*len(lons)*len(depths),1),
                                      grid_vs.reshape(len(lats)*len(lons)*len(depths),1))), 
                           columns=['lon', 'lat', 'depth', 'vp', 'vs'])


    # Check for NAN drama
    assert np.isnan(vmod_interp['vp']).sum() == 0 and np.isnan(vmod_interp['vs']).sum() == 0

    vmod_interp['rho']=vp2rho(vmod_interp['vp'])

    ## open the textfile from before in append mode
    file=open(filename,'w')

    ## Write Header
    # f'whatever str {this is a var:.3f}'
    file.write(name+'\n'+
        str(delta)+'\n'+
        str(len(lats))+' '+f"{np.min(lats):.3f}"+' '+f"{np.max(lats):.3f}"+'\n'+
        str(len(lons))+' '+f"{np.min(lons):.3f}"+' '+f"{np.max(lons):.3f}"+'\n'+
        str(len(depths))+' '+f"{np.min(depths):.3f}"+' '+f"{np.max(depths):.3f}"+'\n'+
        '-99 -99 -99 -99\n'+
        '.FALSE.\n')

    file.close()

    ## Write depth profiles
    # note the sorting looks weird because comparing floats in weird sometimes

    Ndepth=len(depths)
    for lat in lats:
        for lon in lons:
            # print(lat, lon, Ndepth)
            tmp=vmod_interp[(vmod_interp['lat']==lat) & (vmod_interp['lon']==lon)]
            tmp.insert(0,'ind',range(1,len(depths)+1))
            #  print(tmp[['ind','depth','vp','vs','rho']])
            # depthprof=lat_subset.loc[np.abs(lat_subset['lon']-lon)<1e-2,['depth', 'vp','vs','rho']]
            file=open(filename,'a')
            file.write(str(lat)+' '+str(lon)+' '+str(Ndepth)+'\n')
            file.close()
            tmp[['ind','depth','vp','vs','rho']].to_csv(filename, sep=' ', header=False, index=False, mode='a',float_format='%.3f')

