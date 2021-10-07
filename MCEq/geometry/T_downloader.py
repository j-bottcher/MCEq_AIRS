
from __future__ import print_function
from __future__ import division

from datetime import datetime, timedelta, date
import sys
import os
import glob
import numpy as np

# Set working directory

def AIRS_Download(workpath,startDate,endDate):    
    remotedir0 = "https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRS3STD.006/" # Place of AIRS data
    localdir0 = workpath + startDate[0:4] 
    #https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRS3STD.006/2016/AIRS.2016.01.01.L3.RetStd_IR001.v6.0.31.0.G16004140111.hdf
    if not os.path.exists(localdir0):
        os.makedirs(localdir0)

    stdate = datetime.strptime(startDate, "%Y%m%d")
    endate = datetime.strptime(endDate, "%Y%m%d")

    print('Start date:',stdate)
    print('End date:',endate)
    #print('Downloading satellite data ...') 
    d = stdate
    Localfiles = []
    while d <= endate:

        adate = datetime.strftime(d, "%Y%m%d")
        print(adate)
        #d += timedelta(days=1)

        # Remote directory at NASA
        remotedir  = remotedir0 +  adate[0:4] 

        # Local directory (temporary)
        localdir = localdir0  
        localfile  = localdir  + "/" + "AIRS."+adate[0:4]+"."+adate[4:6]+"."+adate[6:8]+ ".*.hdf"
        Localfiles.append(localfile)
        
        if glob.glob(localfile) == []:
            os.chdir(localdir)
            print('Downloading satellite data ...')
            acclist = "AIRS."+adate[0:4]+"."+adate[4:6]+"."+adate[6:8]+ ".*.hdf"
            opts = " --load-cookies ~/.netrc --save-cookies ~/.netrc --keep-session-cookies "
            server = remotedir
            os.system("wget -e robots=off -nd -nc -np -r -l1 --no-parent --reject 'index.html*' -A " + acclist + opts + server) # Download
            print("wget -q -e robots=off -nd -nc -np -r -l1 --no-parent --reject 'index.html*' -A " + acclist + opts + server)
            printfile = glob.glob("AIRS."+adate[0:4]+"."+adate[4:6]+"."+adate[6:8]+ ".*.hdf")
            print(printfile, 'has been created.')
        d += timedelta(days=1)
    print('Satellite data downloaded.')
    return Localfiles

def AIRS_read_in(filename):
    from netCDF4 import Dataset
    hdffile = glob.glob(filename)
    nc = Dataset(hdffile[0])
    nc_var = nc.variables
    lat,long  = nc_var["Latitude"][:,0], nc_var["Longitude"][0,:]
    press_lvl = nc_var['StdPressureLev:ascending_MW_only'][:] #in hPa
    #--- read temperature profile and geopotential height
    #    ascending node
    tprof_all_a = np.array(nc_var[ "Temperature_A"]) # in K ,height,lat,long
    geoph_all_a = np.array(nc_var[ "GPHeight_A"]) #in m ,height,lat,long
    #    descending node
    tprof_all_d = np.array(nc_var[ "Temperature_D"]) #in K , height,lat,long
    geoph_all_d = np.array(nc_var[ "GPHeight_D"]) #in m , height,lat,long
    
    tprof_all_a[tprof_all_a==-9999.] = np.NaN
    geoph_all_a[geoph_all_a==-9999.] = np.NaN
    tprof_all_d[tprof_all_d==-9999.] = np.NaN
    geoph_all_d[geoph_all_d==-9999.] = np.NaN
    
    temp_asc   = tprof_all_a
    geoph_asc  = geoph_all_a
    temp_desc  = tprof_all_d
    geoph_desc = geoph_all_d
    nc.close()
    return press_lvl,lat*np.pi/180. ,long*np.pi/180. ,temp_asc,geoph_asc,temp_desc,geoph_desc
    

if __name__ == '__main__':
    AIRS_Download('/data/user/jbottcher/AIRS_Data/AIRS_Download/','20120125','20120125')
    
