
from datetime import datetime, timedelta, date
import sys
sys.path.insert(0, "/home/jbottcher")
import geopotcalc.gpcalc as gpc
import xarray as xa
import pandas as pd
import os
os.environ["ECCODES_DIR"] = "/home/jbottcher/.local/"
import cfgrib

import glob
import numpy as np
import cdsapi
import pygrib


def Download_ERA5(Outdir, datestart, datestop,time,model_lvl=False):
    if not os.path.exists(Outdir+datestart[0:4]):
        os.makedirs(Outdir+datestart[0:4])
    
    stdate = datetime.strptime(datestart, "%Y%m%d")
    endate = datetime.strptime(datestop, "%Y%m%d")
    #Year,Month, Day = datestart
    cds = cdsapi.Client()
    
    d = stdate
    Localfiles = []
    while d <= endate:
        adate = datetime.strftime(d, "%Y%m%d")
        d     += timedelta(days=1)
        if not model_lvl:
            Filename = Outdir+datestart[0:4]+'/ERA5_'+adate+'_'+time+'.grib'
        else:
            Filename = Outdir+datestart[0:4]+'/ERA5_Model_Lvl_single_day'+adate+'_'+time+'.grib2'
        Year  = adate[0:4]
        Month = adate[4:6]
        Day   = adate[6:8]
        print(Filename)
        if glob.glob(Filename) == []:
            if not model_lvl:
                cds.retrieve(
                    'reanalysis-era5-pressure-levels',
                    {
                        'product_type': 'reanalysis',
                        'format': 'grib',
                        'variable': [
                            'geopotential', 'temperature',
                        ],
                        'pressure_level': [
                            '1', '2', '3',
                            '5', '7', '10',
                            '20', '30', '50',
                            '70', '100', '125',
                            '150', '175', '200',
                            '225', '250', '300',
                            '350', '400', '450',
                            '500', '550', '600',
                            '650', '700', '750',
                            '775', '800', '825',
                            '850', '875', '900',
                            '925', '950', '975',
                            '1000',
                        ],
                        'year': Year,
                        'month': Month,
                        'day': Day,
                        'time': time,
                    },
                    Filename )
            else:
                if glob.glob(Outdir+'ERA5_Model_Lvl_'+Year+Month+'01'+'.grib2') == []:
                    print("Downloading:")
                    cds.retrieve('reanalysis-era5-complete', {
                                    'class': 'ea',
                                    'date': datetime.strftime(d, "%Y-%m-%d"),
                                    'expver': '1',
                                    'levelist': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
                                    'levtype': 'ml',
                                    'param': '129/130/152/133',
                                    'stream': 'oper',
                                    'grid'    : '1.0/1.0',  
                                    'time': '12:00:00',
                                    'type': 'an',
                                }, Filename)
                else:
                        Filename = Outdir+'ERA5_Model_Lvl_'+Year+Month+'01'+'.grib2'
        Localfiles.append(Filename)
    return Localfiles
def read_grib_data(Filename):
    grbs = pygrib.open(Filename)
    pressure_lvls = np.array([
            '1', '2', '3',
            '5', '7', '10',
            '20', '30', '50',
            '70', '100', '125',
            '150', '175', '200',
            '225', '250', '300',
            '350', '400', '450',
            '500', '550', '600',
            '650', '700', '750',
            '775', '800', '825',
            '850', '875', '900',
            '925', '950', '975',
            '1000',
        ]).astype(np.float)
    grbs.seek(0)
    T,H,lat,long = [],[],[],[]
    for i,grb in enumerate(grbs):
        if 'Temperature' in str(grbs.message(i+1)):
            T.append(grb.values) #K
        elif 'Geopotential' in str(grbs.message(i+1)):
            H.append(grb.values) #m^2/s^2
        l = grb.latlons()
        lat = l[0]
        long = l[1]


    T,H  = np.array(T),np.array(H)/9.80665 #K, #m
    return T,H,lat,long,pressure_lvls
def read_gribTwo_Data(Filename, date=None):
    
    
       
    era5_ds = cfgrib.open_datasets(Filename)
    if date is not None:
        for j,x in enumerate(era5_ds):
            era5_ds[j] = x.loc[dict(time=datetime.strftime(date,"%Y-%m-%d"))].squeeze('time')
            
    mldf    = pd.read_csv('~/akbk.csv',index_col='n')

    gpc.set_data(era5_ds, mldf, 137)
    Z = gpc.get_phi().values/9.80665  # calculates the geopential height at all lvls
    P = (1/2.*(gpc.get_p()[:-1]+gpc.get_p()[1:])/1e2)
    T = era5_ds[1]['t'].values
    lat,long = era5_ds[0]['latitude'].values,era5_ds[0]['longitude'].values
    
    return T,Z,lat,long,P
            

    
    
