#!/usr/bin/env python
# -*- coding: utf-8 -*-
# shauser modified version 3

# Imports ---------------------------------------------------------------------

import numpy as np
from scipy.interpolate import RectBivariateSpline
import os
import numpy
import numpy.ma as ma
import glob
import json
from datetime import datetime, timedelta, date
from pyhdf import HDF, SD, VS
import MCEq.geometry.myfunc_py3 as myfunc
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

# Comments --------------------------------------------------------------------

###############################################################################
###############################################################################
###########                     Important Note!                     ###########
###########      since the integration dX has to be conducted       ###########                                                ###########
###########             starting in the outer atmossphere,          ###########
###########            the datasets should be of a format           ###########
########### with height[0], temp[0] being the uppermost data points.###########
###########            Here, data should be an array of form        ###########
###########  Latitudes x Pressure levels x Longitudes (x Energies)  ###########
###############################################################################
###############################################################################

'''
AIRS-INFO
orbital period: 98.8 min
repeat cycle period: 233 orbits =^= 15.986389 days 
ascending orbit: equatorial crossing at 1:30 pm
descending orbit: 1:30 am
'''

# 3d USSA for testing ---------------------------------------------------------

temperature = np.flipud(np.array([287.43,283.20,278.68,268.57,260.81,251.92,241.44,228.58,220.79,216.65,216.65,216.65,216.65,217.23,220.50,223.13,225.02,227.70,232.71,239.22,249.45,257.88,264.03,270.65]))
altitude = np.flipud(np.array([110.9 ,  762.0 ,  1457.3 , 3012.2 , 4206.4 , 5574.4 , 7185.4 , 9164.0 , 10362.9 ,11784.0 ,13608.4 ,16179.7 ,18441.6, 20576.2 ,23848.6 ,26481.2 ,28368.1 ,31054.6 ,33452.6 ,35776.5 ,39429.54 ,   42439.8, 44637.2, 47820.0]))
pressure = np.flipud(np.array([1000.0 , 925.0 ,  850.0 ,  700.0 ,  600.0 ,  500.0  , 400.0 ,  300.0  , 250.0  , 200.0 ,  150.0   ,100.0 ,  70.0 ,   50.0  ,  30.0   , 20.0  ,  15.0  ,  10.0  ,  7.0   ,  5.0  ,   3.0   ,2.0,1.5,1.0]))
Latitude = np.linspace(-np.pi/2,np.pi/2,180)
Longitude = np.linspace(0,2*np.pi,3)
zenith = (Latitude + 270*np.pi/180. )/2.

T = np.zeros((len(zenith),len(altitude),len(Longitude)))
P = np.zeros((len(zenith),len(altitude),len(Longitude)))
H = np.zeros((len(zenith),len(altitude),len(Longitude)))
for i in range(len(zenith)):
    for j in range(len(altitude)):
        for k in range(len(Longitude)):
            T[i,j,k] = temperature[j]
            P[i][j][k] = pressure[j]
            H[i][j][k] = altitude[j]
H = np.array(H)
P = np.array(P)
T = np.array(T)

# Dictionary for further use --------------------------------------------------

USSA = {
"temperature" : T,
"altitude" : H,
"pressure" : pressure,
"latitude" : Latitude,
"longitude": Longitude,
"zenith"   : zenith
}

# Now following are the core calculations for Teff ----------------------------

# This function is used for the readin at the moment --------------------------

def read_in(workpath,startDate="yearmonthday",endDate="yearmonthday"):

	# Make the data easy iterable
    stdate = datetime.strptime(startDate, "%Y%m%d")
    endate = datetime.strptime(endDate, "%Y%m%d")

    print('Current date:',stdate)
    print('Loading satellite data ...')    
    
	# Loop over all dates
    d = stdate
    while d <= endate:
            
        adate = datetime.strftime(d, "%Y%m%d")
        d += timedelta(days=1)
         
        # Set data path
        remotedir0 = "https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRS3STD.006" # Place of AIRS data
        remotedir  = remotedir0 + "/" + adate[0:4]

        # Set working path
        localdir0 = workpath + startDate[0:4] # Directory of downloaded HDF files
        if not os.path.exists(localdir0):
            os.makedirs(localdir0)
        localdir = localdir0
        localfile  = localdir  + "/" + "AIRS."+adate[0:4]+"."+adate[4:6]+"."+adate[6:8]+ ".*.hdf"
        
        if glob.glob(localfile) == []:
            os.chdir(localdir)

            acclist = "AIRS."+adate[0:4]+"."+adate[4:6]+"."+adate[6:8]+ ".*.hdf"
            opts = "--save-cookies ~/.urs_cookies --keep-session-cookies "
            server = remotedir
            os.system("wget -q -nd -nc -np -r -l1 --no-parent -A " + acclist + opts + server) # Download
            
        if glob.glob(localfile) == []: 
            continue
        #--- output file
        hdffile = glob.glob(localfile)[0]
                
                #--- open HDF file for vdata access
        hdf = HDF.HDF(hdffile)
        vs  = hdf.vstart() #initialize VS interface
        #--- read vdata - pressure level, start time and end time
        grid_plvl = myfunc.read_vdata(vs, "StdPressureLev")[0][0]
         #0 - 1000hPa... 23 - 1hPa
        #    ascending node
        sttime_a = myfunc.read_vdata(vs, "AscendingGridStartTimeUTC")[0][0]
        entime_a = myfunc.read_vdata(vs, "AscendingGridEndTimeUTC")[0][0]
        sttime_a = sttime_a[0:10]+" "+sttime_a[11:19]
        entime_a = entime_a[0:10]+" "+entime_a[11:19]
        #    descending node
        sttime_d = myfunc.read_vdata(vs, "DescendingGridStartTimeUTC")[0][0]
        entime_d = myfunc.read_vdata(vs, "DescendingGridEndTimeUTC")[0][0]
        sttime_d = sttime_d[0:10]+" "+sttime_d[11:19]
        entime_d = entime_d[0:10]+" "+entime_d[11:19]
        #--- close HDF file
        vs.end() # terminate VS interface
        hdf.close()
        
        #--- open HDF file for SD access
        sd = SD.SD(hdffile)

        #--   read latitude, longitude grid
        grid_lat = myfunc.read_sd(sd, "Latitude")
        grid_lon = myfunc.read_sd(sd, "Longitude")

        #--- read surface temperature and pressure
        #    ascending node
        surft_all_a = myfunc.read_sd(sd, "SurfAirTemp_A")
        surfp_all_a = myfunc.read_sd(sd, "SurfPres_Forecast_A")
        #    descending node
        surft_all_d = myfunc.read_sd(sd, "SurfAirTemp_D")
        surfp_all_d = myfunc.read_sd(sd, "SurfPres_Forecast_D")
        #--- read temperature profile and geopotential height
        #    ascending node
        tprof_all_a = myfunc.read_sd(sd, "Temperature_A")
        geoph_all_a = myfunc.read_sd(sd, "GPHeight_A")
        #    descending node
        tprof_all_d = myfunc.read_sd(sd, "Temperature_D")
        geoph_all_d = myfunc.read_sd(sd, "GPHeight_D")
        #--- close HDF file
        sd.end()
        
        
        #create datasets fit for my calc_effective_temp
        tprof_all_a=numpy.rot90(tprof_all_a,1) #lat*height*long
        geoph_all_a=numpy.rot90(geoph_all_a,1)
        tprof_all_d=numpy.rot90(tprof_all_d,1) 
        geoph_all_d=numpy.rot90(geoph_all_d,1)
        
        press_lvl = grid_plvl[::-1]
        lat_range = grid_lat[:,0]
        lon_range = grid_lon[0,:] 
        #create masks for invalid data (written as -9999.)          
        temp_asc_mask   = np.ma.masked_equal(tprof_all_a[::-1,:,::-1], -9999. ) 
        geoph_asc_mask  = np.ma.masked_equal(geoph_all_a[::-1,:,::-1], -9999. )         
        temp_desc_mask  = np.ma.masked_equal(tprof_all_d[::-1,:,::-1], -9999. ) 
        geoph_desc_mask = np.ma.masked_equal(geoph_all_d[::-1,:,::-1], -9999. ) 
        
        #sometimes H is invalid, but not T, or otherwise. masks are combined.
        temp_asc = np.ma.array(temp_asc_mask, mask = (temp_asc_mask.mask + geoph_asc_mask.mask))        
        geoph_asc = np.ma.array(geoph_asc_mask, mask = (temp_asc_mask.mask + geoph_asc_mask.mask)) 
        temp_desc = np.ma.array(temp_desc_mask, mask = (temp_desc_mask.mask + geoph_desc_mask.mask))  
        geoph_desc = np.ma.array(geoph_desc_mask, mask = (temp_desc_mask.mask + geoph_desc_mask.mask)) 
        
        '''
        returns:
        1)  press_lvl : 1d array containing pressure levels, starting in space
        2)  zenith /rad in South Pole coords, 180 == North Pole, 90 == Horizon
        3)  Longitude /rad
        4,5,6,7) masked data of shape (zen x presslvl x lon)
        '''
        return press_lvl,lat_range*np.pi/180 ,lon_range* np.pi/180. ,temp_asc[:,::-1,:],geoph_asc[:,::-1,:],temp_desc[:,::-1,:],geoph_desc[:,::-1,:] , localfile

# This function is used for the calculation of Teff ---------------------------

def calc_Teff_3d(radzenith,press_levels,longitude,alti_data,temp_data,detector_datafile,weights=False):
    '''
    zenith in RAD !!!
    zenith, press_levels,longitude: 1d arrays defining the shape (resolution): (lat x height x long) 
    alti_data,temp_data: 3d datasets from AIRS. arrays of form (lat x height x long)
    detector_datafile: string: name of 2d data set as described in Sebastian Schoenens README.
    out: array of shape (Lat x Long), returns effective Temperatures for all given zeniths and longitudes.
    ascdesc = string, either "asc" or "desc"

    if weights = True:
    calculates dW/dH as function of zenith and height, returns 2darray of shape: zenith x presslvls
    input:   dE: 1d array with Energy bin steps (not logarithmic) /GeV
             dX: 3d array Lat x Presslvls x Long (one dX bin for each pressurelvl and datapoint) /kg m**-2
             Aeff: 2d-array of effective detector Area 
             p: 4d array of shape Lat x pressurelvls x Long x energylvls: Production yield after gaisser-formula
             t: temperature data (3d Lat x height x Long)    
    '''
         
    fill_USSA = False
    if fill_USSA == True:
        
        temperature = np.flipud(np.array([287.43,283.20,278.68,268.57,260.81,251.92,241.44,228.58,220.79,216.65,216.65,216.65,216.65,217.23,220.50,223.13,225.02,227.70,232.71,239.22,249.45,257.88,264.03,270.65]))
        altitude = np.flipud(np.array([110.9 ,  762.0 ,  1457.3 , 3012.2 , 4206.4 , 5574.4 , 7185.4 , 9164.0 , 10362.9 ,11784.0 ,13608.4 ,16179.7 ,18441.6, 20576.2 ,23848.6 ,26481.2 ,28368.1 ,31054.6 ,33452.6 ,35776.5 ,39429.54 ,   42439.8, 44637.2, 47820.0]))

        ticker = 0
        for i in range(len(temp_data)):
            for j in range(len(temp_data[0])):
                for k in range(50):
                    if temp_data[i,j,k] <= 0.:                
                        temp_data[i,j,k] = temperature[j]
                        ticker = ticker+1
                    elif temp_data[i,j,k] == 'nan':
                        temp_data[i,j,k] = temperature[j]
                        ticker = ticker+1
                    if alti_data[i,j,k] <= 0.:
                        alti_data[i,j,k] = altitude[j]
                        ticker = ticker+1
                    elif alti_data[i,j,k] =='nan':
                        ticker = ticker+1
        print('filled data ratio:')
        print(float(ticker) / (  float(len(temp_data)) * float(len(temp_data[0])) * float(len(temp_data[0][0]))))
    
    Aeff, E, dE  = Aeff_read_data_sebastian(np.cos(radzenith),name = detector_datafile)
    thetastar = Theta_Star(radzenith,alti_data)
    x,dx = X(alti_data,press_levels,temp_data,thetastar)    
    prod_yield = P_yield(E,thetastar,x,temp_data)
    
    out1,out2 = Teff(dE,dx,Aeff,prod_yield,temp_data)
    print('Teff ist wirklich done.')
    print(np.shape(out2))
    print(out2[0][0])
    
    
    #fill data with closest valid data point (longitudinal closest favoured)
    fill = False
    if fill == True:
        while np.any(out1.mask) == True:
            #for axis in (1,0):#this would move latitudes too
            for shift in (-1,1):                   
                out1_shifted=np.roll(out1,shift=shift,axis=1)
                idx=~out1_shifted.mask * out1.mask
                out1[idx]=out1_shifted[idx]
                
    print('sollte jetzt mit den weights anfangen')
	
    #Weights, optional
    if weights == True:
        print('doing weights')
        num = dE[np.newaxis,np.newaxis,:] * Aeff[:,np.newaxis,:] * prod_yield[:,:,0,:]
        numsum = np.sum(num, axis = 2)
        numsum_dx = dx[:,:,0] * numsum#2d zenith x presslvl
        
        denom = dE[np.newaxis,np.newaxis,:] * Aeff[:,np.newaxis,:] * prod_yield[:,:,0,:]
        denomsum_de = np.sum(denom, axis = 2)
        denomsum_de_dx = np.sum(  (dx[:,:,0] * denomsum_de), axis = 1)#1d x zenith
        weights = numsum_dx / denomsum_de_dx[:,np.newaxis]        
        weights = np.array(weights)
        delta_H = np.zeros(np.shape(weights))
        
        alti_data = alti_data[:,::-1,:]
        for ith in range(len(alti_data)):
            delta_H[ith,:] = np.append([alti_data[ith,0,0]],np.diff(alti_data[ith,:,0]))#unterster Punkt bekommt height above zero zugewiesen
        delta_H = delta_H[:,::-1]#hab nur die schleife falschrum gecoded

        if delta_H.any() == 0.:
            print('error doing weights: divide by zero!')
            
        diffweights = weights / (delta_H/1000.) * 100#nicer in m,percent
        print('done with weights')
        return out1,diffweights
        
    return out1, out2

# calc_Teff_3d uses the following funtions to do so ---------------------------

def Aeff_read_data_sebastian(costheta,name = "IC2012-13-14_Aachen_eff_area.txt"):
    
    '''
    costheta: 1d array of cos(zenith)
    name= filename of det area, format as in Sebastians README, NO # in line 0 (!)
    returns Aeff(E,costheta), E, dE as 1d-arrays of same len()
    '''

    data = np.genfromtxt(name)

    #log_Ebinedge = np.array(data[1:,0])#log(E) GeV
    #Ebinedge = 10** log_Ebinedge
    #Zbinedge = data[0,2;]
    log_Ebins = np.linspace(1.5, 9.25, 32)
    log_Ebincenter = log_Ebins[:-1] + 0.5 * np.diff(log_Ebins)
    Zbins = np.linspace(-1.,0.1,12) # Same for all Aeff .txts
    Zbincenter = Zbins[:-1] + 0.5 * np.diff(Zbins)
 
    #aeff = data[2:,3:] #m**2* 10000.#cm**2
    aeff = data[1:,2:]

    #Emu = np.array((Ebinedge[:-1])+ np.abs(0.5*np.diff(Ebinedge)))    
    ##Zbincenters = Zbinedge[:-1] + np.abs(0.5*np.diff(Zbinedge))
    #dEmu = np.diff(10**log_Ebinedge)
    Emu = 10 ** log_Ebincenter
    dEmu = np.diff(10 ** log_Ebins)  
  
    #out = np.ones((len(aeff),len(costheta)))
    #ratio_nbins_costh = (float(len(costheta))/float(len(aeff[0])))
    #for ith in range(len(costheta)):
    #    if ith == len(costheta)-1:
    #        i_aeff = len(aeff[0]) - 1
    #    else:
    #        i_aeff = int(ith/ratio_nbins_costh)
    #    for iEn in range(len(aeff)):
    #        out[iEn,ith] = aeff[iEn][i_aeff]
    #out = out.T
    #aeff = aeff.T
    #Zbincenter = np.delete(Zbincenter,-1)
    #out = np.ones((len(costheta),len(aeff[0])))
    #for i in range(len(costheta))
    #    index = np.argmin(abs(abs(Zbincenter) - abs(costhet)))
    #    out[i] = aeff[index]
    
    # Non - interpolated
    
    #aeff = aeff.T
    #indices = np.digitize(costheta,Zbins) - 1
    #out = aeff[indices]

    # Interpolated
    
    out = np.zeros((31,180))
    for energy_bin in range(len(aeff[::,0])):
        out[energy_bin] = np.interp(costheta,Zbincenter,aeff[energy_bin])
    out = out.T
 
    return out, Emu, dEmu
    

def Theta_Star(theta,h):
    '''
    theta = zenith in detector coordinates /rad (90..180 degree, 180 == north pole)
    returns Theta_Star /rad as np.array of shape (zenith x height)
    '''   
    r_e = 6371000.#m   
    h_mask = np.ma.masked_equal(h, -9999. ) #refresh some masks just to be sure 

    out = np.arccos(np.sqrt(1-(r_e/(r_e + h_mask))**2 * np.sin(theta[:,np.newaxis,np.newaxis])**2))
    return np.ma.MaskedArray(out,mask = np.ma.getmask(h_mask))

   
def X(height,pressure,temperature,thetaStar):
    '''
    height,temperature,thetaStar = 3d arrays of shape (zenith x pressure levels x long)
    pressure: 1d array of len(pressure_lvls)
    returns 2d (zen x press levels) arrays X, dX in units kg/m**2    
    '''
    print('calculating X,dX')
    r_e = 6371000.#m
    g = 9.81 #m/s**2
    x  = np.zeros(np.shape(thetaStar))#zenith x heights
    dx = np.zeros(np.shape(thetaStar))
    #x ist masked array
    #x bekommt eine mask, wenn keine Daten vorliegen. 
    
        
    for ilong in range(len(height[0][0])):
        for ith in range(len(height)):
            xdepth = 0.
            for ipress in range(len(height[0])):     
                #calculate x = dp/(g dl)
                if np.cos(thetaStar[ith][ipress][ilong]) != 0. and height[ith,ipress,ilong] > 0.:#sonst wird auch in Bergen etc. ein X berechnet obwohl dort keine Luft mehr ist.
                    if ipress == 0:
                        dpress = pressure[ipress]
                    else:
                        dpress = pressure[ipress] - pressure[ipress-1]      
                        
                    dl = np.sqrt( (np.cos(thetaStar[ith][ipress][ilong]))**2 +(2.*height[ith][ipress][ilong] /r_e)*(1.-(np.cos(thetaStar[ith][ipress][ilong]))**2))    
                    xdepth += dpress*100./(g* dl)#muss in pascal sein

                    x[ith,ipress,ilong] = xdepth
    #dx

            
            dx[ith][0][ilong] = np.sqrt(x[ith][1][ilong] * x[ith][0][ilong])
            dx[ith][-1][ilong] = x[ith][-1][ilong] - np.sqrt(x[ith][-1][ilong] * x[ith][-2][ilong])  
                         
            for ipress in range(1,len(pressure)-1):
                if height[ith,ipress+1,ilong] >= 0.:#otherwise, invalid data is taken into account for low heights.
                    dx[ith][ipress][ilong]= np.sqrt(x[ith][ipress+1][ilong] * x[ith][ipress][ilong]) - np.sqrt(x[ith][ipress][ilong] * x[ith][ipress-1][ilong])
           
    print('X,dX --done') 
    return x,dx


def P_yield(E_nu,Theta_star,atm_X,temp):
    '''
    formula see T.Gaisser
    flux model does not account for the Energy-knee (very minor influence here)
    input:
    '''
    print('calculating P')
#_____________parameters_________________
    
    LAMBDA_N = 107.*10#kg/m2
    LAMBDA_Pi = 147.*10
    LAMBDA_K = 152.*10
    r_Pi = 0.573
    r_K = 0.046
    Z_NPi = 0.079
    Z_NK = 0.0118
    BR_Pinu =1.000
    BR_Knu = 0.635
    gamma = 1.64
    Phi_0 = 1.8*10**4#Gev**-gamma / (sr m**2 s)
    #like Sebastian:
    flux = Phi_0 * E_nu**(-(gamma+1.0))
    for e in E_nu:
        if (e > 3e6) :
            flux *= (3e6/e)**0.30
#_____________calculation_________________(following T.Gaisser) 

    A_Pinu = BR_Pinu * Z_NPi / LAMBDA_N * (1-r_Pi)**gamma / (gamma+1) * np.exp(-atm_X / LAMBDA_N)
    A_Knu  = BR_Knu  * Z_NK  / LAMBDA_N * (1-r_K)**gamma  / (gamma+1) * np.exp(-atm_X / LAMBDA_N)
    
    B_Pinu = (gamma+2)/(gamma+1) * (1/(1-r_Pi)) * (LAMBDA_Pi - LAMBDA_N)/(LAMBDA_Pi * LAMBDA_N) * (atm_X * np.exp(-atm_X/LAMBDA_N) / (np.exp(-atm_X/LAMBDA_Pi) - np.exp(-atm_X / LAMBDA_N)))
    B_Knu  = (gamma+2)/(gamma+1) * (1/(1-r_K))  * (LAMBDA_K  - LAMBDA_N)/(LAMBDA_K  * LAMBDA_N) * (atm_X * np.exp(-atm_X/LAMBDA_N) / (np.exp(-atm_X/LAMBDA_K)  - np.exp(-atm_X / LAMBDA_N)))

    epsilon_Pi = 0.524 * temp#T = T(Theta,atm_X)#!!!assuming atmospheric scale height of ideal gas
    epsilon_K  = 3.897 * temp#T=T(Theta,atm_X) 
    
    #array broadcasting for speed
    Theta_star = Theta_star[:,:,:,np.newaxis]
    A_Pinu = A_Pinu[:,:,:,np.newaxis]
    A_Knu = A_Knu[:,:,:,np.newaxis]
    B_Pinu = B_Pinu[:,:,:,np.newaxis]
    B_Knu = B_Knu[:,:,:,np.newaxis]    
    epsilon_Pi = epsilon_Pi[:,:,:,np.newaxis]
    epsilon_K = epsilon_K[:,:,:,np.newaxis]    

    p_output =  flux * (A_Pinu/ (1.+(B_Pinu  * np.cos(Theta_star)/ epsilon_Pi) * E_nu)  +  (A_Knu / (1.+(B_Knu  * np.cos(Theta_star) / epsilon_K) * E_nu))) 
    #Phi_0 * E_nu**-(gamma+1.)
    print('P shape:',np.shape(p_output))
    print('P --done')
    
    return p_output


def Teff(de,dx,aeff,p,t):
    '''
    calculates effective Temp, 
    input:   dE: 1d array with Energy bin steps (not logarithmic) /GeV
             dX: 3d array Lat x Presslvls x Long (one dX bin for each pressurelvl and datapoint) /kg m**-2
             Aeff: 2d-array of effective detector Area (zenith x Energies)
             p: 4d array of shape Lat x pressurelvls x Long x energylvls: Production yield after gaisser-formula
             t: temperature data (3d Lat x height x Long)
    output:  2d array of shape Lat x Long
    '''
    print('calculating Teff\'s ')
    #print 'shapes dX,aeff,p,t,de'
    #print  np.shape(dx),np.shape(aeff),np.shape(p),np.shape(t),np.shape(de)

    de = de[np.newaxis,:]
    a = de * aeff #(l,i)
    #a = a.T #(i,l)
    a = a[:,np.newaxis,np.newaxis,:]
    t = t[:,:,:,np.newaxis]    
    b = p*t #(i,j,k,l)
     
    numerator = a*b
    denominator = a * p    
    #integrate E          
    numerator_E_int = np.sum(numerator,axis = 3)
    denominator_E_int = np.sum(denominator,axis = 3)
    
    numerator_E_int_dx = np.multiply(dx,numerator_E_int)
    denominator_E_int_dx = np.multiply(dx,denominator_E_int)
    #integrate X
    denominator_E_X_int = np.sum(denominator_E_int_dx,axis = 1)    
    numerator_E_X_int = np.sum(numerator_E_int_dx,axis = 1)
                      
    #print 'shape numerator',np.shape(numerator)
    #print 'shape numerator',np.shape(numerator_E_int) 
    #print 'shape numerator',np.shape(numerator_E_X_int) 
    print('Teff --done')
    return np.multiply(numerator_E_X_int,denominator_E_X_int**-1), denominator_E_X_int #,dx,aeff,p,t

# Unused Functions ------------------------------------------------------------
 
def alphasneu(pressure,temperature,height,theta,integrate_E = False):
    print('calculating alphas')

    LAMBDA_N = 107.*10#kg/m2
    LAMBDA_Pi = 147.*10
    LAMBDA_K = 152.*10
    r_Pi = 0.573
    r_K = 0.046
    Z_NPi = 0.079
    Z_NK = 0.0118
    BR_Pinu =1.000
    BR_Knu = 0.635
    gamma = 1.64#2.0 knee !  
          
    Phi_0 = 1.8*10**4#Gev**-gamma / (sr m**2 s)
    
    thetastar = Theta_Star(theta,height)    
    atm_X,dX = X(height,pressure,temperature,thetastar)
    Aeff,E,dE = Aeff_read_data_christian(np.cos(theta),name = "IC2012-13-14_Aachen_eff_area.txt")

    A_pi= BR_Pinu * (Z_NPi * (1-r_Pi)**gamma) / LAMBDA_N  * np.exp(-atm_X / LAMBDA_N)
    A_k = BR_Knu * (Z_NK * (1-r_K)**gamma) / LAMBDA_N  * np.exp(-atm_X / LAMBDA_N)    
    B_Pi= (gamma+2)/(gamma+1) * 1/(1-r_Pi) * (LAMBDA_Pi - LAMBDA_N)/(LAMBDA_Pi * LAMBDA_N) * atm_X * np.exp(-atm_X / LAMBDA_N) / (np.exp(atm_X/LAMBDA_Pi) - np.exp(atm_X/LAMBDA_N))
    B_k = (gamma+2)/(gamma+1) * 1/(1-r_K) *  (LAMBDA_K - LAMBDA_N)/ (LAMBDA_K * LAMBDA_N)  * atm_X * np.exp(-atm_X / LAMBDA_N) / (np.exp(atm_X/LAMBDA_K)  - np.exp(atm_X/LAMBDA_N))   

def alt_calc_Teff_3d(radzenith,press_levels,longitude,alti_data,temp_data,detector_datafile,weights=False):
    '''
    zenith in RAD !!!
    zenith, press_levels,longitude: 1d arrays defining the shape (resolution): (lat x height x long) 
    alti_data,temp_data: 3d datasets from AIRS. arrays of form (lat x height x long)
    detector_datafile: string: name of 2d data set as described in Sebastians README.
    out: array of shape (Lat x Long), returns effective Temperatures for all given zeniths and longitudes.
    ascdesc = string, either "asc" or "desc"

    if weights = True:
    calculates dW/dH as function of zenith and height, returns 2darray of shape: zenith x presslvls
    input:   dE: 1d array with Energy bin steps (not logarithmic) /GeV
             dX: 3d array Lat x Presslvls x Long (one dX bin for each pressurelvl and datapoint) /kg m**-2
             Aeff: 2d-array of effective detector Area 
             p: 4d array of shape Lat x pressurelvls x Long x energylvls: Production yield after gaisser-formula
             t: temperature data (3d Lat x height x Long)    
    '''
         
    fill_USSA = False
    if fill_USSA == True:
        
        temperature = np.flipud(np.array([287.43,283.20,278.68,268.57,260.81,251.92,241.44,228.58,220.79,216.65,216.65,216.65,216.65,217.23,220.50,223.13,225.02,227.70,232.71,239.22,249.45,257.88,264.03,270.65]))
        altitude = np.flipud(np.array([110.9 ,  762.0 ,  1457.3 , 3012.2 , 4206.4 , 5574.4 , 7185.4 , 9164.0 , 10362.9 ,11784.0 ,13608.4 ,16179.7 ,18441.6, 20576.2 ,23848.6 ,26481.2 ,28368.1 ,31054.6 ,33452.6 ,35776.5 ,39429.54 ,   42439.8, 44637.2, 47820.0]))

        ticker = 0
        for i in range(len(temp_data)):
            for j in range(len(temp_data[0])):
                for k in range(50):
                    if temp_data[i,j,k] <= 0.:                
                        temp_data[i,j,k] = temperature[j]
                        ticker = ticker+1
                    elif temp_data[i,j,k] == 'nan':
                        temp_data[i,j,k] = temperature[j]
                        ticker = ticker+1
                    if alti_data[i,j,k] <= 0.:
                        alti_data[i,j,k] = altitude[j]
                        ticker = ticker+1
                    elif alti_data[i,j,k] =='nan':
                        ticker = ticker+1
        print('filled data ratio:')
        print(float(ticker) / (  float(len(temp_data)) * float(len(temp_data[0])) * float(len(temp_data[0][0]))))
    
    Aeff, E, dE  = Aeff_read_data_sebastian(np.cos(radzenith),name = detector_datafile)
    thetastar = Theta_Star(radzenith,alti_data)
    x,dx = alternativeX(alti_data,press_levels,temp_data,thetastar)    
    prod_yield = P_yield(E,thetastar,x,temp_data)
    
    
  
    
    out1 = Teff(dE,dx,Aeff,prod_yield,temp_data)
    
    #fill data with closest valid data point (longitudinal closest favoured)
    fill = True
    if fill == True:
        while np.any(out1.mask) == True:
            #for axis in (1,0):#this would move latitudes too
            for shift in (-1,1):                   
                out1_shifted=np.roll(out1,shift=shift,axis=1)
                idx=~out1_shifted.mask * out1.mask
                out1[idx]=out1_shifted[idx]

    #Weights, optional
    if weights == True:
        print('doing weights')
        num = dE[np.newaxis,np.newaxis,:] * Aeff[:,np.newaxis,:] * prod_yield[:,:,0,:]
        numsum = np.sum(num, axis = 2)
        numsum_dx = dx[:,:,0] * numsum#2d zenith x presslvl
        
        denom = dE[np.newaxis,np.newaxis,:] * Aeff[:,np.newaxis,:] * prod_yield[:,:,0,:]
        denomsum_de = np.sum(denom, axis = 2)
        denomsum_de_dx = np.sum(  (dx[:,:,0] * denomsum_de), axis = 1)#1d x zenith
        weights = numsum_dx / denomsum_de_dx[:,np.newaxis]        
        weights = np.array(weights)
        delta_H = np.zeros(np.shape(weights))
        
        alti_data = alti_data[:,::-1,:]
        for ith in range(len(alti_data)):
            delta_H[ith,:] = np.append([alti_data[ith,0,0]],np.diff(alti_data[ith,:,0]))#unterster Punkt bekommt height above zero zugewiesen
        delta_H = delta_H[:,::-1]#hab nur die schleife falschrum gecoded

        if delta_H.any() == 0.:
            print('error doing weights: divide by zero!')
            
        diffweights = weights / (delta_H/1000.) * 100#nicer in m,percent
        print('done with weights')
        return out1,diffweights
        
    return out1

def calc_theory_alpha_3d(pressure,temperature,height,theta,integrate_E = False):
    print('calculating alphas')

    lambda_N = 75.*10 #kg/m2
    LAMBDA_N = 107.*10
    LAMBDA_Pi = 147.*10
    LAMBDA_K = 152.*10
    r_Pi = 0.573
    r_K = 0.046
    Z_NPi = 0.079
    Z_NK = 0.0118
    BR_Pinu =1.000
    BR_Knu = 0.635
    gamma = 1.64#2.0 knee !        
    Phi_0 = 1.8*10**4#Gev**-gamma / (sr m**2 s)

    
    thetastar = Theta_Star(theta,height)    
    atm_X,dX = X(height,pressure,temperature,thetastar)
    Aeff,E,dE = Aeff_read_data_sebastian(np.cos(theta),name = "IC2012-13-14_Aachen_eff_area.txt")
    
    flux = Phi_0 * E**-(gamma+1)
    print('langer testrun')
    A_Pinu = BR_Pinu * Z_NPi / lambda_N * (1-r_Pi)**gamma / (gamma+1) * np.exp(-atm_X / LAMBDA_N)
    A_Knu  = BR_Knu  * Z_NK  / lambda_N * (1-r_K)**gamma  / (gamma+1) * np.exp(-atm_X / LAMBDA_N)
    
    B_Pinu = (gamma+2)/(gamma+1) * (1/(1-r_Pi)) * (LAMBDA_Pi - LAMBDA_N)/(LAMBDA_Pi * LAMBDA_N) * (atm_X * np.exp(-atm_X/LAMBDA_N) / (np.exp(-atm_X/LAMBDA_Pi) - np.exp(-atm_X / LAMBDA_N)))
    B_Knu  = (gamma+2)/(gamma+1) * (1/(1-r_K))  * (LAMBDA_K  - LAMBDA_N)/(LAMBDA_K  * LAMBDA_N) * (atm_X * np.exp(-atm_X/LAMBDA_N) / (np.exp(-atm_X/LAMBDA_K)  - np.exp(-atm_X / LAMBDA_N)))

    epsilon_Pi = 0.524 * temperature#T = T(Theta,atm_X)#!!!assuming atmospheric scale height of ideal gas
    epsilon_K  = 3.897 * temperature#T=T(Theta,atm_X)    
    print('A_Pinu')
    print( A_Pinu[::5,::5,::5])
    print( 'A_Knu')
    print( A_Knu[::5,::5,::5])  
    print( 'B_Pinu')
    print( B_Pinu[::5,::5,::5])  
    print( 'B_Knu')
    print (B_Knu[::5,::5,::5])   
    
    
    thetastar = thetastar[:,:,:,np.newaxis]
    A_Pinu = A_Pinu[:,:,:,np.newaxis]
    A_Knu = A_Knu[:,:,:,np.newaxis]
    B_Pinu = B_Pinu[:,:,:,np.newaxis]
    B_Knu = B_Knu[:,:,:,np.newaxis]    
    epsilon_Pi = epsilon_Pi[:,:,:,np.newaxis]
    epsilon_K = epsilon_K[:,:,:,np.newaxis] 
    #numerator
    
    numerator = (A_Pinu * B_Pinu * np.cos(thetastar) * E / epsilon_Pi ) / (1. + B_Pinu * np.cos(thetastar) * E / epsilon_Pi)**2+ (A_Knu * B_Knu * np.cos(thetastar) * E / epsilon_K ) / (1. + B_Knu * np.cos(thetastar) * E / epsilon_K)**2
    numerator *= flux * dE
    numerator_x_int = np.sum(numerator * dX[:,:,:,np.newaxis], axis = 1)
    #denominator
    denominator = (A_Pinu ) / (1. + B_Pinu * np.cos(thetastar) * E / epsilon_Pi)**2 +(A_Pinu ) / (1. + B_Pinu * np.cos(thetastar) * E / epsilon_Pi)**2
    denominator *= flux * dE
    denominator_x_int = np.sum(denominator * dX[:,:,:,np.newaxis], axis = 1)
    
    countercheck = 0.
    for check in np.nditer(denominator_x_int):
        if check == 0.:
            countercheck+=1.
    print('denom was zero',countercheck,'times!')
    print(np.size(denominator_x_int))
    
    alphas_th_lon_e = numerator_x_int / denominator_x_int
    
    
    print('hier die alphas dazu')
    acounter = 0.
    for acheck in np.nditer(alphas_th_lon_e):
        if  acheck >= 1. :
            acounter += 1. 
    print('alphas >1 ',acounter)
    print(alphas_th_lon_e[0,0,:])
    if integrate_E == True:
    
        return np.mean(alphas_th_lon_e,axis = 2)
    
    return alphas_th_lon_e

def Aeff_read_data_christian(costheta,name = "IC2012-13-14_Aachen_eff_area.txt"):#here, this is an arbitrary argument
	
    B = np.load("IC2012_13_14_effA_2D.npz")
    coszen_B = B['cos_zen_bins']
    e_B = B['e_bins']
    aeff = B['eff_area']
    Zctr = []
    log_Ectr = []
    
    for i in range(len(coszen_B)-1):
        Zctr.append((coszen_B[i+1]+coszen_B[i])/2.)
    for j in range(len(e_B)-1):
        log_Ectr.append((e_B[j+1]+e_B[j])/2.)
    
    Zctr = np.array(Zctr)
    log_Ectr = np.array(log_Ectr)
    Ectr = 10**(log_Ectr) 
    dE_seb = np.diff(np.exp(e_B))
    out = np.ones((len(aeff),len(costheta)))
    ratio_nbins_costh = (float(len(costheta))/float(len(aeff[0])))

    for ith in range(len(costheta)):
        if ith == len(costheta)-1:

            i_aeff = len(aeff[0]) - 1


        else:
            i_aeff = int(ith/ratio_nbins_costh)

        for iEn in range(len(aeff)):
            out[iEn,ith] = aeff[iEn][i_aeff]
    return out.T,Ectr,dE_seb#[:,::-1]


