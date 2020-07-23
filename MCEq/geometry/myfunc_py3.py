# -*- coding: utf-8 -*-
"""
Created on Mon May 22 12:30:46 2017

@author: Jakob
"""

#!/usr/bin python

import os
import tarfile
import numpy
#from pyhdf import HDF, SD, VS
#import warnings
#warnings.simplefilter('ignore', DeprecationWarning)

#--- function to read vdata
def read_vdata(vs, attribute):
    vd = vs.attach(attribute)
    vdata = vd.read()
    vd.detach()
    return vdata

#--- function to read sd data
def read_sd(sd, attribute):
    sds  = sd.select(attribute)
    sddata = sds.get()
    sds.endaccess()
    return sddata

#--- get mean value over the south pole
def get_mean_2d(data_2d):
    # get mean of 360 longitudinal bin, standard deviation
    # numpy[180,360] -> numpy[180]
    data_2d_ma    = numpy.ma.masked_equal(data_2d, -9999.) # mask no data bin
    data_mean     = data_2d_ma.mean(axis=1)
    data_mean_std = data_2d_ma.std(axis=1)
    data_mean_num = data_2d_ma.count(axis=1)
    data_mean     = numpy.ma.filled(data_mean,     fill_value=-1.)
    data_mean_std = numpy.ma.filled(data_mean_std, fill_value=-1.)
    return data_mean, data_mean_std, data_mean_num

def get_mean_3d(data_3d):
    # get mean of 360 longitudinal bin, standard deviation
    # numpy[24,180,360] -> numpy[24,180]
    data_3d_ma    = numpy.ma.masked_equal(data_3d, -9999.) # mask no data bin
    data_mean     = data_3d_ma.mean(axis=2)
    data_mean_std = data_3d_ma.std(axis=2)
    data_mean_num = data_3d_ma.count(axis=2)
    data_mean     = numpy.ma.filled(data_mean,     fill_value=-1.)
    data_mean_std = numpy.ma.filled(data_mean_std, fill_value=-1.)
    return data_mean, data_mean_std, data_mean_num
#--- write all data
def write_all(
    txtfile, grid_plvl, grid_lat, sttime_a, entime_a, sttime_d, entime_d,
    surft_a, surft_a_std, surft_a_num, surft_d, surft_d_std, surft_d_num,
    surfp_a, surfp_a_std, surfp_a_num, surfp_d, surfp_d_std, surfp_d_num,
    tprof_a, tprof_a_std, tprof_a_num, tprof_d, tprof_d_std, tprof_d_num,
    geoph_a, geoph_a_std, geoph_a_num, geoph_d, geoph_d_std, geoph_d_num):

    fout = open(txtfile, "w")
    #    ascending node
    fout.write( "#------  ---- ascending node\n" )
    fout.write( "# start  %s\n" % sttime_a )
    fout.write( "#   end  %s\n" % entime_a )
    buff = "".join(["  ----------------------------- %6.1f hPa" % elem for elem in grid_plvl])
    fout.write( "".join(["#------  -------------------- surface observation", '%s' % buff, "\n"]) )
    buff = "    geoph    std nave    tprof    std nave"
    fout.write( "".join(["#   lat    surfp    std nave    surft    std nave", buff * 24, "\n"]) )
    for ilat in range(0,180):
        surf_lat_a = '%9.1f%7.1f%5i%9.2f%7.2f%5i' % (
            surfp_a[ilat], surfp_a_std[ilat], surfp_a_num[ilat],
            surft_a[ilat], surft_a_std[ilat], surft_a_num[ilat])
        tprof_lat_a     = [ "%.2f" % elem for elem in tprof_a[0:,ilat] ]
        tprof_lat_a_std = [ "%.2f" % elem for elem in tprof_a_std[0:,ilat] ]
        tprof_lat_a_num = [   "%i" % elem for elem in tprof_a_num[0:,ilat] ]
        geoph_lat_a     = [ "%.1f" % elem for elem in geoph_a[0:,ilat] ]
        geoph_lat_a_std = [ "%.1f" % elem for elem in geoph_a_std[0:,ilat] ]
        geoph_lat_a_num = [   "%i" % elem for elem in geoph_a_num[0:,ilat] ]
        prof_lat_a = tuple(
            ['%9s%7s%5s' * 2 % (geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num)
             for geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num
             in zip(geoph_lat_a, geoph_lat_a_std, geoph_lat_a_num,
                    tprof_lat_a, tprof_lat_a_std, tprof_lat_a_num)])
        fout.write( "".join(['%7.1f' % grid_lat[ilat,0], '%42s' % surf_lat_a, '%42s' * 24 % prof_lat_a, "\n"]) )
    #   descending node
    fout.write( "#------  --- descending node\n" )
    fout.write( "# start  %s\n" % sttime_d )
    fout.write( "#   end  %s\n" % entime_d )
    buff = "".join(["  ----------------------------- %6.1f hPa" % elem for elem in grid_plvl])
    fout.write( "".join(["#------  -------------------- surface observation", '%s' % buff, "\n"]) )
    buff = "    geoph    std nave    tprof    std nave"
    fout.write( "".join(["#   lat    surfp    std nave    surft    std nave", buff * 24, "\n"]) )
    for ilat in range(0,180):
        surf_lat_d = '%9.1f%7.1f%5i%9.2f%7.2f%5i' % (
            surfp_d[ilat], surfp_d_std[ilat], surfp_d_num[ilat],
            surft_d[ilat], surft_d_std[ilat], surft_d_num[ilat])
        tprof_lat_d     = [ "%.2f" % elem for elem in tprof_d[0:,ilat] ]
        tprof_lat_d_std = [ "%.2f" % elem for elem in tprof_d_std[0:,ilat] ]
        tprof_lat_d_num = [   "%i" % elem for elem in tprof_d_num[0:,ilat] ]
        geoph_lat_d     = [ "%.1f" % elem for elem in geoph_d[0:,ilat] ]
        geoph_lat_d_std = [ "%.1f" % elem for elem in geoph_d_std[0:,ilat] ]
        geoph_lat_d_num = [   "%i" % elem for elem in geoph_d_num[0:,ilat] ]
        prof_lat_d = tuple(
            ['%9s%7s%5s' * 2 % (geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num)
             for geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num
             in zip(geoph_lat_d, geoph_lat_d_std, geoph_lat_d_num,
                    tprof_lat_d, tprof_lat_d_std, tprof_lat_d_num)])
        fout.write( "".join(['%7.1f' % grid_lat[ilat,0], '%42s' % surf_lat_d, '%42s' * 24 % prof_lat_d, "\n"]) )
    fout.close()
    #--- compress to tar.gz
    # dirname  = os.path.dirname(txtfile)
    # os.chdir(dirname)
    # basename = os.path.basename(txtfile)
    # tar = os.path.splitext(basename)[0]+".tar.gz"
    # fout = tarfile.open(tar, "w:gz")
    # fout.add(basename)
    # fout.close()
def write_all2(
    txtfile, grid_plvl, grid_lat,grid_lon, sttime_a, entime_a,
        sttime_d, entime_d,tprof_all_a,geoph_all_a,tprof_all_d,geoph_all_d):

    fout = open(txtfile, "w")
    #    ascending node
    fout.write( "#------  ---- ascending node\n" )
    fout.write( "# start  %s\n" % sttime_a )
    fout.write( "#   end  %s\n" % entime_a )
    buff = "".join(["  ----------------------------- %6.1f hPa" % elem for elem in grid_plvl])
    fout.write( "".join(["#------  -------------------- surface observation", '%s' % buff, "\n"]) )
    buff = "    geoph    std nave    tprof    std nave"
    fout.write( "".join(["#   lat    surfp    std nave    surft    std nave", buff * 24, "\n"]) )
    for ilat in range(0,180):
        surf_lat_a = '%9.1f%7.1f%5i%9.2f%7.2f%5i' % (
            surfp_a[ilat], surfp_a_std[ilat], surfp_a_num[ilat],
            surft_a[ilat], surft_a_std[ilat], surft_a_num[ilat])
        tprof_lat_a     = [ "%.2f" % elem for elem in tprof_a[0:,ilat] ]
        tprof_lat_a_std = [ "%.2f" % elem for elem in tprof_a_std[0:,ilat] ]
        tprof_lat_a_num = [   "%i" % elem for elem in tprof_a_num[0:,ilat] ]
        geoph_lat_a     = [ "%.1f" % elem for elem in geoph_a[0:,ilat] ]
        geoph_lat_a_std = [ "%.1f" % elem for elem in geoph_a_std[0:,ilat] ]
        geoph_lat_a_num = [   "%i" % elem for elem in geoph_a_num[0:,ilat] ]
        prof_lat_a = tuple(
            ['%9s%7s%5s' * 2 % (geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num)
             for geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num
             in zip(geoph_lat_a, geoph_lat_a_std, geoph_lat_a_num,
                    tprof_lat_a, tprof_lat_a_std, tprof_lat_a_num)])
        fout.write( "".join(['%7.1f' % grid_lat[ilat,0], '%42s' % surf_lat_a, '%42s' * 24 % prof_lat_a, "\n"]) )
    #   descending node
    fout.write( "#------  --- descending node\n" )
    fout.write( "# start  %s\n" % sttime_d )
    fout.write( "#   end  %s\n" % entime_d )
    buff = "".join(["  ----------------------------- %6.1f hPa" % elem for elem in grid_plvl])
    fout.write( "".join(["#------  -------------------- surface observation", '%s' % buff, "\n"]) )
    buff = "    geoph    std nave    tprof    std nave"
    fout.write( "".join(["#   lat    surfp    std nave    surft    std nave", buff * 24, "\n"]) )
    for ilat in range(0,180):
        surf_lat_d = '%9.1f%7.1f%5i%9.2f%7.2f%5i' % (
            surfp_d[ilat], surfp_d_std[ilat], surfp_d_num[ilat],
            surft_d[ilat], surft_d_std[ilat], surft_d_num[ilat])
        tprof_lat_d     = [ "%.2f" % elem for elem in tprof_d[0:,ilat] ]
        tprof_lat_d_std = [ "%.2f" % elem for elem in tprof_d_std[0:,ilat] ]
        tprof_lat_d_num = [   "%i" % elem for elem in tprof_d_num[0:,ilat] ]
        geoph_lat_d     = [ "%.1f" % elem for elem in geoph_d[0:,ilat] ]
        geoph_lat_d_std = [ "%.1f" % elem for elem in geoph_d_std[0:,ilat] ]
        geoph_lat_d_num = [   "%i" % elem for elem in geoph_d_num[0:,ilat] ]
        prof_lat_d = tuple(
            ['%9s%7s%5s' * 2 % (geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num)
             for geoph, geoph_std, geoph_num, tprof, tprof_std, tprof_num
             in zip(geoph_lat_d, geoph_lat_d_std, geoph_lat_d_num,
                    tprof_lat_d, tprof_lat_d_std, tprof_lat_d_num)])
        fout.write( "".join(['%7.1f' % grid_lat[ilat,0], '%42s' % surf_lat_d, '%42s' * 24 % prof_lat_d, "\n"]) )
    fout.close()
    #--- compress to tar.gz
    # dirname  = os.path.dirname(txtfile)
    # os.chdir(dirname)
    # basename = os.path.basename(txtfile)
    # tar = os.path.splitext(basename)[0]+".tar.gz"
    # fout = tarfile.open(tar, "w:gz")
    # fout.add(basename)
    # fout.close()