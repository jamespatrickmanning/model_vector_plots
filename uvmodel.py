# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:58:18 2013

@author: jmanning
routine to generate a vector plot from user-specified time and place
reads a control file for imput parameters

modified in Dec 2015 to test on both eMOLT and COMET machine
revisited in Dec 2018 to overlay on drifter tracks in CC Bay
"""

from pylab import *
from matplotlib.collections import PolyCollection
import matplotlib.tri as Tri
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4
import sys
import numpy as np
from datetime import timedelta
#from pydap.client import open_url
#lon_bin_size=2.0
#lat_bin_size=2.0
# read HARDCODES listed in control file ######################
urlname=open("ctrl_uvmodel.csv", "r").readlines()[0][27:-1]
depth=int(open("ctrl_uvmodel.csv", "r").readlines()[1][22:-1])
TIME=open("ctrl_uvmodel.csv", "r").readlines()[2][31:-1]
lon_bin_size=float(open("ctrl_uvmodel.csv", "r").readlines()[3][32:-1])
lat_bin_size=float(open("ctrl_uvmodel.csv", "r").readlines()[4][32:-1])
##############################################################
def sh_bindata(x, y, z, xbins, ybins):
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
    return xb,yb,zb_mean

if urlname=="massbay":
    TIME=dt.datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S") 
    now=dt.datetime.now()
    if TIME>now:
         diff=(TIME-now).days
    else:
         diff=(now-TIME).days
    if diff>3:
        print "please check your input start time,within 3 days both side form now on"
        sys.exit(0)   
print urlname
if urlname=="30yr":
      stime=dt.datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S")
      timesnum=stime.year-1981
      standardtime=dt.datetime.strptime(str(stime.year)+'-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
      timedeltaprocess=(stime-standardtime).days
      startrecord=26340+35112*(timesnum/4)+8772*(timesnum%4)+1+timedeltaprocess*24     
      url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?lon,lat,lonc,latc,time,nv,h,siglay,v,u'
else:
    timeperiod=(TIME)-(now-timedelta(days=3))
    startrecord=(timeperiod.seconds)/60/60
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?lon,lat,lonc,latc,time,nv,h,siglay,v,u'
nc = netCDF4.Dataset(url)
print "hold on readin nc"
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
latc = nc.variables['latc'][:]
lonc = nc.variables['lonc'][:]
h = nc.variables['h'][:]
siglay=nc.variables['siglay']
u= nc.variables['u']
v= nc.variables['v']
nv = nc.variables['nv'][:].T - 1
print 'we have the model data now we want to get depth of interest'
print 'where we might be able to speed this up with, for example, utotal=u[startrecord,0,:] for the case of surface?'
utotal=[]
vtotal=[]
if depth==-1: # case of surface flow
  utotal=u[startrecord,0,:]
  vtotal=v[startrecord,0,:]
else:
  for i in range(len(lon)):
    depthtotal=siglay[:,i]*h[i]
    layer=np.argmin(abs(depthtotal+depth))
    utotal.append(u[startrecord,layer,i])
    vtotal.append(v[startrecord,layer,i])
  utotal=np.array(utotal)
  vtotal=np.array(vtotal)

print 'now lets bin the data'
xi = np.arange(min(lon)-0.1,max(lon)+0.1,lon_bin_size)
yi = np.arange(min(lat)-0.1,max(lat)+0.1,lat_bin_size)
xb,yb,ub_mean = sh_bindata(lon[::-1], lat[::-1], utotal, xi, yi)
xb,yb,vb_mean = sh_bindata(lon[::-1], lat[::-1], vtotal, xi, yi)
xxb,yyb = np.meshgrid(xb, yb)
latsize=[min(lat)-0.6,max(lat)+0.6]
lonsize=[min(lon)-0.6,max(lon)+0.6]
print 'and plot'
plt.figure()
m = Basemap(projection='cyl',llcrnrlat=min(latsize),urcrnrlat=max(latsize),\
            llcrnrlon=min(lonsize),urcrnrlon=max(lonsize),resolution='h')#,fix_aspect=False)
m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,5),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,10),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
ub = np.ma.array(ub_mean, mask=np.isnan(ub_mean))
vb = np.ma.array(vb_mean, mask=np.isnan(vb_mean))
Q=m.quiver(xxb,yyb,ub,vb,scale=10)

plt.quiverkey(Q,0.8,0.05,1, '1m/s', labelpos='W')
plt.title(urlname+' Depth:'+str(depth)+' Time:'+str(TIME)) 
plt.savefig(urlname+'_vectors'+str(TIME)+'.png')
