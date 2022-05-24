#!/usr/local/env python


import sys,os,string

import astropy
from astropy.io import fits
import numpy as np
import astropy.io.fits as pyfits
import scipy
from scipy import ndimage
from scipy import ndimage, stats, interpolate, signal
import time


imagen_in=sys.argv[1]
imagen_fin=imagen_in[:-5]
pyfits.info(imagen_in)
#lee headers
header=pyfits.getheader(imagen_in)

#print header['naxis']

#lee la imagen
#imagen=pyfits.getdata(imagen_in,0)
imagen=fits.getdata(imagen_in,0)

# largo de la imagen
Xsize= imagen.size/len(imagen)
Ysize=len(imagen)
print( Xsize,Ysize)


sx = ndimage.sobel(imagen, axis=-1, mode='reflect')
#sy = ndimage.sobel(im, axis=0, mode='constant')
#sob = np.hypot(sx, im)
sob = np.hypot(sx,imagen)


resta=(sob - imagen)  

#pyfits.writeto("RESTA.fits",resta,header)
#pyfits.writeto("RESTA_SX.fits",sx,header)
#pyfits.writeto("RESTA_SOB.fits",sob,header)

#test=ndimage.spline_filter(imagen,order=5)
#pyfits.writeto("RESTA_DATA.fits",test,header)

#crear grid de 100 x 100 e interpolar y restar a ver que queda...



#def median_interpolate(im_data, im_mask, inter_type='cubic', inter_width=64, grid_sampling=10):
def median_interpolate(im_data, inter_type='tps', inter_width=64, grid_sampling=10):
    ny, nx = im_data.shape
    im_data = im_data.copy()
#    im_data[im_mask] = np.nan
  
    y=0
    xx=[]
    yy=[]
    zz=[]
    while y < ny:
        y_width=inter_width
        if y+y_width > ny: y_width = ny-y
        x=0
        while x < nx:
            x_width=inter_width
            if x+x_width > nx: x_width = nx-x
            xx.append(int(x+x_width/2))
            yy.append(int(y+y_width/2))
            gv=np.isfinite(im_data[y:y+y_width,x:x+x_width])
            if np.sum(gv)*1./(x_width*y_width) > 0.7:
                zz.append(np.nanmedian(im_data[y:y+y_width,x:x+x_width],axis=None))
            else:
                zz.append(np.nan)
            x+=x_width
        y+=y_width

    xx = np.array(xx)
    yy = np.array(yy)
    zz = np.array(zz)
    print (zz)

    start_time = time.time()
    gv=np.isfinite(zz)
    if inter_type=='cubic':
        print("GridData cubic interpolation - STARTING")
        grid_x, grid_y = np.meshgrid(np.arange(0,nx,grid_sampling), np.arange(0,ny,grid_sampling))
        im_model_data=interpolate.griddata(np.stack((xx[gv],yy[gv]),axis=-1), zz[gv], (grid_x, grid_y), method='cubic',fill_value=0.) 
    elif inter_type=='inverse':
        print("RBF inverse interpolation - STARTING")
        grid_x, grid_y = np.meshgrid(np.arange(0,nx,grid_sampling), np.arange(0,ny,grid_sampling))
        rbf_model = interpolate.Rbf(xx[gv],yy[gv],zz[gv],function='inverse')
        im_model_data = rbf_model(grid_x,grid_y)
    elif inter_type=='multiquadric':
        print("RBF inverse interpolation - STARTING")
        grid_x, grid_y = np.meshgrid(np.arange(0,nx,grid_sampling), np.arange(0,ny,grid_sampling))
        rbf_model = interpolate.Rbf(xx[gv],yy[gv],zz[gv],function='multiquadric', epsilon=100)
        im_model_data = rbf_model(grid_x,grid_y)
    elif inter_type=='tps':
        print("RBF tps interpolation - STARTING")
        grid_x, grid_y = np.meshgrid(np.arange(0,nx,grid_sampling), np.arange(0,ny,grid_sampling))
        rbf_model = interpolate.Rbf(xx[gv],yy[gv],zz[gv],function='thin_plate')
        im_model_data = rbf_model(grid_x,grid_y)
        
    print("Interpolation done")
    print("Total execution time: %.1f minutes" % ((time.time() - start_time)/60.))

    #    gv_nan=np.isnan(im_data)
    #    if np.sum(gv_nan)>0:
    #        nim_data[gv_nan]=nim_model_data[gv_nan]       
    im_data=np.nan_to_num(im_data)
    
    return im_data, im_model_data

#test2=median_interpolate(imagen, im_mask, inter_type='cubic', inter_width=64, grid_sampling=10)

#ok para el patron de lineas


for i in range(657,665):
    test2,test3=median_interpolate(imagen[i,:,:], inter_type='cubic', inter_width=20,grid_sampling=1)

    #puede ser 8
    #test2,test3=median_interpolate(imagen, inter_type='inverse', inter_width=4,grid_sampling=1) 


    pyfits.writeto("resta_"+imagen_fin+str(i)+"_20test.fits",test2,header)
    pyfits.writeto("model_"+imagen_fin+str(i)+"_20test.fits",test3,header)
