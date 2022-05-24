from time import perf_counter as clock
from os import path

import sys,os,string
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
from numpy import linspace,exp
from numpy.random import randn
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import UnivariateSpline
import time
from  scipy import ndimage
from numpy.polynomial import legendre
from scipy import signal
import multiprocessing
from multiprocessing import Pool,cpu_count,Queue,Manager



start_time1 = time.time()
start_time2 = time.time()


#file='spectra_188_95.fits'
#CUBO='MUSE_NGC1427A_cut.fits'
#CUBO='KARTWHEEL_RESAMPLED_restframe.fits' 
CUBO=sys.argv[1]


hdu = fits.open(CUBO)
t2 = hdu[0].data
CD1=hdu[0].header['CD3_3']
CRVAL=hdu[0].header['CRVAL3']
A1= hdu[0].data
lamb,eje_y,eje_x=A1.shape
z_dim,y_dim,x_dim=A1.shape

A=np.where(np.isnan(A1),0,A1)
LAMBDA=np.arange(1,z_dim+1)*CD1+CRVAL

model='orion_smooth_001.fits'
mod=fits.open(model)
t=mod[1].data



flux2=[]
wave2=[]
z = 0.00

for j in range (len(t['flux'][:])):
    
    flux2.append(t['flux'][j])
    wave2.append(float(t['wavelength'][j])/(1-z))

#global flux_m
#global wave_m

flux_m=np.asarray(flux2)
wave_m=np.asarray(wave2)

#def aplanar(A,t,LAMBDA,flux_m,wave_m,xx,yy):
def aplanar(xx,yy):
    global flux_m
    global wave_m
    global A

    try:

        c = 299792.458  # speed of light in km/s
        #velscale = c*np.log(wave_c[1]/wave_c[0])
        velscale = c*np.log(LAMBDA[1]/LAMBDA[0])
        
        FWHM_gal = 2.1 
        
        mask = (wave_m >= LAMBDA[0]) & (wave_m <= LAMBDA[len(LAMBDA)-1])
        flux_m = flux_m[mask]
        flux_m = flux_m/np.median(flux_m)   # Normalize spectrum to avoid numerical issues
        wave_m = wave_m[mask]
        
        flux_c=A[:,xx,yy]
    
        mediana=np.median(flux_c)
        flux_c=flux_c/np.median(flux_c)
        wave_c=LAMBDA
        regul_err = 0.013
        
        lam_range_gal = np.array([np.min(wave_m), np.max(wave_m)])/(1 + z)
        
        flux_mi=np.interp(wave_c,wave_m,flux_m,left=np.min(wave_c), right=np.max(wave_c))
        flux_mi=np.interp(wave_c,wave_m,flux_m)
        
        template=np.column_stack([wave_c,flux_mi])
        galaxy=np.column_stack([wave_c,flux_c])


        #noise=np.full_like(galaxy[:,1], 0.06635)
        noise=np.full_like(galaxy[:,1], 0.01)
        vel = 0#c*np.log(1 + z)   # eq.(8) of Cappellari (2017)
        vel = c*np.log(1 + z)   # eq.(8) of Cappellari (2017)
        start = [vel, 0.] 
        moments=[2,1]#1,1,1]
        velscale=800#2000
        regul_err = 0.013#0.013

        #    for n in range(len(galaxy[:,1])):
        #        print(galaxy[n,1],flux_c[n],flux_co[n])

        pp = ppxf(template[:,1],galaxy[:,1],noise,velscale,start,vsyst=9,plot=False,
                  clean=True, degree=4, mdegree=1,component=0,linear=False,bias=10)#,reg_dim=reg_dim)

        fit=np.column_stack((wave_c[pp.goodpixels],pp.bestfit[pp.goodpixels]))
        x = np.linspace(-1, 1, len(galaxy))
        mpoly2 = legendre.legval(x, np.append(1, pp.mpolyweights))
        
        s=UnivariateSpline(wave_c[pp.goodpixels],flux_c[pp.goodpixels],k=5,s=90)
        ys=s(wave_c[pp.goodpixels])
    
        #plt.plot(ys)
        flux_mm=np.interp(wave_c,wave_c[pp.goodpixels],ys)
        flux_restado=np.interp(wave_c,wave_c[pp.goodpixels],ys)

        #pp.plot()
        #plt.show()
        #print("Mediana")
        #print(mediana)
    
        #plt.plot(wave_c,flux_c)#*np.median(flux_c))
        #print(np.median(flux_c))
    
        #plt.plot(wave_c,flux_c*mediana,'.')
        #plt.plot(wave_c,(flux_c-flux_restado)*mediana,'red')
        #plt.plot(wave_c[pp.goodpixels],ys,'black')
        #plt.show()
        #A[:,yy,xx]=(flux_c-flux_restado)*mediana

        

        return((flux_c-flux_restado)*mediana)
    except:
        out=np.empty(len(flux_c))
        out[:]=np.nan
        return((out))
   
#return(A)
 
#aplanar(434,89)
#exit()
#aplanar(t,CUBO,hdu,235,235)       
#aplanar(t,CUBO,hdu,200,104)    
#aplanar(t,CUBO,hdu,193,97)
#aplanar(t,CUBO,hdu,191,104)    

tmp= [(x, y) for x in range(y_dim) for y in range(x_dim)]
#tmp= [(x, y) for x in range(200,240) for y in range(50,70)]
#tmp= [(x, y) for x in range(0,80) for y in range(0,80)]
#tmp= [(x, y) for x in range(20,24) for y in range(50,57)]
p,q=np.asarray(tmp).reshape(-1,2).T
#print (p[2],q[2])

#print(p,q)
#print(p.shape,q.shape)

CUBO_extendido=np.asarray(list(map(aplanar,p,q)))

#print(type(CUBO_extendido))
a,b=CUBO_extendido.shape
#print(a,b)
#print(CUBO_extendido)


for i in range (a):
    A1[:,int(p[i]),int(q[i])]=CUBO_extendido[i,:]


#f=open("NO_FIT_flux.txt","w")#




#for j in range(200,240):# (eje_y):
#    try:
#        for i in range(50,70):#(eje_x):
#            try:
#                #A[:,j,i]=aplanar(t,CUBO,hdu,i,j)
#                A1[:,j,i]=aplanar(i,j)
#            except:
#                f.write(str(i)+' '+str(j)+'\n')
#                pass
#    except:
#        f.write(str(i)+' '+str(j)+'\n')
#        pass#

cubo_out=sys.argv[1][:-5]+'_NO_BKGND.fits'
fits.writeto(cubo_out,A1,hdu[0].header)
#f.close()

print("Total execution time: %.2f min" % ((time.time() - start_time1)/60))


#def main():
#    #sys.tracebacklimit = 0 #solo funciona en 3.7?##
#
#    manager = multiprocessing.Manager()
#    return_dict = manager.dict()
#    subprocs = []
#
#    for j in range(90,102):#(eje_y):
#        #try:
#        for i in range(176,182):#(eje_x):
#            try:
#            p=multiprocessing.Process(target=aplanar,args=(t,CUBO,hdu,i,j))
#            A[:,j,i]=aplanar(t,CUBO,hdu,i,j)
#            except:
#                f.write(str(j),str(i))
#                pass
#            except:
#                f.write(str(j),str(i))
#                pass
#            subprocs.append(p)
#        print("Starting all processes now")#
#
#        # now start them all
#        for p in subprocs: p.start()
#    
#        # wait on them all
#        for p in subprocs: p.join()#
#
#        #f.close()    
#        B=return_dicy.values()
#        fits.writeto("TEST_CUBO_paralelo_noNAN.fits",A,hdu[0].header)
#        print("Total execution time: %.2f min" % ((time.time() - start_time1)/60))
#        return True#
#
#main()


#aplanar(t,CUBO,hdu,250,185)
#aplanar(t,CUBO,hdu,250,184)
#aplanar(t,CUBO,hdu,237,108)
#aplanar(t,CUBO,hdu,236,108)
#aplanar(t,CUBO,hdu,238,108)

#aplanar(t,CUBO,hdu,219,246)
#aplanar(t,CUBO,hdu,187,93) 
#aplanar(t,CUBO,hdu,63,57) 
#aplanar(t,CUBO,hdu,20,43)
#aplanar(t,CUBO,hdu,22,196)
#aplanar(t,CUBO,hdu,238,203)
#aplanar(t,CUBO,hdu,236,201)
