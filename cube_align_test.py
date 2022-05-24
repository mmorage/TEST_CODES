#!/usr/bin/env python

import sys,os,string
from pyraf import iraf
from iraf import images,immatch, imalign, imutil,imshift

cubo_base=sys.argv[1]
cubo_desalign=sys.argv[2]
#archivo_ref_points_base=sys.argv[3]
datab=sys.argv[3]


def get_cube_dim(CUBO):
    #"get the dimension of the cube as int 
    #map(int,NNNN) change the list to int"
    head1=iraf.imheader(CUBO,Stdout=1)
    X,Y,Z=map(int,(head1[0].split('[')[1].split("]")[0].split(",")))
    return (X,Y,Z)
 
x_base,y_base,z_base=get_cube_dim(cubo_base)
x_desalign,y_desalign,z_desalign=get_cube_dim(cubo_desalign)

print cubo_base, x_base,y_base,z_base
print cubo_desalign,x_desalign,y_desalign,z_desalign


#       input  = "ANCI_SHIFT1.fits" Input images
#    reference = "ANCI_REF.fits"    Reference image
#       coords = "imalign.txt"      Reference coordinates file
#       output = "TEST_imag"        Output images
#      (shifts = "0.8,1")           Initial shifts file
#     (boxsize = 7)                 Size of the small centering box
#      (bigbox = 11)                Size of the big centering box
#    (negative = no)                Are the features negative ?
#  (background = INDEF)             Reference background level
#       (lower = INDEF)             Lower threshold for data
#       (upper = INDEF)             Upper threshold for data
#    (niterate = 3)                 Maximum number of iterations
#   (tolerance = 0)                 Tolerance for convergence
#    (maxshift = INDEF)             Maximum acceptable pixel shift
# (shiftimages = yes)               Shift the images ?
# (interp_type = "spline3")         Interpolant
#(boundary_typ = "nearest")         Boundary type
#    (constant = 0.0)               Constant for constant boundary extension
#  (trimimages = no)                Trim the shifted images ?
#     (verbose = yes)               Print the centers, shifts, and trim section ?
#        (list = )               
#        (mode = "ql")           



#
iraf.imshift.interp_type="drizzle"
iraf.imshift.interp_type="drizzle"
iraf.imshift.interp_type="drizzle"



#for z in range (1,z_base+1):
#    iraf.imalign(input=cubo_desalign+'[*,*,'+str(z)+']',reference=cubo_base+'[*,*,'+str(z)+']',coords=archivo_ref_points_base,output='align_'+cubo_desalign[:-5]+'_'+str(z)+'.fits',shifts="0.8,1",boxsize=7,bigbox=11,negative="no",background="INDEF",lower="INDEF",upper="INDEF",niterate=3,tolerance=0,maxshift="INDEF",shiftimages="yes",interp_type="spline3",boundary_typ="nearest",constant=0,trimimages="no",verbose="yes",mode="ql")
  

for z in range (1,z_base+1):
#    iraf.imshift(input=cubo_desalign+'[*,*,'+str(z)+']',output='align_'+cubo_desalign[:-5]+'_'+str(z)+'.fits',xshift=X_s,yshift=Y_s,shifts_file= "",interp_type="drizzle",boundary_typ = "nearest",constant = 0.0,mode = "ql")     

    iraf.immatch.gregister(input=cubo_desalign+'[*,*,'+str(z)+']',output='geo_align_'+cubo_desalign[:-5]+'_'+str(z)+'.fits',database=sys.argv[3],transform='INDEF',xmin = 1.0,xmax = 316.0,ymin = 1.0,ymax = 316.0, xscale = 1.0, yscale = 1.0,ncols = 316,nlines = 316,xsample = 1.0,ysample = 1.0)            
    
#creating the aligned datacube
from astropy.io import fits
import numpy as np

#cube_template_nohead=fits.getdata(cubo_desalign,header=False)
#cube_template_head=fits.getheader(cubo_base)
#cube_no_align_head=fits.getheader(cubo_desalign)


cube_template_nohead=fits.getdata(cubo_base,header=False)
cube_template_head=fits.getheader(cubo_base)
cube_no_align_head=fits.getheader(cubo_base)


##Notar que esto es Z,Y,X
print cube_template_nohead.shape
Z=cube_template_nohead.shape[0]
Y=cube_template_nohead.shape[1]
X=cube_template_nohead.shape[2]
print X,Y,Z

#print cube_template_head

cubo_null=cube_template_nohead - cube_template_nohead
print cubo_null.shape

for i in range (1,Z+1):
    slide_cubo_z=fits.getdata('geo_align_'+cubo_desalign[:-5]+'_'+str(i)+'.fits',header=False)
    #print Z,slide_cubo_z.shape
    cubo_null[i-1,:,:]= cubo_null[i-1,:,:] + slide_cubo_z  
    #print cubo_null[Z,50,50],Z,slide_cubo_z[50,50]

#print cubo_null.shape
#print cubo_null[Z,Y,X]

fits.writeto(cubo_desalign[:-5]+'geo_align.fits',cubo_null,cube_no_align_head)

#fits.update(cubo_base,"CUBO_NULL.fits",


#header_template=cube_template[0].header

#data=cube_template[0].data

#Z,Y,X=data.shape

#print Z,Y,X
##n=np.arange(int(Z),int(Y),int(X))
#n=np.arange(3681,325,323)
#hdu=fits.PrimaryHDU(n)
#hdul = fits.HDUList([hdu])
#hdul.writeto('new1.fits')
