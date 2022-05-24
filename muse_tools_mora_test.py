import astropy.io.fits as pyfits
import sys,os,string
import numpy as np


def MUSE_delta_pix(cubo):
    header=pyfits.getheader(cubo)
    return(header['CD3_3'])
##################################
def MUSE_lamda_inicial(cubo):
    header=pyfits.getheader(cubo)
    return(header['CRVAL3'])

################################


def MUSE_SIZE_xyz(cubo):
    """muestra tama~No del CUBO  x,y,z
    
    """
    imagen=pyfits.getdata(cubo,header=False)
    z_dim,y_dim,x_dim=imagen.shape
    return(x_dim,y_dim,z_dim)

##############################################

def MUSE_READ(cubo):
    """ lee MUSE CUBO (solo el DAT) y lo deja en memoria y  header
    """
    imagen=pyfits.getdata(cubo,header=False)
    header=pyfits.getheader(cubo)
    return(imagen,header)
 
#######################################

def MUSE_lee_spectra(cubo,x,y):
    """ LEE un spaxel como spectra
    """
    imagen=pyfits.getdata(cubo,header=False)
    header=pyfits.getheader(cubo)
    z_dim,y_dim,x_dim=imagen.shape

    LAMBDA=np.arange(1,z_dim+1)*header['CD3_3']+header['CRVAL3']
    FLUX=imagen[:,y,x]
    return(LAMBDA,FLUX)

##########################################

def MUSE_lee_slide(cubo,pixel):
    """ LEE un slide del cubo
    """
    imagen=pyfits.getdata(cubo,header=False)
    header=pyfits.getheader(cubo)
    z_dim,y_dim,x_dim=imagen.shape
    data_cube=imagen[pixel-1,:,:]
    return(data_cube)
    
########################################

def MUSE_to_xy_lambda_flux(cubo):
    """Lee cubo (DAT) y entrega array de lambda y cubo (cuentas,eje y, eje x)
    """
    imagen=pyfits.getdata(cubo,header=False)
    header=pyfits.getheader(cubo)

    z_dim,y_dim,x_dim=imagen.shape

    LAMBDA=np.arange(1,z_dim+1)*header['CD3_3']+header['CRVAL3']
    return(LAMBDA,imagen)

#############################################

def MUSE_lambda(cubo):

    imagen=pyfits.getdata(cubo,header=False)
    header=pyfits.getheader(cubo)

    z_dim,y_dim,x_dim=imagen.shape

    LAMBDA=np.arange(1,z_dim+1)*header['CD3_3']+header['CRVAL3']


    return(LAMBDA)

#############################################

def array_map(i,j):
    """ Creates an x,y array for map function and thus avoid loops
    use: 
    p,q=array_map(y_dim,x_dim)
    """
    tmp= [(x, y) for x in range(int(j)) for y in range(int(i)) ]
    p,q=np.asarray(tmp).reshape(-1,2).T
    #p,q=np.asarray(tmp).reshape(-1,2).T
    
    return(p,q)
