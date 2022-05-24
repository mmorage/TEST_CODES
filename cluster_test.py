import sys,os,string

import astropy
from astropy.io import fits as pyfits
import numpy as np 
#import pyfits
import scipy
from scipy import ndimage
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import matplotlib.cm as cm
import pandas as pd
import itertools
import aplpy
from sklearn.preprocessing import normalize




Signal=250 #40 numero ideal 
limit=10000000000000 # absurd number, just in case for fine tunning

#color represnetation of the BPT: Pivotal color points

#x_green,y_green=-0.6,0.3 
#x_rojo,y_rojo  =-0.4,-0.6  
#x_blue,y_blue  =-0.17,0.3


x_green,y_green=-1.7,0.9  
x_rojo,y_rojo  =-0.5,-1 
x_blue,y_blue  =0.5,1





#emission lines


Hbeta_ini=4855    ;Hbeta_end=4866
OIII5007_ini=5001. ;OIII5007_end=5011
Ha_ini=6556       ;Ha_end=6567.
NII6584_ini=6576  ;NII6584_end=6587


#Read MUSE DATA CUBE
imagen_in=sys.argv[1] #TEST_CUBO_fit_continuo_multiplicado.fits #MUSE_NGC1427A_sum.fits  TEST_CUBO_fit_continuo_flux.fits


#extinction law
def K_lambda(l):
    K=-2.156+1.509/l-0.198/(l*l)-0.011/(l*l*l) #Calzetti 1999 
    return(K)

#slice data cube
def IFU_SLICE(imagen3D,lambda1,lambda2):
    global Signal

    header1=pyfits.getheader(imagen3D)
    imagen=(pyfits.getdata(imagen3D))

    print(imagen.shape)

    #read headers
    L_inicial=float(header1['CRVAL3']) 
    #Delta_L=float(header1['CDELT3'])
    Delta_L=float(header1['CD3_3'])
    print( L_inicial)
    print( Delta_L,"\n")
    print(  (lambda1-L_inicial)/Delta_L)
    print(  (lambda2-L_inicial)/Delta_L)
    pix1=int((lambda1-L_inicial)/Delta_L)
    pix2=int((lambda2-L_inicial)/Delta_L)
    print( "pix1=",pix1,"pix2=",pix2)

    
    
    ifu_slice=0
    for z in range(pix1+1,pix2+1):
        
        #if imagen[z,:,:]>10:
        ifu_slice=ifu_slice + imagen[z,:,:]#[0:49][0:49]
        #print(z)
    sn= ifu_slice < Signal
    ifu_slice[sn]=np.nan
    return(ifu_slice*10**(-20))



#Hbeta_ini=4885    ;Hbeta_end=4896
#OIII5007_ini=5031 ;OIII5007_end=5044
#Ha_ini=6594       ;Ha_end=6608
#NII6584_ini=6615  ;NII6584_end=6630



Hbeta=IFU_SLICE(imagen_in,Hbeta_ini,Hbeta_end)
Hbeta_K_L=K_lambda(0.5*(Hbeta_ini+Hbeta_end))
OIII=IFU_SLICE(imagen_in,OIII5007_ini,OIII5007_end)
OIII_K_L= K_lambda(0.5*(OIII5007_ini+OIII5007_end))
Ha=IFU_SLICE(imagen_in,Ha_ini,Ha_end)
Ha_K_L=K_lambda(0.5*(Ha_ini+Ha_end))
NII=IFU_SLICE(imagen_in,NII6584_ini,NII6584_end)
NII_K_L=K_lambda(0.5*(NII6584_ini+NII6584_end))


#Av= 3.1*(np.log10((Ha/Hbeta)) - np.log10(2.86) )/0.232
 
EBV=(np.log10((Ha/Hbeta)) - np.log10(2.86) )/0.232

l_y,l_x=Ha.shape
print(Ha.shape,l_y,l_x)
#exit()
archivo1=open('voro_Hbeta.txt', 'w')
archivo2=open('voro_OIII.txt', 'w')
archivo3=open('voro_Halpha.txt', 'w')
archivo4=open('voro_NII.txt', 'w')

for i in range((l_y)):
    for j in range((l_x)):
        #archivo1.write(str(int(j))+' '+str(int(i))+' '+str(Hbeta[j][i])+' '+str(np.sqrt(Hbeta[j][i]))+'\n')
        #archivo2.write(str(int(j))+' '+str(int(i))+' '+str(OIII[j][i])+' '+str(np.sqrt(OIII[j][i]))+'\n')
        #archivo3.write(str(int(j))+' '+str(int(i))+' '+str(Ha[j][i])+' '+str( np.sqrt(Ha[j][i]))+'\n')
        #archivo4.write(str(int(j))+' '+str(int(i))+' '+str(NII[j][i])+' '+str( np.sqrt(NII[j][i]))+'\n')


        archivo1.write(str(int(j))+' '+str(int(i))+' '+str(Hbeta[i][j]*10**(20))+' '+str(np.sqrt(Hbeta[i][j]*10**(20)))+'\n')
        archivo2.write(str(int(j))+' '+str(int(i))+' '+str(OIII[i][j]*10**(20))+' '+str(np.sqrt(OIII[i][j]*10**(20)))+'\n')
        archivo3.write(str(int(j))+' '+str(int(i))+' '+str(Ha[i][j]*10**(20))+' '+str( np.sqrt(Ha[i][j]*10**(20)))+'\n')
        archivo4.write(str(int(j))+' '+str(int(i))+' '+str(NII[i][j]*10**(-0))+' '+str( np.sqrt(NII[i][j]*10**(20)))+'\n')

archivo1.close()
archivo2.close()
archivo3.close()
archivo4.close()

Ha_Hb=Ha/Hbeta

cut=EBV<0  
cut2=EBV>limit
EBV[cut]=np.nan
EBV[cut2]=np.nan

Ha_corr=Ha*10**(0.4*EBV*Ha_K_L)
Hbeta_corr=Hbeta*10**(0.4*EBV*Hbeta_K_L)
NII_corr=NII*10**(0.4*EBV*NII_K_L)
OIII_corr=OIII*10**(0.4*EBV*OIII_K_L)




#Av/R= (Log(Fha/Fhb) - log 2.86) / 0.232 =E(B-V)

#OIII_BETA=(OIII/Hbeta) #*EBV_nor#*EBV
#NII_ALPHA=(NII/Ha)#*EBV_nor#*EBV,c=np.sqrt((x_bpt+1.0)**2 +(y_bpt+0.5)**2),cmap=cmap,vmin=0,vmax=2)

OIII_BETA=(OIII_corr/Hbeta_corr) 
NII_ALPHA=(NII_corr/Ha_corr)



EBV_nor=EBV/EBV
OB_nor =OIII_BETA/OIII_BETA
NII_nor=NII_ALPHA/NII_ALPHA
Ha_Hb_nor=Ha_Hb/Ha_Hb



NORM_FACTOR=OB_nor*NII_nor*EBV_nor*Ha_Hb_nor

OIII_BETA=np.log10(OIII/Hbeta)#*EBV_nor#*EBV
NII_ALPHA=np.log10(NII/Ha)#*EBV_nor#*EBV

OIII_BETA=OIII_BETA/NORM_FACTOR
NII_ALPHA=NII_ALPHA/NORM_FACTOR
EBV2=EBV/NORM_FACTOR

EBV2[cut]=np.nan
EBV2[cut2]=np.nan






bpt_x=[]
bpt_y=[]
xt=[]
yt=[]

EBV_plot_t=[]


largo_y,ancho_x=OIII_BETA.shape
print(largo_y,ancho_x)
#for i in range((largo_y)):
#    for j in range((296)):
#for i in range(50,largo_y):
#    for j in range(50,ancho_x):

for i in range(largo_y):
    for j in range(ancho_x):
        if np.isnan(NII_ALPHA[i][j]) or np.isnan(OIII_BETA[i][j]) or np.isnan(EBV2[i][j]):
            f=1
        else:
            xt.append(int(j))
            yt.append(int(i))
            bpt_x.append(NII_ALPHA[i][j])
            bpt_y.append(OIII_BETA[i][j])
            EBV_plot_t.append(EBV2[i][j])

x_bpt=np.array(bpt_x)
y_bpt=np.array(bpt_y)
X=np.array(xt)
Y=np.array(yt)
EBV_plot=np.array(EBV_plot_t)

x_bpt_max=np.nanmax(x_bpt)
y_bpt_max=np.nanmax(y_bpt)

x_bpt_min=np.nanmin(x_bpt)
y_bpt_min=np.nanmin(y_bpt)


x=np.arange(1000)*0.01-1.5#[0,1,2,3,4,5,6,7,8,9,10]
y1=0.61/(x-0.5) +1.3 #Kaufmann
y2=0.61/(x-0.47) +1.19 #kewley 2001 THEASTROPHYSICALJOURNAL, 556 : 121,140, 2001 July 20 
y3=(-30.787+1.1358*x+0.27297*x**2)*np.tanh(5.7409*x)-31.093 #Mon. Not. R. Astron. Soc.371,972-982 (2006) 

#cmap=cm.get_cmap("Spectral")
#cmap=cm.get_cmap("twilight")
#cmap=cm.get_cmap("viridis")

X=X*EBV_plot/EBV_plot
Y=Y*EBV_plot/EBV_plot
x_bpt=x_bpt*EBV_plot/EBV_plot
y_bpt=y_bpt*EBV_plot/EBV_plot

#################
#
# Calcular distancia del punto a cada eje del triangulo
# Y esto sera un gradiente de color
#
##################
#######################
#
#   2 
# 1   3
#
######################


def Distance(x_p,y_p,x_v,y_v):
    """ 
    x_p= position x for the point on of the graph
    y_p= position y for the point on of the graph
    x_v= x position of the point to which the distance is  measured.
    y_v= y position of the point to which the distance is  measured.

    """  
    D=(np.sqrt((x_p - x_v)**2 + (y_p - y_v)**2 ))
    return (D)

    
v1max=Distance(x_bpt_max,y_bpt_max,-1.5,-1)
v2max=Distance(x_bpt_max,y_bpt_max,-1.5,1)
v3max=Distance(x_bpt_max,y_bpt_max,1,0.5)

v1min=Distance(x_bpt_min,y_bpt_min,-1.5,-1)
v2min=Distance(x_bpt_min,y_bpt_min,-1.5,1)
v3min=Distance(x_bpt_min,y_bpt_min,1,0.5)


v1=Distance(x_bpt,y_bpt, -1,0.7)

v2=Distance(x_bpt,y_bpt,-1.25,1)

v3=Distance(x_bpt,y_bpt,-1,0)

#rojo=(1,0,0,1)
rojo=(((v1max-v1min))) 
rojo2=((v1/v1min)) 
print(np.nanmax(x_bpt),np.nanmin(y_bpt))
print(np.nanmax(rojo),np.nanmin(rojo))
print(np.nanmax(rojo2),np.nanmin(rojo2))

#rgv
cmap_rg=mpl.colors.LinearSegmentedColormap.from_list('my_rg',['#FF0000','#00FF00'],30000)
cmap_gb=mpl.colors.LinearSegmentedColormap.from_list('my_rv',['#00FF00','#0000FF'],30000)

cmap_rg=mpl.colors.LinearSegmentedColormap.from_list('my_rg',['r','g'],30000)
cmap_bg=mpl.colors.LinearSegmentedColormap.from_list('my_rv',['b','g'],30000)
cmap_rb=mpl.colors.LinearSegmentedColormap.from_list('my_gv',['r','g'],30000)

#1-v1/(np.nanmax(v1)-np.nanmin(v1),
#)

#cmap2=mpl.colors.LinearSegmentedColormap.from_list('my_cmap',[v1+v2+v3],30000)

c=np.sqrt((x_bpt+1.0)**2 +(y_bpt+0.5)**2)

#plt.scatter(x_bpt,y_bpt,marker='.',s=0.5,clim=(-1,1),cmap=cmap)
#plt.plot(x,y1)
#plt.plot(x,y2)
#plt.plot(x,y3)


def point_color(xp,yp,x1,y1,x2,y2,x3,y3):
    #for j in range(len(yp)):
     #   for i in range(len(xp)):
    if math.isnan(xp) or math.isnan(yp):
        color=[1,1,1]
        rt=1
        gt=1
        bt=1
        suma=1.0
    else:    
        rt=Distance(xp,yp,x1,y1)
        gt=Distance(xp,yp,x2,y2)
        bt=Distance(xp,yp,x3,y3)
        suma=rt+gt+bt
        
    rojo=rt/suma
    verde=gt/suma
    azul=bt/suma
    color=np.array([1-rojo,1-verde,1-azul])
    return(color)

#color1=Distance(x_bpt,y_bpt,-1.07,-0.94)
#color2=Distance(x_bpt,y_bpt,-1.2,0.6)
#color3=Distance(x_bpt,y_bpt,-0.05,0.6)

print (Distance(-4,-2,-4,1.5))
print (Distance(-4,-2,-1,1.6))
#exit()

#origen_x
x_o=0.000000001
y_o=0.00001

x_og=3.000001
y_og=0.00001

x_ob=-1.000001
y_ob=0.500001


###finales####
#x_rojo,y_rojo  =-1,-1
#x_green,y_green=1,2
#x_blue,y_blue  =-1,-2.5
##x_blue,y_blue  =-0.2,0.7
############################
#x_rojo,y_rojo  =-0.9,-0.2
##0.8,1.5
##x_green,y_green=-1.3,0.6  #-1.5,1.5
#x_green,y_green=-0.8,0.25 #-0.4,0.7  #-1.5,1.5
##x_blue,y_blue  =-0.7,0.7  #0.8,-1.5
#x_blue,y_blue  =-0.5,-0.1  #0.8,-1.5
# 
##x_blue,y_blue  =-0.2,0.7



Yr=y_o + ((y_rojo-y_o)/(x_rojo-x_o))*(x_bpt-x_o)
Yg=y_og + ((y_green-y_og)/(x_green-x_og))*(x_bpt-x_og)
Yb=y_ob + ((y_blue-y_ob)/(x_blue-x_ob))*(x_bpt-x_ob)
################ FUNCIONA!!!###########

color1=((Distance(x_bpt,y_bpt,x_blue,y_green))/Distance(x_bpt,y_bpt, x_rojo,y_rojo))#*250                           #,Distance(-4,-2,-4,1.5),Distance(-4,-2,-1,1.6)
color2=((Distance(x_bpt,y_bpt,x_blue,y_rojo))/Distance(x_bpt,y_bpt,x_green,y_green))#*250
color3=((Distance(x_bpt,y_bpt,x_green,y_rojo))/Distance(x_bpt,y_bpt,x_blue,y_blue))#*250

#color1=(Distance(x_bpt,y_bpt,-0, 1  )/Distance(x_bpt,y_bpt,x_rojo,y_rojo  ))                           #,Distance(-4,-2,-4,1.5),Distance(-4,-2,-1,1.6)
#color2=np.arcsinh(Distance(x_bpt,y_bpt,-0.9,0.2   )/Distance(x_bpt,y_bpt,x_green,y_green))
#color3=np.arcsinh(Distance(x_bpt,y_bpt,-0.9,0.2  )/Distance(x_bpt,y_bpt,x_blue,y_blue  ))      


#color1=(Distance(x_bpt,y_bpt,x_rojo,x_green)/Distance(x_bpt,y_bpt,x_rojo,y_rojo  ))                           #,Distance(-4,-2,-4,1.5),Distance(-4,-2,-1,1.6)
#color2=(Distance(x_bpt,y_bpt,x_green,y_blue   )/Distance(x_bpt,y_bpt,x_green,y_green))
#color3=(Distance(x_bpt,y_bpt,x_blue,y_rojo   )/Distance(x_bpt,y_bpt,x_blue,y_blue  ))      

#color1=color1**(4/2)
#color2=np.sinh(color2**(1/2))
#color3=np.arcsinh(color3**(4/2))

color1=1.1*color1#**(4/2)
color2=1.05*(color2)#**(1/2))
##color3=1+np.sinh(color3)

suma_color=(color1+color2+color3)#/(color1+color2+color3)



#c=np.append(color1,color2,color3,axis=1) 
col=np.stack([color1/suma_color,color2/suma_color,color3/suma_color,color2/color2]).T
#print =(color1)
#exit()

#hex_col=rgb2hex(int(col[0]*255),int(col[1]*255),int(col[3][2]*255))
#exit()
#hex_col=rgb2hex(20,190,200)
#hex_col=rgb2hex(col[0],col[1],col[2])



#col=np.stack([color1,color2/color2,color3/color3,color2/color2]).T
#bokeh_colors = ["#%02x%02x%02x" % (r, g, b) for r, g, b in A_color[:,0:3]]

#col=np.stack([color1,color2,color3,color2/color2]).T


#col=np.stack([color1,color1/color2,color1/color3,color2/color2]).T

print(c.shape)

#plt.scatter(x_bpt,color1)
#plt.scatter(x_bpt,color2)
#plt.scatter(x_bpt,color3)
#plt.scatter(x_bpt,1-color2/(color3+color2+color1))


xk=np.arange(-10,0.0,0.01)#[0,1,2,3,4,5,6,7,8,9,10]
y1=0.61/(xk-0.05) +1.3 #Kaufmann
x=np.arange(-10,0.4,0.01)#[0,1,2,3,4,5,6,7,8,9,10]
y2=0.61/(x-0.47) +1.19 #kewley 2001 THEASTROPHYSICALJOURNAL, 556 : 121,140, 2001 July 20 
y3=(-30.787+1.1358*x+0.27297*x**2)*np.tanh(5.7409*x)-31.093 #Mon. Not. R. Astron. Soc.371,972-982 (2006) 


y4=np.arange(-10,20,0.01)#[0,1,2,3,4,5,6,7,8,9,10]

x4=-0.596*y4*y4 -0.687*y4-0.655 #Brinchnamm, Pettini and Charlot  MNRAS 385,768 2008

from bokeh.io import output_file, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, widgets
from bokeh.plotting import figure
from bokeh.models.glyphs import Circle
from bokeh.palettes import d3
from bokeh.palettes import brewer
from bokeh.models import Range1d

#import holoviews as hv
#from holoviews import opts




TOOLS="hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"
 #print(colors)
print(len(x_bpt))
print(len(255*color2/suma_color))
c1=np.intc(255*color1/suma_color)
c2=np.intc(255*color2/suma_color)
c3=np.intc(255*color3/suma_color)
print(type(c1))

output_file("INTERACTIVE_BPT_restframe.html")
colorea = [
    "#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b  in zip(c1,c2,c3)
]

 
xk=np.arange(-10,0.0,0.01)#[0,1,2,3,4,5,6,7,8,9,10]
y1=0.61/(xk-0.05) +1.3 #Kaufmann

x=np.arange(-10,0.4,0.01)#[0,1,2,3,4,5,6,7,8,9,10]
y2=0.61/(x-0.47) +1.19 #kewley 2001 THEASTROPHYSICALJOURNAL, 556 : 121,140, 2001 July 20 
y3=(-30.787+1.1358*x+0.27297*x**2)*np.tanh(5.7409*x)-31.093 #Mon. Not. R. Astron. Soc.371,972-982 (2006) 


y4=np.arange(-10,20,0.01)#[0,1,2,3,4,5,6,7,8,9,10]
x4=-0.596*y4*y4 -0.687*y4-0.655 #Brinchnamm, Pettini and Charlot  MNRAS 385,768 2008




data_lineas1={'xk':xk,'y1':y1}
data_lineas2={'x':x,'y2':y2,'y3':y3}
data_lineas3={'x4':x4,'y4':y4}
data= {'x':x_bpt,'y':y_bpt, 'x1':X,'y1':Y,'rc': c1, 'gc':c2,'bc':c3,'colores':colorea }
Data= {'xbpt':x_bpt,'ybpt':y_bpt, 'x1':X,'y1':Y,'rc': c1, 'gc':c2,'bc':c3,'colores':colorea } 



from sklearn.neighbors import KernelDensity
from matplotlib.colors import LogNorm

from astroML.density_estimation import KNeighborsDensity
from sklearn.mixture import GaussianMixture

from scipy import stats

from pandas import DataFrame
from sklearn.cluster import KMeans

import matplotlib.colors as  col_conv

import pyfof



df=DataFrame(Data,columns=['rc','gc','bc'])
dfb = DataFrame(Data,columns=['xbpt','ybpt','zbpt','x1','y1'])



print(pd.DataFrame(df).head())

import hdbscan

#clusterer = hdbscan.HDBSCAN(min_cluster_size=2,cluster_selection_method='leaf').fit(df)
clusterer = hdbscan.HDBSCAN().fit(df)
#a=clusterer.fit(df)

print(clusterer.labels_)
print(clusterer.labels_.max())
print(clusterer.probabilities_)


print(clusterer.labels_.shape)
print(clusterer.labels_.max().shape)
print(clusterer.probabilities_.shape)




#pal = sns.color_palette('deep', 12)
#colors = [sns.desaturate(pal[col], sat) for col, sat in zip(clusterer.labels_,
cm = plt.cm.get_cmap('RdYlBu')
        #                            clusterer.probabilities_)]

results={'xbpt':x_bpt,'ybpt':y_bpt, 'x1':X,'y1':Y,'rc': c1, 'gc':c2,'bc':c3,'colores':colorea,'Nc':clusterer.labels_,'Prob': clusterer.probabilities_} 
import statistics
from statistics import mode


dfrr=DataFrame(results,columns=['xbpt','ybpt', 'x1','y1','rc', 'gc','bc','colores','Nc','Prob'])
dfr=dfrr.loc[dfrr['Nc']>-1]
#print(mode(dfr['Nc']))
print("lala")

#plt.scatter(dfrr['Nc'],dfrr['Prob']) 

#plt.hist(dfrr['Nc'],bins=60)
#plt.show()
#dfrr.sort_values(by=['Nc']) 
#dfr=dfrr.loc[dfrr['Prob']==1]
#pal = sns.color_palette('deep', 12)
#colors = [sns.desaturate(pal[col], sat) for col, sat in zip(clusterer.labels_,
cm = plt.cm.get_cmap('RdYlBu')
        #                            clusterer.probabilities_)]


        
#for index,row in dfr.iterrows():
print(dfr['Nc'])
#plt.scatter(dfr['xbpt'],dfr['ybpt'],c=dfr['Nc'],s=1.5,cmap='hsv')# c=colors, **plot_kwds);
plt.scatter(dfr['xbpt'],dfr['ybpt'],c=dfr['Nc'],s=1.5,cmap='RdYlBu')# c=colors, **plot_kwds);
#plt.scatter(row['xbpt'],row['ybpt'],c=row['Nc'],s=1.5)#,cmap='RdYlBu')# c=colors, **plot_kwds);
#plt.colorbar()
#plt.pause(0.05)
plt.show()

plt.scatter(dfr['x1'].T,dfr['y1'].T,c=dfr['Nc'].T,s=1.5,cmap='RdYlBu')# c=colors, **plot_kwds);
plt.colorbar()
plt.show()


#plt.scatter(dfb['x1'].T,dfb['y1'].T,c=clusterer.probabilities_,s=1.5,cmap='RdYlBu')# c=colors, **plot_kwds);
#plt.colorbar()
#plt.show()
 
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')
ax.scatter(dfr['x1'].T,dfr['y1'].T,dfr['Nc'].T,c=dfr['Nc'],marker='.',cmap='RdYlBu')
plt.show()




fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')
ax.scatter(dfr['xbpt'].T,dfr['ybpt'].T,dfr['Nc'].T,c=dfr['Nc'],marker='.',cmap='RdYlBu')
plt.show()



fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')
ax.scatter(dfr['rc'].T,dfr['bc'].T,dfr['gc'].T,c=dfr['Nc'],marker='.',cmap='RdYlBu')
plt.show()
