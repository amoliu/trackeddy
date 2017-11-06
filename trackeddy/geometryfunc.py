import numpy as np
import pylab as plt
from scipy.interpolate import interp2d,interp1d
from scipy.optimize import curve_fit,leastsq


def fit_ellipse(x,y,diagnostics=False):
    '''
    This function was translated form Matlab to python by Josue Martinez Moreno,
    the original source:
    Copyright (c) 2003, Ohad Gal 
    All rights reserved.

    Redistribution and use in source and binary forms, with or without 
    modification, are permitted provided that the following conditions are 
    met:

    * Redistributions of source code must retain the above copyright 
    notice, this list of conditions and the following disclaimer. 
    * Redistributions in binary form must reproduce the above copyright 
    notice, this list of conditions and the following disclaimer in 
    the documentation and/or other materials provided with the distribution
    For more information go to the main source:
    https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse?requestedDomain=www.mathworks.com
    '''
    orientation_tolerance = 1e-3;
    x=x[:]
    y=y[:]
    mean_x=np.mean(x)
    mean_y=np.mean(y)
    xp=x-mean_x
    yp=y-mean_y
    X=np.array([xp**2,xp*yp,yp**2,xp,yp]).T
    
    a=np.sum(X,axis=0)
    b=np.dot(X.T,X)

    # xb = a: solve b.T x.T = a.T instead 
    x2 = np.linalg.lstsq(b.T, a.T)
    
    x2 = np.dot(a, np.linalg.pinv(b)) 
    a,b,c,d,e=x2
    if ( min(abs(b/a),abs(b/c)) > orientation_tolerance ):
        orientation_rad = 1/2 * np.arctan( b/(c-a) )
        cos_phi = np.cos( orientation_rad )
        sin_phi = np.sin( orientation_rad )
        a,b,c,d,e = [a*cos_phi**2 - b*cos_phi*sin_phi + c*sin_phi**2,0,a*sin_phi**2 + b*cos_phi*sin_phi + \
                     c*cos_phi**2,d*cos_phi - e*sin_phi,d*sin_phi + e*cos_phi]
        mean_x,mean_y=cos_phi*mean_x - sin_phi*mean_y,sin_phi*mean_x + cos_phi*mean_y       
    else:
        orientation_rad = 0;
        cos_phi = np.cos( orientation_rad );
        sin_phi = np.sin( orientation_rad );
    test = a*c;
    if test>0:
        detect='Ellipse'
        status=True
    elif test==0:
        detect='Parabola'
        status=False
    else:
        detect='Hyperbola'
        status=False
    if status==True:
        # make sure coefficients are positive as required
        if (a<0):
            a=-a
            c=-c
            d=-d
            e=-e
        
        # final ellipse parameters
        X0          = mean_x - (d/2)/a;
        Y0          = mean_y - (e/2)/c;
        F           = 1 + (d**2)/(4*a) + (e**2)/(4*c);
        a,b       =  np.sqrt( F/a ),np.sqrt( F/c );    
        long_axis   = 2*max(a,b);
        short_axis  = 2*min(a,b);

        # rotate the axes backwards to find the center point of the original TILTED ellipse
        R           = [[ cos_phi,sin_phi],[-sin_phi,cos_phi ]]
        P_in        = np.dot(R, np.array([X0,Y0]))
        X0_in       = P_in[0]
        Y0_in       = P_in[1]
        
        ver_line        = np.array([ [X0,X0], [Y0-1*b, Y0+1*b]])
        horz_line       = np.array([ [X0-1*a,X0+1*a], [Y0,Y0] ])
        new_ver_line    = np.dot(R,ver_line)
        new_horz_line   = np.dot(R,horz_line)
    
        # the ellipse

        theta_r         = np.linspace(0,2*np.pi);
        ellipse_x_r     = X0 + a*np.cos(theta_r)
        ellipse_y_r     = Y0 + b*np.sin(theta_r)
        rotated_ellipse =  np.dot(R, np.array([ellipse_x_r,ellipse_y_r]))
        
        # pack ellipse into a structure
        ellipse_t = {'a':a,'b':b,'phi':orientation_rad,'X0':X0,'Y0':Y0,\
                     'X0_in':X0_in,'Y0_in':Y0_in,'long_axis':long_axis,\
                     'short_axis':short_axis,'minoraxis':new_horz_line,\
                     'majoraxis':new_ver_line,'ellipse':rotated_ellipse\
                     ,'status':'Cool'}
    else:
        # report an empty structure
        ellipse_t = {'a':'','b':'','phi':'','X0':'','Y0':'',\
                     'X0_in':'','Y0_in':'','long_axis':'',\
                     'short_axis':'','minoraxis':'',\
                     'majoraxis':'','ellipse':'',\
                     'status':detect}
    if diagnostics==True and status==True:
        # draw
        plt.plot( x,y,'b' );
        plt.plot( new_ver_line[0],new_ver_line[1],'r' )
        plt.plot( new_horz_line[0],new_horz_line[1],'r' )
        plt.plot( rotated_ellipse[0],rotated_ellipse[1],'r' )
        plt.show()
    return ellipse_t,status

def PolyArea(x,y):
    '''
    *************** Poligon Area *******************
    This function calculaate the area of a poligon.
    Usage:
    
    Example:
    '''
    area=0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
    return area

def eccentricity(a,b):
    '''
    *************** Ellipse Eccentricity *******************
    This function calculate the eccentricity of a ellipse.
    Usage:
    
    Example:
    a=minoraxis
    b=mayoraxis
    eccen=eccentricity(a,b)
    '''
    if b>a:
        b1=a
        a=b
        b=b1
    eccen=np.sqrt(1-(b**2/a**2))
    return eccen

def find2d(arrayx,arrayy,valuex,valuey):
    '''
    *************** Find 2 Value in 2 list *******************
    This function find values in a list of values.
    Usage:
    
    Example:
    arrayx=list
    arrayx=list
    valuex=value2search
    valuey=value2search
    indexes=find2d(arrayx,arrayy,valuex,valuey)
    '''
    idx=(np.abs(arrayx-valuex)).argmin()
    idy=(np.abs(arrayy-valuey)).argmin()
    return idx,idy

def find(array,value):
    '''
    *************** Find Value in list *******************
    This function find values in a list of values.
    Usage:
    
    Example:
    arrayx=list
    arrayx=list
    valuex=value2search
    valuey=value2search
    indexes=find2d(arrayx,arrayy,valuex,valuey)
    '''
    idx=(np.abs(array-value)).argmin()
    return idx

def contourmaxvalue(contcoordx,contcoordy,var,x,y,levels,date):
    '''
    *************** Find maximum value inside contour *******************
    This function find the maximum value inside an specific contour.
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    idxcheckmax,idycheckmax=find2d(x,y,contcoordx.max(),contcoordy.max())
    idxcheckmin,idycheckmin=find2d(x,y,contcoordx.min(),contcoordy.min())
    #print(idycheckmin,idycheckmax,idxcheckmin,idxcheckmax)
    if len(np.shape(var))==3:
        if levels[0]>0:
            var[var>levels[0]]==np.nan
            sshextrem=np.nanmax(var[date,idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        else:
            var[var<levels[0]]==np.nan
            sshextrem=np.nanmin(var[date,idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        indexes=np.where(var[date,idycheckmin:idycheckmax,idxcheckmin:idxcheckmax]==sshextrem)
    else:
        #print(np.shape(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax]))
        if levels[0]>0:
            var[var>levels[0]]==np.nan
            sshextrem=np.nanmax(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        else:
            var[var<levels[0]]==np.nan
            sshextrem=np.nanmin(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        #print(sshextrem)
        indexes=np.where(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax]==sshextrem)
    coord=[x[idxcheckmin+indexes[1][0]],y[idycheckmin+indexes[0][0]]]
    return coord

def centroidvalue(contcoordx,contcoordy,var,x,y,levels,date):
    '''
    *************** Find centroid inside contour *******************
    This function find the centroid inside an specific contour.
    Usage:
    
    Example:
    
    Author: Josue Martinez Moreno, 2017
    '''
    idxcheckmax,idycheckmax=find2d(x,y,contcoordx.max(),contcoordy.max())
    idxcheckmin,idycheckmin=find2d(x,y,contcoordx.min(),contcoordy.min())
    #print(idycheckmin,idycheckmax,idxcheckmin,idxcheckmax)
    if len(np.shape(var))==3:
        if levels[0]>0:
            var[var>levels[0]]==np.nan
        else:
            var[var<levels[0]]==np.nan
        sum_T=np.sum(var[date,idycheckmin:idycheckmax+1,idxcheckmin:idxcheckmax+1])
        sum_X=np.sum(var[date,idycheckmin:idycheckmax+1,idxcheckmin:idxcheckmax+1],axis=0)
        sum_Y=np.sum(var[date,idycheckmin:idycheckmax+1,idxcheckmin:idxcheckmax+1],axis=1)
        XM=0
        for ii in range(len(sum_X)):
            XM=XM+sum_X[ii]*x[idxcheckmin+ii]
        YM=0
        for ii in range(len(sum_Y)):
            YM=YM+sum_Y[ii]*y[idycheckmin+ii]
        xcpos=XM/sum_T
        ycpos=YM/sum_T
    else:
        if levels[0]>0:
            var[var>levels[0]]==np.nan
        else:
            var[var<levels[0]]==np.nan
        sum_T=np.sum(var[idycheckmin:idycheckmax+1,idxcheckmin:idxcheckmax+1])
        sum_X=np.sum(var[idycheckmin:idycheckmax+1,idxcheckmin:idxcheckmax+1],axis=0)
        sum_Y=np.sum(var[idycheckmin:idycheckmax+1,idxcheckmin:idxcheckmax+1],axis=1)
        XM=0
        for ii in range(len(sum_X)):
            XM=XM+sum_X[ii]*x[idxcheckmin+ii]
        YM=0
        for ii in range(len(sum_Y)):
            YM=YM+sum_Y[ii]*y[idycheckmin+ii]
        xcpos=XM/sum_T
        ycpos=YM/sum_T
    coord=[xcpos,ycpos]
    return coord


def gaus(x,a,x0,sigma):
    '''
    *************** Gaussian function *******************
    This function build a gausian curve.
    Usage:
    
    Example:

    '''
    gauss=a*np.exp(-(x-x0)**2/(2*sigma**2))
    return gauss

def adjust1Gaus(x,y):
    '''
    *************** One Gaussian Fit *******************
    This function fit one gaussian in a curve curve.
    Usage:
    
    Example:

    '''
    gauss_fit = lambda p, x: p[0]*(1/np.sqrt(2*np.pi*(p[2]**2)))*np.exp(-(x-p[1])**2/(2*p[2]**2)) #1d Gaussian func
    e_gauss_fit = lambda p, x, y: (gauss_fit(p,x) -y) #1d Gaussian fit

    v0= [1,10,1,1,30,1] #inital guesses for Gaussian Fit. - just do it around the peaks
    out = leastsq(e_gauss_fit, v0[:], args=(x, y), maxfev=1000, full_output=1) #Gauss Fit
    v = out[0] #fit parameters out
    covar = out[1] #covariance matrix output

    gausfit = gauss_fit(v,x) 
    return gausfit

def adjustMGaus(x,y):
    '''
    *************** Multiple Gaussian Fit *******************
    This function fit multiple gaussian in a curve curve.
    Usage:
    
    Example:

    '''
    gauss_fit = lambda p, x: p[0]*(1/np.sqrt(2*np.pi*(p[2]**2)))*np.exp(-(x-p[1])**2/(2*p[2]**2))+\
            p[3]*(1/np.sqrt(2*np.pi*(p[5]**2)))*np.exp(-(x-p[4])**2/(2*p[5]**2)) #1d Gaussian func
    e_gauss_fit = lambda p, x, y: (gauss_fit(p,x) -y) #1d Gaussian fit

    n=len(x)
    #v0=range(0,n,int(n/10))
    v0=[1,int(n/3),1,1,int(n/2),1,1,2*int(n/3),1]
    #v0= [1,10,1,1,30,1] #inital guesses for Gaussian Fit. - just do it around the peaks
    out = leastsq(e_gauss_fit, v0[:], args=(x, y), maxfev=1000, full_output=1) #Gauss Fit
    v = out[0] #fit parameters out
    covar = out[1] #covariance matrix output

    gausfit = gauss_fit(v,x) # this will only work if the units are pixel and not wavelength
    return gausfit

def rsquard(y,y1):
    '''
    *************** Correlation Coefficient (R^2) *******************
    This function calculate the Pearson Coefficient.
    Usage:
    
    Example:

    '''
    #yhat=y1            # or [p(z) for z in x]
    #ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    #ssreg = np.sum((y-yhat)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    #sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    #R2 = 1 - ssreg / sstot
    R2=np.corrcoef(y, y1)[0,1]
    #print(numpadj[0,1],R2,' geometry file')
    return R2 

def ellipsoidfit(y,yfit,ellipsrsquarefit=0.85,diagnostics=False):
    '''
    *************** Ellipsoid  Fit *******************
    This function check the fitness of an ellipsoid in a curve.
    Usage:
    
    Example:
    
    Author: Josue Martinez Moreno, 2017
    '''
    x=range(0,len(y))
    f=interp1d(x, y)
    xnew=np.linspace(0,len(y)-1,len(yfit))
    eddy2fit=f(xnew)
    
    indxed=find(yfit,yfit.max())
    indxrd=find(eddy2fit,eddy2fit.max())
    
    if yfit[indxed]==yfit[indxed+1] and yfit[indxed]==yfit[indxed-1]:
        indxed=find(yfit,yfit.min())
        indxrd=find(eddy2fit,eddy2fit.min())
    
    eddy2fit=list(eddy2fit)*2
    eddyfitdisplace=np.zeros(len(yfit))
    for ii in range(len(yfit)):
        eddyfitdisplace[ii]=eddy2fit[indxrd-indxed+ii]
    Rsquard=rsquard(eddyfitdisplace,yfit)
    if diagnostics==True:
        #plt.figure()
        #plt.plot(eddy['ellipse'][0][0],eddy['ellipse'][0][1])
        #plt.plot(eddy['contour'][0][0],eddy['contour'][0][1])
        plt.figure()
        #plt.plot(eddy['ellipse'][0][1])
        plt.plot(yfit)
        plt.plot(eddyfitdisplace)
        plt.text(0, np.mean(yfit), str(round(Rsquard,2)))
        plt.show()
    if Rsquard>=ellipsrsquarefit:
        check=True
    else:
        check=False
    return Rsquard,check
    
def extractprofeddy(axis,field,lon,lat,n,gaus='One',kind='linear',gaussrsquarefit=0.65,varname='',diagnostics=False):
    '''
    *************** Extract Profile Eddy *******************
    This function extracts the mayoraxis and minoraxis profile,
    inside an eddy.
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    try:
        fieldnan=field.filled(0)
    except:
        field[~np.isfinite(field)]=0
        fieldnan=field
    #fieldnan=field
    #print(nanmax(field),nanmin(field))
    y=np.linspace(axis[1,0],axis[1,1],n)
    x=np.linspace(axis[0,0],axis[0,1],n)
    
    field2interp=interp2d(lon, lat, fieldnan[:,:], kind=kind)
    field_interp = field2interp(x,y)
    #pcolormesh(lon,lat,field,vmin=-0.000004,vmax=0.000004)

    axisdata=np.zeros([n])
    for ii in range(n):
        axisdata[ii]=field_interp[ii,ii]

    n = len(axisdata)  #the number of data
    x = np.array(range(n))
    y = axisdata
    
    if gaus=='None':
        Rsquared=1
    else:
        if gaus=='One':
            gausfit=adjust1Gaus(x,y)
        elif gaus=='Multiple':
            gausfit=adjustMGaus(x,y)
        else:
            print('Select a gausian method to adjust.')
            return
        
        Rsquared = rsquard(y,gausfit)
        
    if Rsquared >= gaussrsquarefit:
        check=True
    else:
        check=False
        
    if diagnostics==True:
        print('std',varname,' vs fit',Rsquared)
        plt.plot(x,y,'b+:',label='data')
        plt.plot(x,gausfit,'ro:',label='fit')
        plt.legend()
        plt.title('Fit for Time Constant')
        plt.xlabel('Position (n)')
        plt.ylabel(varname)
        plt.show()
    return axisdata,check