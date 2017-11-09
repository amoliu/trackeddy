import numpy as np
import pylab as plt
from scipy.interpolate import interp2d,interp1d
from scipy.optimize import curve_fit,leastsq


def fit_ellipse(x,y,diagnostics=False):
    '''
    **************** fit_ellipse *****************
    Fitting of an ellipse to an array of positions.
    
    Function translated form Matlab to python by Josue Martinez Moreno,
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
    Notes:
    
    Args:
        x,y (array): Coordinates of the datapoints to fit an ellipse.
        diagnostics (boolean): Used to display all the statistics and plots to identify bugs. 
    Returns:
        ellipse_t (dict) - This dictionary contains useful parameters describing completly the ellipsoid ajusted.
        status (boolean) - This value will be true if and only if the the fit corresponds to a ellipse.
    Usage:
    R = np.arange(0,2*pi, 0.01)
    x = 1.5*np.cos(R) + 2 + 0.1*np.random.rand(len(R))
    y = np.sin(R) + 1. + 0.1*np.random.rand(len(R))
    ellipse,status=fit_ellipse(x,y,diagnostics=False)
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
    *************** PolyArea *******************
    Calculate the area of a poligon.
    Notes:
    
    Args:
        x,y (array): Coordinates of the datapoints to fit an ellipse.
    Returns:
        area (float) -  Area contained by the poligon.
    Usage:
        R = np.arange(0,2*pi, 0.01)
        x = 1.5*np.cos(R) + 2 + 0.1*np.random.rand(len(R))
        y = np.sin(R) + 1. + 0.1*np.random.rand(len(R))
        area=PolyArea(x,y)
    '''
    area=0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
    return area

def eccentricity(a,b):
    '''
    *************** eccentricity *******************
    This function calculate the eccentricity of a ellipse.
    Notes:
        
    Args:
        a (float): Mayor axis of an ellipse
        b (float): Minor axis of an ellipse
    Returns:
        eccen (float) - Eccentricity of the ellipsoid with parameters a,b.
    Usage:
        a=0.5
        b=0.3
        eccen=eccentricity(a,b)
    '''
    if b>a:
        b1=a
        a=b
        b=b1
    eccen=np.sqrt(1-(b**2/a**2))
    return eccen

def find2l(arrayx,arrayy,valuex,valuey):
    '''
    *************** find2l *******************
    Find values in two list of values.
    Notes:
        
    Args:
        arrayx (list|array): Array where it will look the closer index to the valuex
        arrayy (list|array): Array where it will look the closer index to the valuey
        valuex (int|float): Value to look for in arrayx.
        valuey (int|float): Value to look for in arrayy.
    Returns:
        idx,idy (int) - Index of closer values.
    Usage:
        arrayx=[0,1,2,3]
        arrayy=[4,5,6,7]
        valuex=2.2
        valuey=6.6
        indexes=find2l(arrayx,arrayy,valuex,valuey)
    '''
    idx=(np.abs(arrayx-valuex)).argmin()
    idy=(np.abs(arrayy-valuey)).argmin()
    return idx,idy

def find(array,value):
    '''
    *************** find *******************
    Find values in a list of values.
    Notes:
        
    Args:
        array (list|array): Array where it will look the closer index to the value
        value (int|float): Value to look for in array.
    Returns:
        idx - Index of closer values.
    Usage:
        array=[0,1,2,3]
        value=2.2
        idx=find(array,value)
    '''
    idx=(np.abs(array-value)).argmin()
    return idx

def find2D(array,value):
    '''
    *************** find2D *******************
    Find values in a 2D array of values.
    Notes:
        
    Args:
        array (list|array): 2D Array where it will look the closer index to the value
        value (int|float): Value to look for in array.
    Returns:
        yp,xp - Index of closer values.
    Usage:
        array=[[0,1,2,3],[0,1,4,3]]
        value=2
        idx,idy=find2D(array,value)
    '''
    yp,xp=np.where(array==value)
    return yp,xp
    
def contourmaxvalue(contcoordx,contcoordy,var,x,y,levels,date=''):
    '''
    *************** contourmaxvalue *******************
    Find the maximum value inside an specific contour.
    Notes:
        
    Args:
        contcoordx (list|array): Contains the coordinates in X of the contour of the field var. 
        contcoordy (list|array): Contains the coordinates in Y of the contour of the field var.
        var (array): 3D Matrix representing a surface (np.shape(var)=(date,len(x),len(y))).
        x (list|arrat): Contains the coordinate X of the grid of var.
        y (list|arrat): Contains the coordinate Y of the grid of var.
        levels (list): Level of the extracted contour.
        date (int): Used if len(var)==3 (i.e. Var contains a time dimension).
    Returns:
        coord (list) - Location of the max value in the grid.
    Usage:
        center_eddy=contourmaxvalue(contcoordx,contcoordx,sshnan,lon,lat,levels,date)
    '''
    idxcheckmax,idycheckmax=find2d(x,y,contcoordx.max(),contcoordy.max())
    idxcheckmin,idycheckmin=find2d(x,y,contcoordx.min(),contcoordy.min())
    #print(idycheckmin,idycheckmax,idxcheckmin,idxcheckmax)
    if len(np.shape(var))==3 or date='':
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
    *************** centroidvalue *******************
    Find the centroid inside an specific contour.
    Notes:
        
    Args:
        contcoordx (list|array): Contains the coordinates in X of the contour of the field var. 
        contcoordy (list|array): Contains the coordinates in Y of the contour of the field var.
        var (array): 3D Matrix representing a surface (np.shape(var)=(date,len(x),len(y))).
        x (list|arrat): Contains the coordinate X of the grid of var.
        y (list|arrat): Contains the coordinate Y of the grid of var.
        levels (list): Level of the extracted contour.
        date (int): Used if len(var)==3 (i.e. Var contains a time dimension).
    Returns:
        coord (list) - Location of the centroid in the grid.
    Usage:
        center_eddy=centroidvalue(contcoordx,contcoordx,sshnan,lon,lat,levels,date)
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
    *************** gaus *******************
    Build a gausian curve.
    Notes:
    
    Args:
        x (list|array): Array of positions.
        a (float): Amplitud of gaussian.
        x0 (float): Center of Gausian.
        sigma (float): Deviation.
    Returns:
        gauss (array) - Array of gaussian values.
    Usage:
        x=np.arange(-5,5,0.1)
        x0=0
        a=3
        sigma=2
        gaussian=gaus(x,a,x0,sigma)
        plot(x,gaussian)
        show()
    '''
    gauss=a*np.exp(-(x-x0)**2/(2*sigma**2))
    return gauss

def adjust1Gaus(x,y):
    '''
    *************** adjust1Gaus *******************
    Fit one gaussian in a curve curve.
    Notes:
    
    Args:
        x(list|array): Coordinates in x of data.
        y(list|array): Data to be ajusted with one gaussian. 
    Returns:
        gausfit(list|array) - Data ajusted.
    Usage:
        x=np.arange(-5,5,0.1)
        x0=0
        a=3
        sigma=2
        gaussian=gaus(x,a,x0,sigma)
        gaussianfit=adjust1Gaus(x,gaussian)
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
    *************** adjustMGaus *******************
    Fit multiple gaussian in a curve curve.
    Notes:
        
    Args:
        x(list|array): Coordinates in x of data.
        y(list|array): Data to be ajusted with multiple gaussians. 
    Returns:
        gausfit(list|array) - Data ajusted.
    Usage:
        x=np.arange(-5,5,0.1)
        x0=0
        a=3
        sigma=2
        gaussian=gaus(x,a,x0,sigma)+gaus(x,a-2,x0+2,sigma-1)
        gaussianfit=adjustMGaus(x,gaussian)
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

def rsquard(y,yfit):
    '''
    *************** rsquard *******************
    Calculate the Pearson Coefficient.
    Notes:
        Make sure the x grid coincide at least with indexes for y and y1.
    Args:
        y(list|array): Original data.
        yfit(list|array): Data ajusted. 
    Returns:
        R2 (float): Pearson Coefficient.
    Usage:
        x=np.arange(-5,5,0.1)
        x0=0
        a=3
        sigma=2
        gaussian=gaus(x,a,x0,sigma)+gaus(x,a-2,x0+2,sigma-1)
        gaussianfit=adjustMGaus(x,gaussian)
        R2=rsquard(gaussian,gaussianfit)
    '''
    #yhat=y1            # or [p(z) for z in x]
    #ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    #ssreg = np.sum((y-yhat)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    #sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    #R2 = 1 - ssreg / sstot
    R2=np.corrcoef(y, yfit)[0,1]
    #print(numpadj[0,1],R2,' geometry file')
    return R2 

def ellipsoidfit(y,yfit,ellipsrsquarefit=0.85,diagnostics=False):
    '''
    *************** ellipsoidfit *******************
    Check the fitness of an ellipsoid in a curve.
    Notes:
        
    Args:
        y(list|array): Original data.
        yfit(list|array): Ellipse ajusted to contour. 
        ellipsrsquarefit (float): [0 > ellipsrsquarefit < 1]
            Pearson Coefficient to validate an ellipse.
        diagnostics (boolean): Used to display all the 
            statistics and plots to identify bugs.
    Returns:
        Rsquard (float) - Fitness of the ellipse.
        check (boolean) - True if gaussian adjust is greater than gaussrsquarefit.
    Usage:
        Check scan_eddym function.
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
    *************** extractprofeddy *******************
    Extracts the profile inside a segment.
    Notes:
        Primary used to extract the mayor or minor axis of an ellipse.
    Args:
        axis (list): Coordinates of an segment.
        field (array): Surface where the profile will be extracted.
        lon,lat (array|list): Coordinates fo the field.
        n (int): Number of desired divisions in the segment.
        gaus (Default:One|None|Multiple): Ajustment to segment.
        kind (Default:linear|cubic|etc): Type of interpolation inside
            segment (For more information check scipy.interpolate interp2d)
        gaussrsquarefit (float): [0 > ellipsrsquarefit < 1]
            Pearson Coefficient to validate an gaussian.
        varname (str): Name of variable, just used for plots.
        diagnostics (boolean): Used to display all the 
            statistics and plots to identify bugs.
    Returns:
        axisdata (array) - Data extracted in the segment.
        check (boolean) - True if gaussian adjust is greater than gaussrsquarefit.
    Usage:
        Check scan_eddym function.
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