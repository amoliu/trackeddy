from __future__ import print_function
import numpy as np
import numpy.ma as ma
import pylab as plt
from trackeddy.datastruct import *
from trackeddy.geometryfunc import *
from trackeddy.init import *
from trackeddy.physics import *
from trackeddy.printfunc import *
from trackeddy.savedata import *
from mpl_toolkits.basemap import Basemap
import sys
import time

def scan_eddym(ssh,lon,lat,levels,date,areamap,mask='',destdir='',physics='',basemap=False,eddycenter='masscenter',ellipsrsquarefit=0.85,eccenfit=0.95,gaussrsquarefit=0.85,checkgauss=True,diagnostics=False,plotdata=False,pprint=True):
    '''
    *************Scan Eddym***********
    Function to identify each eddy using closed contours,
    also this function checks if the elipse adjusted have
    a consistent eccentricity, vorticty and other parameters.
    Usage:
    ssh= Sea Surface Height in cm
    lon,lat=longitude and latitude of your grid.
    levels=where the code will find the closed contours.
    date=date in julian days
    areamap=Section of interest
    mask=Continent mask
    
    Example:
    ssh=Read netCDF4 data with mask or create a mask for your data
    lon=Import your longitude coordinates, if your grid is regular, you can use a linspace instead
    lat=Import your latitude coordinates (same as above).
    levels=List of the levels in which ones you want to find the eddies
    date=Date as Julian Days
    areamap=array([[0,len(lon)],[0,len(lat)]]) Array with the index of your area of interest.
    I used some auxilar functions, each one has his respective author.
    Author: Josue Martinez Moreno, 2017
    '''
    ellipse_path=[]
    contour_path=[]
    mayoraxis_eddy=[]
    minoraxis_eddy=[]
    shapedata=np.shape(ssh)
    if ssh is ma.masked:
        print('Invalid ssh data, must be masked')
        return
    if shapedata == [len(lat), len(lon)]:
        print('Invalid ssh data size, should be [length(lat) length(lon]')
        return
    if np.shape(areamap) == shapedata:
        if np.shape(areamap) == [1, 1] | len(areamap) != len(lat):
            print('Invalid areamap, using NaN for eddy surface area')
        return
    if len(levels)!= 2:
        print('Invalid len of levels, please use the function for multiple levels or use len(levels)==2')
        return
    #Saving mask for future post-processing.  
    
    if mask!='':
        ssh=np.ma.masked_array(ssh, mask)
    sshnan=ssh.filled(np.nan)

    #sshnan=ma.masked_array(okparm, mask=mask[0,:,:])
    #sshnan=sshnan.filled(nan)
    #Obtain the contours of a surface (contourf), this aproach is better than the contour.
    if len(np.shape(lon))== 1 and len(np.shape(lat)) == 1:
        Lon,Lat=np.meshgrid(lon,lat)
    else:
        Lon,Lat=lon,lat
    
    min_x=Lon[0,0]
    min_y=Lat[0,0]
    max_x=Lon[-1,-1]
    max_y=Lat[-1,-1]
    
    if basemap==True:
        fig, ax = plt.subplots(figsize=(10,10))
        m = Basemap(projection='ortho',lat_0=-90,lon_0=-100,resolution='c')
        m.drawcoastlines()
        m.fillcontinents(color='black',lake_color='aqua')
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
        lonm,latm=m(Lon,Lat)
    
        if len(shapedata)==3:
            m.contourf(lonm[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                       latm[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                    sshnan[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],levels=levels)
            plt.show()

        else:
            m.contourf(lonm[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                       latm[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                    sshnan[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],levels=levels)
            plt.show()
            
    if len(shapedata)==3:
        CS=plt.contourf(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                sshnan[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],levels=levels)
    else:
        CS=plt.contourf(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                sshnan[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],levels=levels)
    # Close the contour plot.
    plt.close()
    CONTS=CS.allsegs[:][:]
    #Loop in contours of the levels defined.
    total_contours=0
    eddyn=0
    
    for ii in range(0,np.shape(CONTS)[0]):
        CONTSlvls=CONTS[ii]
        for jj in range(0,np.shape(CONTSlvls)[0]):
            CONTeach=CONTSlvls[jj]
            if (len(CONTeach[:,0]) | len(CONTeach[:,1])) <= 10:
                xx=np.nan
                yy=np.nan
                center=[np.nan,np.nan]
                check=False
            else:
                ellipse,status=fit_ellipse(CONTeach[:,0],CONTeach[:,1],diagnostics=diagnostics)
                if status==True:
                    center = [ellipse['X0_in'],ellipse['Y0_in']]
                    phi = ellipse['phi']
                    axes = [ellipse['a'],ellipse['b']]
                    R = np.arange(0,2.1*np.pi, 0.1)
                    a,b = axes
                    #Ellipse coordinates.
                    xx = ellipse['ellipse'][0]
                    yy = ellipse['ellipse'][1]
                    
                    mayoraxis = ellipse['majoraxis']
                    minoraxis = ellipse['minoraxis']
                    
                    #Area of Contours (contarea) and ellipse (ellipsarea)
                    contarea=PolyArea(CONTeach[:,0],CONTeach[:,1])
                    ellipsarea=PolyArea(xx,yy)
                    
                    # Linear Eccentricity check
                    eccen=eccentricity(a,b)
                    #Record and check how many grid points have land or masked values
                    landcount=0
                    for ii in range(0,len(CONTeach[:,0])):
                        idxcheck,idycheck=find2d(lon,lat,CONTeach[ii,0],CONTeach[ii,1])
                        idxelipcheck,idyelipcheck=find2d(lon,lat,center[0],center[1])
                        if len(shapedata)==3:
                            if sshnan[date,idycheck,idxcheck]==np.nan:
                                landcount=countzeros+1
                        else:
                            if sshnan[idycheck,idxcheck]==np.nan:
                                landcount=countzeros+1
                    if landcount>=len(CONTeach[:,0])/2:
                        #print 'Thisone is land'
                        check=False
                    else:
                        if contarea>ellipsarea:
                            #if contarea/1.5>ellipsarea:
                            #    check=False
                            if eccen<eccenfit and eccen>0:
                                if ellipsarea < 200 and contarea < 200:
                                    if checkgauss==True:
                                        checke=False
                                        checkM=False
                                        checkm=False
                                        
                                        tic=time.time()
                                        ellipseadjust,checke=ellipsoidfit(CONTeach[:,1],ellipse['ellipse']\
                                                                          [1],ellipsrsquarefit=ellipsrsquarefit,\
                                                                          diagnostics=diagnostics)
                                        #print('Time elapsed ellipsoidfit:',str(time.time()-tic))
                                        
                                        if checke==True:
                                            tic=time.time()
                                            if len(shapedata)==3:
                                                profile,checkM=extractprofeddy(mayoraxis,sshnan[date,:,:],lon,lat,50,\
                                                                  gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                  diagnostics=False)
                                                if checkM==True:
                                                    profile,checkm=extractprofeddy(minoraxis,sshnan[date,:,:],lon,lat,50,\
                                                                      gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                      diagnostics=False)
                                            else:
                                                profile,checkM=extractprofeddy(mayoraxis,sshnan[:,:],lon,lat,50,\
                                                                  gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                  diagnostics=False)
                                                if checkM==True:
                                                    profile,checkm=extractprofeddy(minoraxis,sshnan[:,:],lon,lat,50,\
                                                                      gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                      diagnostics=False)
                                        #print('Time elapsed ellipsoidfit:',str(time.time()-tic))
                                        if checkM==True and  checke==True and checkm==True: 
                                            check=True
                                        else:
                                            check=False
                                else:
                                    check=False
                            else:
                                check=False
                        elif contarea<ellipsarea:
                            #if contarea<ellipsarea/1.5:
                                #print 'Removing contour, thisone is really overestimate'
                            #    check=False
                            #elif eccen<0.95 and eccen>0.4:
                            if eccen<eccenfit and eccen>0:
                                if ellipsarea < 200 and contarea<200:
                                    if checkgauss==True:
                                        checke=False
                                        checkM=False
                                        checkm=False
                                        
                                        tic=time.time()
                                        ellipseadjust,checke=ellipsoidfit(CONTeach[:,1],ellipse['ellipse']\
                                                                          [1],ellipsrsquarefit=ellipsrsquarefit,\
                                                                          diagnostics=diagnostics)
                                        #print('Time elapsed ellipsoidfit:',str(time.time()-tic))
                                        
                                        if checke==True:
                                            tic=time.time()
                                            if len(shapedata)==3:
                                                profile,checkM=extractprofeddy(mayoraxis,sshnan[date,:,:],lon,lat,50,\
                                                                  gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                  diagnostics=False)
                                                if checkM==True:
                                                    profile,checkm=extractprofeddy(minoraxis,sshnan[date,:,:],lon,lat,50,\
                                                                      gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                      diagnostics=False)
                                            else:
                                                profile,checkM=extractprofeddy(mayoraxis,sshnan[:,:],lon,lat,50,\
                                                                  gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                  diagnostics=False)
                                                if checkM==True:
                                                    profile,checkm=extractprofeddy(minoraxis,sshnan[:,:],lon,lat,50,\
                                                                      gaus='One',kind='linear',gaussrsquarefit=gaussrsquarefit,\
                                                                      diagnostics=False)
                                        #print('Time elapsed ellipsoidfit:',str(time.time()-tic))
                                        if checkM==True and  checke==True and checkm==True: 
                                            check=True
                                        else:
                                            check=False
                                else:
                                    check=False
                            else:
                                #print 'Removing contour, thisone is really underestimate'
                                check=False
                        else:
                            check=False
                            
                    if diagnostics == True and  check == True:
                        print("Ellipse parameters")
                        print("Ellipse center = ",  center)
                        print("angle of rotation = ",  phi)
                        print("axes (a,b) = ", axes)
                        print("Eccentricity = ",eccen)
                        print("Area (cont,ellips) = ",contarea,ellipsarea)
                        print("Ellipse adjust = ",ellipseadjust)
                    if check==True:# or check==False:
                        ellipse_path.append([xx,yy])
                        contour_path.append([CONTeach[:,0],CONTeach[:,1]])
                        mayoraxis_eddy.append([mayoraxis[0],mayoraxis[1]])
                        minoraxis_eddy.append([minoraxis[0],minoraxis[1]])
                        #Switch from the ellipse center to the position of the maximum value inside de contour
                        if eddycenter == 'maximum':
                            center_eddy=contourmaxvalue(CONTeach[:,0],CONTeach[:,1],sshnan,lon,lat,levels,date)
                        elif eddycenter == 'masscenter':
                            center_eddy=centroidvalue(CONTeach[:,0],CONTeach[:,1],sshnan,lon,lat,levels,date)
                        if eddyn==0:
                            position_eddy=center_eddy
                            position_ellipse=center
                            total_eddy=eddyn
                            area=contarea
                            angle=phi
                            if CS.levels[0] > 0:
                                level=CS.levels[0]
                            else:
                                level=CS.levels[1]
                            levelprnt=level
                        else:
                            position_eddy=np.vstack((position_eddy,center_eddy))
                            position_ellipse=np.vstack((position_ellipse,center))
                            total_eddy=np.vstack((total_eddy,eddyn))
                            area=np.vstack((area,contarea))
                            angle=np.vstack((angle,phi))
                            
                            if CS.levels[0] > 0:
                                levelprnt=CS.levels[0]
                                level=np.vstack((level,levelprnt))
                            else:
                                levelprnt=CS.levels[1]
                                level=np.vstack((level,levelprnt))
                        
                        eddyn=eddyn+1
                    if diagnostics == True and plotdata == True:
                        f, (ax1, ax2) = plt.subplots(1, 2,figsize=(13, 6))
                        if len(shapedata)==3:
                            ax1.contourf(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                                ssh[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]])
                            cc=ax2.pcolormesh(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                                ssh[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                                          vmin=ssh[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].min()\
                                          ,vmax=ssh[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].max())
                            cca=ax2.contour(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                                ssh[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],levels=levels,cmap='jet')
                            ax2.clabel(cca, fontsize=9, inline=1)
                        else:
                            cca=ax1.contourf(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                            sshnan[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],levels=levels)
                            ax1.plot(CONTeach[:,0],CONTeach[:,1],'-r')
                            ax2.plot(CONTeach[:,0],CONTeach[:,1],'-r')
                            ax2.pcolormesh(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                            sshnan[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],vmin=-20,vmax=20)
                            plt.show()
                            f, (ax1, ax2) = plt.subplots(1, 2,figsize=(13, 6))
                            ax1.contourf(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                                    ssh[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]])
                            cc=ax2.pcolormesh(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                                    ssh[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                                              vmin=ssh[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].min()\
                                              ,vmax=ssh[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].max())
                            cca=ax2.contour(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                                    ssh[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],levels=levels,cmap='jet')
                            ax2.clabel(cca, fontsize=9, inline=1)
                        ax1.plot(CONTeach[:,0],CONTeach[:,1],'*r')
                        ax1.plot(xx,yy,'-b')
                        ax1.plot(center[0],center[1],'ob')
                        f.subplots_adjust(right=0.8)
                        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
                        f.colorbar(cc, cax=cbar_ax)
                        ax2.plot(CONTeach[:,0],CONTeach[:,1],'-r')
                        ax2.plot(xx,yy,'-b')
                        ax2.plot(center[0],center[1],'ob')
                        ax2.plot(lon[idxelipcheck],lat[idyelipcheck],'om')
                        ax2.set_ylim([CONTeach[:,1].min(),CONTeach[:,1].max()])
                        ax2.set_xlim([CONTeach[:,0].min(),CONTeach[:,0].max()])
                        plt.show()
                        plt.close()
                total_contours=total_contours+1
            if pprint==True:
                string='Total of contours was: %d - Total of eddies: %d - Level: %.1f   ' % (total_contours,eddyn,levelprnt)
                pt =Printer(); pt.printtextoneline(string)
        position_eddy=np.array(position_eddy)
        position_ellipse=np.array(position_ellipse)
        level=np.array(level)
        mayoraxis_eddy=np.array(mayoraxis_eddy)
        minoraxis_eddy=np.array(minoraxis_eddy)
        eddys=dict_eddym(contour_path,ellipse_path,position_eddy,position_ellipse,mayoraxis_eddy,minoraxis_eddy,\
                         area,angle,total_eddy,level)
        #if destdir!='':
        #    save_data(destdir+'day'+str(date)+'_one_step_cont'+str(total_contours)+'.dat', variable)
        
    return eddys
    
def scan_eddyt(ssh,lat,lon,levels,date,areamap,destdir='',okparm='',diagnostics=False):
    '''
    SCAN_EDDY Scan all of the ssh data passed in (will function correctly if data passed in is a subset)
    ssh: ssh cube with nans for land
    lat: A 1D array of double's that gives the latitude for a given index in ssh data , should be equal to size(ssh, 1)
    lon: A 1D array of double's that gives the longitude for a given index in ssh data, should be equal to size(ssh, 2)
    dates: A 1D array of the dates of ssh data, length should be equal to shape(ssh)[0] 
    destdir: destination directory to save eddies
    '''
    if len(np.shape(ssh))==3:
        if date==0:
            print('Please change the date to the number of iteratios you want')
    else:
        print('Please use the other function scan_eddym')
        return
    for tt in range(0,date):
        print("**********Starting iteration ",tt,"**********")
        eddys=scan_eddym(ssh[tt,:,:],lon,lat,levels,tt,areamap,destdir='',okparm=okparm,diagnostics=diagnostics)
        if tt==0:
            eddytd=dict_eddyt(tt,eddys)
        else:
            eddytd=dict_eddyt(tt,eddys,eddytd) 
        print("**********Finished iteration ",tt,"**********")
    if destdir!='':
        save_data(destdir+str(date),eddies)
    return eddytd

def exeddydt(eddydt,lat,lon,data,threshold,inside='',diagnostics=False):
    '''*************Extract Eddy***********
    Function to extract each eddy in multiple timesteps using closed contours.
    Usage:
    eddydt= Eddy data structure
    lon,lat=longitude and latitude of your grid.
    levels=Level of the contour
    Example:
    Author: Josue Martinez Moreno, 2017
    '''
    justeddy=np.zeros(np.shape(data))
    print('*******Removing of eddies******')
    for key, value in eddydt.items():
        if type(value['time'])==int:
            time=[value['time']]
        else:
            time=[]
            for ii in value['time']:
                time.append(ii[0])
        ct=0 
        for tt in time:
            if len(value['level'])!= 1:
                level=value['level'][ct]
            else:
                level=value['level']
            if type(value['time'])==int:
                lonmi=value['contour'][0][0].min()
                lonma=value['contour'][0][0].max()
                latmi=value['contour'][0][1].min()
                latma=value['contour'][0][1].max()
            else:
                lonmi=value['contour'][ct][0].min()
                lonma=value['contour'][ct][0].max()
                latmi=value['contour'][ct][1].min()
                latma=value['contour'][ct][1].max()
            
            mimcx,mimcy=find2d(lon,lat,lonmi,latmi)
            mamcx,mamcy=find2d(lon,lat,lonma,latma)
            loncm=lon[mimcx-threshold:mamcx+1+threshold]
            latcm=lat[mimcy-threshold:mamcy+1+threshold]
            
            if mimcx==0:
                mimcx=1
            if mimcy==0:
                mimcy=1
            if inside == '':
                datacm=data[tt,mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]-level
                if level > 0:
                    datacm[datacm<=0]=0
                    datacm[datacm>=1000]=0
                elif level < 0:
                    datacm[datacm>=0]=0
                    datacm[datacm<=-1000]=0
            else:
                datacm=data[tt,mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]*1
                insidecm=inside[tt,mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]*1
                if level > 0:
                    insidecm[insidecm<=level]=0
                    insidecm[insidecm>=level]=1
                elif level < 0:
                    insidecm[insidecm>=level]=0
                    insidecm[insidecm<=level]=1
                #if np.shape(insidecm)!=np.shape(datacm):
                #    print('Inside and general field should have the same shape')
                #else:
                datacm=datacm*insidecm
                
            if diagnostics==True:
                plt.figure()
                plt.pcolormesh(lon[mimcx-threshold:mamcx+1+threshold],lat[mimcy-threshold:mamcy+1+threshold],datacm)
                #plt.contourf(lon[mimcx-threshold:mamcx+1+threshold],lat[mimcy-threshold:mamcy+1+threshold],insidecm)
                plt.colorbar()
                cca=plt.contourf(lon[mimcx-threshold:mamcx+1+threshold],lat[mimcy-threshold:mamcy+1+threshold],datacm,alpha=0.5)
                plt.plot(value['contour'][ct][0],value['contour'][ct][1],'-m')
                plt.show()
            justeddy[tt,mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]=datacm
            
            ct=ct+1  
    print('*******End the Removing of eddies******')
    return justeddy

def exeddy(eddydt,lat,lon,data,ct,threshold,inside='',diagnostics=False):
    '''*************Extract Eddy***********
    Function to extract the values of the eddies inside the closed contours.
    Usage:
    eddydt= Eddy data structure
    lon,lat=longitude and latitude of your grid.
    levels=Level of the contour
    Example:
    
    Author: Josue Martinez Moreno, 2017
    '''
    justeddy=np.zeros(np.shape(data))
    print('*******Removing of eddies******')
    for key, value in eddydt.items():
        #print(type(value['level']))
        #print(len(value['level']))
        if len(value['level'])!= 1:
            level=value['level'][ct]
        else:
            level=value['level']
        #print(level)
        rct=value['time']
        #print(len(value['time']))
        if type(value['time'])==int:
            lonmi=np.array(value['contour'][0][0]).min()
            lonma=np.array(value['contour'][0][0]).max()
            latmi=np.array(value['contour'][0][1]).min()
            latma=np.array(value['contour'][0][1]).max()
        else:
            lonmi=value['contour'][ct][0].min()
            lonma=value['contour'][ct][0].max()
            latmi=value['contour'][ct][1].min()
            latma=value['contour'][ct][1].max()
            
        mimcx,mimcy=find2d(lon,lat,lonmi,latmi)
        mamcx,mamcy=find2d(lon,lat,lonma,latma)
        loncm=lon[mimcx-threshold:mamcx+1+threshold]
        latcm=lat[mimcy-threshold:mamcy+1+threshold]

        if mimcx==0:
            mimcx=1
        if mimcy==0:
            mimcy=1
        if inside == '':
            datacm=data[mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]-level
            if level > 0:
                datacm[datacm<=0]=0
                datacm[datacm>=1000]=0
            elif level < 0:
                datacm[datacm>=0]=0
                datacm[datacm<=-1000]=0
        else:
            datacm=data[mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]*1
            insidecm=inside[mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]*1
            if level > 0:
                insidecm[insidecm<=level]=0
                insidecm[insidecm>=level]=1
            elif level < 0:
                insidecm[insidecm>=level]=0
                insidecm[insidecm<=level]=1
            #if np.shape(insidecm)!=np.shape(datacm):
            #    print('Inside and general field should have the same shape')
            #else:
            datacm=datacm*insidecm
            
        if diagnostics==True:
            plt.figure()
            plt.pcolormesh(lon[mimcx-threshold:mamcx+1+threshold],lat[mimcy-threshold:mamcy+1+threshold],datacm)
            plt.contourf(lon[mimcx-threshold:mamcx+1+threshold],lat[mimcy-threshold:mamcy+1+threshold],insidecm)
            plt.colorbar()
            cca=plt.contourf(lon[mimcx-threshold:mamcx+1+threshold],lat[mimcy-threshold:mamcy+1+threshold],datacm,alpha=0.5)
            plt.plot(value['contour'][0],value['contour'][1],'-m')
            plt.show()
            plt.figure()
            plt.pcolormesh(justeddy)
            plt.show()
            plt.close()
            
        justeddy[mimcy-threshold:mamcy+1+threshold,mimcx-threshold:mamcx+1+threshold]=datacm
    print('*******End the Removing of eddies******')
    return justeddy
def analyseddyzt(data,x,y,t0,t1,tstep,maxlevel,minlevel,dzlevel,data_meant='',areamap='',mask='',physics='',eddycenter='masscenter',eccenfit=0.8,ellipsrsquarefit=0.85,gaussrsquarefit=0.65,checkgauss=True,destdir='',saveformat='nc',diagnostics=False,plotdata=False,pprint=False):
    '''
    *************Analys eddy in z and t ***********
    Function to identify each eddy using closed contours, 
    moving in time and contour levels
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    if len(np.shape(data))<3:
        print('If you whant to analyze in time the data need to be 3d [i.e. data(t,x,y)]')
        #return
    if areamap=='':
        areamap=np.array([[0,len(x)],[0,len(y)]])
    if mask == '':
        if len(np.shape(data))<3:
            mask=ma.getmask(data[:,:])
        else:
            mask=ma.getmask(data[0,:,:])
        
    pp =  Printer(); 
    for ii in range(t0,t1,tstep):
        levellist=np.arange(minlevel,maxlevel+dzlevel,dzlevel)
        farlevel=levellist[0]
        if abs(levellist)[0]<abs(levellist)[-1]:
            levellist=np.flipud(levellist)
            farlevel=levellist[0]
        if data_meant=='':
            #print('Be sure the data is an anomaly', end='')
            dataanomaly=data[ii,:,:]
        else:
            dataanomaly=data[ii,:,:]-data_meant
        for ll in levellist:
            if minlevel<0 and maxlevel<0:
                levels=[-500,ll]
            elif minlevel>0 and maxlevel>0:
                levels=[ll,500]
            #tic=time.time()
            eddies=scan_eddym(dataanomaly,x,y,levels,ii,areamap,mask=mask,destdir=destdir\
                          ,physics=physics,eddycenter=eddycenter,checkgauss=checkgauss\
                          ,eccenfit=eccenfit,ellipsrsquarefit=ellipsrsquarefit,gaussrsquarefit=gaussrsquarefit\
                          ,diagnostics=diagnostics,plotdata=plotdata,pprint=pprint)
            #print('ellapse identification:',time.time()-tic)
            if ll == farlevel:
                eddz = dict_eddyz(ii,ll,farlevel,eddies,diagnostics=diagnostics)
            else:
                #tic=time.time()
                eddz = dict_eddyz(ii,ll,farlevel,eddies,eddz,diagnostics=diagnostics)
                #print('ellapse dz:',time.time()-tic)
        if ii==0:
            eddytd=dict_eddyt(ii,eddz)
        else:
            #tic=time.time()
            eddytd=dict_eddyt(ii,eddz,eddytd) 
            #print('ellapse dt:',time.time()-tic)
        pp.timepercentprint(t0,t1,tstep,ii)
    if destdir!='':
        if saveformat=='nc':
            eddync(destdir+str(level)+'.nc',eddytd)
        else:
            np.save(destdir+str(level)+'.npy',eddytd)
    return eddytd

def analyseddyt(data,x,y,level,t0,t1,tstep,data_meant='',areamap='',mask='',physics='',eddycenter='masscenter',ellipsrsquarefit=0.85,eccenfit=0.8,gaussrsquarefit=0.65,checkgauss=True,destdir='',saveformat='nc',diagnostics=False,plotdata=False,pprint=False):
    '''
    *************Analys eddy in t ***********
    Function to identify each eddy using closed contours, 
    moving in time and contour levels
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    if len(np.shape(data))<3:
        print('The data need to have 3d [i.e. data(t,x,y)]')
        #return
    if areamap=='':
        areamap=np.array([[0,len(x)],[0,len(y)]])
    if mask == '':
        if len(np.shape(data))<3:
            mask=ma.getmask(data[:,:])
        else:
            mask=ma.getmask(data[0,:,:])
        
    pp =  Printer(); 
    for ii in range(t0,t1,tstep):
        if data_meant=='':
            #print('Be sure the data is an anomaly', end='')
            dataanomaly=data[ii,:,:]
        else:
            dataanomaly=data[ii,:,:]-data_meant
        #Levels to Analyse, note that one of them is an extreme value,
        #This is because we don't want interference from any other contour.
        if level<0:
            levels=[-500,level]
        elif level>0:
            levels=[level,500]
        eddies=scan_eddym(dataanomaly,x,y,levels,ii,areamap,mask=mask,destdir=destdir\
                      ,physics=physics,eddycenter=eddycenter,checkgauss=checkgauss\
                      ,eccenfit=eccenfit,ellipsrsquarefit=ellipsrsquarefit,gaussrsquarefit=gaussrsquarefit\
                      ,diagnostics=diagnostics,plotdata=plotdata,pprint=pprint)
        if ii==0:
            eddytd=dict_eddyt(ii,eddies)
        else:
            eddytd=dict_eddyt(ii,eddies,eddytd) 
        pp.timepercentprint(t0,t1,tstep,ii)
    if destdir!='':
        if saveformat=='nc':
            eddync(destdir+str(level)+'.nc',eddytd)
        else:
            np.save(destdir+str(level)+'.npy',eddytd)
    return eddytd