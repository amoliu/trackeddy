import numpy as np
import netCDF4 as nc4
from datetime import datetime

def eddync():
    print 'Work in progress'
    
def vargeonc(filename,lat,lon,var,tt,varname,nc_description='',units='',dt='',dim='2D'):
    f = nc4.Dataset(filename,'w', format='NETCDF4')
    f.createDimension('lon', len(lon))
    f.createDimension('lat', len(lat))
    f.createDimension('time', None)
    longitude = f.createVariable('Longitude', 'f4', 'lon')
    latitude = f.createVariable('Latitude', 'f4', 'lat')
    
    time = f.createVariable('Time', 'i4', 'time')
    if dim == '3D':
        f.createDimension('z', len(z))
        levels = f.createVariable('Levels', 'i4', 'z')
        varnc = f.createVariable(varname, 'f4', ('time', 'lat', 'lon', 'z'))
        varnc[tt,:,:,:] = var
        levels[:] = z
        levels.units = 'meters [m]'
    else:
        varnc = f.createVariable(varname, 'f4', ('time', 'lat', 'lon'))
        varnc[tt,:,:] = var
        
    longitude[:] = lon
    latitude[:] = lat
    today = datetime.today()
    time_num = today.toordinal()
    time[tt]=time_num
    
    #Add global attributes
    f.description = nc_description
    f.history = "Created " + today.strftime("%d/%m/%y")
    
    #Add local attributes to variable instances
    longitude.units = 'Longitude [Degrees East]'
    latitude.units = 'Latitude [Degrees North]'
    time.units = 'days since Jan 01, 0001'
    varnc.units = units
    f.close()

def save_data(path, variable):
    with file(path, 'w') as variable_file:
        np.savetxt(variable_file, variable)
    variable_file.close() 

def dict_eddym(contour, ellipse, position,area,number,level):
    contour_data={'Contour':contour,'Ellipse':ellipse,'Position':position,'Area':area,'EddyN':number,'Level':level}
    return contour_data

def dict_eddyt(ts,eddys,eddydt=''):
    '''
    Function to create a dictionary with all the eddies.
    All the keys have the following form:
    {eddyn_0:{neddy,time,position,ellipse,contour},eddyn_1:{neddy,time,position,ellipse,contour}}
    Usage:
    
    '''
    if ts==0 or eddydt=='':
        #Define the variable
        eddydt={'eddyn_'+str(eddys['EddyN'][0][0]):{'neddy':eddys['EddyN'][0],'time':ts,'position':eddys['Position'][0],\
                    'area':eddys['Area'][0],'ellipse':eddys['Ellipse'][0],'contour':eddys['Contour'][0],\
                                                 'level':eddys['Level'][0]}}
        #Append all the new data
        for nn in range(1,len(eddys['EddyN'])):
            eddydt['eddyn_'+str(eddys['EddyN'][nn][0])]={'neddy':eddys['EddyN'][nn],'time':ts,'position':eddys['Position'][nn],\
                    'area':eddys['Area'][nn],'ellipse':eddys['Ellipse'][nn],'contour':eddys['Contour'][nn],\
                                                      'level':eddys['Level'][nn]}
    else:         
        for key, value in eddydt.items():
            #print key
            count_new=0
            #print len(eddys['EddyN'])
            for nn in range(0,len(eddys['EddyN'])):
                maxlon=eddys['Contour'][nn][0].max()
                maxlone=eddys['Ellipse'][nn][0].max()
                maxlat=eddys['Contour'][nn][1].max()
                maxlate=eddys['Ellipse'][nn][1].max()
                minlon=eddys['Contour'][nn][0].min()
                minlone=eddys['Ellipse'][nn][0].min()
                minlat=eddys['Contour'][nn][1].min()
                minlate=eddys['Ellipse'][nn][1].min()
                area=eddys['Area'][nn]
                if len(np.shape(value['position']))<2:
                    eddyxt0=value['position'][0]
                    eddyyt0=value['position'][1]
                    areae=value['area']
                    timee=value['time']
                else:
                    eddyxt0=value['position'][-1,0]
                    eddyyt0=value['position'][-1,1]
                    areae=value['area'][-1]
                    timee=value['time'][-1]
                eddyxt1=eddys['Position'][nn][0]
                eddyyt1=eddys['Position'][nn][1]
                #except ValueError:
                #    eddyxt0=value['position'][0]
                    #print 'lalal'+str(len(eddyxt0))
                #    eddyyt0=value['position'][1]
                #    eddyxt1=eddys['Position'][nn][0]
                #    eddyyt1=eddys['Position'][nn][1]
#                print dir(value)
                #print value['time'],timee,ts
                if (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                    (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                    (areae<=area+area/4 and areae>=area-area/4) and (eddyxt1!=eddyxt0 and eddyyt1!=eddyyt0)\
                    and (int(timee)-ts)>5 and int(timee)!=ts :
                #if (value['neddy']==eddys['EddyN'][nn]):
                    #print 'number',nn,'max',maxlon,'t0',eddyxt0,'t1',eddyxt1,'min',minlon,'area0',areae,'area1',area
                    #print nn,maxlat,eddyyt0,eddyyt1,minlat
                    #print "****Tracking Eddy"+str(nn)+"****"
                    number=value['neddy']
                    time=value['time']
                    #print 'Time '+str(time.append(ts))+str(type(time))
                    #print time,ts
                    position=value['position']
                    ellipse=value['ellipse']
                    contour=value['contour']
                    eddydt['eddyn_'+str(number)]={'neddy':number,'time':np.vstack((time,ts)),\
                                            'position':np.vstack((position,eddys['Position'][nn])),\
                                            'area':np.vstack((areae,area)),\
                                            'ellipse':np.vstack((ellipse,eddys['Ellipse'][nn])),\
                                            'contour':np.vstack((contour,eddys['Contour'][nn])),\
                                            'level':np.vstack((level,eddys['Level'][nn]))}
                else:
                    count_new=count_new+1
                    if count_new==len(eddys['EddyN']):
                        #print '*****New Eddy*****'
                        eddydt['eddyn_'+str(eddys['EddyN'][nn])]={'neddy':eddys['EddyN'][nn],'time':ts,\
                                            'position':eddys['Position'][nn],'area':eddys['Area'][nn],\
                                            'ellipse':eddys['Ellipse'][nn],\
                                            'contour':eddys['Contour'][nn],'level':eddys['Level'][nn]}
        #eddydt={'eddyn'+str(eddys['EddyN']):{'time':ts,'position':eddys['Position'],'ellipse':eddys['Ellipse']\
#                                    ,'contour':eddys['Contour']}}
    #print eddydt
    return eddydt




def dict_eddyz(ts,ll,minlevel,eddys,eddz='',threshold=1.5,diagnostics=False):
    '''
    Function to create a dictionary with all the important eddy in Z.
    All the keys have the following form:
    {eddyn_0:{neddy,time,position,ellipse,contour},eddyn_1:{neddy,time,position,ellipse,contour}}
    Usage:
    
    '''
    if ll==minlevel or eddz=='':
        eddz=eddys
    else:         
        count_new=0
        #Always check in the next one because probable we will have more contours if we get closer to 0.
        checklist=[]
        checklist1=[]
        contour=eddz['Contour']
        ellipse=eddz['Ellipse']
        position=eddz['Position']
        area=eddz['Area']
        number=eddz['EddyN']
        level=eddz['Level']
        for nn1 in range(0,len(eddys['EddyN'])):
            maxlon=eddys['Contour'][nn1,0].max()
            maxlone=eddys['Ellipse'][nn1,0].max()
            maxlat=eddys['Contour'][nn1,1].max()
            maxlate=eddys['Ellipse'][nn1,1].max()
            minlon=eddys['Contour'][nn1,0].min()
            minlone=eddys['Ellipse'][nn1,0].min()
            minlat=eddys['Contour'][nn1,1].min()
            minlate=eddys['Ellipse'][nn1,1].min()
            arealb=eddys['Area'][nn1]
            eddyxlb=eddys['Position'][nn1,0]
            eddyylb=eddys['Position'][nn1,1]
            
            #print eddz['Position'][:,0]
            
            for nn0 in range(0,len(eddz['EddyN'])):
                eddyxlt=eddz['Position'][nn0,0]
                eddyylt=eddz['Position'][nn0,1]
                arealt=eddz['Area'][nn0]
                
                if (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat) and\
                (arealt>=arealb/threshold and arealt<=arealb):
                    contour[nn0]=eddys['Contour'][nn1]
                    ellipse[nn0]=eddys['Ellipse'][nn1]
                    position[nn0]=eddys['Position'][nn1]
                    area[nn0]=eddys['Area'][nn1]
                    level[nn0]=eddys['Level'][nn1]
                    checklist.append(eddys['EddyN'][nn1][0])
                    #print number[nn0][0],eddys['EddyN'][nn1][0]
                    #eddz={'Contour':contour,'Ellipse':ellipse,'Position':position,'Area':area,
                    # 'EddyN':number,'Level':level}
                    if diagnostics==True:
                        print 'Contour is growing'
                        plt.figure()
                        plt.plot(eddz['Contour'][nn0,0],eddz['Contour'][nn0,1],'-r')
                        plt.plot(eddys['Contour'][nn1,0],eddys['Contour'][nn1,1],'-b')
                        plt.show()
                        plt.close()
                        print 'eddyN', nn1,nn0
                        print 'eddyt1pos',eddyxlb,eddyylb
                        print 'eddyt0pos',eddyxlt,eddyylt
                        print 'area',arealb, arealt
                elif (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
                    checklist1.append(eddys['EddyN'][nn1][0])
                
                
                #elif(eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                #    (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
                    # dict_eddym(contour, ellipse, position,area,number,level)
                    # {'Contour':contour,'Ellipse':ellipse,'Position':position,'Area':area,
                    # 'EddyN':number,'Level':level}
#        dict_eddym(contour_path, ellipse_path,position,area,total_eddy,level)
        #print checklist
        #print checklist1
        
        for nn1 in range(0,len(eddys['EddyN'])):
            #print 'tr or fs',(eddys['EddyN'][nn1]!= checklist).all() and (eddys['EddyN'][nn1]!= checklist1).all()
            #print 'eddyn',eddys['EddyN'][nn1]
            if (eddys['EddyN'][nn1]!= checklist).all() and (eddys['EddyN'][nn1]!= checklist1).all():
                #print 'test'
                contour=list(contour)
                contour.append(eddys['Contour'][nn1])
                ellipse=list(ellipse)
                ellipse.append(eddys['Ellipse'][nn1])
                position=np.vstack((position,eddys['Position'][nn1]))
                area=np.vstack((area,eddys['Area'][nn1]))
                level=np.vstack((level,eddys['Level'][nn1]))
                number=np.vstack((number,len(eddz['EddyN'])+1))

        eddz={'Contour':contour,'Ellipse':ellipse,'Position':position,'Area':area,\
              'EddyN':number,'Level':level}
    return eddz