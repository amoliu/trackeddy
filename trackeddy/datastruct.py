import numpy as np
import netCDF4 as nc4
from datetime import datetime
import pylab as plt
import warnings
warnings.filterwarnings("ignore")

def eddync():
    print('Work in progress')
    
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

def dict_eddym(contour, ellipse, position_eddy,position_ellipse,majoraxis_eddy,minoraxis_eddy,area,angle,number,level):
    contour_data={'Contour':contour,'Ellipse':ellipse,'Position':position_eddy,'PositionEllipse':position_ellipse,'MajorAxis':majoraxis_eddy,'MinorAxis':minoraxis_eddy,'Area':area,'Angle':angle,'EddyN':number,'Level':level}
    return contour_data

def dict_eddyt(ts,eddys,eddydt=''):
    '''
    Function to create a dictionary with all the eddies.
    All the keys have the following form:
    {eddyn_0:{neddy,time,position,ellipse,contour},eddyn_1:{neddy,time,position,ellipse,contour}}
    Usage:
    
    '''

    #print('eddys',eddys['MinorAxis'])
    if ts==0 or eddydt=='':
        #Define the variable
        eddydt={'eddyn_'+str(eddys['EddyN'][0][0]):{'neddy':eddys['EddyN'][0],'time':ts,'position':eddys['Position'][0],\
                'area':eddys['Area'][0],'ellipse':[eddys['Ellipse'][0]],'contour':[eddys['Contour'][0]],\
                'angle':eddys['Angle'][0],'position_eddy':eddys['PositionEllipse'][0],'level':eddys['Level'][0],\
                'majoraxis':eddys['MajorAxis'][0],'minoraxis':eddys['MinorAxis'][0]}}
        #Append all the new data
        for nn in range(1,len(eddys['EddyN'])):
            eddydt['eddyn_'+str(eddys['EddyN'][nn][0])]={'neddy':eddys['EddyN'][nn],'time':ts,'position':eddys['Position'][nn],\
                    'area':eddys['Area'][nn],'ellipse':[eddys['Ellipse'][nn]],'contour':[eddys['Contour'][nn]],\
                    'angle':eddys['Angle'][nn],'position_eddy':eddys['PositionEllipse'][nn],'level':eddys['Level'][nn],\
                    'majoraxis':eddys['MajorAxis'][nn],'minoraxis':eddys['MinorAxis'][nn]}
    else: 
        checklist=[]
        ceddies=[]
        #print(eddydt)
        for key, value in list(eddydt.items()):
            #print key
            #print len(eddys['EddyN'])
            
            #number=value['neddy']
            
            #if isinstance(number, np.float64) or isinstance(number, np.int64):
            #    number = number
            #else: 
            #    number = number[0]
                
            #level=value['level']
            #print 'Time '+str(time.append(ts))+str(type(time))
            #print time,ts
            #position=value['position']
            #position_eddy=value['position_eddy']
            #ellipse=value['ellipse']
            #print(key,value['neddy'],value['ellipse'])
            #contour=value['contour']
            #print(key)
            #print(type(contour),contour)
            #print(type(ellipse),ellipse)
            minoraxis= value['minoraxis']
            majoraxis= value['majoraxis']
            number=value['neddy']
            time=value['time']
            position=value['position']
            position_eddy=value['position_eddy']
            ellipse=value['ellipse']
            contour=value['contour']
            level= value['level']
            if isinstance(number, np.float64) or isinstance(number, np.int64):
                number = number
            else: 
                number = number[0]
            
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
                angle=eddys['Angle'][nn]                
                if len(np.shape(value['position']))<2:
                    eddyxt0=value['position'][0]
                    eddyyt0=value['position'][1]
                    areae=value['area']
                    anglee=value['angle']
                    timee=value['time']
                else:
                    eddyxt0=value['position'][-1,0]
                    eddyyt0=value['position'][-1,1]
                    areae=value['area'][-1]
                    anglee=value['angle'][-1]
                    timee=value['time'][-1]
                eddyxt1=eddys['Position'][nn][0]
                eddyyt1=eddys['Position'][nn][1]
                
                #ellipse=value['ellipse'].append(eddys['Ellipse'][nn])
                #contour=value['contour'].append(eddys['Contour'][nn])
                #print(ellipse)
                #print([value['ellipse']]+[eddys['Ellipse'][nn]])
                #Check this code again
                if (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                    (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                    int(timee)!=ts and (ts-int(timee))<5 and (areae<=area+area/4 and areae>=area-area/4):
                #   
                #    if ellipse==[None]:
                #        print('stop')
                #        return
                #    print(ellipse)
                    print(np.shape(minoraxis),np.shape(eddys['MinorAxis'][nn]))
                    eddydt['eddyn_'+str(number)]={'neddy':number,'time':np.vstack((value['time'],ts)),\
                                            'position':np.vstack((position,eddys['Position'][nn])),\
                                            'area':np.vstack((areae,area)),'angle':np.vstack((value['angle'],angle)),\
                                            'ellipse':ellipse+[eddys['Ellipse'][nn]],\
                                            'contour':contour+[eddys['Contour'][nn]],\
                                            'position_eddy':np.vstack((position_eddy,eddys['PositionEllipse'][nn])),\
                                            'level':np.vstack((level,eddys['Level'][nn])),\
                                            'minoraxis':np.vstack((minoraxis,eddys['MinorAxis'][nn])),\
                                            'majoraxis':np.vstack((majoraxis,eddys['MajorAxis'][nn]))}
                    print(np.shape(np.vstack((minoraxis,eddys['MinorAxis'][nn]))),np.shape(eddys['MinorAxis'][nn]))
                #print(ellipse.append(eddys['Ellipse'][nn]))
                #checklist.append(number)
                
                #if (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                #    (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                #    (areae<=area+area/4 and areae>=area-area/4) and int(timee)!=ts and (ts-int(timee))<5:# 
                        
                #    eddydt['eddyn_'+str(number)]={'neddy':number,'time':np.vstack((value['time'],ts)),\
                #                            'position':np.vstack((position,eddys['Position'][nn])),\
                #                            'area':np.vstack((areae,area)),'angle':np.vstack((value['angle'],angle)),\
                #                            'ellipse':np.vstack((ellipse,eddys['Ellipse'][nn])),\
                #                            'contour':np.vstack((contour,eddys['Contour'][nn])),\
                #                            'position_eddy':np.vstack((position_eddy,eddys['PositionEllipse'][nn])),\
                #                            'level':np.vstack((level,eddys['Level'][nn])),\
                #                            'minoraxis':np.vstack((minoraxis,eddys['MinorAxis'][nn])),\
                #                            'majoraxis':np.vstack((majoraxis,eddys['MajorAxis'][nn]))}

                    checklist.append(number)
                
                
            ceddies.append(number)

        for nn in range(0,len(eddys['EddyN'])):
            if (eddys['EddyN'][nn]!=checklist).all():
                #print '*****New Eddy*****'
                #print('max value number',np.max(ceddies))
                number=np.max(ceddies)+nn+1
                if isinstance(number, np.float64) or isinstance(number, np.int64):
                    number = number
                    #print('numberint',str(number),eddys['EddyN'][nn][0])
                else: 
                    number = eddys['EddyN'][nn][0]
                #   print('numberlist',str(number))
                
                #print(eddydt[:]['neddy'])
                #print(number)
                eddydt['eddyn_'+str(number)]={'neddy':number,'time':ts,\
                                    'position':eddys['Position'][nn],'area':eddys['Area'][nn],'angle':eddys['Angle'][nn],\
                                    'position_eddy':eddys['PositionEllipse'][nn],'ellipse':[eddys['Ellipse'][nn]],\
                                    'contour':[eddys['Contour'][nn]],'level':eddys['Level'][nn],\
                                    'minoraxis':eddys['MinorAxis'][nn],'majoraxis':eddys['MajorAxis'][nn]}
                #print('New?',number)
        #eddydt={'eddyn'+str(eddys['EddyN']):{'time':ts,'position':eddys['Position'],'ellipse':eddys['Ellipse']\
#                                    ,'contour':eddys['Contour']}}
    #print eddydt
    return eddydt


def dict_eddyz(ts,ll,maxlevel,eddys,eddz='',threshold=1.5,diagnostics=False):
    '''
    Function to create a dictionary with all the important eddy in Z.
    All the keys have the following form:
    {eddyn_0:{neddy,time,position,ellipse,contour},eddyn_1:{neddy,time,position,ellipse,contour}}
    Usage:
    
    '''
    if eddz=='' and maxlevel==ll:
        eddz=eddys
    else:         
        count_new=0
        #Always check in the next one because probable we will have more contours if we get closer to 0.
        checklist=[]
        checklist1=[]
        contour=eddz['Contour']
        ellipse=eddz['Ellipse']
        position=eddz['Position']
        position_ellipse=eddz['PositionEllipse']
        area=eddz['Area']
        angle=eddz['Angle']
        majoraxis=eddz['MajorAxis']
        #print('eddz',majoraxis)
        minoraxis=eddz['MinorAxis']
        number=eddz['EddyN']
        level=eddz['Level']
        #print eddys['EddyN']
        #print np.shape(eddys['Contour'])
        if type(eddys['EddyN'])==int:
            maxlon=np.squeeze(eddys['Contour'])[0][0].max()
            maxlone=np.squeeze(eddys['Ellipse'])[0][0].max()
            maxlat=np.squeeze(eddys['Contour'])[0][1].max()
            maxlate=np.squeeze(eddys['Ellipse'])[0][1].max()
            minlon=np.squeeze(eddys['Contour'])[0][0].min()
            minlone=np.squeeze(eddys['Ellipse'])[0][0].min()
            minlat=np.squeeze(eddys['Contour'])[0][1].min()
            minlate=np.squeeze(eddys['Ellipse'])[0][1].min()
            arealb=eddys['Area']
            anglelb=eddys['Angle']
            eddyxlb=np.squeeze(eddys['Position'])[0]
            eddyylb=np.squeeze(eddys['Position'])[1]
            
            #print eddz['Position'][:,0]
            if type(eddz['EddyN'])==int:
                eddzlist=[eddz['EddyN']]
            else:
                eddzlist=eddz['EddyN']
            for nn0 in range(0,len(eddzlist)):
                eddyxlt=eddz['Position'][nn0,0]
                eddyylt=eddz['Position'][nn0,1]
                arealt=eddz['Area'][nn0]
                
                if (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat) and\
                (arealt>=arealb/threshold and arealt<=arealb):
                    contour[nn0]=np.squeeze(eddys['Contour'])
                    ellipse[nn0]=np.squeeze(eddys['Ellipse'])
                    position[nn0]=np.squeeze(eddys['Position'])
                    area[nn0]=eddys['Area']
                    angle[nn0]=eddys['Angle']
                    majoraxis[nn0]=eddys['MajorAxis']
                    minoraxis[nn0]=eddys['MinorAxis']
                    level[nn0]=eddys['Level']
                    checklist.append(eddys['EddyN'])
                    #print number[nn0][0],eddys['EddyN'][nn1][0]
                    #eddz={'Contour':contour,'Ellipse':ellipse,'Position':position,'Area':area,
                    # 'EddyN':number,'Level':level}
                    if diagnostics==True:
                        print('Contour is growing')
                        plt.figure()
                        plt.plot(eddz['Contour'][nn0,0],eddz['Contour'][nn0,1],'-r')
                        plt.plot(np.squeeze(eddys['Contour'])[0],np.squeeze(eddys['Contour'])[1],'-b')
                        plt.show()
                        plt.close()
                        print('eddyN',nn0)
                        print('eddyt1pos',eddyxlb,eddyylb)
                        print('eddyt0pos',eddyxlt,eddyylt)
                        print('area',arealb, arealt)
                elif (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
                    checklist1.append(eddys['EddyN'])

            if (eddys['EddyN']!= checklist):# and (eddys['EddyN']!= checklist1):
                #print 'test'
                contour=list(contour)
                contour.append(np.squeeze(eddys['Contour']))
                ellipse=list(ellipse)
                ellipse.append(np.squeeze(eddys['Ellipse']))
                position=np.vstack((position,np.squeeze(eddys['Position'])))
                position_ellipse=np.vstack((position_ellipse,np.squeeze(eddys['PositionEllipse'])))
                area=np.vstack((area,eddys['Area']))
                angle=np.vstack((angle,eddys['Angle']))
                majoraxis=list(majoraxis)
                majoraxis.append(eddys['MajorAxis'])
                minoraxis=list(minoraxis)
                minoraxis.append(eddys['MinorAxis'])
                level=np.vstack((level,eddys['Level']))
                number=np.vstack((number,len(number)+1))
                    
        else:
            for nn1 in range(0,len(eddys['EddyN'])):
                maxlon=eddys['Contour'][nn1][0].max()
                maxlone=eddys['Ellipse'][nn1][0].max()
                maxlat=eddys['Contour'][nn1][1].max()
                maxlate=eddys['Ellipse'][nn1][1].max()
                minlon=eddys['Contour'][nn1][0].min()
                minlone=eddys['Ellipse'][nn1][0].min()
                minlat=eddys['Contour'][nn1][1].min()
                minlate=eddys['Ellipse'][nn1][1].min()
                arealb=eddys['Area'][nn1]
                eddyxlb=eddys['Position'][nn1][0]
                eddyylb=eddys['Position'][nn1][1]
                
                #print eddz['Position'][:,0]
                #print(type(eddz['EddyN']))
                #print(eddz['EddyN'])
                if type(eddz['EddyN'])==int:
                    eddzlist=[eddz['EddyN']]
                else:
                    eddzlist=eddz['EddyN']
                for nn0 in range(0,len(eddzlist)):
                    eddyxlt=eddz['Position'][nn0,0]
                    eddyylt=eddz['Position'][nn0,1]
                    arealt=eddz['Area'][nn0]
                    
                    if (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                    (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat) and\
                    (arealt>=arealb/threshold and arealt<=arealb):
                        contour[nn0]=eddys['Contour'][nn1]
                        ellipse[nn0]=eddys['Ellipse'][nn1]
                        position[nn0]=eddys['Position'][nn1]
                        position_ellipse[nn0]=eddys['PositionEllipse'][nn1]
                        area[nn0]=eddys['Area'][nn1]
                        angle[nn0]=eddys['Angle'][nn1]
                        majoraxis[nn0]=eddys['MajorAxis'][nn1]
                        minoraxis[nn0]=eddys['MinorAxis'][nn1]
                        level[nn0]=eddys['Level'][nn1]
                        checklist.append(eddys['EddyN'][nn1][0])
                        #print number[nn0][0],eddys['EddyN'][nn1][0]
                        #eddz={'Contour':contour,'Ellipse':ellipse,'Position':position,'Area':area,
                        # 'EddyN':number,'Level':level}
                        if diagnostics==True:
                            print('Contour is growing')
                            #plt.figure()
                            #plt.plot(eddz['Contour'][nn0,0],eddz['Contour'][nn0,1],'-r')
                            #plt.plot(eddys['Contour'][nn1,0],eddys['Contour'][nn1,1],'-b')
                            #plt.show()
                            #plt.close()
                            print('eddyN', nn1,nn0)
                            print('eddyt1pos',eddyxlb,eddyylb)
                            print('eddyt0pos',eddyxlt,eddyylt)
                            print('area',arealb, arealt)
                    elif (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                    (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
                        checklist1.append(eddys['EddyN'][nn1][0])
                    
                 
            for nn1 in range(0,len(eddys['EddyN'])):
                #print 'tr or fs',(eddys['EddyN'][nn1]!= checklist).all() and (eddys['EddyN'][nn1]!= checklist1).all()
                #print 'eddyn',eddys['EddyN'][nn1]
                if (eddys['EddyN'][nn1]!= checklist).all():# and (eddys['EddyN'][nn1]!= checklist1).all():
                    #print 'test'
                    contour=list(contour)
                    contour.append(eddys['Contour'][nn1])
                    ellipse=list(ellipse)
                    ellipse.append(eddys['Ellipse'][nn1])
                    majoraxis=list(majoraxis)
                    majoraxis.append(eddys['MajorAxis'][nn1])
                    minoraxis=list(minoraxis)
                    minoraxis.append(eddys['MinorAxis'][nn1])
                    position=np.vstack((position,eddys['Position'][nn1]))
                    position_ellipse=np.vstack((position_ellipse,eddys['PositionEllipse'][nn1]))
                    area=np.vstack((area,eddys['Area'][nn1]))
                    angle=np.vstack((angle,eddys['Angle'][nn1]))
                    
                    level=np.vstack((level,eddys['Level'][nn1]))
                    number=np.vstack((number,len(number)+1))
            #print(checklist)   
            #print(number)
        eddz={'Contour':contour,'Ellipse':ellipse,'Position':position,'PositionEllipse':position_ellipse,'Area':area,\
              'MajorAxis':majoraxis,'MinorAxis':minoraxis,'Angle':angle,'EddyN':number,'Level':level}
    return eddz