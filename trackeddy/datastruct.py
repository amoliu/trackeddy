import numpy as np

def save_data(path, variable):
    with file(path, 'w') as variable_file:
        np.savetxt(variable_file, variable)
    variable_file.close() 

def dict_eddym(contour, ellipse, position,area,number):
    contour_data={'Contour':contour,'Ellipse':ellipse,'Position':position,'Area':area,'EddyN':number}
    return contour_data

def dict_eddyt(ts,eddys,eddydt=''):
    '''
    Function to create a dictionary with al the important
    of each eddy.
    All the keys have the following form:
    {eddyn_0:{neddy,time,position,ellipse,contour},eddyn_1:{neddy,time,position,ellipse,contour}}
    Usage:
    
    '''
    print('ts',ts)
    if ts==0 or eddydt=='':
        eddydt={'eddyn_'+str(eddys['EddyN'][0]):{'neddy':eddys['EddyN'][0],'time':ts,'position':eddys['Position'][0],\
                    'area':eddys['Area'][0],'ellipse':eddys['Ellipse'][0],'contour':eddys['Contour'][0]}}
        for nn in range(1,len(eddys['EddyN'])):
            eddydt['eddyn_'+str(eddys['EddyN'][nn])]={'neddy':eddys['EddyN'][nn],'time':ts,'position':eddys['Position'][nn],\
                    'area':eddys['Area'][nn],'ellipse':eddys['Ellipse'][nn],'contour':eddys['Contour'][nn]}
    else:         
        for key, value in eddydt.items():
            print key
            count_new=0
            for nn in range(0,len(eddys['EddyN'])):
                maxlon=eddys['Contour'][nn,0].max()
                maxlone=eddys['Ellipse'][nn,0].max()
                maxlat=eddys['Contour'][nn,1].max()
                maxlate=eddys['Ellipse'][nn,1].max()
                minlon=eddys['Contour'][nn,0].min()
                minlone=eddys['Ellipse'][nn,0].min()
                minlat=eddys['Contour'][nn,1].min()
                minlate=eddys['Ellipse'][nn,1].min()
                area=eddys['Area'][nn]
                if len(np.shape(value['position']))<2:
                    eddyxt0=value['position'][0]
                    eddyyt0=value['position'][1]
                    areae=value['area']
                else:
                    #print np.shape(value['position'])
                    eddyxt0=value['position'][-1,0]
                    eddyyt0=value['position'][-1,1]
                    areae=value['area'][-1]
                eddyxt1=eddys['Position'][nn][0]
                eddyyt1=eddys['Position'][nn][1]
                #except ValueError:
                #    eddyxt0=value['position'][0]
                    #print 'lalal'+str(len(eddyxt0))
                #    eddyyt0=value['position'][1]
                #    eddyxt1=eddys['Position'][nn][0]
                #    eddyyt1=eddys['Position'][nn][1]
                if (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                    (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                    (areae<=area+area/4 and areae>=area-area/4) and (eddyxt1!=eddyxt0 and eddyyt1!=eddyyt0):
                #if (value['neddy']==eddys['EddyN'][nn]):
                    print 'number',nn,'max',maxlon,'t0',eddyxt0,'t1',eddyxt1,'min',minlon,'area0',areae,'area1',area
                    print nn,maxlat,eddyyt0,eddyyt1,minlat
                    print "****Tracking Eddy"+str(nn)+"****"
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
                                            'contour':np.vstack((contour,eddys['Contour'][nn]))}
                else:
                    count_new=count_new+1
                    if count_new==len(eddys['EddyN']):
                        print '*****New Eddy*****'
                        eddydt['eddyn_'+str(eddys['EddyN'][nn])]={'neddy':eddys['EddyN'][nn],'time':[ts],\
                                            'position':eddys['Position'][nn],'area':eddys['Area'][nn],\
                                            'ellipse':eddys['Ellipse'][nn],\
                                            'contour':eddys['Contour'][nn]}
        #eddydt={'eddyn'+str(eddys['EddyN']):{'time':ts,'position':eddys['Position'],'ellipse':eddys['Ellipse']\
#                                    ,'contour':eddys['Contour']}}
    #print eddydt
    return eddydt