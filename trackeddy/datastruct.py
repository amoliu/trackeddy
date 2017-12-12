import numpy as np
import netCDF4 as nc4
from datetime import datetime
import pylab as plt
import warnings
from trackeddy.savedata import *
warnings.filterwarnings("ignore")

def dict_eddym(contour, ellipse, position_eddy,position_ellipse,majoraxis_eddy,minoraxis_eddy,area,angle,number,level,gaussianfitdict):
    '''
    ********************** dict_eddym **********************
    Create a dictionary with a basic structure of the eddies. Can contain the structure of one eddy or multple eddies.
    Notes: 
        
    Args:
        contour (array): Close contour used to track any gaussian (Eddy).
        ellipse (array): Ellipse ajusted to identify a gaussian shape. 
        position_eddy (array): Positions (X, Y) describing the location of the feature. 
        position_ellipse (array): Positions (X, Y) describing the location of the ajusted ellipse. 
        majoraxis_eddy (array): Positions (X, Y) describing the mayor axis of the ajusted ellipse.
        minoraxis_eddy (array): Positions (X, Y) describing the minor axis of the ajusted ellipse.
        area (array): Area of the feature.
        angle (array): Angle of the ajusted ellipse starting on the 0 degree convention (counterclockwise).
        number (list): Number of the feature. Flag used to track features in time.
        level (array): Contour level used during the identification of the feature.
        GaussianFitdict(dict): Syntetic field of a 2D gaussian adjust.
    Returns:
        contour_data - Dictionary containing all the relevant information of the features.
        The new dictionary containing all the information has the following form:        
        {'Contour':contour,'Ellipse':ellipse,'Position':position_eddy,\
        'PositionEllipse':position_ellipse,'MajorAxis':majoraxis_eddy,\
        'MinorAxis':minoraxis_eddy,'Area':area,'Angle':angle,'EddyN':number,'Level':level,'2DGaussianFit':gausssianfitp}
    Usage:
        eddys=dict_eddym(contour_path,ellipse_path,position_eddy,position_ellipse,mayoraxis_eddy,minoraxis_eddy,\
                     area,angle,total_eddy,level)
    '''
    contour_data={'Contour':contour,'Ellipse':ellipse,'Position':position_eddy,'PositionEllipse':position_ellipse,'MajorAxis':majoraxis_eddy,'MinorAxis':minoraxis_eddy,'Area':area,'Angle':angle,'EddyN':number,'Level':level,'2DGaussianFit':gaussianfitdict}
    return contour_data

def dict_eddyt(ts,eddys,eddydt=''):
    '''
    ********************** dict_eddyt **********************
    Create a dictionary with all the eddies and it's track on time. When 'ts==0' or eddydt is not defined it start creating a new dictionary otherwise it grows the dictionary.
    Notes:
        Check for some possible bugs where it's saving more than once some features.
    Args:
        ts (int): Step or counter used to move in time.
        eddys (dict): Dictionary containing all the relevant information of the features. Check "dict_eddym".
        eddydt(dict): Output dictionary of this function used to grow the features. Check "dict_eddyt".
    Returns:
        eddydt - Dictionary containig each eddy and it's track on time.
        All the keys have the following form:
        {eddyn_0:{neddy,time,position,ellipse,contour,...},eddyn_1:{neddy,time,position,ellipse,contour,...},...}
    Usage:
        if tt==0:
            eddytd=dict_eddyt(tt,eddys)
        else:
            eddytd=dict_eddyt(tt,eddys,eddytd) 
        
    '''
    if ts==0 or eddydt=='':
        eddydt={'eddyn_'+str(eddys['EddyN'][0][0]):{'neddy':eddys['EddyN'][0],'time':np.array([ts]),'position':eddys['Position'][0],\
                'area':eddys['Area'][0],'ellipse':[eddys['Ellipse'][0]],'contour':[eddys['Contour'][0]],\
                'angle':eddys['Angle'][0],'position_eddy':eddys['PositionEllipse'][0],'level':eddys['Level'][0],\
                'majoraxis':eddys['MajorAxis'][0],'minoraxis':eddys['MinorAxis'][0],\
                '2dgaussianfit':eddys['2DGaussianFit'][0],'timetracking':True}}
        for nn in range(1,len(eddys['EddyN'])):
            eddydt['eddyn_'+str(eddys['EddyN'][nn][0])]={'neddy':eddys['EddyN'][nn],'time':np.array([ts]),\
                    'position':eddys['Position'][nn],\
                    'area':eddys['Area'][nn],'ellipse':[eddys['Ellipse'][nn]],'contour':[eddys['Contour'][nn]],\
                    'angle':eddys['Angle'][nn],'position_eddy':eddys['PositionEllipse'][nn],'level':eddys['Level'][nn],\
                    'majoraxis':eddys['MajorAxis'][nn],'minoraxis':eddys['MinorAxis'][nn],\
                    '2dgaussianfit':eddys['2DGaussianFit'][nn],'timetracking':True}
    else: 
        checklist=[]
        checklist1=[]
        ceddies=[]

        for key, value in list(eddydt.items()):
            if value['timetracking']==False:
                status='No Tracking'
            else:
                minoraxis= value['minoraxis']
                majoraxis= value['majoraxis']
                number=value['neddy']
                times=value['time']
                position=value['position']
                position_eddy=value['position_eddy']
                ellipse=value['ellipse']
                contour=value['contour']
                level= value['level']
                gaussianfit=value['2dgaussianfit']
                if isinstance(number, np.float64) or isinstance(number, np.int64):
                    number = number
                else: 
                    number = number[0]
                    
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
                
                    if len(np.shape(minoraxis))==3:
                        minoraxis=np.squeeze(minoraxis)
                    if len(np.shape(majoraxis))==3:
                        majoraxis=np.squeeze(majoraxis)
                        
                if type(eddys['EddyN'])== int:
                    maxlon=eddys['Contour'][0][0].max()
                    maxlone=eddys['Ellipse'][0][0].max()
                    maxlat=eddys['Contour'][0][1].max()
                    maxlate=eddys['Ellipse'][0][1].max()
                    minlon=eddys['Contour'][0][0].min()
                    minlone=eddys['Ellipse'][0][0].min()
                    minlat=eddys['Contour'][0][1].min()
                    minlate=eddys['Ellipse'][0][1].min()
                    area=eddys['Area']
                    angle=eddys['Angle']            
                    eddyxt1=eddys['Position'][0]
                    eddyyt1=eddys['Position'][1]
                    
                    if (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                        (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                        int(timee)!=ts and (ts-int(timee))<5 and (checklist1!=eddys['EddyN']):
                        ## Last condition added 27 of sepbermber to stop repeating eddies trackings. ##
                        eddydt['eddyn_'+str(number)]={'neddy':number,'time':np.vstack((value['time'],ts)),\
                                            'position':np.vstack((position,eddys['Position'])),\
                                            'area':np.vstack((value['area'],area)),'angle':np.vstack((value['angle'],angle)),\
                                            'ellipse':ellipse+[eddys['Ellipse'][0]],\
                                            'contour':contour+[eddys['Contour'][0]],\
                                            'position_eddy':np.vstack((position_eddy,eddys['PositionEllipse'])),\
                                            'level':np.vstack((level,eddys['Level'])),\
                                            'minoraxis':np.vstack((minoraxis,np.squeeze(eddys['MinorAxis'][0]))),\
                                            'majoraxis':np.vstack((majoraxis,np.squeeze(eddys['MajorAxis'][0]))),\
                                            '2dgaussianfit':np.vstack((gaussianfit,eddys['2DGaussianFit'][0])),\
                                            'timetracking':True}
                        
                        checklist1.append(eddys['EddyN'])
                        checklist.append(number)
                    elif (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                           (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                        int(timee)!=ts and (ts-int(timee))>=5:
                        #print('Add Removal')
                        eddydt['eddyn_'+str(number)]['timetracking']=False    
                else: 
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
                        eddyxt1=eddys['Position'][nn][0]
                        eddyyt1=eddys['Position'][nn][1]
                        #print('time_step',np.shape(gaussianfit),np.shape(eddys['2DGaussianFit'][nn]))
                        #print('time',gaussianfit,eddys['2DGaussianFit'][nn])
                        #print(np.shape(minoraxis),np.shape(np.squeeze(eddys['MinorAxis'][nn])))
                        if (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                            (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                            int(timee)!=ts and (ts-int(timee))<5 and (checklist1!=eddys['EddyN'][nn][0]).all():
                            ## Last condition added 27 of sepbermber to stop repeating eddies trackings. ##
                            eddydt['eddyn_'+str(number)]={'neddy':number,'time':np.vstack((value['time'],ts)),\
                                            'position':np.vstack((position,eddys['Position'][nn])),\
                                            'area':np.vstack((value['area'],area)),'angle':np.vstack((value['angle'],angle)),\
                                            'ellipse':ellipse+[eddys['Ellipse'][nn]],\
                                            'contour':contour+[eddys['Contour'][nn]],\
                                            'position_eddy':np.vstack((position_eddy,eddys['PositionEllipse'][nn])),\
                                            'level':np.vstack((level,eddys['Level'][nn])),\
                                            'minoraxis':np.vstack((np.squeeze(minoraxis),np.squeeze(eddys['MinorAxis'][nn]))),\
                                            'majoraxis':np.vstack((np.squeeze(majoraxis),np.squeeze(eddys['MajorAxis'][nn]))),\
                                            '2dgaussianfit':np.vstack((gaussianfit,eddys['2DGaussianFit'][nn])),\
                                            'timetracking':True}
                        
                            checklist1.append(eddys['EddyN'][nn][0])
                            checklist.append(number)
                        elif (eddyxt1<=maxlon and eddyxt1>=minlon and eddyyt1<=maxlat and eddyyt1>=minlat) and\
                            (eddyxt0<=maxlon and eddyxt0>=minlon and eddyyt0<=maxlat and eddyyt0>=minlat) and\
                            int(timee)!=ts and (ts-int(timee))>=5:
                            #print('Add Removal')
                            eddydt['eddyn_'+str(number)]['timetracking']=False
                    
                ceddies.append(number)
                
        counter=0
        checklist2=[]
        checklist3=[]
        for key, value in list(eddydt.items()):
            if value['timetracking']==True:
                number=value['neddy']
                if isinstance(number, np.float64) or isinstance(number, np.int64):
                    number = number
                else: 
                    number = number[0]
                if type(eddys['EddyN'])== int:
                    number1=eddys['EddyN']
                    
                    if (number!=checklist) and (number1!=checklist1) and (number!=checklist2)\
                        and (number1!=checklist3):
                        checklist2.append(number)
                        checklist3.append(number1)
                        counter=counter+1
                        number=np.max(ceddies)+counter
                        if isinstance(number, np.float64) or isinstance(number, np.int64):
                            number = number
                        else: 
                            number = eddys['EddyN'][nn][0]
                        eddydt['eddyn_'+str(number)]={'neddy':number,'time':np.array([ts]),\
                                        'position':eddys['Position'][nn],'area':eddys['Area'],'angle':eddys['Angle'],\
                                        'position_eddy':eddys['PositionEllipse'],'ellipse':[eddys['Ellipse'][0]],\
                                        'contour':[eddys['Contour'][0]],'level':eddys['Level'],\
                                        'minoraxis':eddys['MinorAxis'][0],'majoraxis':eddys['MajorAxis'][0],\
                                        '2dgaussianfit':eddys['2DGaussianFit'][0],\
                                        'timetracking':True}
                else:
                    for nn in range(0,len(eddys['EddyN'])):
                        number1=eddys['EddyN'][nn]
                    
                        if (number!=checklist).all() and (number1!=checklist1).all() and (number!=checklist2).all()\
                            and (number1!=checklist3).all():
                            checklist2.append(number)
                            checklist3.append(number1)
                            counter=counter+1
                            number=np.max(ceddies)+counter
                            if isinstance(number, np.float64) or isinstance(number, np.int64):
                                number = number
                            else: 
                                number = eddys['EddyN'][nn][0]
                            
                            eddydt['eddyn_'+str(number)]={'neddy':number,'time':np.array([ts]),\
                                        'position':eddys['Position'][nn],'area':eddys['Area'][nn],'angle':eddys['Angle'][nn],\
                                        'position_eddy':eddys['PositionEllipse'][nn],'ellipse':[eddys['Ellipse'][nn]],\
                                        'contour':[eddys['Contour'][nn]],'level':eddys['Level'][nn],\
                                        'minoraxis':eddys['MinorAxis'][nn],'majoraxis':eddys['MajorAxis'][nn],\
                                        '2dgaussianfit':eddys['2DGaussianFit'][nn],\
                                        'timetracking':True}
    return eddydt


def dict_eddyz(ts,ll,maxlevel,eddys,eddz='',threshold=1.5,diagnostics=False):
    '''
    ********************** dict_eddyz **********************
    Create a dictionary with all the eddies and it's develop in delta eta, where the biggest contour is assigned as the contour that will be tracked.
    Notes:
        Check for some possible bugs where it's saving more than once some features.
    Args:
        ts (int): Step or counter used to move in time.
        ll (int): Level analysed used to move in delta eta, where delta eta is the changes in elevation of the surface.
        maxlevel (int): Max level defined by the user. 
        eddys (dict): Dictionary containing all the relevant information of the features. Check "dict_eddym".
        eddz(dict): Output dictionary of this function used to grow the features. Check "dict_eddyz".
        threshold (float): Used to grow the detection; R=R0+threshold, where R0 is the radius of the contour.
        diagnostics (boolean): Used to display all the statistics and plots to identify bugs.
    Returns:
        eddz - Dictionary containing the largest contour and parameters of each features.
        All the keys have the following form:
        {eddyn_0:{neddy,time,position,ellipse,contour},eddyn_1:{neddy,time,position,ellipse,contour}}
    Usage:
        if ll == 0:
            eddz = dict_eddyz(ii,ll,farlevel,eddies,diagnostics=diagnostics)
        else:
            eddz = dict_eddyz(ii,ll,farlevel,eddies,eddz,diagnostics=diagnostics)
    '''
    checklist=[]
    checklist1=[]
    if eddz=='' or maxlevel==ll:
        eddz=eddys
    else:         
        count_new=0
        #Always check in the next one because probable we will have more contours if we get closer to 0.
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
        gauss2d=np.array(eddz['2DGaussianFit'])
        #print('goodone:',np.shape(gauss2d))
        #print eddys['EddyN']
        #print np.shape(eddys['Contour'])
        if type(eddys['EddyN'])==int and type(eddz['EddyN'])==int:
            #print('int - int')
            maxlon=np.squeeze(eddys['Contour'])[0].max()
            maxlone=np.squeeze(eddys['Ellipse'])[0].max()
            maxlat=np.squeeze(eddys['Contour'])[1].max()
            maxlate=np.squeeze(eddys['Ellipse'])[1].max()
            minlon=np.squeeze(eddys['Contour'])[0].min()
            minlone=np.squeeze(eddys['Ellipse'])[0].min()
            minlat=np.squeeze(eddys['Contour'])[1].min()
            minlate=np.squeeze(eddys['Ellipse'])[1].min()
            arealb=eddys['Area']
            anglelb=eddys['Angle']
            eddyxlb=np.squeeze(eddys['Position'])[0]
            eddyylb=np.squeeze(eddys['Position'])[1]
            eddyxlt=eddz['Position'][0]
            eddyylt=eddz['Position'][1]
            arealt=eddz['Area']
                
            if (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
            (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):# and\
            #(arealt>=arealb/threshold and arealt<=arealb):
                contour=eddys['Contour']
                ellipse=eddys['Ellipse']
                position=eddys['Position']
                area=eddys['Area']
                angle=eddys['Angle']
                majoraxis=eddys['MajorAxis']
                minoraxis=eddys['MinorAxis']
                level=eddys['Level']
                print('1.1',np.shape(gauss2d))
                gauss2d=eddys['2DGaussianFit']
                print('1.2',np.shape(gauss2d))
                checklist=eddys['EddyN']
                if diagnostics==True:
                    print('Contour is growing')
                    plt.figure()
                    plt.plot(eddz['Contour'][0],eddz['Contour'][1],'-r')
                    plt.plot(np.squeeze(eddys['Contour'])[0],np.squeeze(eddys['Contour'])[1],'-b')
                    plt.show()
                    plt.close()
                    print('eddyN',eddz['EddyN'])
                    print('eddyt1pos',eddyxlb,eddyylb)
                    print('eddyt0pos',eddyxlt,eddyylt)
                    print('area',arealb, arealt)
            #elif (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
            #    (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
            #        checklist1.append(eddys['EddyN'])

            if (eddys['EddyN']!= checklist):# and (eddys['EddyN']!= checklist1):
                #print 'test'
                #print('1')
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
                print('1.3',np.shape(gauss2d),np.shape(eddys['2DGaussianFit']))
                gauss2d=np.vstack((gauss2d,eddys['2DGaussianFit']))
                print('1.4',np.shape(gauss2d))
                number=np.vstack((number,number+1))
            print('1',np.shape(gauss2d))
        elif type(eddys['EddyN'])==int and type(eddz['EddyN'])!=int:
            #print('int - no int')
            maxlon=np.squeeze(eddys['Contour'])[0].max()
            maxlone=np.squeeze(eddys['Ellipse'])[0].max()
            maxlat=np.squeeze(eddys['Contour'])[1].max()
            maxlate=np.squeeze(eddys['Ellipse'])[1].max()
            minlon=np.squeeze(eddys['Contour'])[0].min()
            minlone=np.squeeze(eddys['Ellipse'])[0].min()
            minlat=np.squeeze(eddys['Contour'])[1].min()
            minlate=np.squeeze(eddys['Ellipse'])[1].min()
            arealb=eddys['Area']
            anglelb=eddys['Angle']
            eddyxlb=np.squeeze(eddys['Position'])[0]
            eddyylb=np.squeeze(eddys['Position'])[1]
            #print(eddyxlb,eddyylb)
            #print(maxlon,minlon,maxlat,minlat)
            #print(eddz['EddyN'])
            for nn0 in range(0,len(eddz['EddyN'])):
                eddyxlt=eddz['Position'][nn0,0]
                eddyylt=eddz['Position'][nn0,1]
                arealt=eddz['Area'][nn0]
                #print(nn0)
                if (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):# and\
                #(arealt>=arealb/threshold and arealt<=arealb):
                    #print(nn0,'lalala',eddys['EddyN'])
                    contour[nn0]=np.squeeze(eddys['Contour'])
                    ellipse[nn0]=np.squeeze(eddys['Ellipse'])
                    position[nn0]=np.squeeze(eddys['Position'])
                    area[nn0]=eddys['Area']
                    angle[nn0]=eddys['Angle']
                    majoraxis[nn0]=eddys['MajorAxis']
                    minoraxis[nn0]=eddys['MinorAxis']
                    level[nn0]=eddys['Level']
                    gauss2d[nn0]=eddys['2DGaussianFit']
                    checklist=eddys['EddyN']
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
                #elif (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                #(eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
                #    checklist1.append(eddys['EddyN'])
            #print(eddys['EddyN'], checklist)
            #print(eddys['EddyN']!= checklist)
            if (eddys['EddyN']!= checklist):# and (eddys['EddyN']!= checklist1):
                #print 'test'
                #print('2',eddys['EddyN'],checklist)
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
                #print(np.shape(gauss2d),gauss2d)
                #print(np.shape(eddys['2DGaussianFit']),eddys['2DGaussianFit'])
                gauss2d=np.vstack((gauss2d,eddys['2DGaussianFit']))
                number=np.vstack((number,len(number)+1))
            #print('2',np.shape(gauss2d))
        elif type(eddys['EddyN'])!=int and type(eddz['EddyN'])==int:
            #print('no int - int')
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
                eddyxlt=eddz['Position'][0]
                eddyylt=eddz['Position'][1]
                arealt=eddz['Area']
                    
                if (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):# and\
                #(arealt>=arealb/threshold and arealt<=arealb):
                    contour=eddys['Contour'][nn1]
                    ellipse=eddys['Ellipse'][nn1]
                    position=eddys['Position'][nn1]
                    position_ellipse=eddys['PositionEllipse'][nn1]
                    area=eddys['Area'][nn1]
                    angle=eddys['Angle'][nn1]
                    majoraxis=eddys['MajorAxis'][nn1]
                    minoraxis=eddys['MinorAxis'][nn1]
                    level=eddys['Level'][nn1]
                    #print('3.0',np.shape(gauss2d))
                    gauss2d=[eddys['2DGaussianFit'][nn1]]
                    #print('3.1',np.shape(gauss2d))
                    checklist.append(eddys['EddyN'][nn1][0])
                    if diagnostics==True:
                        print('Contour is growing')
                        print('eddyN', nn1)
                        print('eddyt1pos',eddyxlb,eddyylb)
                        print('eddyt0pos',eddyxlt,eddyylt)
                        print('area',arealb, arealt)
                #elif (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                #(eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
                #    checklist1.append(eddys['EddyN'][nn1][0])
                    
                #print(eddz['EddyN'])
            for nn1 in range(0,len(eddys['EddyN'])):
                #print 'tr or fs',(eddys['EddyN'][nn1]!= checklist).all() and (eddys['EddyN'][nn1]!= checklist1).all()
                #print 'eddyn',eddys['EddyN'][nn1]
                if (eddys['EddyN'][nn1]!= checklist).all():# and (eddys['EddyN'][nn1]!= checklist1).all():
                    #print 'test'
                    #print('3')
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
                    #print('3.2',np.shape(gauss2d))
                    gauss2d=np.vstack((gauss2d,[eddys['2DGaussianFit'][nn1]]))
                    #print('3.3',np.shape(gauss2d))
                    #print(number)
                    if type(number)==int:
                        numz=number+1
                    else:
                        numz=len(number)+1
                    number=np.vstack((number,numz))
            #print('3',np.shape(gauss2d))       
        else:
            #print('no int - no int')
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
                for nn0 in range(0,len(eddz['EddyN'])):
                    eddyxlt=eddz['Position'][nn0,0]
                    eddyylt=eddz['Position'][nn0,1]
                    arealt=eddz['Area'][nn0]
                    
                    if (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                    (eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):# and\
                    #(arealt>=arealb/threshold and arealt<=arealb):
                        contour[nn0]=eddys['Contour'][nn1]
                        ellipse[nn0]=eddys['Ellipse'][nn1]
                        position[nn0]=eddys['Position'][nn1]
                        position_ellipse[nn0]=eddys['PositionEllipse'][nn1]
                        area[nn0]=eddys['Area'][nn1]
                        angle[nn0]=eddys['Angle'][nn1]
                        majoraxis[nn0]=eddys['MajorAxis'][nn1]
                        minoraxis[nn0]=eddys['MinorAxis'][nn1]
                        level[nn0]=eddys['Level'][nn1]
                        gauss2d[nn0]=np.array(eddys['2DGaussianFit'][nn1])
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
                    #elif (eddyxlb<=maxlon and eddyxlb>=minlon and eddyylb<=maxlat and eddyylb>=minlat) and\
                    #(eddyxlt<=maxlon and eddyxlt>=minlon and eddyylt<=maxlat and eddyylt>=minlat):
                    #    checklist1.append(eddys['EddyN'][nn1][0])
                    
                 
            for nn1 in range(0,len(eddys['EddyN'])):
                #print 'tr or fs',(eddys['EddyN'][nn1]!= checklist).all() and (eddys['EddyN'][nn1]!= checklist1).all()
                #print 'eddyn',eddys['EddyN'][nn1]
                if (eddys['EddyN'][nn1]!= checklist).all():# and (eddys['EddyN'][nn1]!= checklist1).all():
                    #print 'test'
                    #print('4')
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
                    #print(np.shape(gauss2d),np.shape(eddys['2DGaussianFit'][nn1]))
                    gauss2d=np.vstack((gauss2d,[eddys['2DGaussianFit'][nn1]]))
                    number=np.vstack((number,len(number)+1))
            #print(checklist)   
            #print(number)
        #print('z_step',np.shape(gauss2d))
            #print('4',np.shape(gauss2d))
        eddz={'Contour':contour,'Ellipse':ellipse,'Position':position,'PositionEllipse':position_ellipse,'Area':area,\
              'MajorAxis':majoraxis,'MinorAxis':minoraxis,'Angle':angle,'EddyN':number,'Level':level,'2DGaussianFit':gauss2d}
    return eddz

def joindict(dict1,dict2):
    '''
    ********************** joindict **********************
    Join two dictionaries containing the features parameters in time an delta eta. Useful when datasets are in multiple files.
    Notes:
        Check for some possible bugs where it's saving more than once some features.
    Args:
        dict1 (dict): Dictionary containing all the relevant information of the features. Check "dict_eddym".
        dict2 (dict): Dictionary containing all the relevant information of the features. Check "dict_eddym".
    Returns:
        dict1 - Dictionary containing both input dictionaries.
        Same keys that the initial dicionaries.
    Usage:
        dictjoin=joindict(eddydt1,eddydt2)
    '''
    checklist=[]
    checklist1=[]
    checklist2=[]
    for key, value in list(dict1.items()):
        check=False
        if type(value['time'])==int or len(value['time'])==1:
            eddyxt0=value['position'][0]
            eddyyt0=value['position'][1]
            if value['time']>=89-5:
                check=True
                timee=value['time']
        else:
            eddyxt0=value['position'][-1][0]
            eddyyt0=value['position'][-1][1]
            if value['time'][-1]>=89-5:
                check=True
                timee=value['time'][-1]
        if check==True:
            if type(value['time'])==int:
                lonmi0=value['contour'][0][0].min()
                lonma0=value['contour'][0][0].max()
                latmi0=value['contour'][0][1].min()
                latma0=value['contour'][0][1].max()
            else:
                lonmi0=value['contour'][-1][0].min()
                lonma0=value['contour'][-1][0].max()
                latmi0=value['contour'][-1][1].min()
                latma0=value['contour'][-1][1].max()
                
            for key1, value1 in list(dict2.items()):
                if type(value1['time'])==int:
                    ts=value1['time']
                    eddyxt1=value1['position'][0]
                    eddyyt1=value1['position'][1]
                else:
                    ts=int(value1['time'][0])
                    eddyxt1=value1['position'][0][0]
                    eddyyt1=value1['position'][0][1]
                if ts<=5:
                    if (eddyxt1<=lonma0 and eddyxt1>=lonmi0 and eddyyt1<=latma0 and eddyyt1>=latmi0) and\
                        (eddyxt0<=lonma0 and eddyxt0>=lonmi0 and eddyyt0<=latma0 and eddyyt0>=latmi0):
                        dict1[key]={'neddy':int(value['neddy']),'time':np.vstack((value['time'],value1['time']+timee)),\
                                            'position':np.vstack((value['position'],value1['position'])),\
                                            'area':np.vstack((value['area'],value1['area'])),\
                                            'angle':np.vstack((value['angle'],value1['angle'])),\
                                            'ellipse':value['ellipse']+value1['ellipse'],\
                                            'contour':value['contour']+value1['contour'],\
                                            'position_eddy':np.vstack((value['position_eddy'],value1['position_eddy'])),\
                                            'level':np.vstack((value['level'],value1['level'])),\
                                            'minoraxis':np.vstack((value['minoraxis'],value1['minoraxis'])),\
                                            'majoraxis':np.vstack((value['majoraxis'],value1['majoraxis']))}
                        checklist.append(int(value['neddy']))
                        checklist1.append(int(value1['neddy']))
                        checklist2.append(key)
    neweddycount=1
    for key1, value1 in list(dict2.items()):
        if type(value1['neddy']!=checklist1) is np.ndarray:
            check=(value1['neddy']!=checklist1).all()
        else:
            check=(value1['neddy']!=checklist1)
        if check:
            dict1['eddyn_'+str(len(dict1)+neweddycount)]={'neddy':len(dict1)+neweddycount,\
                                            'time':value1['time']+timee,\
                                            'position':value1['position'],\
                                            'area':value1['area'],\
                                            'angle':value1['angle'],\
                                            'ellipse':value1['ellipse'],\
                                            'contour':value1['contour'],\
                                            'position_eddy':value1['position_eddy'],\
                                            'level':value1['level'],\
                                            'minoraxis':value1['minoraxis'],\
                                            'majoraxis':value1['majoraxis']}
            neweddycount=neweddycount+1
    return dict1