import numpy as np

### NEED 2 FIX THIS FUNCTION #########
def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.ones([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    #print S
    E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a
######################################

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.ones([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    #print S\n",
    try:
        E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
    except:
        aa=[[  2.46360347e+10,  8.86065412e+09,   3.18685354e+09,  -1.36178743e+08, -4.89784392e+07,   7.52745780e+05],
         [  8.86065412e+09,   3.18685354e+09,   1.14619832e+09,  -4.89784392e+07, -1.76157812e+07,   2.70735119e+05],
         [  3.18685354e+09,   1.14619832e+09,   4.12248210e+08, -1.76157812e+07, -6.33578186e+06,   9.73738194e+04],
         [ -1.36178743e+08,  -4.89784392e+07,  -1.76157812e+07,   7.52745780e+05, 2.70735119e+05,  -4.16090533e+03],
         [ -4.89784392e+07,  -1.76157812e+07,  -6.33578186e+06,   2.70735119e+05, 9.73738194e+04,  -1.49652765e+03],
         [  7.52745780e+05,   2.70735119e+05,   9.73738194e+04,  -4.16090533e+03, -1.49652765e+03,   2.30000000e+01]]
        E, V =  np.linalg.eig(np.dot(np.linalg.inv(aa), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    #print up, down1,down2
    res1=np.sqrt(abs(up)/abs(down1))
    res2=np.sqrt(abs(up)/abs(down2))
    return np.array([res1, res2])

def PolyArea(x,y):
    area=0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
    return area

def eccentricity(a,b):
    if b>a:
        b1=a
        a=b
        b=b1
    eccen=np.sqrt(1-(b**2/a**2))
    #print sqrt(a**2-b**2)
    return eccen

def find2d(arrayx,arrayy,valuex,valuey):
    idx=(np.abs(arrayx-valuex)).argmin()
    idy=(np.abs(arrayy-valuey)).argmin()
    return idx,idy

def find(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx

def contourmaxvalue(contcoordx,contcoordy,var,x,y,levels,date):
    idxcheckmax,idycheckmax=find2d(x,y,contcoordx.max(),contcoordy.max())
    idxcheckmin,idycheckmin=find2d(x,y,contcoordx.min(),contcoordy.min())
    #print(idycheckmin,idycheckmax,idxcheckmin,idxcheckmax)
    if len(np.shape(var))==3:
        if levels[0]>0:
            sshextrem=np.nanmax(var[date,idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        else:
            sshextrem=np.nanmin(var[date,idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        indexes=np.where(var[date,idycheckmin:idycheckmax,idxcheckmin:idxcheckmax]==sshextrem)
    else:
        #print(np.shape(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax]))
        if levels[0]>0:
            sshextrem=np.nanmax(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        else:
            sshextrem=np.nanmin(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax])
        #print(sshextrem)
        indexes=np.where(var[idycheckmin:idycheckmax,idxcheckmin:idxcheckmax]==sshextrem)
    coord=[x[idxcheckmin+indexes[1][0]],y[idycheckmin+indexes[0][0]]]
    return coord