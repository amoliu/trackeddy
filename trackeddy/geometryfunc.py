import numpy as np

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