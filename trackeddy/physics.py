import numpy as np
import seawater as sw
from sympy.physics.vector import curl

def okuboweissparm(u,v,z):
    du_x=np.gradient(u[z,:,:],axis=1)
    dv_x=np.gradient(v[z,:,:],axis=1)
    du_y=np.gradient(u[z,:,:],axis=0)
    dv_y=np.gradient(v[z,:,:],axis=0)
    sn=du_x-dv_y
    ss=dv_x+du_y
    w=vorticity(u,v,z)
    owparm=sn**2+ss**2-w**2
    return owparm
    
def vorticity2D(u,v):
    dv_x=np.gradient(v[:,:],axis=1)
    du_y=np.gradient(u[:,:],axis=0)
    w=dv_x-du_y
    return w

    
#def vorticity3D(u,v,w,z):
#    dv_x=np.gradient(v[z,:,:],axis=1)
#    du_y=np.gradient(u[z,:,:],axis=0)
#    w=dv_x-du_y
#    return w

def geostrophicssh(eta,lat,lon):
    print('Work in progress')

def EkE(eta,u,v):
    print('Work in progress')
    
def PVort(S,T,P,U,V):
    theta=np.zeros(np.shape(S))
    for ii in range(0,np.shape(S)[1]):
        for jj in range(0,np.shape(S)[2]):
            theta[:,ii,jj]=sw.ptmp(S[:,ii,jj],T[:,ii,jj],P)
    rho=sw.dens(S,T,P)
    omega=7.2921159*10**(-5)
    zeta=2*omega+curl(U,V)
    Q=(1/rho)*(zeta)*np.gradient(theta)