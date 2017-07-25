def okuboweissparm(u,v,z):
    du_x=gradient(u[z,:,:],axis=1)
    dv_x=gradient(v[z,:,:],axis=1)
    du_y=gradient(u[z,:,:],axis=0)
    dv_y=gradient(v[z,:,:],axis=0)
    sn=du_x-dv_y
    ss=dv_x+du_y
    w=vorticity(u,v,z)
    owparm=sn**2+ss**2-w**2
    return owparm
    
def vorticity(u,v,z):
    dv_x=gradient(v[z,:,:],axis=1)
    du_y=gradient(u[z,:,:],axis=0)
    w=dv_x-du_y
    return w

def geostrophicssh(eta,lat,lon):
    print 'Work in progress'

def EkE(eta,u,v):
    print 'Work in progress'