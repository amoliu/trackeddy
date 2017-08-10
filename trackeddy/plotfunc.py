import matplotlib as plt
import cmocean as cm

def basemap_mplot(x,y,data,scale='Lin',vmin='',vmax='',cmap=cm.cm.thermal,xan=1,yan=1,figsize=(20,20)):
    fig, ax = plt.subplots(xan, yan, figsize=figsize)

    X,Y=np.meshgrid(x,y)
        
    for ii in range(0,xan):
        for jj in range(0,yan):
            if xan==1 and yan==1:
                ttl=plt.title("SSHa", fontsize=20)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection='ortho',lat_0=-90,lon_0=-100,resolution='c')
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=plt
            elif xan!=1 and yan==1:
                ttl=ax[ii].set_title("SSHa", fontsize=20)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection='ortho',lat_0=-90,lon_0=-100,resolution='c',ax=axes[ii])
                lonm,latm=map(Lon,Lat)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[ii]
            elif xan==1 and yan!=1:
                ttl=ax[jj].set_title("SSHa", fontsize=20)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection='ortho',lat_0=-90,lon_0=-100,resolution='c',ax=axes[jj])
                lonm,latm=map(Lon,Lat)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[jj]
            else:
                ttl=ax[ii,jj].set_title("SSHa", fontsize=20)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection='ortho',lat_0=-90,lon_0=-100,resolution='c',ax=axes[ii,jj])
                lonm,latm=map(Lon,Lat)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                axes[0,0].pcolormesh(lonm,latm,ssha,vmin=-20,vmax=20)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[ii,jj]
                
            lonm,latm=map(X,Y)
            if scale =='Lin' and vmin == '' and vmax == '':
                vmin=data_masked.min()
                vmax=data_masked.max()
                im=m.pcolormesh(x,y,data_masked,cmap=cmap,vmin=vmin,vmax=vmax)
            elif scale == 'SymLog' and vmin == '' and vmax == '':
                vmin=data_masked.min()
                vmax=data_masked.max()
                im=m.pcolormesh(x,y,data_masked,cmap=cmap,norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=vmin, vmax=vmax))
            elif scale == 'Lin' and vmin != '' and vmax != '':
                im=m.pcolormesh(x,y,data_masked,cmap=cmap,vmin=vmin, vmax=vmax)
            elif scale == 'SymLog' and vmin != '' and vmax != '':
                im=m.pcolormesh(x,y,data_masked,cmap=cmap,norm=colors.SymLogNorm(linthresh=0.01, linscale=0.01,vmin=vmin, vmax=vmax))
            elif scale == 'Log' and vmin != '' and vmax != '':
                im=m.pcolormesh(x,y,abs(data_masked),cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax))
            elif scale == 'Log2+' and vmin != '' and vmax != '':
                lev_exp = np.arange(np.floor(np.log10(0.001)-1),np.ceil(np.log10(0.3)+1))
                levs = np.power(10, lev_exp)
                im=m.contourf(x,y,abs(data_masked),levs,cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax),alpha=0.5)
    return fig,im,cax,ax   