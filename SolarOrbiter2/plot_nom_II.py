import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])
cmap = mpl.cm.jet
cmap.set_under('w')
hmap = mpl.cm.seismic

def cut_data(dat,time,t0,t1):
    cdat = {}
    ctime = {}
    for p in dat.keys():
        m = (time[p]>=t0) * (time[p]<t1)
        ctime[p] = time[p][m] 
        cdat[p] = dat[p][m] 
    return ctime, cdat
        
def plot_ts(idat,itime,ebins=ebins,res = '1min', head = 0, period = None, save = False, norm = False, overflow = True, esquare = False):
    if period:
        time,dat = cut_data(idat,itime,period[0]-dt.timedelta(seconds=59),period[1])
    else:
        time,dat = itime, idat
    pldat = []
    pltime = []
    if type(res) == dt.timedelta:
        for i in range(16):
            if head in [0,1]:
                tmptime = [time['STEP_C_%i_%.2i'%(head,i)][0]]
                while tmptime[-1]<time['STEP_C_%i_%.2i'%(head,i)][-1]:
                    tmptime.append(tmptime[-1]+res)
            if head in [0,1]:
                tdat = 1.*dat['STEP_C_%i_%.2i'%(head,i)]
                ttime = time['STEP_C_%i_%.2i'%(head,i)]
                tmpdat = []
                for en in range(tdat.shape[1]):
                    H,t = np.histogram(ttime,bins = tmptime,weights = tdat[:,en])
                    tmpdat.append(H)
                tmpdat = np.array(tmpdat).T
                pldat.append(tmpdat)
                pltime.append(np.array(tmptime[:-1]))
            else:
                pldat.append(1.*dat['STEP_C_0_%.2i'%(i)]-dat['STEP_C_1_%.2i'%(i)])
    elif res == '1min':
        for i in range(16):
            if head in [0,1]:
                pldat.append(1.*dat['STEP_C_%i_%.2i'%(head,i)])
                pltime.append(time['STEP_C_%i_%.2i'%(head,i)])
            else:
                pldat.append(1.*dat['STEP_C_0_%.2i'%(i)]-dat['STEP_C_1_%.2i'%(i)])
                pltime.append(time['STEP_C_0_%.2i'%(i)])
            if esquare:
                pldat[-1]= pldat[-1]/(np.diff(ebins)**2)[np.newaxis,:]
                #pldat[-1]= pldat[-1]/ np.diff(ebins)[np.newaxis,:]
    elif res == '1s':
        for i in range(16):
            if head in [0,1]:
                pldat.append(dat['STEP_C_%i_%.2i'%(head,i)]/60.)
                pltime.append(time['STEP_C_%i_%.2i'%(head,i)])
            else:
                pldat.append(dat['STEP_C_0_%.2i'%(i)]/60.-dat['STEP_C_1_%.2i'%(i)]/60.)
                pltime.append(time['STEP_C_0_%.2i'%(i)])
            if esquare:
                pldat[-1]= pldat[-1]/(np.diff(ebins)**2)[np.newaxis,:]
        for i in range(16):
            if head in [0,1]:
                pldat.append(dat['STEP_M_%i_%.2i'%(head,i)])
                pltime.append(time['STEP_M_%i_%.2i'%(head,i)])
            else:
                pldat.append(dat['STEP_M_0_%.2i'%(i)]-dat['STEP_M_1_%.2i'%(i)])
                pltime.append(time['STEP_M_0_%.2i'%(i)])
            if esquare:
                pldat[-1]= pldat[-1]/(np.diff(ebins)[8:40]**2)[np.newaxis,:]
    if not overflow:
        for i in range(len(pldat)):
            pldat[i][:,-1] = 0.
    if norm == 'tmax':
        vmax = np.zeros(pltime[0].shape)
        for i in range(len(pldat)):
            vpmax = np.amax(pldat[i],axis=1)
            vmax[vpmax>vmax] = vpmax[vpmax>vmax]
    elif norm == 'ptmax' or norm == 'ptemax':
        vmax = np.zeros((16,pltime[0].shape[0]))
        for i in range(16):
            vpmax = np.amax(pldat[i],axis=1)
            vmax[i][vpmax>vmax[i]] = vpmax[vpmax>vmax[i]]
    elif norm == 'max' or norm == 'logmax':
        vmax = 0.
        for i in range(len(pldat)):
            vpmax = np.amax(pldat[i])
            vmax = max(vpmax,vmax)
    elif norm == 'pemax':
        vmax = np.zeros((16,pldat[i].shape[1]))
        for i in range(16):
            vpmax = np.amax(pldat[i],axis=0)
            vmax[i] = vpmax
    elif norm == 'emax':
        vmax = np.zeros(pldat[i].shape[1])
        for i in range(16):
            vpmax = np.amax(pldat[i],axis=0)
            vmax[vpmax>vmax] = vpmax[vpmax>vmax]
    ### Plotting
    fig = plt.figure(figsize = (15,10))
    fig.subplots_adjust(wspace = 0, hspace = 0)
    ax = []
    for i in range(16):
        pdat = pldat[i]
        ptime = np.append(pltime[i],pltime[i][-1]+dt.timedelta(seconds=60))
        if i == 0:
            ax.append(fig.add_subplot(4,5,3))
        else:
            ax.append(fig.add_subplot(4,5,5+i,sharex = ax[0], sharey = ax[0]))

        ax[i].set_yscale('log')
        if type(res) == dt.timedelta:
            if not norm:
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.1)
            elif norm == 'tmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
            elif norm == 'ptmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
            elif norm == 'ptemax':
                pdat2 = (pdat.T/vmax[i]).T
                vemax = np.amax(pdat2,axis=0)
                print(pdat.shape)
                print(vemax.shape)
                pdat2 = (pdat2/vemax)
                pdat2[np.isinf(pdat2)] = np.NaN
                tmp = ax[i].pcolormesh(ptime,ebins,pdat2.T, cmap = cmap, vmin = np.nanmin(pdat2),vmax = 1.)
            elif norm == 'max':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
            elif norm == 'logmax':
                tmp = ax[i].pcolormesh(ptime,ebins,np.log10(pdat).T, cmap = cmap, vmin = 0.,vmax = np.log10(vmax))
            elif norm == 'pemax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
            elif norm == 'emax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)

        elif res == '1min':
            if not norm:
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.1)
            elif norm == 'tmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
            elif norm == 'ptmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
            elif norm == 'ptemax':
                pdat2 = (pdat.T/vmax[i]).T
                vemax = np.amax(pdat2,axis=0)
                print(pdat.shape)
                print(vemax.shape)
                pdat2 = (pdat2/vemax)
                pdat2[np.isinf(pdat2)] = np.NaN
                tmp = ax[i].pcolormesh(ptime,ebins,pdat2.T, cmap = cmap, vmin = np.nanmin(pdat2),vmax = 1.)
            elif norm == 'max':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
            elif norm == 'logmax':
                tmp = ax[i].pcolormesh(ptime,ebins,np.log10(pdat).T, cmap = cmap, vmin = 0.,vmax = np.log10(vmax))
            elif norm == 'pemax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
            elif norm == 'emax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)

        elif res == '1s':
            if not norm:
                tmp = ax[i].pcolormesh(ptime,ebins[:9],pdat[:,:8].T, cmap = cmap, vmin = 0.1)
                tmp = ax[i].pcolormesh(ptime,ebins[40:],pdat[:,40:].T, cmap = cmap, vmin = 0.1)
            elif norm == 'tmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
            elif norm == 'ptmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
            elif norm == 'max':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
            elif norm == 'logmax':
                tmp = ax[i].pcolormesh(ptime,ebins[:9],np.log10(pdat[:,:8]).T, cmap = cmap, vmin = np.log10(1./60.),vmax = np.log10(vmax))
                tmp = ax[i].pcolormesh(ptime,ebins[40:],np.log10(pdat[:,40:]).T, cmap = cmap, vmin = np.log10(1./60.),vmax = np.log10(vmax))
            elif norm == 'pemax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
            elif norm == 'emax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
            pdat = pldat[i+16]
            ptime = np.append(pltime[i+16],pltime[i+16][-1]+dt.timedelta(seconds=1))
            if not norm:
                tmp = ax[i].pcolormesh(ptime,ebins[8:41],pdat.T, cmap = cmap, vmin = 0.1)
            elif norm == 'tmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
            elif norm == 'ptmax':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
            elif norm == 'max':
                tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
            elif norm == 'logmax':
                tmp = ax[i].pcolormesh(ptime,ebins[8:41],np.log10(pdat).T, cmap = cmap, vmin = np.log10(1./60.),vmax = np.log10(vmax))
            elif norm == 'pemax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
            elif norm == 'emax':
                tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
            
            
        if (norm and i == 0): # or not norm:
            tax = fig.add_subplot(4,50,31)
            if norm == 'tmax':
                plt.colorbar(tmp,cax = tax, label = 'Counts(t)/max(Counts(t)')
            elif norm == 'ptmax':
                plt.colorbar(tmp,cax = tax, label = 'Counts(t)/max(Counts(t)')
            elif norm == 'max':
                plt.colorbar(tmp,cax = tax, label = 'Counts')
            elif norm == 'logmax':
                plt.colorbar(tmp,cax = tax, label = 'Log10(C)')
            elif norm == 'pemax':
                plt.colorbar(tmp,cax = tax, label = 'C(E,pixel)/max(C(E,pixel)')
            elif norm == 'emax':
                plt.colorbar(tmp,cax = tax, label = 'C(E)/max(C(E)')

        ax[i].hlines(ebins[8],ptime[0],ptime[-1])
        ax[i].hlines(ebins[40],ptime[0],ptime[-1])
        if i == 0:
            tl = [l.get_text() for l in ax[i].get_xticklabels()]
            ax[i].set_title('Head %i'%head)
        ax[i].text(0.05,0.05,'%i'%(i),transform=ax[i].transAxes,color = 'k', backgroundcolor = 'w', fontsize = 6)
        if i not in [0,1,6,11]:
            for t in ax[i].get_yticklabels():
                t.set_visible(False)
            #ax[i].set_yticklabels([])
        if i < 11:
            for t in ax[i].get_xticklabels():
                t.set_visible(False)
            #ax[i].set_xticklabels([])
        if i > 10:
            ax[i].tick_params(labelrotation=90)
        if i == 13:
            ax[i].set_xlabel('Date')
        if i == 6:
            ax[i].set_ylabel('Energy [keV]')
    if period:
        ax[0].set_xlim(period[0],period[1])
    if save:
        if type(save) == str:
            plt.savefig(save + 'head%i/'%head + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
            #plt.savefig(save + 'head%i/'%head + 'TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
        else:
            plt.savefig('TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
            #plt.savefig('TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))


def plot_ts_nativ(dat,time,ebins=ebins,head = 0, save = False, norm = False, single = True):
    if norm == 'tmax':
        vmax = np.zeros(time['STEP_C_%i_00'%(head)].shape)
        for i in range(16):
            vpmax = np.amax(dat['STEP_C_%i_%.2i'%(head,i)],axis=1)
            vmax[vpmax>vmax] = vpmax[vpmax>vmax]
    elif norm == 'max':
        vmax = 0.
        for i in range(16):
            vpmax = np.amax(dat['STEP_C_%i_%.2i'%(head,i)])
            vmax = max(vpmax,vmax)
    if single:
        for i in range(16):
            fig = plt.figure(figsize = (15,10))
            ax = fig.gca()
            ax.set_yscale('log')
            if not norm:
                tmp = ax.pcolormesh(np.append(time['STEP_A_%i_%.2i'%(head,i)],time['STEP_A_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60)),ebins[:9],dat['STEP_A_%i_%.2i'%(head,i)][:,:8].T, cmap = cmap, vmin = 0.1)
                plt.colorbar(tmp, label = 'low')
                tmp = ax.pcolormesh(np.append(time['STEP_M_%i_%.2i'%(head,i)],time['STEP_M_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=1)),ebins[8:41],dat['STEP_M_%i_%.2i'%(head,i)].T, cmap = cmap, vmin = 0.1)
                plt.colorbar(tmp, label = 'mid')
                tmp = ax.pcolormesh(np.append(time['STEP_A_%i_%.2i'%(head,i)],time['STEP_A_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60)),ebins[40:],dat['STEP_A_%i_%.2i'%(head,i)][:,8:].T, cmap = cmap, vmin = 0.1)
                plt.colorbar(tmp, label = 'high')
            elif norm == 'tmax':
                tmp = ax.pcolormesh(np.append(time['STEP_C_%i_%.2i'%(head,i)],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60)),ebins,dat['STEP_C_%i_%.2i'%(head,i)].T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
                plt.colorbar(tmp)

            ax.hlines(ebins[8],time['STEP_C_%i_%.2i'%(head,i)][0],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60))
            ax.hlines(ebins[40],time['STEP_C_%i_%.2i'%(head,i)][0],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60))
            ax.set_title('Head %i, Pixel %i'%(head,i))
            if save:
                plt.savefig('TS_%i-%i_H%i_P%i.png'%(time['STEP_C_%i_%.2i'%(head,i)][0].day, time['STEP_C_%i_%.2i'%(head,i)][-1].day,head,i))
    else:
        fig = plt.figure(figsize = (15,10))
        fig.subplots_adjust(wspace = 0, hspace = 0)
        for i in range(16):
            if i == 0:
                ax = fig.add_subplot(4,5,3)
            else:
                ax = fig.add_subplot(4,5,5+i)
                
            ax.set_yscale('log')
            if not norm:
                tmp = ax.pcolormesh(np.append(time['STEP_C_%i_%.2i'%(head,i)],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60)),ebins,dat['STEP_C_%i_%.2i'%(head,i)].T, cmap = cmap, vmin = 0.1)
            elif norm == 'tmax':
                tmp = ax.pcolormesh(np.append(time['STEP_C_%i_%.2i'%(head,i)],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60)),ebins,dat['STEP_C_%i_%.2i'%(head,i)].T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
            elif norm == 'max':
                tmp = ax.pcolormesh(np.append(time['STEP_C_%i_%.2i'%(head,i)],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60)),ebins,dat['STEP_C_%i_%.2i'%(head,i)].T, cmap = cmap, vmin = 0.9,vmax = vmax)

            if (norm and i == 0) or not norm:
                tax = fig.add_subplot(4,50,31)
                plt.colorbar(tmp,cax = tax)
            ax.hlines(ebins[8],time['STEP_C_%i_%.2i'%(head,i)][0],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60))
            ax.hlines(ebins[40],time['STEP_C_%i_%.2i'%(head,i)][0],time['STEP_C_%i_%.2i'%(head,i)][-1]+dt.timedelta(seconds=60))
            if i == 0:
                ax.set_title('Head %i'%head)
            ax.text(0.05,0.05,'%i'%(i),transform=ax.transAxes,color = 'k', backgroundcolor = 'w', fontsize = 6)
            if i not in [0,1,6,11]:
                ax.set_yticklabels([])
            if i < 11:
                ax.set_xticklabels([])
            if i == 13:
                ax.set_xlabel('Date')
            if i == 6:
                ax.set_ylabel('Energy [keV]')
        if save:
            plt.savefig('TS_%i-%i_H%i_%s.png'%(time['STEP_C_%i_%.2i'%(head,i)][0].day, time['STEP_C_%i_%.2i'%(head,i)][-1].day,head,norm))
            

def plot_ts_diff(dat,time,ebins=ebins, save = False, single = True, bint = False):
    vmax = 0
    vmin = 0
    diffdat = np.zeros((16,dat['STEP_C_0_00'].shape[0],dat['STEP_C_0_00'].shape[1]))
    difft = np.zeros((16,dat['STEP_C_0_00'].shape[0]),dtype = dt.datetime)
    for i in range(16):
        diffdat[i] = (dat['STEP_C_0_%.2i'%i]-dat['STEP_C_1_%.2i'%i])
        difft[i] = time['STEP_C_0_%.2i'%i]
    if bint:
        print(bint)
        tb = [time['STEP_C_0_00'][0]]
        while tb[-1] < difft[0][-1]:
            tb.append(tb[-1]+dt.timedelta(seconds=bint))
        tb = np.array(tb)
        tdat = np.zeros((16,tb.shape[0]-1,diffdat.shape[2]))
        for i in range(16):
            for j in range(diffdat.shape[2]):
                H,x = np.histogram(difft[i],tb,weights = diffdat[i,:,j])
                tdat[i,:,j] = H
        diffdat = tdat
        difft = np.zeros((16,tb.shape[0]-1),dtype = dt.datetime)
        for i in range(16):
            difft[i] = tb[:-1]
    for i in range(16):
        vpmax = np.amax(diffdat[i])
        vmax = max(vpmax,vmax)
        vpmin = np.amin(diffdat[i])
        vmin = min(vpmin,vmin)
    vlim = max(abs(vmin),abs(vmax))
    if single:
        for i in range(16):
            vmax = np.amax(diffdat[i])
            vpmin = np.amin(diffdat[i])
            vlim = max(abs(vmax),abs(vmin))
            fig = plt.figure(figsize = (15,10))
            ax = fig.gca()
            ax.set_yscale('log')
            tmp = ax.pcolormesh(np.append(difft[i],difft[i][-1]+dt.timedelta(seconds=max(60,bint))),ebins,diffdat[i].T, cmap = hmap, vmin = -vlim, vmax = vlim)
            plt.colorbar(tmp)
            ax.hlines(ebins[8],difft[i][0],difft[i][-1]+dt.timedelta(seconds=max(60,bint)))
            ax.hlines(ebins[40],difft[i][0],difft[i][-1]+dt.timedelta(seconds=max(60,bint)))
            ax.set_title('Pixel %i'%i)
            if type(save) == str:
                plt.savefig(save + 'DIFFTS_%i-%i_P%i.png'%(difft[i][0].day, difft[i][-1].day,i))
            else:
                plt.savefig('DIFFTS_%i-%i_P%i.png'%(difft[i][0].day, difft[i][-1].day,i))
    else:
        fig = plt.figure(figsize = (15,10))
        fig.subplots_adjust(wspace = 0, hspace = 0)
        for i in range(16):
            if i == 0:
                ax = fig.add_subplot(4,5,3)
            else:
                ax = fig.add_subplot(4,5,5+i)
                
            ax.set_yscale('log')
            print(difft.shape,diffdat.shape)
            print(difft[i])
            tmp = ax.pcolormesh(np.append(difft[i],difft[i][-1]+dt.timedelta(seconds=max(60,bint))),ebins,diffdat[i].T, cmap = hmap, vmin = -vlim, vmax = vlim)

            if (i == 0):
                tax = fig.add_subplot(4,50,31)
                plt.colorbar(tmp,cax = tax)
            ax.hlines(ebins[8],difft[i][0],difft[i][-1]+dt.timedelta(seconds=max(60,bint)))
            ax.hlines(ebins[40],difft[i][0],difft[i][-1]+dt.timedelta(seconds=max(60,bint)))
            ax.text(0.05,0.05,'%i'%(i),transform=ax.transAxes,color = 'k', backgroundcolor = 'w', fontsize = 6)
            if i not in [0,1,6,11]:
                ax.set_yticklabels([])
            if i < 11:
                ax.set_xticklabels([])
            if i == 13:
                ax.set_xlabel('Date')
            if i == 6:
                ax.set_ylabel('Energy [keV]')
        if type(save) == str:
            plt.savefig(save + 'diff/DIFFTS_%i-%i_all.png'%(difft[i][0].day, difft[i][-1].day))
        else:
            plt.savefig('DIFFTS_%i-%i_all.png'%(difft[i][0].day, difft[i][-1].day))


            
def plot_sumspec(dat,time,ebins=ebins,head = 0, save = False,single = False):
    maxv = 0 
    for i in range(16):
        if np.amax(dat['STEP_C_%i_%.2i'%(head,i)].sum(axis=0)) > maxv:
            maxv = np.amax(dat['STEP_C_%i_%.2i'%(head,i)].sum(axis=0))

    if single:
        for i in range(16):
            fig = plt.figure(figsize = (15,10))                    
            ax = fig.gca()                                     
            ax.set_yscale('log')
            ax.set_xscale('log')                                                                                                                                  
            ax.bar(ebins[:-1]+np.diff(ebins)/2. , dat['STEP_C_%i_%.2i'%(head,i)].sum(axis=0), width = np.diff(ebins))
            ax.vlines(ebins,0.5,maxv*1.1,colors ='k')                                                       
            ax.set_ylim(0.5,maxv*1.1)                                                                        
            ax.set_title('Head %i, Pixel %i'%(head,i))
            if save:
                plt.savefig('HIST_%i-%i_H%i_P%i.png'%(time['STEP_C_%i_%.2i'%(head,i)][0].day, time['STEP_C_%i_%.2i'%(head,i)][-1].day,head,i))
    else:
        fig = plt.figure(figsize = (15,10))
        fig.subplots_adjust(wspace = 0, hspace = 0)
        ax = []
        for i in range(16):
            if i == 0:
                ax.append(fig.add_subplot(4,5,3))
                ax[i].set_title('Head %i'%head)
            else:
                ax.append(fig.add_subplot(4,5,5+i))
                
            ax[i].set_yscale('log')
            ax[i].set_xscale('log')                                                                                                                                  
            ax[i].bar(ebins[:-1]+np.diff(ebins)/2. , dat['STEP_C_%i_%.2i'%(head,i)].sum(axis=0), width = np.diff(ebins))
            ax[i].vlines(ebins,0.5,maxv*1.1,colors ='k', linewidths = 0.25)                                                       
            ax[i].set_ylim(0.5,maxv*1.1)                                                                        
            ax[i].text(0.1,0.9,'%i'%(i),transform=ax[i].transAxes,color = 'k', backgroundcolor = 'w', fontsize = 6)
            if i not in [0,1,6,11]:
                ax[i].set_yticklabels([])
            if i < 11:
                ax[i].set_xticklabels([])
            if i == 13:
                ax[i].set_xlabel('Energy [keV]')
            if i == 6:
                ax[i].set_ylabel('Counts')
        if save:
            plt.savefig('HIST_%i-%i_H%i_all.png'%(time['STEP_C_%i_%.2i'%(head,i)][0].day, time['STEP_C_%i_%.2i'%(head,i)][-1].day,head))
