import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_hannes as ld
import plot_nom_II as pt
import math
import mag

from PIL import Image
import glob

from scipy.optimize import curve_fit

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


class STEP():
    def __init__(self,year,month,day,rpath = '/media/tomruge/Toms_SSD/data/step_v0008/',rpath_mag = None,magnet_default_path=False,lastofmonth=False):
        '''Magnetfeld wird gleich mitgeladen, wenn rpath_mag übergeben wird. Wird magnet_default_path auf True gesetzt wird automatisch der 
        korrekte Dateipfad für die Uni-Rechner genutzt.'''
        self.ebins = ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])
        
        # Fehlt da aber nicht last of day?
        if lastofmonth:
            if month!=12:
                self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year,month+1,1)), products=('M','A'))
            else:
                self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year+1,1,1)), products=('M','A'))
        else:
            self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)), products=('M','A'))
        print('STEP-Data loaded successfully.')
        
        # Combining data (Main and Auxiliary Product), how does this product look like?
        ld.combine_data(self.itime, self.idata)
        print(self.idata)
        #print(self.itime)
        print('STEP-Data combined successfully.')
        
        # Loading MAG-Data
        if magnet_default_path == True:
            self.magnet_default_path = '/data/projects/solo/mag/l2_soar/rtn_1minute'
            rpath_mag = self.magnet_default_path

        if type(rpath_mag) == str:
            if lastofmonth:
                if month!=12:
                    self.mag = mag.MAGdata(path = rpath_mag, period = (dt.datetime(year,month,day),dt.datetime(year,month+1,1)))
                else:
                    self.mag = mag.MAGdata(path = rpath_mag, period = (dt.datetime(year,month,day),dt.datetime(year+1,1,1)))
            else:
                self.mag = mag.MAGdata(path = rpath_mag, period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)))

    def cut_data(self,t0,t1):
        # Was macht das?
        cdat = {}
        ctime = {}
        for p in self.idata.keys():
            m = (self.itime[p]>=t0) * (self.itime[p]<t1)
            ctime[p] = self.itime[p][m] 
            cdat[p] = self.idata[p][m] 
        return ctime, cdat

    def data_prep(self,ebins=ebins,res = '1min', head = 0, period = None, norm = False, overflow = True, esquare = False):
        '''Vorbereitung der STEP-Daten basierend auf Lars Skripten.'''
        vmax = None # Falls keine Normierung angegben wird, gebe ich None zurück. Was soll vmax eigentlich darstellen?

        if period:
            time,dat = self.cut_data(period[0]-dt.timedelta(seconds=59),period[1])
        else:
            time,dat = self.itime, self.idata
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

        return pldat, pltime, vmax
    
    def data_boxes(self,pldat,box_list):
        '''Nehme pldat und, um die Daten auf die Boxen einzuschränken. Dabei werden die Werte, die ignoriert werden sollen durch 0.0 ersetzt.
        box_data wird zurückgegeben.'''
        # Erster Index von pdat müsste die Zeitreihe sein, der zweite der Energie-Bin.
        # [[[timelow,timeup],[energylow,energyup]],[[timelow,timeup],[energylow,energyup]]]
        box_data = []
        for pdat in pldat:
            ayuda = np.zeros(shape=(len(pdat),len(pdat.T)),dtype='float')
            for box in box_list:
                tlow = box[0][0]
                tup = box[0][1]
                elow = box[1][0]
                eup = box[1][1]

                # Beschleunigung durch Nutzung von Masken???
                # Genauer ansehen, wenn ich in Laufzeitprobleme laufe...
                # Code auf jeden Fall schon so angepasst, dass data_boxes time bekommt
                # mask = ()

                # Loope erst durch die entsprechende Zeitreihe und dann jeweils durch die Energiebins
                for t in range(tlow,tup+1):
                    for e in range(elow,eup+1):
                        ayuda[t][e] = pdat[t][e]
            box_data.append(ayuda)
        return box_data
    
    def cut_straight_line(self,pldat,time_ind = (0,10),e_ind = (0,10)):
        '''Funktion, um Daten unterhalb einer Geraden in den 2d-Arrays der STEP-Daten zu entfernen. Ähnliche Funktionsweise, wie data_boxes().
        Werte ohne counts werden mit 0.0  gefüllt, um bei der Mittelwertsbildung egal zu sein. Vorsichtig sein!!! Funktion wird für Spezialfälle fehleranfällig sein!'''
        modified_data = []
        # Loope  durch die Pixel
        for pdat in pldat:
            ayuda = np.zeros(shape=(len(pdat),len(pdat.T)),dtype='float')
            # Hier aus den Start- und Endindizes für time und energy eine lineare Funktion basteln...
        
        modified_data.append(ayuda)
        return modified_data


    def step_plot(self,xlabel,ylabel,title):
        '''Erzeugt den bekannten STEP-Plot mit allen 16 Pixeln. Gibt Liste mit axes zurück.'''
        fig = plt.figure(figsize = (15,10))
        fig.subplots_adjust(wspace = 0, hspace = 0)
        ax = []        
        for i in range(16):
            if i == 0:
                ax.append(fig.add_subplot(4,5,3))
            else:
                if i == 1:
                    ax.append(fig.add_subplot(4,5,5+i))
                else:
                    ax.append(fig.add_subplot(4,5,5+i,sharex = ax[1], sharey = ax[1]))              
                if i == 6:
                    ax[i].set_ylabel(ylabel)
                if i == 13:
                    ax[i].set_xlabel(xlabel)
            if i not in [0,1,6,11]:
                for t in ax[i].get_yticklabels():
                    t.set_visible(False)
            if i < 11:
                for t in ax[i].get_xticklabels():
                    t.set_visible(False)
        ax[0].set_title(title)
        return fig, ax
    
    def plot_ts(self,ebins=ebins,res = '1min', head = 0, period = None, save = False, norm = False, overflow = True, esquare = False, box_list= False):
        '''Plottet die Zeitreihen der Energien. Übergebe verschachtelte Liste mit Grenzen für Boxen, die eingezeichnet werden. Grenzen als Indizes der Listen übergeben.'''
        pldat, pltime, vmax = self.data_prep(ebins,res,head,period,norm,overflow,esquare)
        fig, ax = self.step_plot('Date', 'Energy [keV]', 'Head %i'%head)
        for i in range(16):
            pdat = pldat[i]
            ptime = np.append(pltime[i],pltime[i][-1]+dt.timedelta(seconds=60))

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
                
            # PLots der Boxen, Grenzen als  Indizes übergeben.
            # Erster Index von pdat müsste die Zeitreihe sein, der zweite der Energie-Bin.
            # [[[timelow,timeup],[energylow,energyup]],[[timelow,timeup],[energylow,energyup]]]
            if type(box_list) == list:
                for box in box_list:
                    tlow = ptime[box[0][0]]
                    tup = ptime[box[0][1]]
                    elow = ebins[box[1][0]]
                    eup = ebins[box[1][1]]
                    ax[i].vlines([tlow,tup],elow,eup,color='firebrick')
                    ax[i].hlines([elow,eup],tlow,tup,color='firebrick')
                
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
            if i > 10:
                ax[i].tick_params(labelrotation=90)
            if period:
                ax[0].set_xlim(period[0],period[1])
        if save:
            if type(save) == str:
                plt.savefig(save + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig(save + 'head%i/'%head + 'TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
            else:
                plt.savefig('TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig('TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
        print('Plotted TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res) + ' successfully.')
        
    def calc_energy_means(self,ebins=ebins,res = '1min', head = -1, period = (dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)), box_list=None, window_width=5, pixel_list=[i for i in range(1,16)], norm = 'tmax', overflow = True, esquare = False):
        '''Wichtig: Die Energie ist in keV angegeben!!!
        Falls andere Zeitauflösung als Minuten gewählt wird, kann es Probleme mit window_width und der Berechnung der Zeit geben.'''
        i = 0
        pixel_means = [[] for i in range(16)]     # Liste mit Listen der Mittelwerte der einzelnen Pixel. Die erste Liste enthält die Zeitstempel (jeweils Mitte der Zeitfenster)
        
        if type(box_list) == list:
            little_helper_dat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]
            pldat = self.data_boxes(little_helper_dat,box_list)
        else:
            pldat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]

        while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
            pixel_means[0].append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
            # Berechnung der Mittelwerte:
            for k in pixel_list:
                pdat = pldat[k][i*window_width:(i+1)*window_width]
                integral = np.sum(pdat,axis=0)
                # calculating mean
                mean = 0.0
                for j in range(len(integral)):
                    # Für Bestimmung des Mittelwertes der Verteilung wird Mitte der Bins gewählt
                    mean += integral[j]*0.5*(ebins[j+1]+ebins[j])
                mean = mean/np.sum(integral)
                pixel_means[k].append(mean)
            i +=1
        
        for i in range(len(pixel_means)):
            # Array erstellen
            pixel_means[i] = np.array(pixel_means[i])
            
        return pixel_means
    
    def calc_sliding_energy_means(self,ebins=ebins,res = '1min', head = -1, period = (dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)), box_list=None, window_width=5, sliding_width=4, pixel_list=[i for i in range(1,16)], norm = 'tmax', overflow = True, esquare = False):
        '''Wichtig: Die Energie ist in keV angegeben!!!
        Falls andere Zeitauflösung als Minuten gewählt wird, kann es Probleme mit window_width und der Berechnung der Zeit geben.
        Nehme immer vier benachbarte energy-bins und berechne den Mittelwert (Sowohl der normierten counts und der Energie). Dieses Fenster gleitet über alle bins. 
        Nehme dann aus diesen Werten den maximalen Wert. Versuche so die Effekte der Daten-Boxen zu minimieren, weil die Tail-Länge an der Seite den Mittelwert beeinflusst.'''
        i = 0
        pixel_means = [[] for i in range(16)]     # Liste mit Listen der Mittelwerte der einzelnen Pixel. Die erste Liste enthält die Zeitstempel (jeweils Mitte der Zeitfenster)
        
        if type(box_list) == list:
            little_helper_dat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]
            pldat = self.data_boxes(little_helper_dat,box_list)
        else:
            pldat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]

        while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
            pixel_means[0].append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
            # Berechnung der Mittelwerte:
            for k in pixel_list:
                pdat = pldat[k][i*window_width:(i+1)*window_width]
                integral = np.sum(pdat,axis=0)
                # calculating mean
                # Breite des Sliding-windows sollte natürlicher Teiler der 56 energy-bins sein und wird in Zahl der Indizes angegeben.
                j_indizes=[i*sliding_width for i in range(int(len(integral)/sliding_width))]
                j_indizes=[i for i in range(0,len(integral)-sliding_width)]
                means_energy = [] # mittlere Energien
                means_normedcounts = [] # Mittelung der normierten Counts im Window
                for j in j_indizes:
                    # Für Bestimmung des Mittelwertes der Verteilung wird Mitte der Bins gewählt
                    mean = 0.0
                    count =  0.0
                    for l in range(sliding_width):
                        mean += integral[j+l]*0.5*(ebins[j+l+1]+ebins[j+l])
                        count += integral[j+l]*(ebins[j+l+1]-ebins[j+l])      # Multipliziere hier nochmal mit der Bin-Breite; Keine Normierung auf eins, damit ich nachher max, count wählen kann
                    mean = mean/np.sum(integral[j:j+sliding_width])
                    means_energy.append(mean)
                    means_normedcounts.append(count)
                total_mean = means_energy[means_normedcounts.index(max(means_normedcounts))]  # liefert mir Energiemittelwert für maximalen normierten Count beim Sliding-Window
                pixel_means[k].append(total_mean)
            i +=1
        
        for i in range(len(pixel_means)):
            # Array erstellen
            pixel_means[i] = np.array(pixel_means[i])
            
        return pixel_means
    
    def plot_energy_means(self,pixel_means, rpath, pixel_list=[i for i in range(1,16)]):
        '''Plottet die übergebenen Mittelwerte für alle Pixel.'''
        fig, ax = plt.subplots(figsize=(10,6))
        for i in pixel_list:
            ax.scatter(pixel_means[0],pixel_means[i],marker='x',label=f'pixel {i}')
        ax.set_ylabel('mean of energy [keV]')
        ax.set_xlabel('time')
        ax.tick_params(axis='x',labelrotation=45)
        plt.title('energy means')
        plt.legend()
        plt.savefig(rpath)
