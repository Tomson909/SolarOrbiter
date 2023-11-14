######################################################################
# Author : Lars Berger                                               #
# Date : 23.09.2020                                                  #
#                                                                    #
# Load SOLO/MAG data currently from CDF files provided by Tim H.     #
# directly, later CSV format might/should follow to allow working    #
# without use of a cdf interface                                     #
######################################################################

import numpy as np
import datetime as dt
import os
from cdflib import CDF,cdfepoch
import matplotlib.pyplot as plt

class MAGdata(object):
    def __init__(self, period = (dt.datetime(2022,1,1),dt.datetime(2022,1,2)), mag_path = 'default', frame = 'srf'):
        '''Nutze automatisch Level 2 Daten. Wenn ich einen Tag laden will, muss ich bei period den nächsten Tag als obere Grenze angeben.
        Ich kann angeben, ob ich die Daten im SRF oder in RTN (minütlich oder 8 Hz) laden will ('srf','rtn_min','rtn_8hz').
        Im SRF wird die x-Koordinate unter R, die y-Koordinate unter T und die z-Koordinate unter N gespeichert. Übergebe die Winkel immer
        in Radians an die Funktionen.
        Schreibe Loader so, dass nur einzelne Tage geladen werden können. Das macht das Laden der Flow-Vektoren deutlich einfacher.'''
    
        if mag_path == 'default':
            if frame == 'srf':
                self.path = '/data/projects/solo/mag/l2_soar/srf'
            if frame == 'rtn_min':
                self.path = '/data/projects/solo/mag/l2_soar/rtn_1minute'
            if frame == 'rtn_8hz':
                self.path = '/data/projects/solo/mag/l2_soar/rtn_high_time_res'
        else:
            self.path = mag_path
        self.period = period
        self.frame = frame          
        self.B_R = np.zeros((0))
        self.B_T = np.zeros((0))
        self.B_N = np.zeros((0))
        self.time = []
        self.flow_vector = np.array([[-0.8412,  0.4396,  0.3149], [-0.8743,  0.457,   0.1635], [-0.8862,  0.4632, -0.    ], [-0.8743,  0.457,  -0.1635], [-0.8412,  0.4396, -0.315 ], [-0.7775,  0.5444,  0.3149], [-0.8082,  0.5658,  0.1635], [-0.8191,  0.5736,  0.    ], [-0.8082,  0.5659, -0.1634], [-0.7775,  0.5444, -0.3149], [-0.7008,  0.6401,  0.3149], [-0.7284,  0.6653,  0.1634], [-0.7384,  0.6744, -0.    ], [-0.7285,  0.6653, -0.1635], [-0.7008,  0.6401, -0.315 ]])
        self.pitchangles = []
        self.pitchangles_err = []
        # Blickrichtung der einzelnen Pixel (Phi, Theta) (Wichtig es wird von der Sonne zum Spacecraft geschaut!!!):
        # from where does he have this information?
        self.fov = np.array([[155,24],[155,12],[155,0],[155,-12],[155,-24],[145,24],[145,12],[145,0],[145,-12],[145,-24],[135,24],[135,12],[135,0],[135,-12],[135,-24]])
        # Hier nochmal Blickrichtungen vom Spacecraft zur Sonne:
        fov_rtn = np.array([[-25,24],[-25,12],[-25,0],[-25,-12],[-25,-24],[-35,24],[-35,12],[-35,0],[-35,-12],[-35,-24],[-45,24],[-45,12],[-45,0],[-45,-12],[-45,-24]])
        
        self._load()
        # self._calc_angles()
        self._calc_pw()
        # self._calc_pw_err()
        
    def _load(self):
        d = self.period[0]
        while d < self.period[1]:
            y = d.year
            m = d.month
            dd = d.day
            if self.frame == 'rtn_8hz':
                l = os.listdir(self.path + '/%.4i/'%y)
                l.sort()
                ll = ['solo_L2_mag-rtn-normal_%.4i%.2i%.2i'%(y,m,dd) in s for s in l]
                fname =''
                for i,n in enumerate(ll):
                    if n:
                        fname = self.path + '/%.4i/'%(y) + l[i]

                print(d,fname)
            if self.frame == 'rtn_min':
                l = os.listdir(self.path + '/%.4i/'%y)
                l.sort()
                ll = ['solo_L2_mag-rtn-normal-1-minute_%.4i%.2i%.2i'%(y,m,dd) in s for s in l]
                fname =''
                for i,n in enumerate(ll):
                    if n:
                        fname = self.path + '/%.4i/'%(y) + l[i]

                print(d,fname)
            if self.frame == 'srf':
                l = os.listdir(self.path + '/%.4i/'%y)
                l.sort()
                ll = ['solo_L2_mag-srf-normal_%.4i%.2i%.2i'%(y,m,dd) in s for s in l]
                fname =''
                for i,n in enumerate(ll):
                    if n:
                        fname = self.path + '/%.4i/'%(y) + l[i]

                print(d,fname)
            try:
                dat = CDF(fname)
                t = dat.varget('EPOCH')
                t = cdfepoch.encode_tt2000(t)
                for tt in t:
                    self.time.append(dt.datetime.fromisoformat(str(tt)[:26]))
                    self.time[-1] = self.time[-1].replace(tzinfo = None)
                if self.frame == 'srf':
                    b = dat.varget('B_SRF')
                if self. frame == 'rtn_min' or self.frame == 'rtn_8hz':
                    b = dat.varget('B_RTN')
                self.B_R = np.append(self.B_R,np.array(b[:,0]))
                self.B_T = np.append(self.B_T,np.array(b[:,1]))
                self.B_N = np.append(self.B_N,np.array(b[:,2]))
            except:
                print('No MAG data : ', d)
            d = d + dt.timedelta(days = 1)
        self.time = np.array(self.time).flatten()
        self.B_R = np.array(self.B_R).flatten()
        self.B_T = np.array(self.B_T).flatten()
        self.B_N = np.array(self.B_N).flatten()

    def _calc_angles(self):
        '''Winkel Phi und Theta der Magnetfeldvektoren im Spacecraft-Frame'''
        self.phi = np.arctan2(self.B_T,self.B_R)
        self.theta = np.arctan2(self.B_N,np.sqrt(self.B_T**2+self.B_R**2))
        
    def pw(self,flow,B):
        '''Übergebe den particle flow-Vektor als Geschwindigkeit und den Magnetfeldvektor (am besten in SRF) und berechne die Pitchwinkel über das Skalarprodukt.'''
        len_flow = np.sqrt(flow[0]**2 + flow[1]**2 + flow[2]**2)
        len_B = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
        argument = (flow[0]*B[0] + flow[1]*B[1] + flow[2]*B[2])/len_flow/len_B
        result = np.arccos(argument)
        return result
        
    def _calc_pw(self):
        '''Berechne die Pitchwinkel für die Elektronen, welche auf STEP treffen in erster Näherung.
        Dafür wird der Winkel zwischen dem particle flow vector der Pixel und dem Magnetfeld herangezogen.'''
        for i in range(15):
            self.pitchangles.append(self.pw(self.flow_vector[i],np.array([self.B_R,self.B_T,self.B_N])))
        self.pitchangles = np.array(self.pitchangles)
        
        
    def pw_err(self,v_phi,v_theta):
        ''' Implementierung der Fehler der Pitchwinkel. Nehme erstmal an, dass der Fehler des absoluten Feldes von 0.2 nT
        nur in in B_N-Richtung auftritt. Die Magnetfelder sind in nT gegeben.'''
        B_N_err = 0.2
        B_theta_err = 1/(1+self.B_N*self.B_N/(self.B_T*self.B_T+self.B_R*self.B_R)) * 1/np.sqrt(self.B_T*self.B_T+self.B_R*self.B_R) * B_N_err
        v_theta = np.radians(v_theta)
        v_phi = np.radians(v_phi)
        B_theta = np.radians(self.theta)
        B_phi = np.radians(self.phi)
        
        pw_err = np.abs((np.sin(B_phi)*np.sin(v_phi)+np.cos(B_phi)*np.cos(v_phi))*np.cos(v_theta)*np.sin(B_theta)-np.sin(v_theta)*np.cos(B_theta)/np.sqrt(1-(np.sin(v_theta)*np.sin(B_theta)+np.sin(B_phi)*np.sin(v_phi)*np.cos(v_theta)*np.cos(B_theta)+np.cos(B_phi)*np.cos(v_phi)*np.cos(v_theta)*np.cos(B_theta))**2))*B_theta_err
        return pw_err
        
        
    def _calc_pw_err(self):
        '''Berechne die Fehler der Pitchwinkel mittels der Implementierung unter pw_err.''' 
        for i in range(15):
            v_phi = self.fov[i][0]
            v_theta = self.fov[i][1]
            self.pitchangles_err.append(self.pw_err(v_phi,v_theta))
        self.pitchangles_err = np.array(self.pitchangles_err)

    def average_pw(self,period,window_width=5):
        '''Berechnung der gemittelten Pitchwinkel'''
        # Maske, da ich nur die Magnetfelddaten innerhalb von period brauche:
        mask = (self.time > period[0]) * (self.time <= period[1])
        pw = [[] for i in range(15)]
        pw_time = []
        
        i = 0
        while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
            pw_time.append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
            for k in [i for i in range(1,16)]:
                # Mittelung der Pitchwinkel (k-1, da ich keine Zeit im array stehen habe)
                pw_data = self.pitchangles[k-1][mask]
                new_pw = np.sum(pw_data[i*window_width:(i+1)*window_width])/window_width
                pw[k-1].append(new_pw)
            i +=1
        return pw, pw_time
        
        
        
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
    
    
    def pw_ts_step_plot(self,err=False):
        '''Bei err == True werden Fehler mit dargestellt.'''

        fig, ax = self.step_plot(r'time', r'$\varphi$ [°]', 'Time series of pitch angle from ' + str(self.period[0]) + ' to ' + str(self.period[1]))

        for i in range(1,16):
            if err == True:
                ax[i].errorbar(self.time,np.degrees(self.pitchangles[i-1]),yerr=np.degrees(self.pitchangles_err[-1]),marker='x')
            else:
                ax[i].scatter(self.time,np.degrees(self.pitchangles[i-1]),marker='x')
            ax[i].tick_params(axis='x', labelrotation=90)
        plt.savefig('pw_ts_step_plot.png')
        plt.show()
        
        
    def pw_ts(self,rpath = '',pixel_list = range(16),period=None,window_width=None):
        '''Plottet die Pitchwinkel für alle angegebenen Pixel.'''
        fig, ax = plt.subplots(figsize=(10,6))

        if period != None and window_width == None:
            mask = (self.time >= period[0]) * (self.time < period[1])
            for i in pixel_list:
                ax.plot(self.time[mask],np.degrees(self.pitchangles[i-1][mask]),marker='x',label=f'pixel {i}')
        elif period != None and window_width != None:
            pw, pw_time =self.average_pw(period=period,window_width=window_width)
            for i in pixel_list:
                ax.plot(pw_time,np.degrees(pw[i-1]),marker='x',label=f'pixel {i}')
        else:
            for i in pixel_list:
                ax.plot(self.time,np.degrees(self.pitchangles[i-1]),marker='x',label=f'pixel {i}')

        ax.set_ylabel('pitch angle [°]')
        ax.set_xlabel('time')
        ax.tick_params(axis='x',labelrotation=45)
        plt.title('pitch angle')
        plt.legend()
        plt.savefig('baum.png')


    def mag_ts(self,rpath,period=None):
        '''Plottet das Magnetfeld'''
        fig, ax = plt.subplots(figsize=(10,6))

        if period != None:
            mask = (self.time >= period[0]) * (self.time < period[1])
            ax.plot(self.time[mask],self.B_R[mask],marker='x',label=r'$B_\mathrm{R}$')
            ax.plot(self.time[mask],self.B_T[mask],marker='x',label=r'$B_\mathrm{T}$')
            ax.plot(self.time[mask],self.B_N[mask],marker='x',label=r'$B_\mathrm{N}$')
        else:
            ax.plot(self.time,self.B_R,marker='x',label=r'$B_\mathrm{R}$')
            ax.plot(self.time,self.B_T,marker='x',label=r'$B_\mathrm{T}$')
            ax.plot(self.time,self.B_N,marker='x',label=r'$B_\mathrm{N}$')

        ax.set_ylabel('magnetic field component [nT]')
        ax.set_xlabel('time')
        ax.tick_params(axis='x',labelrotation=45)
        if self.frame == 'srf':
            ayuda = 'SRF'
        if self.frame == 'rtn_min' or self.frame == 'rtn_8hz':
            ayuda = 'RTN'
        plt.title(f'magnetic field in {ayuda}-coordinates')
        plt.legend()
        plt.savefig(rpath)
        
    def srf_ts(self,period,savepath):
        '''time series of magnetic field and pitch angles in srf'''
        mask = (period[0] <= self.time) * (period[1] > self.time)
        
        fig, (ax1, ax2) = plt.subplots(2,1,figsize=(10,10),sharex=True)
        ax1.plot(self.time[mask],self.B_R[mask],label='X')
        ax1.plot(self.time[mask],self.B_T[mask],label='Y')
        ax1.plot(self.time[mask],self.B_N[mask],label='Z')
        ax1.legend()
        ax1.set_ylabel('magnetic field component [nT]')
        ax1.set_title(f'Investigation of Magnetic Field: {period[0]} to {period[1]}')
    
        for i in [1,2,3,4,5]:
            ax2.plot(self.time[mask],np.degrees(self.pitchangles[i-1][mask]),label=f'pixel {i}')
        ax2.legend()
        ax2.set_xlabel('time')
        ax2.set_ylabel('pitch angle [°]')
        ax2.tick_params(axis='x', labelrotation=45)
    
        plt.subplots_adjust(hspace=0.001)
        plt.savefig(savepath)
        plt.close('all')