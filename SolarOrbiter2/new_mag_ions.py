import numpy as np
import matplotlib.pyplot as plt
from cdflib import *
from cdflib import cdfepoch
import datetime as dt
import os
import pandas as pd
# use seaborn for plots
import seaborn as sns
from fix_download_swa import _load_swa_from_SolOArchive
# for loading data from SolO Archive
import sunpy_soar
from sunpy.net import Fido
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Product,Identifier
# exclude warnings, maybe I should not do this, but I don't know how to fix it
import warnings

#install astropy
from astropy.time import Time
warnings.filterwarnings("ignore",category=FutureWarning)
warnings.filterwarnings("ignore",category=RuntimeWarning)
#import constants
from scipy.constants import e, m_e, m_p
# import data import
from load_nom_II import load_nom, combine_data

class CORRECT_PIXEL_AGING:
    def __init__(self,reference_pixel = 1,period=(dt.datetime(2022, 11, 18), dt.datetime(2022, 11, 19)),pixels=range(15),heads = (0,1)) -> None:
        self.period = period
        self.pixels = pixels
        self.heads = heads
        self.reference_pixel = reference_pixel
        #load all necessary data
        self.dat_SWA = SWA(period=period,pixels=pixels)
        self.dat_MAG = MAG(period=period,pixels=pixels)
        self.dat_STEP = STEP(period=period,pixels=pixels,head=heads)
        # Loading Step data is not necessary because 
        # correct the velocity vector
        self.correct_velocity()
        self.calc_pw()

    # could do this with one function but i am too lazy
    def average_sunwind_velocity_for_vel_vec(self):
        # the timestamps of the velocity vector are not the same as the timestamps of the SWA data
        # therefore I have to average the SWA data for the timestamps of the velocity vector
        # the fastest would be with a histogram
        v_x,_ = np.histogram(self.dat_SWA.V_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'],weights=self.dat_SWA.V_SRF['Vx'])
        v_y,_ = np.histogram(self.dat_SWA.V_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'],weights=self.dat_SWA.V_SRF['Vy'])
        v_z,_ = np.histogram(self.dat_SWA.V_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'],weights=self.dat_SWA.V_SRF['Vz'])
        n_x,_ = np.histogram(self.dat_SWA.V_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'])
        n_y,_ = np.histogram(self.dat_SWA.V_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'])
        n_z,_ = np.histogram(self.dat_SWA.V_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'])

        self.dat_SWA.V_SRF = pd.DataFrame(columns=['time','Vx','Vy','Vz'])
        self.dat_SWA.V_SRF['time'] = self.dat_STEP.velocity_vector[f'{0}_time'][:-1]
        self.dat_SWA.V_SRF['Vx'] = v_x/n_x
        self.dat_SWA.V_SRF['Vy'] = v_y/n_y
        self.dat_SWA.V_SRF['Vz'] = v_z/n_z
        #the averaging works perfectly 


    def average_magnetic_field_for_vel_vec(self):
        print('mag data before averaging:',self.dat_MAG.B_SRF)
        # the timestamps of the velocity vector are not the same as the timestamps of the SWA data
        B_x,_ = np.histogram(self.dat_MAG.B_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'],weights=self.dat_MAG.B_SRF['Bx'])
        B_y,_ = np.histogram(self.dat_MAG.B_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'],weights=self.dat_MAG.B_SRF['By'])
        B_z,_ = np.histogram(self.dat_MAG.B_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'],weights=self.dat_MAG.B_SRF['Bz'])
        n_x,_ = np.histogram(self.dat_MAG.B_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'])
        n_y,_ = np.histogram(self.dat_MAG.B_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'])
        n_z,_ = np.histogram(self.dat_MAG.B_SRF['time'],bins=self.dat_STEP.velocity_vector[f'{0}_time'])

        self.dat_MAG.B_SRF = pd.DataFrame(columns=['time','Bx','By','Bz'])
        self.dat_MAG.B_SRF['time'] = self.dat_STEP.velocity_vector[f'{0}_time'][:-1]
        self.dat_MAG.B_SRF['Bx'] = B_x/n_x
        self.dat_MAG.B_SRF['By'] = B_y/n_y
        self.dat_MAG.B_SRF['Bz'] = B_z/n_z
        #the averaging works perfectly
    
    def correct_velocity(self):
        # subtract SWA wind from velocity vector
        print('correct velocity vector...')

        self.bin_vel_vector = pd.DataFrame()
        for i in self.pixels:
            # but use only the first SWA timestamp just to see whether it wokrs as expected
            self.bin_vel_vector[f'{i}_x'] = self.dat_STEP.bin_vel_vector[f'{i}_x'] - self.dat_SWA.V_SRF['Vx'][100]
            self.bin_vel_vector[f'{i}_y'] = self.dat_STEP.bin_vel_vector[f'{i}_y'] - self.dat_SWA.V_SRF['Vy'][100]
            self.bin_vel_vector[f'{i}_z'] = self.dat_STEP.bin_vel_vector[f'{i}_z'] - self.dat_SWA.V_SRF['Vz'][100]
        # plot the corrected velocity vectors with scatterplot in 3d space
        '''fig = plt.figure(figsize=(20, 10), dpi=300)
        ax = fig.add_subplot(111, projection="3d")
        for i in self.pixels:
            if i == 1:
                ax.scatter(
                    self.bin_vel_vector[f'{i}_x'],
                    self.bin_vel_vector[f'{i}_y'],
                    self.bin_vel_vector[f'{i}_z'],
                    label=f"Pixel {i} corrected",
                )
                # also plot the bin velocity vector from the step data
                ax.scatter(
                    self.dat_STEP.bin_vel_vector[f'{i}_x'],
                    self.dat_STEP.bin_vel_vector[f'{i}_y'],
                    self.dat_STEP.bin_vel_vector[f'{i}_z'],
                    label=f"Pixel {i}",
                )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        # also turn opn legend but tiny
        ax.legend(self.pixels,fontsize='xx-small')
        plt.title("Velocity vectors for each energy bin")
        plt.show()
        '''




    def angle_between_vectors(self,v_vel,v_mag):
        # calculate the angle between two vectors
        len_v_vel = np.sqrt(v_vel[0]**2 + v_vel[1]**2 + v_vel[2]**2)
        len_v_mag = np.sqrt(v_mag[0]**2 + v_mag[1]**2 + v_mag[2]**2)
        # calculate the dot product
        dot = abs(v_vel[0]*v_mag[0] + v_vel[1]*v_mag[1] + v_vel[2]*v_mag[2])
        # calculate the angle
        angle = np.arccos(dot/(len_v_vel*len_v_mag))
        return angle

    def calc_pw(self):
        # calculate the angle between the velocity vector and the magnetic field
        # the angle between v and B is the pitch angle
        print('Calculating pitch angle...')
        print('Sunwind velocity vector:',self.dat_SWA.V_SRF)
        print('Magnetic field vector:',self.dat_MAG.B_SRF)
        self.pitchangles = pd.DataFrame()

        
        for i in self.pixels:
        # calculate the angle between the velocity vector and the magnetic field
            v_vel = np.array([self.bin_vel_vector[f'{i}_x'],self.bin_vel_vector[f'{i}_y'],self.bin_vel_vector[f'{i}_z']])
            v_mag = np.array([self.dat_MAG.B_SRF['Bx'][0],self.dat_MAG.B_SRF['By'][0],self.dat_MAG.B_SRF['Bz'][0]])
            self.pitchangles[f'{i}_data'] = np.degrees(self.angle_between_vectors(v_vel[:,:-1] ,v_mag))
        # print energy bins velocity before correction
        print('Energy bins velocity before correction:',self.dat_STEP.bin_vel_vector)
        # print energy bins velocity after correction
        print('Energy bins velocity after correction:',self.bin_vel_vector)
        # plotte die pitchwinkel, auf der x den Pitchwinkel und auf der y die Energie
        fig, ax = plt.subplots(figsize=(20, 10), dpi=300)
        for i in self.pixels:
            ax.plot(self.pitchangles[f'{i}_data'],self.dat_STEP.primary_energy[:-1]*1000)
        ax.set_xlabel('Energy in keV')
        ax.set_ylabel('Pitch angle')
        # tiny legend
        ax.legend(self.pixels,fontsize='xx-small')
        plt.title('Pitch angle')
        plt.show()
        # everything seems to work, so apply it correctly to all timeintervalls and all pixels
  
    def plot_pitchangle(self):
        # plot the pitch angle
        fig, ax = plt.subplots(figsize=(20, 10), dpi=300)
        for i in self.pixels:
            ax.plot(self.pitchangles[f'{i}_time'],self.pitchangles[f'{i}_data'])
        ax.set_xlabel('Time')
        ax.set_ylabel('Pitch angle')
        ax.legend(self.pixels)
        plt.title('Pitch angle')
        plt.savefig(f'pitchangle_{self.period[0].strftime("%Y-%m-%d")}_{self.period[1].strftime("%Y-%m-%d")}.png')

        

class MAG(CORRECT_PIXEL_AGING):
    def __init__(
            self, period=(dt.datetime(2022, 12, 12), dt.datetime(2022, 12, 13)),pixels=range(15)
        ) -> None:
            # load only 1 day for now
            self.period = period
            self.mag_path = "/media/tomruge/Toms_SSD/data/Mag"
            self.frame = "SRF"
            # specify pixels which will be considered in plots
            self.pixels = pixels
            # init array of magnetic field vectors in spacecraft reference frame
            self.B_SRF = pd.DataFrame()
            # init array of pitch angles
            self.pitchangles = []         
            # init array of flow vectors, the flow vectors are the velocity vectors of the particles hiting each ppixel, which is beacause of geometry the case
            # the flow vectors are in SRF?
            self.flow_vector = np.array(
                [
                    [-0.8412, 0.4396, 0.3149],
                    [-0.8743, 0.457, 0.1635],
                    [-0.8862, 0.4632, -0.0],
                    [-0.8743, 0.457, -0.1635],
                    [-0.8412, 0.4396, -0.315],
                    [-0.7775, 0.5444, 0.3149],
                    [-0.8082, 0.5658, 0.1635],
                    [-0.8191, 0.5736, 0.0],
                    [-0.8082, 0.5659, -0.1634],
                    [-0.7775, 0.5444, -0.3149],
                    [-0.7008, 0.6401, 0.3149],
                    [-0.7284, 0.6653, 0.1634],
                    [-0.7384, 0.6744, -0.0],
                    [-0.7285, 0.6653, -0.1635],
                    [-0.7008, 0.6401, -0.315],
                ]
            )
            # load magnetic fiield data
            self.load_mag()
                      
    def load_mag(self):
        print("Loading Mag data...")
        # load the timeinterval
        # create a array of dates for each day in period, period is a tuple of datetime objects
        dates = np.array([])
        current_date = self.period[0]
        while current_date < self.period[1]:
            dates = np.append(dates,current_date) 
            current_date += dt.timedelta(days=1)
        files = ["solo_L2_mag-srf-normal_"+i.strftime("%Y%m%d")+"_V01.cdf" for i in dates]
        #load the data
        self.B_SRF = pd.DataFrame(columns=['time','Bx','By','Bz'])
        for file in files:
            dat = CDF(self.mag_path + "/" + file)
            # load timestamps intp array
            dummy = pd.DataFrame(columns=['time','Bx','By','Bz'])
            #dummy['time'] = np.array([cdfepoch.to_datetime(ts) for ts in dat.varget("Epoch")])
            dummy['time'] = cdfepoch.to_datetime(np.array(dat.varget("Epoch")))
            # Space craft reference frame coordinates
            magnetic_field_srf = dat.varget("B_SRF")
            dummy['Bx'] = magnetic_field_srf[:, 0]
            dummy['By'] = magnetic_field_srf[:, 1]
            dummy['Bz'] = magnetic_field_srf[:, 2]
            dummy.fillna(0, inplace=True)
            self.B_SRF = pd.concat([self.B_SRF,dummy])        
        
class SWA(CORRECT_PIXEL_AGING):
    def __init__(self, period=(dt.datetime(2022, 12, 12), dt.datetime(2022, 12, 13)),pixels=range(15)) -> None:
        # load only 1 day for now
        self.period = period
        self.swa_path = "/media/tomruge/Toms_SSD/data/SWA_Data/"
        self.frame = "SRF"
        # specify pixels which will be considered in plots
        self.pixels = pixels
        # init array of pitch angles
        self.pitchangles = []
        # init pd array for velocity vectors and time

        # init array of flow vectors, the flow vectors are the velocity vectors of the particles hiting each ppixel, which is beacause of geometry the case
        # the flow vectors are in SRF?
        self.flow_vector = np.array(
                [
                    [-0.8412, 0.4396, 0.3149],
                    [-0.8743, 0.457, 0.1635],
                    [-0.8862, 0.4632, -0.0],
                    [-0.8743, 0.457, -0.1635],
                    [-0.8412, 0.4396, -0.315],
                    [-0.7775, 0.5444, 0.3149],
                    [-0.8082, 0.5658, 0.1635],
                    [-0.8191, 0.5736, 0.0],
                    [-0.8082, 0.5659, -0.1634],
                    [-0.7775, 0.5444, -0.3149],
                    [-0.7008, 0.6401, 0.3149],
                    [-0.7284, 0.6653, 0.1634],
                    [-0.7384, 0.6744, -0.0],
                    [-0.7285, 0.6653, -0.1635],
                    [-0.7008, 0.6401, -0.315],
                ]
            )
            # load solar wind velocity data
        self._load_swa()

    def time_interval_list(self, period):
        dummy = []  
        for i in range(int(np.ceil((period[1]-period[0])/dt.timedelta(days=1)))):
            dummy.append(period[0]+dt.timedelta(days=i))
        return dummy 

    def _load_swa(self):
        """Check if the data from SolO Archive is already downloaded."""
        print("Loading SWA data...")
        dates = self.time_interval_list(self.period)
        # check if file already exists
        self.V_SRF = pd.DataFrame(columns=['time','Vx','Vy','Vz'])
        for date in dates:
            path = (
                self.swa_path
                + "/"
                + "/solo_L2_swa-pas-grnd-mom_"
                + date.strftime("%Y%m%d")
                + "_V02"
                + ".cdf"
            )
            if not os.path.isfile(path):
                _load_swa_from_SolOArchive(date)
            dat = CDF(path)
            # load velocity data and save it in arrays
            velocity_data = dat.varget("V_SRF")
            #load timestamps
            dummy = pd.DataFrame(columns=['time','Vx','Vy','Vz'])
            dummy['time'] = cdfepoch.to_datetime(np.array(dat.varget("Epoch")))
            dummy['Vx'] = velocity_data[:, 0]*1000
            dummy['Vy'] = velocity_data[:, 1]*1000
            dummy['Vz'] = velocity_data[:, 2]*1000
            dummy.fillna(0, inplace=True)
            self.V_SRF = pd.concat([self.V_SRF,dummy])
            #solo_L2_swa-pas-grnd-mom_20211208_V02.cdf
            
class STEP(CORRECT_PIXEL_AGING):
    def __init__(
        self,
        pixels=range(15),
        head=(0, 1),
        period=(dt.datetime(2021, 10, 22), dt.datetime(2021, 10, 23)),
        rpath="/media/tomruge/Toms_SSD/data/step_v0008/",
    ):
        self.ebins = np.array(
            [
                0.98,
                2.144,
                2.336,
                2.544,
                2.784,
                3.04,
                3.312,
                3.6,
                3.92,
                4.288,
                4.672,
                5.088,
                5.568,
                6.08,
                6.624,
                7.2,
                7.84,
                8.576,
                9.344,
                10.176,
                11.136,
                12.16,
                13.248,
                14.4,
                15.68,
                17.152,
                18.688,
                20.352,
                22.272,
                24.32,
                26.496,
                28.8,
                31.36,
                34.304,
                37.376,
                40.704,
                44.544,
                48.64,
                52.992,
                57.6,
                62.72,
                68.608,
                74.752,
                81.408,
                89.088,
                97.28,
                105.984,
                115.2,
                125.44,
                137.216,
                149.504,
                162.816,
                178.176,
                194.56,
                211.968,
                230.4,
                372.736,
            ]
        )
        # these are the primary energy bins for the main channel
        self.primary_energy = np.array([0.005725, 0.006103, 0.006521, 0.006998, 0.007548, 0.008154,
                                        0.008802, 0.009491, 0.01026 , 0.01109 , 0.01198 , 0.01297 ,
                                        0.01405 , 0.01516 , 0.01655 , 0.01786 , 0.01933 , 0.02092 ,
                                        0.02261 , 0.02446 , 0.02652 , 0.0287  , 0.03101 , 0.03346 ,
                                        0.03618 , 0.0394  , 0.04254 , 0.04589 , 0.04977 , 0.0539  ,
                                        0.05826 , 0.06288 ])
        self.primary_energy_bin_width = np.array([0.003233, 0.003021, 0.002852, 0.002754, 0.002682, 0.002626,
                                                0.002598, 0.002599, 0.00261 , 0.00261 , 0.0026  , 0.00269 ,
                                                    0.00278 , 0.00288 , 0.00272 , 0.00272 , 0.00278 , 0.00279 ,
                                                    0.00281 , 0.0029  , 0.00295 , 0.003   , 0.00305 , 0.00319 ,
                                                    0.00322 , 0.00314 , 0.00335 , 0.00388 , 0.00413 , 0.00436 ,
                                                    0.00462 , 0.00512 ])
        self.flow_vector = np.array(
                [
                    [-0.8412, 0.4396, 0.3149],
                    [-0.8743, 0.457, 0.1635],
                    [-0.8862, 0.4632, -0.0],
                    [-0.8743, 0.457, -0.1635],
                    [-0.8412, 0.4396, -0.315],
                    [-0.7775, 0.5444, 0.3149],
                    [-0.8082, 0.5658, 0.1635],
                    [-0.8191, 0.5736, 0.0],
                    [-0.8082, 0.5659, -0.1634],
                    [-0.7775, 0.5444, -0.3149],
                    [-0.7008, 0.6401, 0.3149],
                    [-0.7284, 0.6653, 0.1634],
                    [-0.7384, 0.6744, -0.0],
                    [-0.7285, 0.6653, -0.1635],
                    [-0.7008, 0.6401, -0.315],
                ]
            )
        self.pixel = pixels
        self.head = head
        self.period = period
        self.rpath = rpath
        self.sliding_time_int = dt.timedelta(hours=1)
        '''self.products = [key for key in self.dat.keys() if 'STEP_C' in key]
        self.flow_vector = np.array(
                [
                    [-0.8412, 0.4396, 0.3149],
                    [-0.8743, 0.457, 0.1635],
                    [-0.8862, 0.4632, -0.0],
                    [-0.8743, 0.457, -0.1635],
                    [-0.8412, 0.4396, -0.315],
                    [-0.7775, 0.5444, 0.3149],
                    [-0.8082, 0.5658, 0.1635],
                    [-0.8191, 0.5736, 0.0],
                    [-0.8082, 0.5659, -0.1634],
                    [-0.7775, 0.5444, -0.3149],
                    [-0.7008, 0.6401, 0.3149],
                    [-0.7284, 0.6653, 0.1634],
                    [-0.7384, 0.6744, -0.0],
                    [-0.7285, 0.6653, -0.1635],
                    [-0.7008, 0.6401, -0.315],
                ]
            )'''
        self.velocity_for_each_en_bin()

    def velocity_for_each_en_bin(self):
        '''based on the primary bin calculate the velocity'''
        self.bin_vel = np.sqrt(2*self.primary_energy*e*1000000/m_p)
        #Now for all Pixels
        self.bin_vel_vector = pd.DataFrame()
        for index, pixel in enumerate(self.pixel):
            dummy = np.outer(self.bin_vel,self.flow_vector[index])
            self.bin_vel_vector[f'{pixel}_x'] = dummy[:,0]
            self.bin_vel_vector[f'{pixel}_y'] = dummy[:,1]
            self.bin_vel_vector[f'{pixel}_z'] = dummy[:,2]

        '''# plot all points now in 3d space
        fig = plt.figure(figsize=(20, 10), dpi=300)
        ax = fig.add_subplot(111, projection="3d")
        for i in self.pixel:
            ax.scatter(
                self.bin_vel_vector[f'{i}_x'],
                self.bin_vel_vector[f'{i}_y'],
                self.bin_vel_vector[f'{i}_z'],
                label=f"Pixel {i}",
            )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        plt.title("Velocity vectors for each energy bin")
        plt.show()
        '''
       


    def load(self):
        """Lade die Daten aus dem angegebenen Zeitraum und speichere sie in self.dat und self.time."""
        print('Loading STEP data...')
        self.time, self.dat = load_nom(period=self.period, rpath=self.rpath, pixels=self.pixel, heads=self.head) 
        print('Combining data...')
        combine_data(time=self.time, dat=self.dat,pixels=self.pixel, heads=self.head)
    
    def get_velocity(self,energies = None):
        # 1000 beacause energies of primary are in MeV,so m/s is returned 
        return np.sqrt(2*energies*1000000*e/m_p)
    
    def proton_ebins(self):
        lowest_energy = 5.7
        highest_energy = 62.7
        #self.ebins = self.ebins[(self.ebins >= lowest_energy) & (self.ebins <= highest_energy)] 
        # mask self.ebins but keep shape
        self.ebins = np.ma.masked_where((self.ebins < lowest_energy) | (self.ebins > highest_energy), self.ebins)

    def calc_velocity_vetor(self):
        print('Calculating velocity vector...')
        self.sliding_mean()
        self.velocity_vector = pd.DataFrame()
        for i, key in enumerate(self.pixel):
            # ermittle die mittlere Geschindigkeit für die Energie
            dummy = self.get_velocity(self.sliding_en[f'{key}_data'])
            # jetzt müssen noch die Geschwindigkeiten mit den flow vectors vektoriesiert werden
            dummy  = np.outer(dummy,self.flow_vector[key])
            self.velocity_vector[f'{key}_x'] = dummy[:,0]
            self.velocity_vector[f'{key}_y'] = dummy[:,1]
            self.velocity_vector[f'{key}_z'] = dummy[:,2]
            self.velocity_vector[f'{key}_time'] = self.sliding_en[f'{key}_time']
        print('self.velocity_vector',self.velocity_vector)
    
    def sliding_mean(self):
        '''Calculate sliding mean of data.'''
        print('Calculating sliding mean...')    
        # init dataframe for sliding mean
        self.sliding_en = pd.DataFrame()
        # calculate the sliding mean with the given time interval in sliding_time_int
        for k,product in enumerate(self.products):
            dummy = np.array([])
            for i in range(len(self.time[product])):
                # create timeintervall
                time_int = (self.time[product][i], self.time[product][i] + self.sliding_time_int)
                # create mask for timeintervall
                mask = np.array((self.time[product] >= time_int[0]) & (self.time[product] < time_int[1]))
                # calculate sliding mean energy deposition and also normalize it
                dummy = np.append(dummy, np.sum(self.dat[product][mask]*self.ebins[1:])/np.sum(self.dat[product][mask]))
            self.sliding_en[f'{k}_time'] = self.time[product]
            self.sliding_en[f'{k}_data'] = dummy

    def better_sliding_mean(self):
        # i want to calculate the sliding mean with numpy arrays because it is faster
        print('Calculating sliding mean...')
        # init dataframe for sliding meas
        self.sliding_en = pd.DataFrame()
        #init array for timeinterval
        time_int = []
        for i,timestamp in enumerate(self.time[self.products[0]]):
            time_int.append((timestamp,timestamp+self.sliding_time_int))
        # init array for sliding mean
        self.sliding_en = pd.DataFrame()
        for pixels in self.products:
            pass
        
    def plot_uncorrected_sliding_mean(self):
        """Plotte die sliding mean für die unkorrigierten Daten."""
        print('Plotting uncorrected sliding mean...')
        fig, ax = plt.subplots(figsize=(20, 10), dpi=300)
        for i in self.pixel:
            ax.semilogy(self.sliding_en[f'{i}_time'], self.sliding_en[f'{i}_data'])
        ax.set_xlabel('Time')
        ax.set_ylabel('Sliding mean energy deposition in keV')
        ax.legend(self.products)
        plt.title('Sliding mean energy deposition for uncorrected data')
        plt.savefig('sliding_mean_uncorrected.png')


# list of electron events
list_electron_ev = [(dt.datetime(2022, 3, 18), dt.datetime(2022, 3, 19))]
pixels = range(2)
period = (dt.datetime(2022, 8, 1), dt.datetime(2022, 8, 2))


#CORRECT_PIXEL_AGING(period=period,pixels=pixels)

# now correct somehow the swa loading issue
CORRECT_PIXEL_AGING()
#STEP()
#STEP(period=period,pixels=pixels)