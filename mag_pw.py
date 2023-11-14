import numpy as np
import matplotlib.pyplot as plt
from cdflib import *
from cdflib import cdfepoch
import datetime as dt
import os
import pandas as pd

# for loading data from SolO Archive
import sunpy_soar
from sunpy.net import Fido
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Product
# exclude warnings, maybe I should not do this, but I don't know how to fix it
import warnings

#install astropy
from astropy.time import Time
warnings.filterwarnings("ignore",category=FutureWarning)
warnings.filterwarnings("ignore",category=RuntimeWarning)

# import data import
from load_nom_II import load_nom, combine_data

class CORRECT_PIXEL_AGING:
    def __init__(self,reference_pixel = 1,period=(dt.datetime(2022, 11, 18), dt.datetime(2022, 11, 19)),pixels=range(15),heads = (0,1)) -> None:
        self.period = period
        self.pixels = pixels
        self.heads = heads
        self.reference_pixel = reference_pixel
        #self.dat_SWA = SWA(period=period,pixels=pixels)
        self.dat_MAG = MAG(period=period,pixels=pixels)
        self.dat_STEP = STEP(pixels = pixels,head=heads,period=period)
        self.correct_energy_deposit()
        self.plot_corrected_sliding_mean()
        self.dat_STEP.plot_uncorrected_sliding_mean()
    
    def formula(self,reference_pitchangle = None,pitchangle= None):
        # problem: the timestamps are probably not the same, so I have to interpolate the data
        # but they seem to have the same shape, which is good and makes my life better
        return np.cos(pitchangle)**2/np.cos(reference_pitchangle)**2   
    
    def average_correction_factors(self,step_time = None, pw_time = None, product_step = None,correction_factors = None):
        '''Average the correction factors to the timestamps of the STEP data.'''
        #product needs to be 
        print('Average correction factors...')
        # init array for correction factors
        '''print('pw time',pw_time)
        print('step time',step_time)'''
        dummy = []
        for i in range(len(step_time)-1):
            # create timeintervall
            time_int = (step_time[i], step_time[i+1])
            # create mask for timeintervall
            mask = np.array((pw_time >= time_int[0]) & (pw_time < time_int[1]))
            # calculate average correction factor
            dummy.append(np.mean(correction_factors[mask]))
        '''print('mask',mask)
        print('shape of mask',mask.shape)
        print('masked correction factors',correction_factors[mask])
        print('mean of masked correction factors',np.mean(correction_factors[mask]))
        print('dummy',dummy)'''
        self.correction_factors= np.array(dummy)

        print('Average correction factors...done')

    def correct_energy_deposit(self):
        """Korrigiere die Energieabscheidung der Elektronen für die Pixel.
        Es ist im Prinzip eine VErschiebung der Zeitachse.
        Ich brauche:
                    - die Pitchwinkel der Pixel
                    - für jeden Zeitstempel der STEP DATEN eine korrigierte Energie,
                      multiplikation der deponierten Energie mit Korrekturfaktor"""
        """Formel für die Korrektur der Energieabscheidung der Teilchen."""
        print('Correcting energy deposit...')
        # WOOOOW es funktioniert, ich habe die korrigierten Energien
        self.corrected_mean_energy = pd.DataFrame()
        for i in self.pixels:
            product = self.dat_STEP.products[i]
            #correction factor
            self.correction_factors = self.formula(reference_pitchangle=self.dat_MAG.pitchangles[f'PW_{int(self.reference_pixel)}'], pitchangle=self.dat_MAG.pitchangles[f'PW_{i}'])
            # bevause correctionfactors have way more timestamps than the STEP data
            #print arguments of self.average_correction_factors to see what is going on
            '''print('step time',self.dat_STEP.sliding_en[f'{product}_time'])
            print('pw time',self.dat_MAG.pitchangles['time'])
            print('correction factors',self.correction_factors)'''
            self.average_correction_factors(product_step = f'STEP_C{i}',step_time = self.dat_STEP.sliding_en[f'{i}_time'], pw_time = self.dat_MAG.pitchangles['time'],correction_factors = self.correction_factors)
            # in average_correction_factos is something wrong
            self.corrected_mean_energy[i] = self.dat_STEP.sliding_en[f'{i}_data'][:-1]*self.correction_factors   
            self.corrected_mean_energy[f'{i}_time'] = self.dat_STEP.sliding_en[f'{i}_time'][:-1]  

    def plot_corrected_sliding_mean(self):
        """Plotte die sliding mean für die korrigierten Daten."""
        print('Plotting corrected sliding mean...')
        fig, ax = plt.subplots(figsize=(20, 10), dpi=300)
        for i in self.pixels:
            ax.semilogy(self.corrected_mean_energy[f'{i}_time'], self.corrected_mean_energy[i])
        ax.set_xlabel('Time')
        ax.set_ylabel('Sliding mean energy deposition in keV')
        ax.legend(self.dat_STEP.products)
        plt.title('Sliding mean energy deposition for corrected data')
        plt.savefig(f'sliding_mean_corrected.png')
        plt.close('all')    

class MAG(CORRECT_PIXEL_AGING):
    def __init__(
            self, period=(dt.datetime(2022, 12, 12), dt.datetime(2022, 12, 13)),pixels=range(15)
        ) -> None:
            # load only 1 day for now
            self.period = period
            self.mag_path = "/data/projects/solo/mag/l2_soar/srf/2022"
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
            # calculate pitchangles
            self.calculate_pitchangle()  
                      
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
        
    def pw(self, B = None, flow_vector = None):
        """Übergebe den particle flow-Vektor als Geschwindigkeit und den Magnetfeldvektor (am besten in SRF) und berechne die Pitchwinkel über das Skalarprodukt."""
        len_flow = np.sqrt(flow_vector[0] ** 2 + flow_vector[1] ** 2 + flow_vector[2] ** 2)
        len_B = np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
        argument = (flow_vector[0] * B[0] + flow_vector[1] * B[1] + flow_vector[1] * B[2]) / len_flow / len_B
        result = np.arccos(argument)
        return result

    def calculate_pitchangle(self):
        """Berechne die Pitchwinkel für die Elektronen, welche auf STEP treffen in erster Näherung.
        Dafür wird der Winkel zwischen dem particle flow vector der Pixel und dem Magnetfeld herangezogen.
        """
        print("Calculating pitchangles...")
        # I think here is something wrong because my pitchangles have weird values and are from 0 to 180deg.
        self.pitchangles = pd.DataFrame()
        for i in self.pixels:
            # using SRF coordinates
            self.pitchangles[f'PW_{i}'] = self.pw(B=np.array([self.B_SRF['Bx'],self.B_SRF['By'],self.B_SRF['Bz']]),flow_vector=self.flow_vector[i])
        self.pitchangles['time'] = np.array(self.B_SRF['time'])
                        
    def plot_pitchangles(self):
        """Plotte die Pitchwinkel für die Elektronen, welche auf STEP treffen in erster Näherung."""
        print("Plotting pitchangles...")
        fig2, ax = plt.subplots(figsize=(20, 10), dpi=300)
        for i in range(len(self.pixels)):
            ax.plot(self.B_SRF['time'], self.pitchangles[i])
        ax.set_xlabel("Time")
        ax.legend([f'Pixel {i}' for i in self.pixels])
        ax.set_ylabel("Pitchangle in deg")
        plt.title("Pitchangle of electrons hitting STEP")
        plt.savefig("pitchangle.png")
        plt.close("all")        
        
class SWA(CORRECT_PIXEL_AGING):
    def __init__(self, period=(dt.datetime(2022, 12, 12), dt.datetime(2022, 12, 13)),pixels=range(15)) -> None:
        # load only 1 day for now
        self.period = period
        self.swa_path = "SWA_Data"
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
    
    def _load_swa(self):
        """Check if the data from SolO Archive is already downloaded."""
        print("Loading SWA data...")
        period_diff_format = (
            self.period[0].strftime("%Y-%m-%d"),
            self.period[1].strftime("%Y-%m-%d"),
        )
        # create array of dates for each day in period, period is a tuple of datetime objects
        dates = np.arange(
            np.datetime64(period_diff_format[0]),
            np.datetime64(period_diff_format[1]) + np.timedelta64(1, "D"),
            np.timedelta64(1, "D"),
            dtype="datetime64",
        )
        dates = [dt.datetime.strptime(str(x), "%Y-%m-%d") for x in dates]
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
                # if not, download it
                self._load_swa_from_SolOArchive(date)
            dat = CDF(path)
            # load velocity data and save it in arrays
            velocity_data = dat.varget("V_SRF")
            #load timestamps
            dummy = pd.DataFrame(columns=['time','Vx','Vy','Vz'])
            dummy['time'] = cdfepoch.to_datetime(np.array(dat.varget("Epoch")))
            dummy['Vx'] = velocity_data[:, 0]
            dummy['Vy'] = velocity_data[:, 1]
            dummy['Vz'] = velocity_data[:, 2]
            dummy.fillna(0, inplace=True)
            self.V_SRF = pd.concat([self.V_SRF,dummy])

    def _load_swa_from_SolOArchive(self, date):
        print("Lade Daten von SolO Archive...")
        # Create search attributes
        instrument = Instrument("SWA")
        time = Time(
            self.period[0].strftime("%Y-%m-%d"),
            (self.period[0] + dt.timedelta(days=1)).strftime("%Y-%m-%d"),
        )
        level = Level(2)
        identifier = Product("swa-pas-grnd-mom")
        # Do search
        result = Fido.search(instrument, time, level, identifier)
        # Download files
        Fido.fetch(result, path="SWA_Data")
       
class STEP(CORRECT_PIXEL_AGING):
    def __init__(
        self,
        pixels=range(16),
        head=(0, 1),
        period=(dt.datetime(2021, 10, 22), dt.datetime(2021, 10, 23)),
        rpath="/data/projects/solo/step_v0008/",
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
        self.pixel = pixels
        self.head = head
        self.period = period
        self.rpath = rpath
        self.sliding_time_int = dt.timedelta(hours=1)
        self.load()
        self.sliding_mean()
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
        
        
    def load(self):
        """Lade die Daten aus dem angegebenen Zeitraum und speichere sie in self.dat und self.time."""
        print('Loading STEP data...')
        self.time, self.dat = load_nom(period=self.period, rpath=self.rpath, pixels=self.pixel, heads=self.head)
        print('Combining STEP data...')
        combine_data(time=self.time, dat=self.dat,pixels=self.pixel, heads=self.head)
    
    def sliding_mean(self):
        '''Calculate sliding mean of data.'''
        print('Calculating sliding mean...')    
        self.products = [key for key in self.dat.keys() if 'STEP_C' in key]
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

    def correct_energy_electrons(self):
        pass

# list of electron events
list_electron_ev = [(dt.datetime(2022, 3, 18), dt.datetime(2022, 3, 19))]
pixels = range(15)
period = list_electron_ev[0]
dat_SWA = SWA(period=period,pixels=pixels)
#dat_MAG = MAG(period=period,pixels=pixels)
#dat_STEP = STEP(pixels = pixels[1:],head=[0],period=list_electron_ev[0],rpath='/data/projects/solo/step_v0008/')
#CORRECT_PIXEL_AGING(reference_pixel=1,period=list_electron_ev[0],pixels=pixels,heads=[0])