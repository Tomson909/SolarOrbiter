import sunpy_soar
from sunpy.net import Fido
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Identifier
import numpy as np
import os


#load data from SolO Archive

def load_SWA(period):
    print('Lade Daten von SolO Archive')
    # Create search attributes
    instrument = Instrument('SWA')
    time = Time('2022-11-19', '2022-11-20')
    level = Level(2)
    identifier = Identifier('swa-pas-grnd-mom')
    
    # Do search
    result = Fido.search(instrument, time, level, identifier)
    
    # Download files
    file = Fido.fetch(result)
    # not sure whether its a cdf file or not but something like this should work
    return file

def write_SWA(file):
    print('Schreibe Daten in txt-Datei')
    # load timestamps
    timestamps = np.array(file.varget('EPOCH'))
    # convert timestamps to datetime objects
    timestamps = np.array([file.cdfepoch.to_datetime(x) for x in timestamps])
    # load data in rtn coordinate system
    swgeschw = file.varget('V_SRF')
    pressure = file.varget('P_SRF')
    temperature = file.varget('TXTYTZ_SRF')
    data = np.column_stack((timestamps, swgeschw[:,0], swgeschw[:,1], swgeschw[:,2], pressure[:,0], pressure[:,1], pressure[:,2], temperature[:,0], temperature[:,1], temperature[:,2]))
    column_names = ['Zeit', 'V_x', 'V_y', 'V_z', 'P_x', 'P_y', 'P_z', 'T_x', 'T_y', 'T_z']
    np.savetxt('sonnenwind.txt', data, delimiter=',', header=','.join(column_names), comments='')
    return 0

def main(period):
    # transform datetime objects to strings 
    period = (period[0].strftime('%Y-%m-%d'), period[1].strftime('%Y-%m-%d'))
    # load data
    file = load_SWA(period=period)
    # write data to txt file, each day 1 txt file
    write_SWA(file)    

    
def check_file_already_load_writtem(period):
    period = (period[0].strftime('%Y-%m-%d'), period[1].strftime('%Y-%m-%d'))
    # create array of dates for each day in period, period is a tuple of datetime objects
    dates = np.arange(np.datetime64(period[0]), np.datetime64(period[1])+np.timedelta64(1,'D'), np.timedelta64(1,'D'), dtype='datetime64')
    # check if file already exists
    for date in dates:
        if os.path.isfile(date.strftime('%Y-%m-%d') + '.txt'):
            print('Datei existiert bereits')
            return True
    
