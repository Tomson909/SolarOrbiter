import sunpy_soar
from sunpy.net import Fido
import datetime as dt
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Product,Identifier
import numpy as np

def _load_swa_from_SolOArchive(date):
    '''
    If level 2 data not available, it downloads level 3 without notification,
    which is shitty. So we have to check if level 2 data is available.
    '''
    print("Lade Daten von SolO Archive...")
    # Create search attributes
    instrument = Instrument("SWA")
    time = Time(date, date + dt.timedelta(days=1))
    level = Level(2)
    identifier = Product("swa-pas-grnd-mom")
    # Do search
    print(f'Idetifier: {identifier} \n Instrument: {instrument} \n Level: {level} \n Time: {time}')
    result = Fido.search(instrument, time, level, identifier)
    # Download files
    Fido.fetch(result, path = "/media/tomruge/Toms_SSD/data/SWA_Data")
