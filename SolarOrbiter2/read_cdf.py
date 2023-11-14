import sunpy_soar
from sunpy.net import Fido
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Identifier

from 

def mache():
    print('hallo')
    # Create search attributes
    instrument = Instrument('SWA')
    time = Time('2022-11-19', '2022-11-20')
    level = Level(2)
    identifier = Identifier('swa-pas-grnd-mom')

    # Do search
    result = Fido.search(instrument, time, level, identifier)
    print(result)

    # Download files
    files = Fido.fetch(result,path='')
    print(files)
    
if __name__ == '__main__':
    mache()