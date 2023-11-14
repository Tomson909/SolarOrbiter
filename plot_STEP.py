import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
from Hausarest.STEP import STEP

def plot_STEP():
    first_date = dt.datetime(2020, 1, 1)
    for i in range(0, 1500):
        try:
            date = first_date + dt.timedelta(days=i)
            print(date)
            year = date.strftime("%Y")
            month = date.strftime("%m")
            day = date.strftime("%d")
            data = STEP(year=int(year), month=int(month), day=int(day))
            # save and plot data
            data.plot_ts(save='/media/tomruge/Toms_SSD/data/STEP_Plots/2020',  head=0)
            data.plot_ts(save='/media/tomruge/Toms_SSD/data/STEP_Plots/2020', head=1)
            data.plot_ts(save='/media/tomruge/Toms_SSD/data/STEP_Plots/2020',  head=-1)
        except Exception as e:
            print(e)
            continue    
        plt.close('all')
plot_STEP()