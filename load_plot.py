#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 12:22:19 2023

@author: ruge
"""

import cdflib
from cdflib import CDF
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def load_cdf(filename):
    data = CDF(filename)
    # load timestamps
    timestamps = np.array(data.varget('EPOCH'))
    # convert timestamps to datetime objects
    timestamps = np.array([cdflib.cdfepoch.to_datetime(x) for x in timestamps])
    # load data in rtn coordinate system
    swgeschw = data.varget('V_SRF')
    pressure = data.varget('P_SRF')
    temperature = data.varget('TXTYTZ_SRF')
    
    return timestamps, swgeschw, pressure, temperature


def plot_data(velocity, pressure, temperature, timestamps):
    #close all open figures
    plt.close('all')
    #initialize figure
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(20,10))
    
    #filter out fill values, which are -1e31
    logic = velocity != -1e31
    #plot velocity
    ax1.plot(timestamps[logic[:,0]], velocity[:,0][logic[:,0]])
    ax1.plot(timestamps[logic[:,1]], velocity[:,1][logic[:,1]])
    ax1.plot(timestamps[logic[:,2]], velocity[:,2][logic[:,2]])
    ax1.set_title('Sonnenwindgeschwindigkeit')
    ax1.set_xlabel('Zeit')
    ax1.set_ylabel('Geschwindigkeit [km/s]')
    ax1.legend(['x','y','z'])
    
    #filter out fill values, which are -1e31
    logic = pressure != -1e31
    #plot pressure
    ax2.plot(timestamps[logic[:,0]], pressure[:,0][logic[:,0]])
    ax2.plot(timestamps[logic[:,1]], pressure[:,1][logic[:,1]])
    ax2.plot(timestamps[logic[:,2]], pressure[:,2][logic[:,2]])
    ax2.set_title('Sonnenwinddruck')
    ax2.set_xlabel('Zeit')
    ax2.set_ylabel('Druck [nPa]')
    ax2.legend(['x','y','z'])
    
    #filter out fill values, which are -1e31
    logic = temperature != -1e31
    #plot temperature
    ax3.plot(timestamps[logic[:,0]], temperature[:,0][logic[:,0]])
    ax3.plot(timestamps[logic[:,1]], temperature[:,1][logic[:,1]])
    ax3.plot(timestamps[logic[:,2]], temperature[:,2][logic[:,2]])
    ax3.set_title('Sonnenwindtemperatur')
    ax3.set_xlabel('Zeit')
    ax3.set_ylabel('Temperatur [K]')
    ax3.legend(['x','y','z'])
    plt.tight_layout()


    #save figure
    plt.savefig('sonnenwind.pdf')

def main():
    timestamps, swgeschw, pressure, temperature = load_cdf('solo_L2_swa-pas-grnd-mom_20221119_V02.cdf')
    plot_data(swgeschw, pressure, temperature, timestamps)
    
if __name__ == '__main__':
    main()