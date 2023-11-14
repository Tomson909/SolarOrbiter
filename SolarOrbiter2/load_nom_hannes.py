import os
import datetime as dt
import numpy as np

ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])


def load_nom(rpath = '/data/projects/solo/step_v0008/',period = (dt.datetime(2021,10,22),dt.datetime(2021,10,23)), products = ('M','A'), heads = (0,1), pixels = range(16)):
    dat = {}
    time = {}
    prods = []
    for pt in products:
        for h in heads:
            for pix in pixels:
                prods.append('STEP_%s_%i_%.2i'%(pt,h,pix))
    for p in prods:
        dat[p] = []
        time[p] = []
    
    flist = np.sort(os.listdir(rpath))
    for p in prods:
        pflist = []
        for f in flist:
            tdate = dt.datetime.fromisoformat(f[:10])
            if p in f and tdate>=period[0] and tdate<period[1]:
                fin = open(rpath + f,'r')
                cols = len(fin.readline().split())
                for s in fin:
                    k = s.split()
                    if '_A_' in p:
                        time[p].append(dt.datetime.fromisoformat(k[0])-dt.timedelta(seconds = 60))
                    else:
                        time[p].append(dt.datetime.fromisoformat(k[0])-dt.timedelta(seconds = 1))
                    dat[p].append([int(v) for v in k[1:]])
    for p in prods:
        time[p] = np.array(time[p])
        dat[p] = np.array(dat[p])
    return time,dat

def combine_data(time,dat):
    """
    combines a (full) set of AUX and MAIN STEP SCI data to 1 min res products with 56 energy-channel and adds it to the dictionary
    """
    # combine AUX + MAIN products to 1 minute COMBINED product   -> Timestamps are assumed to be at the begining of the aquisition period (i.e. the begin of the minute)
    # First find unique timestamps of all 16 pixels (to have the final products at the same shape <-> deal with data gaps ...)
    tb = np.array([],dtype = dt.datetime)
    ts = np.array([],dtype = dt.datetime)
    for p in range(16):
        for h in range(2):
            k = '%i_%.2i'%(h,p)
            t = time['STEP_A_' + k]
            tb = np.unique(np.append(t,tb))
            ts = np.unique(np.append(t,ts))
    tb = np.append(tb,tb[-1] + dt.timedelta(seconds = 60))

    for p in range(16):
        for h in range(2):
            k = '%i_%.2i'%(h,p)
            t = time['STEP_A_' + k]
            td = np.zeros((tb.shape[0]-1,56))
            for i in range(8):
                H,x = np.histogram(t,bins = tb,weights=dat['STEP_A_' + k][:,i])
                td[:,i] = H
            for i in range(16):
                H,x = np.histogram(t,bins = tb,weights=dat['STEP_A_' + k][:,i+8])
                td[:,i+40] = H
            for i in range(32):
                H,x = np.histogram(time['STEP_M_' + k],bins = tb,weights=dat['STEP_M_' + k][:,i])
                td[:,8+i] = H
            dat['STEP_C_' + k] = td
            time['STEP_C_' + k] = tb[:-1]
            
            

