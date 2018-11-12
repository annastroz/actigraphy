
# coding: utf-8

# ### Actigraphy module
# 
# #### Opening .csv files, filtering and convolution, splitting to days with curve fitting before. 
# 
# Main result: estimates of fitting parameters with their errors.

# In[4]:

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import islice
import dateutil
import matplotlib.dates as md
import datetime
import pandas as pd
from scipy.optimize import curve_fit
import scipy.signal as ss
#import seaborn as sns


# In[137]:

def csv_file_opener(filename):

    '''
    Opens GENEActiv .csv file of a given 'filename' parameter. 
    Function returns four lists: 
    x_axis – values from x axis
    y_axis – values from y axis
    z_axis – values from z axis
    mvm – sum of vector magnitudes/6000 (one value given for each epoch)
    dates – datetime objects for each record respectively
    
    '''
    # Sampling freq is the same in all of the data files, [= 100 Hz] but samples which I used are averaged for 60-s epochs.
    
    with open(filename, 'r') as csvfile:
        dataset  = csv.reader(x.replace('\0', '') for x in csvfile) 
        counter  = 0
        raw_data = []
        for i, row in enumerate(dataset):
            if i >= 100:
                counter += 1
                raw_data.append(row)
    datestrings       = []
    x_axis            = np.array([]) 
    y_axis            = np.array([])
    z_axis            = np.array([])
    lux               = np.array([])
    
    btn               = np.array([])
    temp              = np.array([])
    svm               = np.array([])
    x_stdev           = np.array([])
    y_stdev           = np.array([])
    z_stdev           = np.array([])
    time              = np.arange(counter)

    for row in raw_data:
        datestrings.append(row[0][:-4])
        x_axis            = np.append(x_axis, float(row[1]))
        y_axis            = np.append(y_axis, float(row[2]))
        z_axis            = np.append(z_axis, float(row[3]))
        lux               = np.append(lux, float(row[4]))
        btn               = np.append(btn, bool(row[5]))
        temp              = np.append(temp, float(row[6]))
        svm               = np.append(svm, float(row[7]))
        x_stdev           = np.append(x_stdev, float(row[8]))
        y_stdev           = np.append(y_stdev, float(row[9]))
        z_stdev           = np.append(z_stdev, float(row[10]))

    dates = [dateutil.parser.parse(s) for s in datestrings]

    mvm = svm/6000
    
    return(x_axis, y_axis, z_axis, mvm, dates)

    
def mvm_plot(dates, mvm, figsize=[10, 8], xticks=True, ylim=False):
    
    """Plots the results of sum of vector magnitude in time. SVM is given as a sum over samples in one epoch, and for every sample 
        vector magnitude is calculated: (x**2 + y**2 + z**2)**(1/2) - 1g. 
        Therefore, if Fs = 100 Hz, and epoch length was 60 secs, then the SVM was calculated from 6000 samples during the bin-csv conversion. 
    
        Parameters:
        mvm – from the csv_reader
        dates – from the csv_reader        
    """
    
    dates_for_axis = []
    running_data   = dates[0]
    if xticks:
        while running_data < dates[-1]:
            dates_for_axis.append(running_data)
            running_data = running_data + datetime.timedelta(hours=6)
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(dates, mvm, linewidth=.6)
    plt.legend([r'$MVM$'], loc=1)
    plt.title(r'$Mean$ $Vector$ $Magnitude$ $(MVM)$', fontsize=18)
    plt.xlabel(r'$Date$', fontsize=14)
    plt.ylabel(r'$MVM$ $value$ $[m/s^2]$', fontsize=14)
         
    ax          = plt.gca()
    ax.set_xticks(dates_for_axis)
    xfmt        = md.DateFormatter('%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    ax.tick_params(direction='in')
    ax.grid(False)
    plt.xticks(rotation=90)
    fig_size    = plt.rcParams["figure.figsize"]
    fig_size[0] = figsize[0]
    fig_size[1] = figsize[1]
    plt.xticks(rotation=90)
    if ylim:
        plt.ylim((0, 500))
    plt.show()
    

def sixAMindex(dat):
    
    '''
    Function which helps to start days splitting from the exact hour of the first day, ie. 6:00 AM.
    '''
    for i in range(len(dat)):
        if dat[i].hour == 6 and dat[i].minute == 0:
            return i
        
def filt_n_convolve(mvm, windowlength):

    '''
    Using bandpass Butterworth filter to MVM data. Then convolution with rectangular normalized window of given *windowlength* parameter. 
    '''    
#     wp     = [0.3, 0.5]
#     ws     = [0.2, 0.6]
#     rp     = 1
#     rs     = 60
#     Fs     = 1/60                        
#     f      = np.arange(0.01,Fs/2,0.01) 
#     t      = np.arange(len(mvm))
#     [n,Wn] = ss.buttord(wp, ws, rp, rs)
#     [b,a]  = ss.butter(n,Wn, btype='bandpass')

#     mvm_filt = ss.filtfilt(b, a, mvm)

    w = np.ones(windowlength)
    avgd = np.convolve(abs(mvm), w/w.sum(), mode='valid')

    return avgd

def mvm_n_fit_plot(dates, mvm, day, fitted_curve, figsize=[10, 8], xticks=True, ylim=False):
   
    '''Same as mvm_plot, but appends the fitted curve calculated with curvefit.'''

    dates_for_axis = []
    running_data   = dates[0]
    if xticks:
        while running_data < dates[-1]:
            dates_for_axis.append(running_data)
            running_data = running_data + datetime.timedelta(hours=6)
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(dates[:len(dates)-59], mvm, 'darkblue', linewidth=1.)
    plt.plot(dates[:len(dates)-59], fitted_curve, 'orange', linewidth=1.)
 #   plt.plot(dates[:len(dates)-59], day[:len(dates)-59], 'gray', linewidth=.4)
    plt.legend(['$MVM$', '$Fit$'], loc=1, fontsize=9)
        
    ax          = plt.gca()
    ax.set_xticks(dates_for_axis)
    xfmt        = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    ax.tick_params(direction='in')
    ax.grid(False)
    plt.xticks(rotation=90)
    fig_size    = plt.rcParams["figure.figsize"]
    fig_size[0] = figsize[0]
    fig_size[1] = figsize[1]
    if ylim:
        plt.ylim((0, 500))
    plt.show()

def f( x, p0, p1, p2):
    T = 24 * 60
    return p0*np.sin(2*np.pi*x/T + p1) + p2

def ff(x, p):
    return f(x, *p)

def fit_curvefit(p0, datax, datay, function, yerr=0):
    
    '''
    Function which calculates curve fit using least squares algorithm. 
    
    p0 – vector of initial guesses for amplitude, phase and mesor.
    datax, datay – vector of time and MVM.
    function – model function, here sinusoidal. (defined as the function f)
    
    '''

    pfit, pcov = curve_fit(f,datax,datay,p0=p0, epsfcn=0.0001)
    
    error      = [] 
    
    for i in range(len(pfit)):
        try:
            error.append(np.absolute(pcov[i][i])**0.5)
        except:
            error.append( 0.00 )
    pfit_curvefit = pfit
    perr_curvefit = np.array(error)
    fitted_curve  = ff(datax, pfit_curvefit)

    return fitted_curve, pfit_curvefit, perr_curvefit 

def recording2days(dat, mvm, windowlength, Ndays=4):
    
    ''' Major function which calculates everything by referring to functions above. 
    It takes Ndays from the signal, starting from 6:00 AM of the first day. Filtering,
    convolution and curve fitting are implemented. 
    
    Returns estimates and their errors for every day of recording respectively.
    
    '''
    
    index_to_start = sixAMindex(dat)
    T = 24 * 60
    
    day_by_day_estimates = np.array([])
    day_by_day_errors = np.array([])
    plt.suptitle(r'$Mean$ $Vector$ $Magnitude$ $(MVM)$', fontsize=18)
    for i in range(1, Ndays+1):
        plt.subplots_adjust(hspace=1.7)
        dat_day = dat[index_to_start + (i-1)*T : index_to_start + i*T]
        day     = mvm[index_to_start + (i-1)*T : index_to_start + i*T]
        plt.subplot(Ndays, 1, i)
        plt.title('$Day$ ${}$'.format(i), fontsize=16)
        plt.xlabel('$Hour$', fontsize=14)
        plt.ylabel('$MVM$ $value$ $[m/s^2]$', fontsize=14)
        avgd = filt_n_convolve(day, windowlength)
        guess_mean = np.mean(avgd)
        guess_std = 3*np.std(avgd)/(2**0.5)
        guess_phase = 0
        err_stdev = 0
        t = np.arange(len(avgd))
        pstart = [guess_std, guess_phase, guess_mean]
        fitted_curve, estimates, perr = fit_curvefit(pstart, t, avgd, ff)
        mvm_n_fit_plot(dat_day, avgd, day, fitted_curve, figsize=[10,4], xticks=True)
        
        
        day_by_day_estimates = np.append(day_by_day_estimates, np.array([estimates]))
        day_by_day_errors = np.append(day_by_day_errors, np.array([perr]))
    return(day_by_day_estimates, day_by_day_errors)    

