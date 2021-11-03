#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DEFINING EVENT CLASS
Created on Tue Oct 19 11:33:01 2021

@author: daniel
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq
colors = ['b','g','r','c','m','y','k','tab:gray']

class Event:
    
    def __init__(self, ADC_counts, channel=0,nbits=14,baseline=8192,dt=4,Vpp=2,stime=None):
        """

        Parameters
        ----------
       ADC_counts : LIST
            List of all the adc counts.
       channel : INTEGER, optional
            Channel from which the data was taken. The default is 0.
        nbits : INTEGER, optional
            Number of bits of ADC_counts. The default is 14.
        baseline : INTEGER, optional
            ADC_count of baseline. The default is 8192.
        dt : INTEGER, optional
            Separation between successive samples in nanoseconds. The default is 4.
        Vpp : FLOAT, optional
            Full scale of voltages. The default is 2.
        stime : INTEGER, optional
            Time at which the event starts. The default is None

        Returns
        -------
        None.

        """
        self.bulk      = ADC_counts
        self.ch        = channel
        self.starttime = stime
        self.adcrange  = 2**nbits
        self.basel     = baseline
        self.vpp       = Vpp
        self.Dt        = dt
    
    
    def ADC_counts(self):
        return self.bulk
    
    
    def channel(self):
        return self.ch
    
    
    def stime(self):
        return self.starttime
    
    
    def ADC_range(self):
        return self.adcrange
    
    
    def Vpp(self):
        return self.vpp
    
    def nsamples(self):
        """
        

        Returns
        -------
        INTEGER
            Number of samples of ADC_counts.

        """
        return len(self.bulk)
    
    
    def smp_times(self):
        """

        Returns
        -------
            Array with time corresponding to each ADC_counts.

        """
        return [self.Dt*i for i in range(len(self.bulk))]
    
    
    def ADC_rebase(self):
        """
        
        Returns
        -------
        List of ADC counts with baseline as reference
        """
        ADCreb = [item-self.basel for item in self.bulk]
        return ADCreb
    
    
    def evn_selection(self,ini,fin):
        """

        Parameters
        ----------
        ini : INTEGER
            Index of first count to select.
        fin : INTEGER
            Index of last count to select.

        Returns
        -------
        A list of selected items between index ini and fin.

        """
        return self.bulk[ini:fin]
    
    
    def fourier(self):
        """

        Returns
        -------
        Fourier frequencies (MHz) and Fourier transform. Returns only positive part of the spectrum.

        """
        n = self.nsamples()
        yC = self.vpp / self.adcrange
        xC = 1.e-6
        yF = yC * (2.0/n) * np.abs(fft(self.bulk)[:n//2])
        xF = xC * fftfreq(n,(self.Dt)*1.e-9)[:n//2] 
        return xF,yF
    
    
    def evn_display(self,ax):
        """

        Parameters
        ----------
        ax : 
            axis where to plot ADC_counts.

        Returns
        -------
        Display of the event.

        """
        y = self.bulk
        #times = [self.Dt*i for i in range(len(y))]
        return ax.plot(self.smp_times(),y,linewidth=0.5,color=colors[self.ch],label='Ch'+str(self.ch))
    
    
    def evn_display_fourier(self,ax):
        """

        Parameters
        ----------
        ax : MATPLOTLIB.PYPLOT AXIS
            Axis where to plot Fourier transform of event.

        Returns
        -------
        Display of the event.

        """
        x,y = self.fourier()
        return ax.plot(x,y, linewidth=0.5,color=colors[self.ch],label='Ch'+str(self.ch))

"""
#Naive test
evento = Event([8456,3500,11242,6532,4551],0,14,8192,4,0.5)

print(evento.ADC_counts())
print(evento.channel())
print(evento.stime())
print(evento.ADC_range())
print(evento.Vpp())
print(evento.nsamples())
print(evento.ADC_rebase())
print(evento.evn_selection(2,4))

fig = plt.figure(0)
ax1 = fig.add_subplot()
evento.evn_display(ax1)

fig1 = plt.figure(1)
ax1 = fig1.add_subplot()
evento.evn_display_fourier(ax1)
"""