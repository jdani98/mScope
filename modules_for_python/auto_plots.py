#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLASSES TO PLOT
Created on Wed Oct 20 11:26:21 2021

@author: daniel
"""

import numpy as np
import scipy.stats as sct
import matplotlib.pyplot as plt

list_of_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown',
                  'tab:pink','tab:grey','tab:olive','tab:cyan']

class histo_stats:
    
    def __init__(self, data):
        """

        Parameters
        ----------
        data : TUPLE
            It must contain one or more lists with data to plot its histogram.

        Returns
        -------
        None.

        """
        self.dat = data
    
    def leng(self):
        """

        Returns
        -------
        Tuple with the number of data of each array.

        """
        leng = ()
        for ristra in self.dat:
            leng += (len(ristra),)
        return leng
    
    
    def mean(self):
        """

        Returns
        -------
        Tuple with mean of each array of data. If it is only one array, the tip for reading it to add a comma, i.e.
        x, = self.mean()

        """
        mean = ()
        for ristra in self.dat:
            mean += (np.mean(ristra),)
        return mean
    
    
    def stdev(self):
        """

        Returns
        -------
        Tuple with standard deviation of each array of data.

        """
        stdev = ()
        for ristra in self.dat:
            stdev += (np.std(ristra),)
        return stdev
    
    
    def skew(self):
        """

        Returns
        -------
        Tuple with skewness of each array of data.

        """
        skwn = ()
        for ristra in self.dat:
            skwn += (sct.skew(ristra),)
        return skwn
    
    
    def kurt(self):
        """

        Returns
        -------
        Tuple with kurtosis of each array of data.

        """
        kurt = ()
        for ristra in self.dat:
            kurt += (sct.kurtosis(ristra),)
        return kurt
    
    
    def first_mm(self,n):
        """

        Parameters
        ----------
        n : INTEGER
            Number of maxima to return.

        Returns
        -------
        Tuple with lists of events with maximum deviation from the mean and its value.

        """
        firstmm = ()
        for ristra in self.dat:
            firstmmi = []
            media = np.mean(ristra)
            desv = np.abs(ristra-media)
            for i in range(n):
                imax = np.argmax(desv)
                desv[imax] = 0
                firstmmi += [imax,ristra[imax]]
            firstmm += (firstmmi,)
        return firstmm
    
    
    def autocolors(self):
        """

        Returns
        -------
        Tuple of automatic colors for plotting histogram.

        """
        autocol = ()
        for i in range(len(self.dat)):
            autocol += (list_of_colors[i],)
        return autocol
    
    def autostyle(self):
        """

        Returns
        -------
        Tuple of automatic styles of histogram.

        """
        autosty = (True)
        for i in range(1,len(self.dat)):
            autosty += (False,)
        return autosty
    
    
    def autobins(self):
        """

        Returns
        -------
        Tuple of automatic bins for histogram.

        """
        autob = ()
        for i in range(len(self.dat)):
            autob += (None,)
        return autob
    
    
#    def plot_histo(self, ax, labels, show_stats=True, colors=autocolors(self), style=autostyle(), bins=autobins(),
    def plot_histo(self, ax, labels=None, show_stats=True, pos_stats=None,
                   colors=None, style=None, bins=None, rang=None,density=False,
                   weights=None,cumulative=False, bottom=None,
                   histtype='bar',align='mid',log=False):
        """

        Parameters
        ----------
        ax: MATPLOTLIB.AXIS
            Axis where to plot the histogram.
        colors : TUPLE(STRING), optional
            Colors for each histogram. The default is list_of_colors.
        labels: TUPLE(STRING), optional
            Labels for each histogram. The default is None.
        show_stats: BOOL, optional
            If true, statistics are shown in a panel. The default is True.
        pos_stats: ARRAY-LIKE or TUPLE, optional
            Must contain pairs of values: the first is the horizontal position of the panel,
            the second is the vertical position. Values relative to the axes of the plot.
            The default is None and this means an automatic set.
        style : TUPLE(BOOL), optional
            Fill or not each histogram. The default is first filled and rest non-filled.
        bins : TUPLE(INTEGER), optional
            Contains numbers of bins for ach histogram. The default is None.
        rang : TUPLE, optional
            The lower and upper range of the bins. The default is None.
        density : BOOL, optional
            If True, draw and return a probability density:. The default is False.
        weights : ARRAY-LIKE OR NONE, optional
            An array of weights, of the same shape as data. The default is None.
        cumulative : BOOL OR -1, optional
            See matplotlib.pyplot.hist manual. The default is False.
        bottom : ARRAY-LIKE, SCALAR OR NONE, optional
            See matplotlib.pyplot.hist manual. The default is None.
        histtype : STRING, optional
            See matplotlib.pyplot.hist manual. The default is 'bar'.
        align : STRING, optional
            See matplotlib.pyplot.hist manual. The default is 'mid'.
        log : BOOL, optional
            If True, the histogram axis will be set to a log scale. The default is False.

        Returns
        -------
        Axis.

        """
        #####
        if labels == None:
            labels = ()
            for i in range(len(self.dat)): labels += (None,)
        if pos_stats == None:
            nhistos = len(self.dat)
            half = 0.25*nhistos/2
            pos_stats = []
            for i in range(nhistos): pos_stats.append([0.9,0.5-half+0.25*i])
        if colors == None:
            colors = ()
            for i in range(len(self.dat)): colors += (list_of_colors[i],)
        if style == None:
            style = ()
            for i in range(len(self.dat)): style += (False,)
        if bins == None:
            bins = ()
            for i in range(len(self.dat)): bins += (None,)
        if rang == None:
            rang = ()
            for i in range(len(self.dat)): rang += (None,)
        #####
        tN = (); tbines = (); tpatches = ()
        for i,ristra in enumerate(self.dat):
            N, bines, patches = ax.hist(ristra,bins=bins[i],range=rang[i],density=density,weights=weights,
                    cumulative=cumulative,bottom=bottom,histtype=histtype,align=align,
                    log=log,color=colors[i],edgecolor=colors[i], fill=style[i],label=labels[i])
            if show_stats:
                texto = ('n: %i \nmean: %6.2f \nstd: %6.2f \nskwn: %6.2f \nkurt: %6.2f' 
                               %(self.leng()[i],self.mean()[i],self.stdev()[i],self.skew()[i],self.kurt()[i]))
                props = dict(boxstyle='round',color=colors[i],fill=True,alpha=0.2)
                ax.text(pos_stats[i][0], pos_stats[i][1], texto, horizontalalignment='left', verticalalignment='center', 
                        bbox=props, transform=ax.transAxes)
            tN += (N,); tbines += (bines,); tpatches += (patches,)
                
        return tN, tbines, tpatches

"""
### Naive test
datos1 = [0,3,6]
datos2 = [0,4,8]
datos3 = [-2,-2,0]

histograma = histo_stats((datos1,datos2,datos3))
#print(histograma.mean())

fig1 = plt.figure(0)
ax1 = fig1.add_subplot()
#histograma.plot_histo(ax1,('1','2','3'),show_stats=True,pos_stats=((1.0,0.4),(1.0,0.6),(1.0,0.8)),colors=('r','b','k'),style=(True,False,False),bins=(None,None,None))
histograma.plot_histo(ax1)
ax1.set_xlabel('x')
ax1.set_ylabel('counts')
ax1.legend()
"""