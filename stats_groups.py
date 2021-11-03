#!usr/bin/python3
# Author: Jose Daniel Viqueira Cao
#         josedaniel.viqueira@rai.usc.es
# Project: ANALYSING DATA FROM AN ARRAY OF ANTENNAS

"""
INSTRUCTIONS:
Program prepared to read files with name formatted as 'filter[f]_r[r]_(no)sparks.txt', where [f]
is the number of filter (1,2,3 or 4) and [r] the number of ring 1 or 2).
If name format is different, change lines > nfilter = int(fname[6]) >nring   = int(fname[9])
> Ch = 2*(nfilter-1) + (nring-1)

It divides all the samples in several groups with <Nevents> each one, being each event a set of
<nsamp> samples. For each group, it returns the mean, std, max and min (simple stats) of the
samples and the gaussian parameters from a fit to the points of the histogram (gaussian).

=> In CONFIGURATION (variable 'options'):
   * fname: Name of file to read the data
   * outname: Name of file to write the outputs: simple statistics and gaussian parameters. New 
     data are appended to previous. Two files are created: one .txt and one .csv
   * options:
       - in_dir: directory to work, where input file is read and outputs are written
       - ADC_ch: number of ADC channels of the digitizer
       - deltaT: time between successive samples in nanoseconds, see specifications of digitizer
       - Vpp: voltage peak-peak set in the digitizer
       - ngroups: number of groups of events to generate.Can be string to stay as free parameter
         provided that <Nevents> is a integer.
       - Nevents: number of events per group. Can be string to stay as free parameter provided
         that <ngroups> is a integer.
       - Show_plots: True to show plots of histograms
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sco
import os
from collections import Counter
import csv

#import functions as fc

### CONFIGURATION #############################################################
fname = 'filter4_r2_nosparks.txt'  # FILENAME
outname = 'summary210921'

options = {	'in_dir':		'Datos_dipolo_210910',
			'ADC_ch':		int(2**14),
			'nsamp':		1030,
            'deltaT':       4,
            'Vpp':          0.5,
            'ngroups':      'free',
            'Nevents':      10000,
            'Show_plots':   False
			}
###############################################################################


cwd = os.getcwd()  # main current working directory
os.chdir(cwd+'/'+options['in_dir'])
list_dir = os.listdir()

"""
file_list = []
for item in list_dir:
#	if item[0] == 'f' and item[3] == 'r' and item[7:9] == '10': # -!- impose the condition
#	if item[0:4] == 'wave':
    if item[0:6] == 'filter' and item[11:13]=='no':
        print(item)
        file_list.append(item)
file_list = sorted(file_list)
"""

nfilter = int(fname[6])
nring   = int(fname[9])
Ch = 2*(nfilter-1) + (nring-1)

g = 0
i = 0
ng = options['ngroups']
Ne = options['Nevents']

if type(ng) == str:
    file = open(fname)
    i=0
    print(' Pre-reading...')
    for line in file:
        i+=1
    Ns = Ne*options['nsamp']
    ng=i//Ns
    print(' Number of lines in file:', i)
    print(' Number of groups automatically set to '+str(ng))
    file.close()

if type(Ne) == str:
    file = open(fname)
    i=0
    print(' Pre-reading...')
    for line in file:
        i+=1
    Ne=i//(ng*options['nsamp'])
    print(' Number of lines in file:', i)
    print(' Number of events per group automatically set to '+str(Ne))
    file.close()

Ns = Ne*options['nsamp']

#####################################################
DATA = []
file = open(fname)
i=0
for line in file:
    if i%Ns == 0:
        g+=1
        if g>ng: break
        print('Group '+str(g)+' created')
        DATA.append([])
    DATA[g-1].append(int(line[:-1]))
    i+=1  
file.close()
#####################################################


colors = ['b','g','r','c','m','y','k','tab:gray']

def gaussian(x,A,x0,sgm):
	return A*np.exp(-(x-x0)**2/(2*sgm**2))


file2 = open(outname+'.txt',"a")
header = ('\n \n# STATISTICS OF FILE '+str(fname))
parheader = ('\n# ADC_ch: %i, nsamp: %i, deltaT: %.2f, Vpp: %.2f' 
             %(options['ADC_ch'],options['nsamp'],options['deltaT'],options['Vpp']))
tag =       ('\n#--------------Simple statistics-------------------|----------Gaussian----------')
subheader = ('\n#    Gr  Ch     N     MEAN     STD     MAX    MIN  |   A         X0     SIGMA   ')
file2.write(header+parheader+tag+subheader)

file3 = open(outname+'.csv',"a")
writer = csv.writer(file3)


plt.close('all')
for i,group in enumerate(DATA):
    counts = Counter(group)
    rango  = np.arange(0,options['ADC_ch'])
    histo  = np.array([0]*options['ADC_ch'])
    for j in range(len(rango)):
        histo[j] = counts[j]
	
    media = np.mean(group)
    desvs = np.std(group)
    maxim = np.max(group)
    minim = np.min(group)
    popt, pcov = sco.curve_fit(gaussian, rango,histo, p0 = [max(histo),media,desvs])
	
    A,x0,sgm = popt
	
    line = ('\n  %5i  %2i  %5i  %7.1f  %7.1f  %5i  %5i  %7.1f   %7.1f  %7.1f ' 
    %(i+1,Ch,Ne,media,desvs,maxim,minim,A,x0,sgm))
    file2.write(line)
    
    writer.writerow([i+1,Ch,Ne,media,desvs,maxim,minim,A,x0,sgm])
	
    if options['Show_plots']:
        #fig = plt.figure('histo_Ch'+str(Ch)+'_Gr'+str(i+1),figsize=(10,8))
        fig, ax = plt.subplots()
        plt.title('Histogram of Ch'+str(Ch)+' g'+str(i+1))
        plt.plot(rango,histo,'.',color=colors[Ch],label='Ch'+str(Ch))
        plt.plot(rango,gaussian(rango,A,x0,sgm),'--',color=colors[Ch-1])
        plt.xlabel('ADC channels')
    
        header1 = ('SIMPLE STATS \n')
        texto1 = ((' N. events: %i \n' %Ne)+
                  (' Mean: %.1f \n' %media)+
                  (' Std.:  %.1f \n' %desvs)+
                  (' Max.: %i \n' %maxim)+
                  (' Min.: %i \n' %minim))
        header2 = ('\nFIT PARAMS \n')
        texto2 = ((' $x_0=$ %.1f \n' %x0)+
                  (' $\sigma=$ %.1f \n' %sgm)+
                  (' $A=$  %.2f') %A)
        texto = header1+texto1+header2+texto2
    
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(0.78, 0.40, texto, transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=props)
        plt.legend()
        plt.grid()
        plt.show()
        #plt.savefig('histo_'+fname[-4]+str(i+1)+'.png')

file2.close()
file3.close()
