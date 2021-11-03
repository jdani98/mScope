#!usr/bin/python3
# Author: Jose Daniel Viqueira Cao
#         josedaniel.viqueira@rai.usc.es
# Project: ANALYSING DATA FROM AN ARRAY OF ANTENNAS

"""
INSTRUCTIONS
Program prepared to read summary of stats per groups. Input file must contain the following
columns: 'Gr','Ch','N','MEAN','STD','MAX','MIN','A','X0','SIGMA'
It returns mean, std, relative std in % and kurtosis for each of the following quantities:
MEAN, STD, MAX, MIN

=> In CONFIGURATION:
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
#import scipy.optimize as sco
import os
from collections import Counter
import csv
import pandas as pd
import scipy.stats as sct

#import functions as fc

### CONFIGURATION #############################################################
fname = 'summary210921.csv'  # FILENAME
outname = 'summary210921_stats'

options = {	'in_dir':		'Datos_dipolo_210910',
			'ADC_ch':		int(2**14),
			'nsamp':		1030,
            'deltaT':       4,
            'Vpp':          0.5,
            'ngroups':      'free',
            'Nevents':      10000,
            'Show_plots':   True
			}
###############################################################################


cwd = os.getcwd()  # main current working directory
os.chdir(cwd+'/'+options['in_dir'])
list_dir = os.listdir()

DATAf = pd.read_csv(fname, sep =',', header = None)
DATAf.columns = ['Gr','Ch','N','MEAN','STD','MAX','MIN','A','X0','SIGMA']

def STATS(i):
	Ch = DATAf[DATAf['Ch'] == i]
	#
	mean_mean   = np.mean(Ch['MEAN'])
	std_mean    = np.std(Ch['MEAN'])
	stdrel_mean = std_mean/mean_mean*100
	kurt_mean   = sct.kurtosis(Ch['MEAN'])
	#
	mean_std    = np.mean(Ch['STD'])
	std_std     = np.mean(Ch['STD'])
	stdrel_std  = std_std/mean_std*100
	kurt_std    = sct.kurtosis(Ch['STD'])
	#
	mean_max    = np.mean(Ch['MAX'])
	std_max     = np.std(Ch['MAX'])
	stdrel_max  = std_max/mean_max*100
	kurt_max    = sct.kurtosis(Ch['MAX'])
	#
	mean_min    = np.mean(Ch['MIN'])
	std_min     = np.std(Ch['MIN'])
	stdrel_min  = std_min/mean_min*100
	kurt_min    = sct.kurtosis(Ch['MIN'])
	
	return Ch, np.array([[mean_mean  , mean_std  , mean_max  , mean_min  ],
	                     [std_mean   , std_std   , std_max   , std_min   ],
	                     [stdrel_mean, stdrel_std, stdrel_max, stdrel_min],
	                     [kurt_mean  , kurt_std  , kurt_max  , kurt_min  ]])

def text_box(sm):
	header =  '             MEAN      STD     MAX      MIN '
	line1  = (' mean     %8.2f %8.2f %8.2f %8.2f' %(sm[0][0],sm[0][1],sm[0][2],sm[0][3]))
	line2  = (' std      %8.2f %8.2f %8.2f %8.2f' %(sm[1][0],sm[1][1],sm[1][2],sm[1][3]))
	line3  = (' stdR       %6.2f   %6.2f   %6.2f   %6.2f' %(sm[2][0],sm[2][1],sm[2][2],sm[2][3]))
	line4  = (' kurtosis   %6.2f   %6.2f   %6.2f   %6.2f' %(sm[3][0],sm[3][1],sm[3][2],sm[3][3]))
	return '\n'+header+'\n'+line1+'\n'+line2+'\n'+line3+'\n'+line4

file2 = open(outname+'.txt','w')
file3 = open(outname+'.csv','w')
writer = csv.writer(file3)
for channel in range(8):
	Ch_data, Ch_stats = STATS(channel)
	if options['Show_plots']:
		fig, ax = plt.subplots()
		plt.hist(Ch_data['STD'], rwidth=0.8)
		plt.xlabel('STD')
		plt.ylabel('counts')
		plt.title('Stats of std for groups of channel '+str(channel))
		plt.savefig(outname+'_'+str(channel)+'.png')
		plt.show()
	#plt.figtext(0.0,0.9,Ch_data['STD'].describe())
	#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	#plt.text(0.78, 0.40, text_box(Ch_stats), transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=props)
	writer.writerow([])
	writer.writerow(['Ch',channel])
	writer.writerow(['MEAN','STD','MAX','MIN'])
	writer.writerow(Ch_stats[0])
	writer.writerow(Ch_stats[1])
	writer.writerow(Ch_stats[2])
	writer.writerow(Ch_stats[3])
	#print('Channel ', channel)
	#print('MEAN  STD   Max   Min')
	#print(Ch_stats)
	file2.write(('\n\nChannel %i' %channel))
	file2.write(text_box(Ch_stats))
file2.close()
file3.close()
