#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDReader Class
Created on Tue Oct 19 11:31:56 2021

@author: daniel
"""

from os import path
import numpy
import struct
import modules_for_python.events as evm
#import events as evm
import os
#import matplotlib.pyplot as plt

# ==================================
# Class to read WaveDump data files
# ==================================

class WDReader:
    
    def __init__(self, fileName, channel=0,nbits=14,length=1030,dt=4,Vpp=2):
        """

        Parameters
        ----------
        FileName : STRING
            Name of file to read. Include extension .txt or .dat.
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
            Time at which the event starts. The default is None.

        Returns
        -------
        None.

        """
        self.FileName = path.abspath(fileName)
        if fileName[-4:]=='.txt': self.filetype = 'r'  # read ASCII
        if fileName[-4:]=='.dat': self.filetype = 'rb' # read binary
        #self.open      = open(self.FileName, self.filetype)
        self.ch        = channel
        self.adcrange  = 2**nbits
        self.len       = length
        self.vpp       = Vpp
        self.Dt        = dt
        
        
    def read_flat(self):
        """

        Returns
        -------
        Flat list with all read ADC_counts.

        """
            
        if self.filetype == 'r':
            flatdata = numpy.loadtxt(self.FileName,dtype='i')
        
        if self.filetype == 'rb':
            flatdata = numpy.fromfile(self.FileName,dtype='<i2')
        return flatdata
    
    
    def read_whead(self):
        """
        
        Returns
        -------
        Channel and dictionary with tuples of time + ADC_counts of each event.

        """
        with open(self.FileName,self.filetype) as f:
            data0 = f.read()
        
        if self.filetype == 'r':
            lenheader = 7
            data1 = data0.split('\n')[:-1]
            channel = int(data1[2][9:])
            length1 = int(data1[0][15:])
            if length1 != self.len: print('WARNING. Lenghts of event do not match')
            total_length_event = self.len + lenheader
            nevents = len(data1)//total_length_event
            resto   = len(data1)%total_length_event
            if resto !=0: print('WARNING. Number of data does not match with event length')
            data = {}
            for i in range(nevents):
                iev = int(data1[i*total_length_event+3][14:])
                time = int(data1[i*total_length_event+5][20:])
                bulk = data1[i*total_length_event+7:(i+1)*total_length_event]
                bulk = [int(item) for item in bulk]
                data[iev] = (time,bulk)
        
        if self.filetype == 'rb':
            nbytes = struct.unpack('I',data0[0:4])[0]
            #print(nbytes)
            channel = struct.unpack('I',data0[12:16])[0]
            nevents = len(data0)//nbytes
            resto   = len(data0)%nbytes
            if resto != 0: print('WARNING. Lenghts of event do not match')
            data = {}
            for i in range(nevents):
                time = struct.unpack('I', data0[i*nbytes+20:i*nbytes+24])[0]
                bulk = struct.unpack('1030h', data0[i*nbytes+24:(i+1)*nbytes])
                bulk = [item for item in bulk]
                data[i] = (time,bulk)
        
        return channel,data
    
    
    def autoread(self):
        """
        Selects reading functions automatically and
        Returns
        -------
        the same format as read_whead()

        """
        try:
            flatdata = self.read_flat()
            # See if number of points matches with event length:
            nevents_nohead = len(flatdata)//self.len
            resto_nohead   = len(flatdata)%self.len
            # See if number of bytes matches with event byte-length:
            bytes_per_event = 4*6 + self.len*2
            nbytes = os.path.getsize(self.FileName)
            nevents_head    = nbytes//bytes_per_event
            resto_head      = nbytes%bytes_per_event
            
            if resto_nohead == 0:
                data = {}
                for i in range(nevents_nohead):
                    time = None
                    bulk = [flatdata[i] for i in range(i*self.len,(i+1)*self.len)]
                    data[i] = (time,bulk)
                channel = self.ch
                print('  INFO. Reading without header')
            elif resto_head == 0:
                channel, data = self.read_whead()
                print('  INFO. Reading with header')
            else:
                print('  ERROR. Length of file does not match with events length')
        
        except:
            #channel, data = self.read_whead()
            print('  ERROR. Something was wrong when reading file')
        
        return channel, data
    
    
    def make_events(self,nevent):
        """
        
        Parameters
        ----------
        nevent : INTEGER
            Number of event to make its class.

        Returns
        -------
        An event as object.

        """
        channel, data = self.autoread()
        event = data[nevent]
        time = event[0]
        bulk = event[1]
        create_event = evm.Event(ADC_counts=bulk,channel=self.ch,dt=self.Dt,Vpp=self.vpp,stime=time)
        return create_event



class readChannels:
    
    def __init__(self, directory,preffix='wave', nbits=14,length=1030,dt=4,Vpp=2):
        """

        Parameters
        ----------
        directory : STRING
            Directory where the files with data from each channel are.
        
        preffix: STRING
            Fixed characters in the name of each file, before the number of channel. The default is 'wave'.
        
        Returns
        -------
        None.

        """
        
        cwd = os.getcwd()  # main current working directory
        if cwd in directory: os.chdir(directory)
        else: os.chdir(cwd+'/'+directory)
        list_dir = os.listdir()

        file_list = []
        for item in list_dir:
            if item[0:4] == 'wave':
                size = os.path.getsize(item)
                file_list.append((item,size/1000))
        file_list = sorted(file_list)
        for item in file_list: print('file %20s   size %12.1f KB' %item)
        
        DATA  = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[]}
        nch_list = []
        lch = []
        eventos = {}
        for fname0 in file_list:
            fname = fname0[0]
            nch = int(fname[4]) # for format wave[ch].dat or wave[ch].txt -!- change if filename format is different
            nch_list.append(nch)
            print(' Reading file ' + fname + '...')
            data0 = WDReader(fname,nch,nbits,length,dt,Vpp)
            channel, data = data0.autoread()
            eventos[channel] = data0 # saving WDReader object into dictionary
            lch.append(len(data))
            DATA[nch] = data
        os.chdir(cwd)
        print(' Working directory:', os.getcwd())
        
        self.all_data = DATA
        self.file_list = file_list
        self.leng_ch = lch
        self.eventsobj = eventos
        #return DATA
    
    def data(self):
        """

        Returns
        -------
        Dictionary with data.

        """
        datos = self.all_data
        return datos
    
    
    def nevents(self):
        """

        Returns
        -------
        Number of events in each channel.

        """
        assertion = True
        for lchi in self.leng_ch:
            if lchi != self.leng_ch[0]:
                assertion = False
                break
        if not assertion: print('Any length of channel does not match with the others. Check if files are strictly from different channels')
        return self.leng_ch[0]
    
    
    def files_list(self):
        """

        Returns
        -------
        List of files which fulfill the preffix condition.

        """
        lista = [item[0] for item in self.file_list]
        return lista
    
    
    def make_events(self,nevent):
        """

        Parameters
        ----------
        nevent : INTEGER
            Number of event.

        Returns
        -------
        Dictionary with result of WDReader.make_events.

        """
        eventos = {}
        for iev,evento in (self.eventsobj).items():
            eventos[iev] = evento.make_events(nevent)
        return eventos
    
        

"""
### Naive test
#fichero = WDReader('../DATA/antenas_210909/wave0.txt')
#fichero = WDReader('../DATA/ant_header_binascii_compar/wave0.dat')
fichero = WDReader('../DATA/waves_211018_1828_2h52/wave2.dat')
#datos   = fichero.read_flat()
#channel, datos = fichero.autoread()

#print(type(datos))
#print(datos)
#print(datos[0])

evento = fichero.make_events(0)

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
"""

"""
### Naive test 2
datos = readChannels('../DATA/waves_211018_1828_2h52')
DATA = datos.data()
nev  = datos.nevents()
lista= datos.files_list()
print(nev,lista)

datos.make_events(0)
"""