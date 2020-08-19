# Author: Robert Beda
# Last Updated: August 18, 2020
# Script for downloading Time Series into local-directory files.
# Scroll to below function definitions to see user instructions.

# Necessary other files: Obsrun_Endtimes.txt, in this 
# script's directory.

import numpy as np
from gwpy.timeseries import TimeSeries
from gwpy.time import to_gps
from gwpy.segments import Segment
import os

#### Functions in this Code:
#    List_Runs(start, end)
#    Pair_Up(l, delind)
#    Seg_Split(start, end, frame)
#    Grab_Series(start, end, channel, frame, nproc)
#    Make_Dir(foldl, end)

def List_Runs(start, end):
    """Return the list of observing periods overlapping with 
    time segment (start, end).
    Arguments"""
    periods = []
    runends = open('Obsrun_Endtimes.txt', 'r').readlines()
    for i in range(0, 3):
        obsstart = int(runends[i][1:11])
        obsend = int(runends[i][16:26])
        if obsstart<int(end) and int(start)<obsend:
                periods.append(i+1)
    return periods


def Pair_Up(l, delind):
    """Return a list of tuples formed from adjacent items in list.
    Arguments:
    l -- List
    delind -- Index to be removed if l has an odd number of elements
    Returns:
    List of tuples formed from adjacent items in list, excluding l[delind]
    in the case of an odd number of elements"""
    length = len(l)
    tupled = []
    if length%2 == 0:
        for i in range(0, length//2):
            tupled.append((l[i*2], l[i*2+1])) #easy if the list is even
    else:
        true_delind = (delind+1)%length - 1 #allows for delind >= length
        l = l[0:true_delind] + l[true_delind+1:] #remove desired list entry
        for i in range(0, length//2):
            tupled.append((l[i*2], l[i*2+1])) #easy if the list is even
    return tupled


def Seg_Split(start, end, frame):
    """Return a list of time segments that combine to form the period between
    start and end, with information gaps excluded.
    Arguments:
    start -- All times in all output segments are at or after this time
    end -- All times in all output segments are at or before this time
    frame -- String such as 'L1_R' informed by desired observatory and 
             type of frame..
    
    Returns:
    segs -- a list of tuples, each representing a usable time segment."""
    empties = []
    obsruns = List_Runs(start, end)
    runends = open('Obsrun_Endtimes.txt', 'r').readlines()
    if len(obsruns)==0:
        return []
    if len(obsruns)>1:
        for i in range(len(obsruns)-1):
            empties.append(Segment(int(runends[obsruns[i]-1][16:26]), int(runends[obsruns[i+1]-1][1:11])))
    for run in obsruns:
        segstrings = open('Gaps/O{}-{}_Gaps.txt'.format(str(run), frame), 'r').readlines()
        for i in range(len(segstrings)):
            segstrings[i] = Segment(int(segstrings[i][1:11]), int(segstrings[i][16:26]))
        empties = empties+segstrings
    starts = [start]+[seg[1] for seg in empties]
    ends = [end]+[seg[0] for seg in empties]
    bookends = starts + ends
    bookends.sort()
    if bookends[-1] == end:
        bookends = bookends[bookends.index(start):]
    else:
        bookends = bookends[bookends.index(start):bookends.index(end)+1]
    if bookends[1] in starts:
        bookends = bookends[1:]
    if bookends[-2] in ends:
        bookends = bookends[:-1]
    return Pair_Up(bookends, 0)


def Grab_Series(start, end, channel, frame, nproc):
    """Return a list of Time Series representing the known Segments from within [start ... end)
    Arguments:
    start -- GPS time. All times in all output segments are at or after this time.
    end -- GPS time. All times in all output segments are at or before this time.
    channel -- String for the channel from which all Time Series in (start, end) will be obtained.
    frame -- String such as 'L1_R' informed by desired observatory and 
             type of frame.
    nproc -- Integer for the number of processors to use in reading each time series.
    
    
    Returns:
    The list of Time Series
    """
    modesl = []
    for ends in Seg_Split(start, end, frame):
        modes = TimeSeries.get(channel, ends[0], ends[1], frametype=frame, verbose=True, nproc=nproc)
        modesl.append(modes)
    return modesl


def Make_Dir(foldl, end):
    """Make a directory from the script's directory following the 
    order of entries in a list of folder names. Return to directory 'end'.
    If any directory names in foldl are already used, use those instead.
    Arguments:
    foldl -- A list of folder name strings whose order defines the filepath.
    end -- The directory to which the function brings the user.
    
    Returns:
    None"""
    for folder in foldl:
        entries = os.listdir(os.getcwd())
        if not folder in entries:
            os.system('mkdir {}'.format(folder))
        os.chdir(folder)
    if not end=='':
        os.chdir(end)
    return None 

# This script will save Time Series data within a specified time
# period in a directory named '{dirname}/{observing run}/{channel}'
# for a each of a set of given channel names. To use it, simply 
# tweak the values laid out in all lines through 163 as is 
# appropriate for your purposes.

# Input start and end times of desired Time Series files:
start = to_gps('Nov 1 2019 01:00:00')
end = to_gps('Mar 20 2020 00:00:00')

# Input name for directory to store files,
# or '' for storage in this script's directory:
dirname = 'Local_Data' 
# Input string for name of observing run containing 'start':
obsrun = 'O3b'
# Input frametype for desired Time Series data:
frame = 'L1_R'
# Input the desired number of parallel processors to assign to data retrieval:
procs = 10

# Input list of channels whose Time Series are to be saved.
chan1 = 'L1:ISI-HAM6_SENSCOR_X_FADE_TIME_LEFT_MON'
chan2 = 'L1:ISI-HAM6_SENSCOR_X_FADE_CUR_CHAN_MON'
chan3 = 'L1:ISI-GND_STS_ITMY_Z_BLRMS_3_10'
chan4 = 'L1:ISI-GND_STS_ITMY_Z_BLRMS_100M_300M'
chan5 = 'L1:PEM-EY_WIND_WEATHER_MPS'

chans = [chan1, chan2, chan3, chan4, chan5]

# No need to edit beyond here unless debugging or tinkering.
# Do so at your own risk!

for chan in chans:
    print(chan+'...')
    if dirname=='':
        Make_Dir([obsrun, chan], os.path.dirname(os.path.abspath(__file__)))
        for series in Grab_Series(start, end, chan, frame, procs):
            series.write('{}/{}/{}.hdf5'.format(obsrun, chan, str(series.t0.value)))
    else:
        Make_Dir([dirname, obsrun, chan], os.path.dirname(os.path.abspath(__file__)))
        for series in Grab_Series(start, end, chan, frame, procs):
            series.write('{}/{}/{}/{}.hdf5'.format(dirname, obsrun, chan, str(series.t0.value)))