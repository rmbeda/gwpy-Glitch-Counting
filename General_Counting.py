# Author: Robert Beda
# Last Updated: August 18, 2020
# Script for plotting glitch rates in various IFO configurations,
# with thresholds applied based on environmental trends.
# Scroll to below function definitions to see user instructions.

# Necessary other files: Tseries.py and consequently 
# Obsrun_Endtimes.txt, both in this script's directory.
import numpy as np
from gwpy.timeseries import TimeSeries
from scipy import stats
from gwpy.segments import DataQualityFlag
from gwpy.segments import Segment
from gwpy.segments import SegmentList
from gwtrigfind import find_trigger_files
from gwpy.table import EventTable
from gwpy.table.filters import in_segmentlist
from gwpy.plot import Plot
from gwpy.time import to_gps
from gwpy.time import from_gps 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import math

#### Functions in this Code:
#    List_Runs(start, end)
#    Pair_Up(l, delind)
#    Seg_Split(start, end, frame)
#    Grab_Series(start, end, channel, frame, nproc)
#    Grab_Sfiles(start, end, run, channel, frame)
#    Flag(modes, strict)
#    Flags(modes, lostrict)
#    Make_Flags(modes, strictls)
#    Group_Flags(modesl, strictls)
#    Net_Flag(modesll, strictlsl)
#    Net_Ends(seglist)
#    Next_End(start, segl)
#    Cut_Time_2(end, jumps, catsecs, cut=0)
#    Get_Rates_3(chunks, segs, verbose = False)
#    Sort_Vals(ts, source, catsecs)
#    Datetimegen(t)
#    Title(t0, tf, descript)
#    Make_Dir(foldl, end)
#    File_Name(t0, tf, info, filetype, extras = [])
#    Time_Scatter(times, deps, xname, yname, title, fname, flags, flagnames, 
#                 s=10, c='b', marker='o', figsize=[12,4])
#    Save_Hist(source, binsize, fname, xlabel='', title='')
#    Save_Hist_2(sources, binsizes, labels, fname, xlabel='', title='')
#    Save_Data(lists, t0, tf, info, extras=[])


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
    The list of Time Series."""
    modesl = []
    for ends in Seg_Split(start, end, frame):
        modes = TimeSeries.get(channel, ends[0], ends[1], frametype=frame, verbose=True, nproc=nproc)
        modesl.append(modes)
    return modesl


def Grab_Sfiles(start, end, run, channel, frame):
    """Return a list of Time Series pulled from files stored in this script's directory.
    Arguments:
    run -- String representing an observing run (ie O1, O2, O3a, etc).
    Otherwise as for Grab_Series().
    
    Returns:
    As for Grab_Series()."""
    modesl = []
    entries = os.listdir('Local_Data/{}/{}'.format(run, channel))
    for ends in Seg_Split(start, end, frame):
        for entry in entries:
            file = 'Local_Data/{}/{}/{}'.format(run, channel, entry)
            series = TimeSeries.read(file, channel)
            if series.times[0].value<=ends[0]<series.times[-1].value:
                modes = TimeSeries.read(file, channel, ends[0], ends[-1], verbose=True)
                modesl.append(modes)
                break
    return modesl


def Flag(modes, strict):
    """Return the DQ flags for a particular Times Series, 
    reflecting when a particular condition has been met.
    Arguments:
    modes -- Time Series to be investigated.
    strict -- A condition of form (mode, '{'max'/'min'/'eq'}{'eq'/''}') 
              that can be met for the corresponding time to be in an active segment in 'modes'.
    
    Returns:
    flag -- A DQFlag meeting the requested requirements."""
    del modes.unit
    if strict[1][0:3] == 'max':
        if len(strict[1])==3:
            return (modes<strict[0]).to_dqflag(name='')
        else:
            return (modes<=strict[0]).to_dqflag(name='')
    elif strict[1][0:3] == 'min':
        if len(strict[1])==3:
            return (modes>strict[0]).to_dqflag(name='')
        else:
            return (modes>=strict[0]).to_dqflag(name='')
    else:
        return (modes==strict[0]).to_dqflag(name='')

    
def Flags(modes, lostrict):
    """Return the DQ flags for a particular Times Series, 
    reflecting when a particular set of conditions has been met.
    Arguments:
    modes -- Time Series to be investigated.
    lostrict -- Contains a set of conditions that can be met
                for the corresponding time to be in an active segment in 'modes'.
    
    Returns:
    flag -- A DQFlag meeting the requested requirements."""
    del modes.unit
    if len(lostrict)==1:
        return Flag(modes, lostrict[0])
    else:
        return Flag(modes, lostrict[0]) & Flags(modes, lostrict[1:])

        
def Make_Flags(modes, strictls):
    """Return the DQ flags for a particular Time Series, for times
    when particular sets of modes are, or are not, active.
    Arguments:
    modes -- Time Series to be investigated.
    strictls -- Contains a list for each set of conditions that can be met
                for the corresponding time to be in an active segment in 'modes'.
    
    Returns:
    flag -- A DQFlag meeting the requested requirements."""
    del modes.unit # Prepare for conversion to dqflag
    if len(strictls)==1:
        return Flags(modes, strictls[0])
    else:
        return Flags(modes, strictls[0])|Make_Flags(modes, strictls[1:])

    
def Group_Flags(modesl, strictls):
    """Return the DQ flags for a set of Time Series from a common channel, 
    for times when particular sets of modes are, or are not, active.
    Arguments:
    modesl -- List of Time Series from a single channel to be investigated.
    strictls -- Contains a list for each set of conditions that can be met
                for the corresponding time to be in an active segment.
    
    Returns:
    flag -- A DQFlag meeting the requested requirements."""
    if len(modesl)==1:
        return Make_Flags(modesl[0], strictls)
    else:
        return Make_Flags(modesl[0], strictls)|Group_Flags(modesl[1:], strictls)

    
def Net_Flag(modesll, strictlsl):
    """Produces the DQ Flag containing the times within (start, end)
    during which each TimeSeries meets its restrictions with corresponding
    index in strictlsl.
    Arguments:
    seriesl -- List of TimeSeries lists, each with a channel common to all 
               entries, that must all meet their associated criteria in 
               strictlsl to establish an 'active' time. Use >=1 entries.
    strictlsl -- See Make_Flags to understand the structure of entries in this list. 
                 Use exactly one entry per channel.
    
    Returns:
    flag -- DQFlag representing the times when all named channels
            meet their associated criteria in 'strictlsl'"""
    if len(modesll) == 1:
        return Group_Flags(modesll[0], strictlsl[0])
    else:
        return Group_Flags(modesll[0], strictlsl[0])&Net_Flag(modesll[1:], strictlsl[1:]) 
    

def Net_Ends(seglist):
    """Returns the effective 'start' and 'end' times of a SegmentList
    with at least one entry.
    Arguments:
    seglist -- A SegmentList.
    
    Returns:
    beg -- The earliest time represented in 'seglist'.
    end -- All times represented in 'seglist' predate this."""
    beg = seglist[0][0]
    end = seglist[0][1]
    for seg in seglist:
        if seg[0]<beg:
            beg = seg[0]
        if seg[1]>end:
            end = seg[1]
    return beg, end


def Next_End(start, segl):
    """Find the segment end in a list of segments that is soonest
    after given time 'start'.
    Arguments:
    start -- Numerical value for time after which we want the next segment end.
    segl -- List of time segments to search.
    
    Returns:
    nextend -- The numerical value of the next segment end after 'start'"""
    bestdiff = False
    nextend = None
    for seg in segl:
        beg = float(seg[0])
        end = float(seg[1])
        if beg>start:
            if bestdiff==False or beg-start<bestdiff:
                nextend=beg
                bestdiff=nextend-start
            else:
                break
        elif end>start:
            if bestdiff==False or end-start<bestdiff:
                nextend=end
                bestdiff=nextend-start
            else:
                break
    return nextend


def Cut_Time_2(jumps, segs, cut=0):
    """Produces a list of start times for time segments represented in 'segs'
    that default to length 'jumps', but that are cut short by switches between
    Segments 'segs'.
    Arguments:
    jumps -- Default separation (in seconds) of entries in output time list.
    segs -- SegmentList containing all time segments to cut up.
    cut -- Minimum length of output segments, as a fraction of 'jumps'.
    
    Returns:
    partstarts -- A list of start times (of the LIGOTimeGPS datatype)
                  for time chunks that try to have 
                  length 'jump', but that are sometimes shorter or longer 
                  near the ends of Segments in the SegmentLists in 'catsecs'.
    """
    partstarts = []
    for seg in segs:
        start_i = seg[0]
        end_i = seg[-1]
        k = (end_i - start_i)/jumps
        l = math.floor(k)
        r = k - l 
        if k < cut:
            continue
        elif r >= cut:
            list_o_jumps = [start_i + n*jumps for n in range(l+1)]
            # print('R >= Cut last =', list_o_jumps[-1])
            # print('segend:', seg[-1])
            partstarts = partstarts + list_o_jumps
        elif r < cut: 
            list_o_jumps = [start_i + n*jumps for n in range(l)]
            # print('R < Cut last =', list_o_jumps[-1])
            # print('segend:', seg[-1])
            partstarts = partstarts + list_o_jumps
    partstarts.sort()
    return partstarts


def Get_Rates_3(chunks, segs, verbose = False):
    """Returns the glitch rates for a given set of time chunks
    defined by a list of start times, with an end time at the last entry.
    
    Arguments:
    chunks -- Sorted list of times representing the beginnings of the 
              time periods for which rate is to be calculated, with 'end' 
              tacked on.
    segs -- Ordered and non-overlapping SegmentList such that every 
            element in 'chunks' (except the last one) is in an entry in 
            'segs'.
    verbose -- Set to 'True' if you want to see the ends of each chunk in
               'chunks' printed as it is processed.
    
    Returns:
    normcounts -- A list of glitch rates (Hz) associated with each time
                  period represented in 'chunks'."""
    traced = False
    normcounts = []
    j = 0
    for i in range(len(chunks)-1):
        while not chunks[i] in segs[j]:
            j = j+1
        segend = segs[j][1]
        if chunks[i+1]>segend:
            chunkend = segend
        else:
            chunkend = chunks[i+1]
        if verbose:
            print(from_gps(chunks[i]), from_gps(chunkend))
        files = find_trigger_files('L1:GDS-CALIB_STRAIN', 'Omicron', chunks[i], chunkend)
        if len(files)>0:
            events = EventTable.read(files, format='ligolw', tablename='sngl_burst', 
                                     columns=['peak','peak_time_ns', 'peak_frequency', 'snr'])
            events = events[(events['peak']>=chunks[i]) & (events['peak']<chunkend)]  
            counts = len(events['peak'])
            length = chunkend - chunks[i]
            normcount = counts/(length)
            normcounts.append(normcount)
        else:
            normcounts.append(0)
        
    return normcounts


def Sort_Vals(ts, source, catsecs):
    """Sort values from a time-corresponding list into categories
    based on a given set of time segments.
    Arguments:
    ts -- List of times corresponding to entries in 'source'.
    source -- List of values to be sorted into approriate
              time-informed categories.
    catsecs -- List of SegmentLists denoting the different
               categories of times.
    
    Returns:
    bigbins -- List of lists with categorised items from 'source'
               in an order of categories according to the order 
               of lists in 'catsecs'."""
    bigbins=[[] for cat in catsecs]
    for i in range(len(ts)):
         for j in range(len(bigbins)):
            if ts[i] in catsecs[j]:
                bigbins[j].append(source[i])
    return bigbins


def Datetimegen(t):
    """Generate a date and a time in my own format.
    Arguments:
    t -- A GPS time (unit removed).
    
    Returns:
    A date of form DDMMYYYY and a time of form HH:MM:SS."""
    gpst = round(t, 2)
    gpst = from_gps(gpst)
    gpst = str(gpst)[0:22] #get a string to slice
    date = gpst[8:10]+gpst[5:7]+gpst[2:4] #so that now gpst[0:10] has form DDMMYYYY
    time = gpst[11:] #the time is everything after the date
    return date, time


def Title(t0, tf, descript):
    """Produce a plot title for a plot of data from
    within a certain period  of time.
    Arguments:
    t0, tf -- Times defining the segment [t0, tf) encompassed
              by plotted data.
    descript -- String describing the plotted informationn"""
    sd = Datetimegen(float(t0))[0]
    st = Datetimegen(float(t0))[1]
    ed = Datetimegen(float(tf))[0]
    et = Datetimegen(float(tf))[1]
    if sd==ed:
        return'{}_{}-{} {}'.format(sd, st, et, descript)
    else:
        return '{}-{} {}'.format(sd, ed, descript)


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
    

def Parse_End(date1, time1, date2, time2):
    """Formats the 'end date' portion of a file name for the 
    'File_Name()' function."""
    if date2[5:]==date1[5:]:
        date2 = date2[0:5]
        if date2[3:5]==date1[3:5]:
            date2 = date2[0:2]
            if date2[0:2]==date1[0:2]:
                return '{}'.format(time2)
    return '{}_{}'.format(date2, time2)


def File_Name(t0, tf, info, filetype, extras = []):
    """Generate a name and target directory for a file containing
    time information relating to a specific time period.
    Prepares filepath for storage.
    
    The name is a string of form 'Time1-DDMMYYYY2_Time2{figtype}.{filetype}',
    stored in directory LLO_EQ/YY/MM/DD/{extras}. 
    Portions of DDMMYYYY2 are omitted as makes sense based on shared 
    times between t0 and tf. Intended file format should match 'filetype'.
    Arguments:
    t0, tf -- Start and end times (respectively) of period represented
              by the file to be named.
    info -- String describing the format/type of information in file
            (eg 'Hist' for a histogram).
    filetype -- String describing the format of the file to be named.
    extras -- List of strings denoting  any extra desired file path steps.
    
    Returns:
    A directory and file name for a file representing data for
    the time period of Segment(t0, tf)."""
    date1, time1 = Datetimegen(float(t0))
    date2, time2 = Datetimegen(float(tf))
    filepath = [str(date1[4:]), str(date1[2:4]), str(date1[0:2])]+extras
    Make_Dir(filepath, os.path.dirname(os.path.abspath(__file__)))
    exstring = ''
    if extras != []:
        for file in extras:
            exstring = exstring+file+'/'
    return '{}/{}/{}/{}{}-{}{}.{}'.format(date1[4:], date1[2:4], date1[0:2], exstring, time1, Parse_End(date1, time1, date2, time2), info, filetype)


def Time_Scatter(times, deps, xname, yname, title, fname, flags, flagnames, 
            s=10, c='b', marker='o', figsize=[12,4]):
    """Make a scatter plot of some time-dependent data in lists, with arbitrary segment bars
    underneath.
    Arguments:
    times -- List of times to plot.
    deps -- List of time-dependent data values to plot
    xname -- String describing x-axis information
    yname -- String describing y-axis information
    flags -- A list of DQflags to plot in segment bars under the plot
    flagnames -- List of strings with which to label the segment bars.
                 Use an empty string to 'skip' an  entry.
    s -- Integer denoting plot marker size.
    c -- String for plot marker colour.
    marker -- String for plot marker shape.
    figsize -- Standard figsize keyword argument.
    
    Returns:
    None"""
    toplot = EventTable(data=[times, deps], masked=None, names=[xname, yname])
    plot = toplot.scatter(xname, yname, s=1, c='r', marker='s', figsize=[16,6], ylabel=yname, title=title)
    for i in range(len(flags)):
        plot.add_segments_bar(flags[i], label=flagnames[i])
    plot.savefig(fname, dpi=300, format='png', transparent=True)
    plot.close()

    
def Save_Hist(source, binsize, fname, xlabel='', title=''):
    """Save a histogram with a given single size of numerical bins for denoting
    categories of values in 'source'.
    Arguments:
    source -- List of information to be categorised.
    binsize -- Desired constant size of each category of data from 'source'.
    fname -- File Name under which the resulting figure is to be saved.
    xlabel -- Label for the categories defined on the horizontal axis.
    title -- Desired title of the saved output figure.
    
    Returns:
    None"""
    mini = min(source)
    mini = mini-mini%binsize
    maxi = max(source)
    maxi = maxi+binsize-maxi%binsize
    n, bins, patches = plt.hist(source, bins=int((maxi-mini)/binsize), 
                                range=(mini, maxi), facecolor='pink')
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title(title)
    plt.savefig(fname, dpi=300, format='png', transparent=True)
    plt.close()
    return (n, bins, patches)


def Save_Hist_2(sources, binsizes, labels, fname, xlabel='', title=''):
    """Save a histogram depicting a bunch of normalized distributions on the 
    same set of axes.
    Arguments:
    sources -- List of lists of information to be categorised.
    binsizes -- Desired constant sizes of each category of data for each entry
                in 'sources'.
    labels -- List of strings identifying the distributions in 'sources'.
    fname -- File Name under which the resulting figure is to be saved.
    xlabel -- Label for the categories defined on the horizontal axis.
    title -- Desired title of the saved output figure.
    
    Returns:
    None"""
    plt.figure()
    for source in sources:
        binsize = binsizes[sources.index(source)]
        mini = min(source)
        mini = mini-mini%binsize
        maxi = max(source)
        maxi = maxi+binsize-maxi%binsize
        n, bins, patches = plt.hist(source, bins=int((maxi-mini)/binsize), alpha=0.5, 
                                    density = True, label=labels[sources.index(source)])
    plt.xlabel(xlabel)
    plt.ylabel('Normalized Counts')
    plt.title(title)
    plt.legend()
    plt.savefig(fname, dpi=300, format='png', transparent=True)
    plt.close()
    return (n, bins, patches)


def Save_Data(lists, t0, tf, info, extras=[]):
    """Save each list in a tuple of time-related numerical lists
    to a .txt file using np.savetxt.
    
    Arguments:
    lists -- Tuple containing lists of numerical data.
    t0 -- Start of time period represented by data.
    tf -- end of time period represented by data.
    info -- Tuple of strings, each with some information about 
            the data in the list at corresponding index in
            'lists'.
    extras -- Argument for the 'File_Name()' function.
    
    Returns:
    None"""
    for i in range(len(lists)):
        array = np.array(lists[i])
        np.savetxt(File_Name(t0, tf, info[i], 'txt', extras=extras), array)
    return None
    
    
# This script will create a directory 
# '{year in 3rd millenium}/{month #}/{day #}/{skip_load_direct}'
# based on the 'start' time it is given. This directory will store resulting
# data files, plots, and statistical results.

# To use the script, You have three options:
# 1) Uncomment lines 679, 681, and 693, and comment out lines 680, 682, and 694, 
#    then run the script (slow).
# 2) Save Time Series for all relevant channels using the 'Tseries.py' script 
#    in this script's directory, and then run this script using
#    'start' and 'end' values appropriate to available Time Series files.
# 3) If you would like to make plots and run KS analysis on data already
#    obtained using this script, simply set 'skip_load' to 'True' in line 655
#    before running.

# Additionally, edit code from here through line 663, and within the range
# of lines 677-707, to tweak values such as relevant start and end times, 
# frametypes to be used for Time Series, the name of a directory to house output 
# files, names and channels for desired detector state categories, and relevant 
# channels for applying environmental thresholds.

# Input start and end times for desired glitch rate plot/distribution:
start = to_gps('Nov 1 2019 01:00:00')
end = to_gps('Mar 20 2020 00:00:00')
# Name the observing run to which search is limited:
obsrun = 'O3b'
# Name the frametype for all Time Series. 
# Note that the script has not been tested with minute trends.
frame = 'L1_R'
# Input the desired number of parallel processors to assign to data retrieval:
procs = 10
# Default length of time over which glitch rates are to be averaged
# (should be more than one minute if using minute trends):
jumps = 30
# Minimum fraction of 'jumps' between a time chunk's start and the end of
# the its 'parent' segment.
cutoff = 0.5
showprog = True # Set to 'True' if you want to see the times whose associated rate
                 # calculations are in progress.
skip_load = True # Set to 'True' if you want to use the secondary intermediate files
                  # that will let you skip glitch rate generation. Use for repeated runs.
# Define the name of the directory within the directory
# '{year in 3rd millenium}/{month #}/{day #}' that will store script-generated files.
skip_load_direct = 'Wind_Thresh_5'
tablet = ['Observing', 'Transition', 'EQ'] # The list of names for configurations
                                                # defined by their respective entries in 'configsecs'
                                                # below.
histbinwidths = [0.1, 0.1, 0.1] # Desired bin width for glitch rate categories in 
                                # each configuration's histogram plot.

if skip_load:
    partstarts = list(np.loadtxt(File_Name(start, end, 'partstarts', 'txt', extras = [skip_load_direct])))
    rates = list(np.loadtxt(File_Name(start, end, 'glitchrates', 'txt', extras = [skip_load_direct])))
    nomflags = DataQualityFlag.read(File_Name(start, end, '{}_dqflag'.format(tablet[0]), 'hdf5', extras = [skip_load_direct]))
    transflags = DataQualityFlag.read(File_Name(start, end, '{}_dqflag'.format(tablet[1]), 'hdf5', extras = [skip_load_direct]))
    EQflags = DataQualityFlag.read(File_Name(start, end, '{}_dqflag'.format(tablet[2]), 'hdf5', extras = [skip_load_direct]))
    configs = [nomflags, transflags, EQflags]
    configsecs = [config.active for config in configs]
else:
    # Category-specific cutoffs
    print('Grabbinng cat-cutoff info...')
    indchan1 = 'L1:ISI-HAM6_SENSCOR_X_FADE_TIME_LEFT_MON'
    indchan2 = 'L1:ISI-HAM6_SENSCOR_X_FADE_CUR_CHAN_MON'
    # transmodes = Grab_Series(start, end, indchan1, frame, procs)
    transmodes = Grab_Sfiles(start, end, obsrun, indchan1, frame)
    # EQmodes = Grab_Series(start, end, indchan2, frame, procs)
    EQmodes = Grab_Sfiles(start, end, obsrun, indchan2, frame)

    transtricts = [[(0, 'min')]]
    EQstricts = [[(5, 'eqeq')],[(6,'eqeq')],[(7, 'eqeq')]]

    # Widely-applied cutoffs
    print('Grabbinng env-cutoff info...')
    envchan1 = 'L1:ISI-GND_STS_ITMY_Z_BLRMS_3_10' # anthro
    envchan2 = 'L1:ISI-GND_STS_ITMY_Z_BLRMS_100M_300M' # micro
    envchan4 = 'L1:PEM-EY_WIND_WEATHER_MPS' # wind
    envchans = [envchan1, envchan2, envchan4]
    # env_series = [Grab_Series(start, end, chan, frame, procs) for chan in envchans]
    env_series = [Grab_Sfiles(start, end, obsrun, chan, frame) for chan in envchans]

    envstricts1 = [[(500, 'maxeq')]]
    envstricts2 = [[(1000, 'maxeq')]]
    envstricts4 = [[(5, 'maxeq')]]
    env_strictsl = [envstricts1, envstricts2, envstricts4]

    # Observing times
    observing = DataQualityFlag.query('L1:DMT-ANALYSIS_READY:1', start, end)
    base_flags = observing&Net_Flag(env_series, env_strictsl)
    transflags = base_flags&Group_Flags(transmodes, transtricts)
    EQflags = base_flags&Group_Flags(EQmodes, EQstricts)-transflags
    nomflags = base_flags-(EQflags|transflags)
    configs = [nomflags, transflags, EQflags]
    
# No need to edit beyond here unless debugging or tinkering.
# Do so at your own risk!

    for i in range(len(configs)):
        configs[i].write(File_Name(start, end, '{}_dqflag'.format(tablet[i]), 'hdf5', extras = [skip_load_direct]))
    
    configsecs = [config.active for config in configs]
    
    # Conglomerate configsecs to define an all-inclusive SegmentList:
    allsegs = []
    for i in range(len(configsecs)):
        allsegs = allsegs+[seg for seg in configsecs[i]]
    allsegs.sort()
    
    # Get the times that start the 'time chunks', with the end tacked on:
    partstarts = Cut_Time_2(jumps, allsegs, cut=cutoff)
    partstarts.append(end)
    partstarts = [float(partstart) for partstart in partstarts]

    # Get glitch rates with indices associated with times in partstarts:
    rates = Get_Rates_3(partstarts, allsegs, verbose = showprog)

    # Save the data from the analysis:
    Save_Data((partstarts, rates), start, end, ('partstarts', 'glitchrates'), extras=[skip_load_direct])

# Make a scatter plot:
Time_Scatter(partstarts[:-1], rates, 'Glitch Times', 'Glitch Rates', 
             Title(start, end, 'Glitch Rate'), File_Name(start, end, 'Glitches', 'png', extras=[skip_load_direct]),
             configs, tablet, 
             s=1, c='r', marker='s', figsize=[16,6])

# Bin glitch rates according to categories of their associated times:
bins = Sort_Vals(partstarts[:-1], rates, configsecs)

# Make a histogram using 'bins', save histogram data:
for i in range(len(bins)):
    if bins[i] == []:
        print('No {} Mode Glitches!'.format(tablet[i]))
        hist_vals.append(None)
        continue
    hist_contents = Save_Hist(bins[i], histbinwidths[i],
                              File_Name(start, end, '{}_Hist'.format(tablet[i]), 'png', extras = [skip_load_direct]),
                              xlabel='Rates',
                              title='{} Mode Distribution of Glitch Rates Averaged over ~{}s intervals'.format(tablet[i], jumps))

Save_Hist_2(bins, histbinwidths, tablet, File_Name(start, end, 'Comp_Hist', 'png', extras = [skip_load_direct]), xlabel='Rates', title='Distributions of Glitch Rates Averaged over ~{}s intervals'.format(jumps))

# run K-S tests comparing the first configuration category with each of the others:
KS_val1, pval1 = stats.ks_2samp(bins[0], bins[1])
KS_val2, pval2 = stats.ks_2samp(bins[0], bins[2])

file = open(File_Name(start, end, 'Stat_Tests', 'txt', extras = [skip_load_direct]), 'w')
file.write('{}/{} KS:{} \n'.format(tablet[0], tablet[1], KS_val1))
file.write('{}/{} pval:{} \n'.format(tablet[0], tablet[1], pval1))
file.write('{}/{} KS:{} \n'.format(tablet[0], tablet[2], KS_val2))
file.write('{}/{} pval:{} \n'.format(tablet[0], tablet[2], pval2))
file.close()
