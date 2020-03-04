r"""
Based on code from util.py osv_1.3\lib\util\

"""
import os
import sys
import datetime
import logging
import glob

import numpy as np
from astropy.io import fits

import clock
import path_utils

logging.getLogger('').handlers = []
logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG, filename = "gbm.log")

from gbm_download_daily import data_date

def hhmmss_to_sod(str_hhmmss):
    m = str_hhmmss.split(':')
    sod = int(m[0])*3600.0 + int(m[1])*60.0 + float(m[2])
    return sod 

def read_poshist(pos_file, verbose = True):
    '''
    Extract Quaternions, Position, Time & Geo Coordinates from file.
    Poshist files for days prior to March 2009 either have the spacecraft lat &
    lon set to zero or the fields are missing altogheter. This should be caught
    by the try except block in place.
    '''
    dtorad=180./math.acos(-1.)
    data=fits.getdata(pos_file,ext=1)
    nt=np.size(data)
    sc_time=data.SCLK_UTC
    sc_quat=np.zeros((nt,4),float)
    sc_pos=np.zeros((nt,3),float)
    sc_coords=np.zeros((nt,2),float)
    try:
        sc_coords[:,0]=data.SC_LON
        sc_coords[:,1]=data.SC_LAT
    except:
        if verbose:
            mes = ''
            mes += '*** No geographical coordinates available '
            mes += 'for this file: %s' %pos_file
            print(mes)
            
    sc_quat[:,0]=data.QSJ_1
    sc_quat[:,1]=data.QSJ_2
    sc_quat[:,2]=data.QSJ_3
    sc_quat[:,3]=data.QSJ_4
    sc_pos[:,0]=data.POS_X
    sc_pos[:,1]=data.POS_Y
    sc_pos[:,2]=data.POS_Z
    return sc_time,sc_pos,sc_quat,sc_coords

def read_pha(pha_file, gti = False, qualMask = True, tOffset = True,):
    """
    Extract Counts & Time From a GBM Pha File
    If gti is true then only data corresponding to gti is returned
    if qualMask is true then apply quality mask
    if tOffset is true then add tzero to time arrays
    
    13.01.10: Fixed a bug where gtis where calculated between time rather
    than time & endtime.
    """
    data = fits.open(pha_file)
    gtis = np.array((data[3].data.START,data[3].data.STOP)) 
    qual = data[2].data['QUALITY']
    if qualMask:
        qual = (qual == 0)
    else:
        qual = (qual != 99)
    time = data[2].data.TIME [qual]
    endtime = data[2].data.ENDTIME [qual]
    if gti:
        #User wants gtis
        if gtis.size:
            #gtis are valid
            for i in range(0,len(gtis[0])):
                print (gtis[0,i], gtis[1,i])
                print ("gti {:d}: {:f} {:f}".format(i, gtis[0,i], gtis[1,i]))
                temp_index=np.where((time >= gtis[0,i]) & (endtime <= gtis[1,i]))
                if i == 0: 
                    gti_index = temp_index[0]
                else: 
                    gti_index = np.concatenate((gti_index,temp_index[0]))
            gti_index = (gti_index,)
        else:       
            #gtis are not valid
            gti_index = (np.empty(0),) #np.where(time == -1)
    else:       
        #user does not want gtis, return all data
        gti_index = np.where(time != -1)
    #Temporal offset (Tzero4) may not exist - check
    tzero = 0
    if tOffset:
        if data[2].header.__contains__('TZERO4'):
            tzero = data[2].header['TZERO4']
        print ("tOffset: ", tzero)
    pha_counts = data[2].data['COUNTS'][qual]
    t_start    = data[2].data.TIME[qual]  + tzero
    t_end      = data[2].data.ENDTIME[qual] + tzero
    t_exposure = data[2].data.EXPOSURE[qual]
    eMin = data[1].data.E_MIN
    eMax = data[1].data.E_MAX
    data.close()
    if gti:
        return t_start, t_end, t_exposure, pha_counts, gti_index, eMin, eMax
    else:
        return  t_start, t_end, t_exposure, pha_counts,  eMin, eMax

def rebin_gbm(x, y, exp, err = [], resolution = [], trange = []):
    '''
    Rebin GBM CSPEC or CTIME data. Takes in x,y,exp, where x is the bin
    centre or bin edges, y is the counts array (bins*chan), err is an array
    of errors (same shape as counts array) and exp is the exposure per bin.
    
    Returns x1, y1, exp1, err, the rebinned data. 
    
    If resolution is not passed then the coarsest resolution from the data is 
    selected. If x1 is a 2xnbins array (i.e. the bin edges) then the same shape 
    array will be returned. An optional parameter trange can also be passed - 
    this is a 2x1 list which contains the edges of the data to be binned.
    
    Errors: If no error is passed then the statistical error is assumed to arise
    from counting error and is given by N^1/2 where N is the number of counts in 
    a bin. If an array of errors is passed then the error on a bin is found by
    summing in quadrature the errors of each bin contributing to it. This will
    fail if the desired resolution is lower than the native resolution of the 
    input data -> but you probably shouldn't be trying to resample the data in
    this case anyway.
    
    '''
    if resolution == []:
        resolution = exp.max().round(3)
    
    #Is x the bin centres or bin edges?
    #Irregardless we need to have edges for error calculation 
    #if user passed an err array
    if len(x.shape)!=1:
        #we have bin edges - get bin centres
        x_edges = x
        bin_edge = True
        x=(x[:,1]-x[:,0])/2 +x[:,0]
    else:
        x_edges=np.column_stack((x-exp/2 ,x+exp/2 ))
        bin_edge = False
    if trange == []:
        x1=np.arange(x[0],x[-1],resolution)
        #x1=np.arange(x[0,0],x[-1,1],resolution)
    else:
        x1=np.arange(trange[0],trange[1],resolution)
    
    #Define new arrays to house rebinned data
    nchan=y[0,:].size
    nbin=x1.size
    y1=np.zeros((nbin,nchan))
    
    # first interpolate exposure
    if bin_edge:
        binWidth = x_edges[:,1] - x_edges[:,0]        
        binWidth1 = np.ones(nbin)*resolution
        exp1 = np.interp(x1,x,exp/binWidth)*binWidth1
    else:
        # no quite correct if there is deadtime, but the 
        # best that can be done without bin edges
        exp1 = np.ones(nbin)*resolution

    for i in range(0,nchan):
        y1[:,i]=np.interp(x1,x,y[:,i]/exp)*exp1
    if err==[]:
        #Statistical error
        err1=np.sqrt(y1)
    else:
        err1=np.zeros((nbin,nchan))
        for i in range(0,nchan):
            err1[:,i]=np.interp(x1,x,err[:,i]/exp)*exp1
        
    if bin_edge:
        #If user passed bin edges return bin edges, otherwise return bin centres
        xi=x1 - resolution/2. #(exp1/2)
        xj=x1 + resolution/2. # (exp1/2)
        x1=np.column_stack((xi,xj))

    return x1,y1,exp1,err1

def get_interval(t_start, t_end, t_exposure, pha_counts, T_0, T_min, T_max):

    arr_bool = np.logical_and(t_start > T_0+T_min, t_end < T_0+T_max)
    return t_start[arr_bool]-T_0, t_end[arr_bool]-T_0, t_exposure[arr_bool], pha_counts[arr_bool]

def rebin(t_start, t_end, t_exposure, pha_counts, n_sum):

    n_max = t_start.size//n_sum * n_sum
    t_start_new = t_start[:n_max:n_sum]
    t_end_new = t_end[n_sum-1:n_max:n_sum]
    shape = t_exposure.shape[0]
    t_exposure_new = t_exposure[:n_max].reshape([shape//n_sum, n_sum]).sum(1)

    shape_cnts = pha_counts.shape 
    pha_counts_new = pha_counts[:shape_cnts[0]//n_sum*n_sum,:].reshape([shape_cnts[0]//n_sum, n_sum, shape_cnts[1]]).sum(1)

    return t_start_new, t_end_new, t_exposure_new, pha_counts_new

def write_acsii(ascii_file, channel_min, channel_max, t_start, t_end, t_exposure, pha_counts):
    """
    "{:s} {:f}\n".format(clock.fermi2utc(row[0]).strftime("%Y-%m-%d %H:%M:%S")
    """

    #t_start = np.array([hhmmss_to_sod(clock.fermi2utc(t).strftime("%H:%M:%S.%f")) for t in t_start])
    #t_end =  np.array([hhmmss_to_sod(clock.fermi2utc(t).strftime("%H:%M:%S.%f")) for t in t_end])
   
    f = open(ascii_file, 'w')
    f.write("Ti    Tf   dT Exposure ")
    for k in range(channel_min, channel_max+1):
        f.write("Ch{:d} ".format(k))
    f.write("\n")

    for i in range(t_start.size):
        f.write("{:10f} {:10f} {:8.3f} {:8.3f} ".format(t_start[i], t_end[i],t_end[i]-t_start[i], t_exposure[i]))
        for k in range(pha_counts[i].size):
           f.write("{:.0f} ".format(pha_counts[i,k]))
        f.write("\n")

    f.close()

def det2ascii(path, str_date, str_det, channel_min, channel_max, T0, resolution):

    #pha_file = 'glg_ctime_n{:s}_{:s}'.format(str_det, str_date)
    #pha_file = path_utils.get_files(path, pattern=pha_file, prefix=True, all=False)

    pha_file = 'glg_ctime_n{:s}_{:s}*pha'.format(str_det, str_date)
    pha_file = glob.glob(os.path.join(path,pha_file))[0]
    #print(pha_file)
    #return
    head, tail  = os.path.split(pha_file)    
    out_file = os.path.join(path, tail.split('.')[0])

    #pha_file = os.path.join(path, pha_file)
    #out_file = os.path.join(path, out_file)

    res_ctime = 0.256 #s
    n_sum = int(resolution/res_ctime)

    if resolution >=1:
        str_res = "{:d}s".format(int(resolution))
    else:
        str_res = "{:d}ms".format(int(resolution*1000))

    ascii_file = "{:s}_{:s}.txt".format(out_file, str_res)
    test_ascii_file = "{:s}_test.txt".format(out_file)

    t_start, t_end, t_exposure, pha_counts, gti_index, eMin, eMax = read_pha(pha_file, gti=True, qualMask=True, tOffset=False)
    print("Emin: ", eMin)
    print("Emax: ", eMax)
    pha_counts = pha_counts[:,channel_min-1:channel_max]

    #print "File {:s} contains data\nfrom {:s} to {:s}".format(pha_file, clock.fermi2utc(t_start[0]).strftime("%Y-%m-%d %H:%M:%S"), 
    #                clock.fermi2utc(t_end[-1]).strftime("%Y-%m-%d %H:%M:%S"))

    #write_acsii(ascii_file_256ms, t_start, t_end, t_exposure, pha_counts, T0=T0_GBM)

    t_start, t_end, t_exposure, pha_counts = get_interval(t_start, t_end, t_exposure, pha_counts, T0, -1000.0, 1000.0 )
    #write_acsii(test_ascii_file, channel_min, channel_max, t_start, t_end, t_exposure, pha_counts)

    t_start_new, t_end_new, t_exposure_new, pha_counts_new = rebin(t_start, t_end, t_exposure, pha_counts, n_sum)
    write_acsii(ascii_file, channel_min, channel_max, t_start_new, t_end_new, t_exposure_new, pha_counts_new)

def qw(s):
    return s.split()

def main():
    """
    Date: 180307

    Channel  1  2  3   4   5   6    7    8
    Emin:  4.5 12 27  50 102 294  538  983
    Emax:   12 27 50 102 294 538  983  2000
    """
    date = '20190902'
    T0_GBM = 11370.000 # s UT 

    path = './'+date

    lst_dets = qw('0 1 2 3 4 5 6 7 8 9 a b')
    #lst_dets = ['0',]

    #data = data_date(date[2:])
    #data.download_data(lst_dets, False, True, False)

    channel_min = 1
    channel_max = 7
    resolution = 30 #s
    
    for det in lst_dets:
        det2ascii(path, date[2:], det, channel_min, channel_max, T0_GBM, resolution)


main()