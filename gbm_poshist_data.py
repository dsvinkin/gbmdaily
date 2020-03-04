r"""
Based on code from util.py osv_1.3\lib\util\

"""

import os
import datetime
import time

import numpy as np
from astropy.io import fits

import clock

def calc_angles(sc_time, sc_pos, sc_quat, src_ra, src_dec):
    """
    Calculate GBM source angles & pointing
    
    Modified 28.11.11 to also calculate the angles for BGO
    
    Note: 06.12.11 I think the loop where the angles are calculated may be 
    superfluous. It this is the case, then getting rid of this loop should
    speed up this function significantly.
    """
    dtorad=180./np.arccos(-1.)
    nt=np.size(sc_time)
    #Calculate Direction Cosines
    scx=np.zeros((nt,3),float)
    scx[:,0]=(sc_quat[:,0]**2-sc_quat[:,1]**2-sc_quat[:,2]**2+sc_quat[:,3]**2)
    scx[:,1]=2.*(sc_quat[:,0]*sc_quat[:,1] + sc_quat[:,3]*sc_quat[:,2])
    scx[:,2]=2.*(sc_quat[:,0]*sc_quat[:,2] - sc_quat[:,3]*sc_quat[:,1])
    scy=np.zeros((nt,3),float)
    scy[:,0]=2.*(sc_quat[:,0]*sc_quat[:,1] - sc_quat[:,3]*sc_quat[:,2])
    scy[:,1]=(-sc_quat[:,0]**2+sc_quat[:,1]**2-sc_quat[:,2]**2+sc_quat[:,3]**2)
    scy[:,2]=2.*(sc_quat[:,1]*sc_quat[:,2] + sc_quat[:,3]*sc_quat[:,0])
    scz=np.zeros((nt,3),float)
    scz[:,0]=2.*(sc_quat[:,0]*sc_quat[:,2] + sc_quat[:,3]*sc_quat[:,1])
    scz[:,1]=2.*(sc_quat[:,1]*sc_quat[:,2] - sc_quat[:,3]*sc_quat[:,0])
    scz[:,2]=(-sc_quat[:,0]**2-sc_quat[:,1]**2+sc_quat[:,2]**2+sc_quat[:,3]**2)

    #Calculate Source coordinates
    source_pos=np.zeros((3),float)

    fra=src_ra/dtorad
    fdec=src_dec/dtorad
    source_pos[0]=np.cos(fdec)*np.cos(fra)
    source_pos[1]=np.cos(fdec)*np.sin(fra)
    source_pos[2]=np.sin(fdec)
    sdotprod=source_pos / np.sqrt(source_pos[0]*source_pos[0]+source_pos[1]*source_pos[1]+source_pos[2]*source_pos[2])
    sc_source_pos=np.zeros((nt,3),float)
    sc_source_pos[:,0]=scx[:,0]*source_pos[0]+scx[:,1]*source_pos[1]+scx[:,2]*source_pos[2]
    sc_source_pos[:,1]=scy[:,0]*source_pos[0]+scy[:,1]*source_pos[1]+scy[:,2]*source_pos[2]
    sc_source_pos[:,2]=scz[:,0]*source_pos[0]+scz[:,1]*source_pos[1]+scz[:,2]*source_pos[2]

    #Define Fermi/GBM detector Geometries
    #The following coordinates are taken from Meegan et al., 2009
    det_zen=[20.58, 45.31, 90.21, 45.24, 90.27, 89.79, 20.43, 
             46.18, 89.97, 45.55, 90.42, 90.32, 90, 90]
    det_az=[45.89, 45.11, 58.44, 314.87, 303.15,  3.35, 224.93,
             224.62, 236.61, 135.19, 123.73, 183.74, 0, 180]
    dets=np.array(['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb','b0','b1'])
    ndet = 14
    det_index = np.arange(0,ndet)
    dtorad_arr = np.ones((ndet))*dtorad
    det_zen = det_zen / dtorad_arr
    det_az = det_az / dtorad_arr
    
    #Calculate Detector unit vectors
    det_unit = np.zeros((ndet,3),float)
    det_unit[:,0] = np.sin(det_zen[:])*np.cos(det_az[:])
    det_unit[:,1] = np.sin(det_zen[:])*np.sin(det_az[:])
    det_unit[:,2] = np.cos(det_zen[:])
    
    distfromz = np.zeros((nt),float)
    distfromgeo = np.zeros((nt),float)    
    distfromdet = np.zeros((nt,ndet),float)
    for i in range(0, nt):
        dotprod = np.zeros((3),float)
        dotprod[0] = -sc_pos[i,0]/ np.sqrt(sc_pos[i,0]*sc_pos[i,0]+sc_pos[i,1]*sc_pos[i,1]+sc_pos[i,2]*sc_pos[i,2])
        dotprod[1] = -sc_pos[i,1]/ np.sqrt(sc_pos[i,0]*sc_pos[i,0]+sc_pos[i,1]*sc_pos[i,1]+sc_pos[i,2]*sc_pos[i,2])   
        dotprod[2] = -sc_pos[i,2]/ np.sqrt(sc_pos[i,0]*sc_pos[i,0]+sc_pos[i,1]*sc_pos[i,1]+sc_pos[i,2]*sc_pos[i,2])

        zdotprod = np.zeros((3),float)
        zdotprod[0] = scz[i,0]/ np.sqrt(scz[i,0]*scz[i,0]+scz[i,1]*scz[i,1]+scz[i,2]*scz[i,2])
        zdotprod[1] = scz[i,1]/ np.sqrt(scz[i,0]*scz[i,0]+scz[i,1]*scz[i,1]+scz[i,2]*scz[i,2])
        zdotprod[2] = scz[i,2]/ np.sqrt(scz[i,0]*scz[i,0]+scz[i,1]*scz[i,1]+scz[i,2]*scz[i,2])

        distfromgeo[i]= dtorad* np.arccos(dotprod[0]*sdotprod[0]+dotprod[1]*sdotprod[1]+dotprod[2]*sdotprod[2])
        distfromz[i] =  dtorad* np.arccos(sdotprod[0]*zdotprod[0]+sdotprod[1]*zdotprod[1]+sdotprod[2]*zdotprod[2])      
        distfromdet[i,:] = dtorad*np.arccos(det_unit[:,0]*sc_source_pos[i,0]+det_unit[:,1]*sc_source_pos[i,1]+det_unit[:,2]*sc_source_pos[i,2])        

    return distfromz, distfromgeo, distfromdet

def read_poshist(pos_file, verbose = True):
    '''
    Extract Quaternions, Position, Time & Geo Coordinates from file.
    Poshist files for days prior to March 2009 either have the spacecraft lat &
    lon set to zero or the fields are missing altogheter. This should be caught
    by the try except block in place.
    '''

    dtorad=180./np.arccos(-1.)
    data=fits.getdata(pos_file, ext=1)

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

class Regions:
    '''
    Determine the time of selection regions
    '''
    def __init__(self, zero_met, tmin, tmax, offset, 
                 orbit_period = 5737.70910239):

        self.zero_met = zero_met
        self.period = orbit_period 
        ranges = {}
        shifts = {}
        for i in offset:
            if i == 'src':
                tempShift = 0
            else:
                tempShift = orbit_period * float(i)
            shifts.update({i: tempShift})

        signs = {'pre': -1, 'pos': +1, 'src': 0}   
     
        for i in shifts.keys():
            if i != 'src':
                loop = ['pre','pos']
            else:
                loop = ['src']
            for range in loop:
                if range == 'src':
                    index = range
                else:
                    index = range+i
                temp = np.array([zero_met + tmin + signs[range] * shifts[i], 
                      zero_met + tmax + signs[range] * shifts[i]])
                ranges.update({index: temp}) 
        self.ranges = ranges
        self.offset = offset

    def __str__(self):
        mes = "Offsets for selection Periods\n" 
#        mes = mes + "Zero offset: "\n"
#        mes = mes + "Negative offset (14 orbits): " + str(self.pre_14) + "\n"
#        mes = mes + "Positive offset (14 orbits): " + str(self.pos_14) + "\n"
#        mes = mes + "Negative offset (16 orbits): " + str(self.pre_16) + "\n"
#        mes = mes + "Positive offset (16 orbits): " + str(self.pos_16) + "\n"
#        mes = mes + "Negative offset (30 orbits): " + str(self.pre_30) + "\n"
#        mes = mes + "Positive offset (30 orbits): " + str(self.pos_30) + "\n"
        return mes

class Poshist_data:
    def __init__(self,pos_files):
        '''
        Read in data for a list of POSHIST Files 
        '''
        for i in pos_files:
            poshist_data = read_poshist(i, verbose = False)
            if i == pos_files[0]:
                self.sc_time = poshist_data[0]
                self.sc_pos = poshist_data[1]
                self.sc_quat = poshist_data[2]
                self.sc_coords = poshist_data[3]
            else:
                self.sc_time = np.concatenate((self.sc_time, poshist_data[0]))
                self.sc_pos = np.concatenate((self.sc_pos, poshist_data[1]))
                self.sc_quat = np.concatenate((self.sc_quat, poshist_data[2]))
                self.sc_coords = np.concatenate((self.sc_coords,
                                                 poshist_data[3]))
        self.period = None
        self.rises = None
        self.sets = None
        self.occTI = None
        
    def calculate_angles(self, regions, ra, dec):
        '''
        Calculate detector angles for region of interest. Result is stored in 
        a dictionary, the keys of which are the detectors
        '''
        self.ang_region = regions
        ranges = regions.ranges
        print ("Time range: ", ranges['src'][0], ranges['src'][1])

        offset = regions.offset
        print ("offset: ", offset)
        # First calculate angles during ROI
        dur_indices = ((self.sc_time >ranges['src'][0]) & 
            (self.sc_time < ranges['src'][1]))
        self.dur_pointing, pointingGeo, distfromdet = calc_angles(
                                    self.sc_time[dur_indices],
                                    self.sc_pos[dur_indices],
                                    self.sc_quat[dur_indices],
                                    ra, dec)
        self.dur_t = self.sc_time[dur_indices]
        pointing = { 'src': self.dur_pointing}
        times = {'src': self.dur_t}
        for i in offset:
            if i == 'src':
                continue
            for j in ['pre', 'pos']:
                tRange = ranges[j + i]
                boolIndex = ((self.sc_time > tRange[0]) &
                              (self.sc_time < tRange[1]))
                t = self.sc_time[boolIndex]
                pointTemp,pointingGeo, detAngles = calc_angles(self.sc_time[boolIndex],
                                                        self.sc_pos[boolIndex],
                                                        self.sc_quat[boolIndex],
                                                        ra, dec)
                times.update({j + i: t})
                pointing.update({j + i: pointTemp})
        self.pointing = pointing
        self.times = times
        dets = np.array(['n0','n1','n2','n3','n4','n5','n6','n7','n8',
                         'n9','na','nb','b0','b1'])
        det_angles = {}
        for ang, det in zip(distfromdet.transpose(), dets):
            det_angles.update({det: ang})

        self.det_angles = det_angles
        self

    def calc_period(self):
        '''
        Calculate orbital period of Fermi. Assumes circular motion, not quite
        correct.
        '''
        G = 6.67428e-11    # m^3 kg^-1 s^-2
        M = 5.9722e24      # kg Mass Earth  
        r = (np.sum(self.sc_pos**2., 1))**(1/2.)
        r_avg = np.average(r)
        r_cubed = (r_avg)**3.
        factor = r_cubed/(G*M)
        period = 2. * np.pi * np.sqrt(factor)
        self.period = period

    def get_gti(self):
        '''
        Determine what NaI, BGO detectors have angles <60, <90 respectively.
        '''
        gtis = {}
        for det in self.det_angles:
            ang = self.det_angles[det]
            #BGO & NaI have different criteria for good detector selecitons
            if det =='b0' or det == 'b1':
                good_ang = 90
            else:
                good_ang = 60
            bool_list = ang < good_ang
            gti = util.make_gti(self.dur_t, bool_list)
            if gti == ([],[]):
                gti = None
            gtis.update({det: gti})
        self.gti = gtis

    def get_steps(self, ra, dec):
        ''' 
        Get Occultation Step Times and determine corresponding time intervals
        '''
        rises, sets = util.calc_occ_steps(ra, dec, self.sc_time, self.sc_pos)
        self.rises, self.sets = rises, sets
        
        # Now lets make the time intervals (TIs) corresponding to the times when 
        # the source is occulted. There are several possibilities
        # 1) There are an equal number of rises and sets
        # 2) There are more sets than rises
        # 3) There are more rises than sets
        # If we have 1), and the first step occurs before the first rise then
        # we can simply fill the TIs with the sets & rises. If however, the
        # first rise occurs before the first step, we have to stick the earliest
        # possible time onto the start of the set array.
        # If we have 2), then the first set should be less than the first rise,
        # and we can simply add the last possible time to the rise list
        # If we have 3), then the first rise should be less than the first set,
        # and we can simply add the first possible time to the set list.

        if rises.size == sets.size:
            if sets[0] < rises[0]:
                occI = list(sets)
                occJ = list(rises)
            else:
                occI = [self.sc_time[0]]
                occI.extend(list(sets))
                occJ = list(rises)
                occJ.extend([self.sc_time[-1]])
        elif sets.size > rises.size:
            if sets[0] < rises[0]:
                occI = list(sets)[:-1]
                occJ = list(rises)                
            else:
                mes = '*** More sets than rises, yet first rise occurs before \
                    first set - shouldn\'t be possible. Likely to crash soon!'
                print (mes)
        elif sets.size < rises.size:
            if rises[0] < sets[0]:
                occI = list(sets)
                # Ignore the first rise
                occJ = list(rises)[1:]
            else:
                mes = '*** More rises than sets, yet first set occurs before \
                    first rise - shouldn\'t be possible. Likely to crash soon!'
                print (mes)

        self.occTI = occI, occJ

def print_angles(data, angle_file):

    dets = ['n0','n1','n2','n3','n4','n5','n6','n7','n8',
                         'n9','na','nb']

    f = open(angle_file,'w')
    f.write('T ')
    for det in dets: 
        f.write("{:s} ".format(det))
    f.write("\n")

    T0 = data.ang_region.zero_met

    for i in range(data.dur_t.size):
        f.write("{:f} ".format(data.dur_t[i]-T0))
        for det in dets:
            f.write("{:8.2f} ".format(data.det_angles[det][i]))
        f.write("\n")

    f.close()        

def main():
   
   date = '20190902'
   T0 = 11370.000 # s UT 
   ra, dec =  279.472820, +61.497984 # Transient RA Dec, deg

   path = './'+date

   poshist_file = '{:s}/glg_poshist_all_{:s}_v00.fit'.format(path, date[2:])
   angle_file = '{:s}/poshist_angles_{:s}_{:05d}.txt'.format(path, date[2:], int(T0))  

   data = Poshist_data((poshist_file,))

   timestruct = clock.parsetime(date) + datetime.timedelta(seconds=T0)
   print(clock.parsetime(date))
   print("Date and time: ", timestruct)
   region = Regions(clock.utc2fermi(timestruct), -500, 500, ['src',])

   data.calculate_angles(region, ra, dec)
   
   print_angles(data, angle_file)

if __name__ == '__main__':
    main()
