#!/usr/bin/env python

'''
Script to download continuous data from Fermi GBM FTP

Downloads CTIME/CSPEC/POSHIST Files 
'''
import os
import sys
import argparse
import logging
import ftplib
from ftplib import FTP, FTP_TLS


def parse_args():
    # define parent parser which contains common arguments
    parentParser = argparse.ArgumentParser(add_help = False)
    parentParser.add_argument('date',type = str,  
                              help = 'Date to download - should be in YYMMDD format',
                              action = 'store',)
    parentParser.add_argument('--ctime', help = 'Download CTIME data',
                              action = 'store_true', default = False)
    parentParser.add_argument('--cspec', help = 'Download CSPEC data',
                              action = 'store_true', default = False)
    parentParser.add_argument('--poshist', help = 'Download POSHIST data',
                              action = 'store_true', default = False)
    parentParser.add_argument('--dets', help = "Detectors", type =str, nargs='*')

    args = parentParser.parse_args()
    return args

class data_date:
    def __init__(self, date,):
        '''Define date & ftp'''
        yr = '20' + date[0:2]
        mt = date[2:4]
        dy = date[4:6]
        self.ftp_dir = "fermi/data/gbm/daily/{:s}/{:s}/{:s}/current".format(yr, mt, dy)

    def nlst(self, ftp, str_pattern):
        files = []
        try:
            files = ftp.nlst(str_pattern)
        except ftplib.error_temp as resp:
            if str(resp) == "450 No files found":
                logging.warning("No {:s} files in this directory".format(str_pattern))
            else:
                logging.error("FTP error: {:s}".format(resp))
        return files

    def download_data(self, dets, cspec=False, ctime=False, poshist=False, path=''):
        '''Login to FTP & Download Data'''
        server = 'heasarc.gsfc.nasa.gov' #'legacy.gsfc.nasa.gov'

        if not (cspec | ctime | poshist):
            logging.error("No arguments passed - returning ...")
            return None, None, None

        logging.info("Connecting ...")
        ftp = FTP_TLS(server)
        ftp.login()
        ftp.prot_p()
        logging.info("Connected to {:s}".format(server))

        try:
            ftp.cwd(self.ftp_dir)
        except ftplib.error_perm:
            logging.error("Cannot cd to folder {:s}! Date may be incorrect.".format(self.ftp_dir))
            return None, None, None

        ctime_files = self.nlst(ftp, 'glg_*ctime*pha')
        cspec_files = self.nlst(ftp,'glg_*cspec*pha')
        poshist_files = self.nlst(ftp,'glg_*poshist*fit')

        if ctime:
            logging.info("Downloading CTIME")
            for name in ctime_files:
                det = name.split("_")[2]
                if (dets == None) or det in (dets):
                    logging.info("Downloading {:s}".format(name))
                    ftp.retrbinary('RETR ' + name, open(os.path.join(path, name),'wb').write)
        if cspec:
            logging.info("Downloading CSPEC")
            for name in cspec_files:
                det = name.split("_")[2]
                if (dets == None) or det in (dets):
                    logging.info("Downloading {:s}".format(name))
                    ftp.retrbinary('RETR ' + name, open(os.path.join(path, name),'wb').write)
        if poshist:
            logging.info("Downloading POSHIST")
            for name in poshist_files:
                logging.info("Downloading {:s}".format(name))
                ftp.retrbinary('RETR ' + name, open(os.path.join(path, name),'wb').write)
        ftp.quit()

        return ctime_files, cspec_files, poshist_files

def main():
    
    date = '20200318'

    path = './'+date
    
    data = data_date(date[2:])
    data.download_data(None, cspec=False, ctime=True, poshist=True, path=path)

if __name__ == '__main__':
    main()
