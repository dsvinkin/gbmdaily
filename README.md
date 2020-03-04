# gbmdaily
Scripts for downloading and plotting Fermi-GBM ctime data and
angles between a sky location and detector axes.

## Usage
Download ctime and poshist data using `gbm_download_daily.py`. 
Convert ctime to ASCII with `gbm_ctime2ascii.py` and extract incident angles with `gbm_poshist_data.py`
Plot the data using `plot_gbm_ctime_daily.py` and `plot_gbm_angles.py`

## Depends on
astropy

## Acknowledgments
Code is based on Fermi GBM Orbital Background Subtraction Tool (osv_1.3)
https://fermi.gsfc.nasa.gov/ssc/data/p7rep/analysis/user/
