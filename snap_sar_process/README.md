## Download and Processing SAR data using SNAP 
This folder provides scripts for downloading and pre-processing SAR data


### Install 
```
   # pre-processing (GDAL will automatically be installed)
   conda create -n sar python=3.9
   conda activate sar
   conda install -c conda-forge asf_search
   conda install -c conda-forge geopandas
   conda install -c conda-forge rasterio
   conda install -c conda-forge psutil
   conda install -c conda-forge scikit-image
   conda install -c conda-forge sentineleof  # download orbit files https://github.com/scottstanie/sentineleof
   

   
```
    
    