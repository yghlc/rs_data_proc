## Downloading and pre-processing of remote sensing data 
This repo provides scriptos for downloading and pre-processing different sources of remote sensing data

### Downloading and pre-processing planet images
    
    planetScripts/download_planet_img.py: given an extent and date range, 
    it will download Planet images. For the large area (such as Tibet), 
    you should divie the entire region to many samll grids (50 by 50 km).
    
    planetScripts/get_planet_image_list.py: the the metadata and coverage of
    downloade Planet images.     

### Install notes
Some functions in this repo depend on [DeeplabforRS](https://github.com/yghlc/DeeplabforRS.git),
please clone them:

    git clone https://github.com/yghlc/DeeplabforRS.git ~/codes/PycharmProjects/DeeplabforRS

    
    