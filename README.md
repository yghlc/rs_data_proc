## Downloading and pre-processing of remote sensing data 
This repo provides scripts for downloading and pre-processing different sources of remote sensing data

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


## Citation
If this repository is useful for your projects, please consider citing our papers:

```
@article{huang2023identifying,
  title={Identifying active retrogressive thaw slumps from ArcticDEM},
  author={Huang, Lingcao and Willis, Michael J and Li, Guiye and Lantz, Trevor C and Schaefer, Kevin and Wig, Elizabeth and Cao, Guofeng and Tiampo, Kristy F},
  journal={ISPRS Journal of Photogrammetry and Remote Sensing},
  volume={205},
  pages={301--316},
  year={2023},
  publisher={Elsevier}
}
```
