## Purpose
The objective of this code is to calculate the light time (time of flight) between any given DSN ground station and a queried point on the Moon. This should also return the uncertainy in that delay. 

## Setup
Python3.10  
Ensure that the following folder path exists:   
lunar-delay/  
├── resources/  
│   ├── dem_to_spice.py  
│   ├── LDEM_80S_20M.tif  
│   └── kernels/  
│       ├── de440.bsp  
│       ├── dss_17_prelim_itrf93_190814.bsp  
│       ├── earth_200101_990827_predict.bpc  
│       ├── earth_assoc_itrf93.tf  
│       ├── earth_topo_20103.tf  
│       ├── earthstns_itrf93_201023.bsp  
│       ├── metakernel.tm  
│       ├── moon_de440_250416.tf  
│       ├── moon_pa_de440_200625.bpc  
│       ├── naif0012.tls  
│       └── pck00011.tpc  
├── lunar-delay.py/  
└── main.py/  

Due to file size restrictions on GitHub, two files will need to be downloaded seperately:  
LDEM_80S_20m.tif  
de440.bsp (https://naif.jpl.nasa.gov/pub/naif/pds/wgc/kernels/spk/de440.bsp)  


## Kernel Rationale
de440.bsp - SPK: Solar system model  
dss_17_prelim_itrf93_190814.bsp - PCK: DSS17 constants  
earth_200101_990827_predict.bpc - PCK: Earth's orientation data + constants  
earth_assoc_itrf93.tf - FK: Earth's frame kernel  
earthstns_itrf93_201023.bsp - PCK: Orientation data + constants for DSN stations  
earth_topo_201023.tf - FK: Frame kernel for the DSN stations  
moon_de440_220930.tf - FK: Lunar frame kernel  
moon_pa_de440_200625.bpc - PCK: Lunar orientation data + constants using the DE440 model  
naif0012.tls - LSK: Leapseconds (as of July 1, 2026)  
pck00011.tpc - PCK: Moon radii and fallback IAU orientation  

## Why not just use spice.spkpos(..., abcorr="LT", ...)?
SPICE has a built in "LT" aberration correction, but this has certain limitations. The main limitation is that this works between two vodies that have their own SPK segments, such as the Moon's center and the Earth's center. Since a surface point is being queried, it needs to be reconstructed everytime since is on a rotating and orbiting body. 
Functionally, this does the same thing as SPICE, but applied to the reconstructed point. 



