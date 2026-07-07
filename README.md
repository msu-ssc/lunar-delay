## Purpose
The objective of this code is to calculate the light time (time of flight) between any given DSN ground station and a queried point on the Moon. This should also return the uncertainy in that delay. 

## Limitations
### LDEM Limitations
Since the LDEM's used are only for points below 80S, it should only be used for this points. The LDEMs also only have a 20m resolution. 
### DSS Limitations
DSS Stations not included in 810-005, 301, Rev. M, Table 7 must be updated in  DSN_Station_Location_Uncertainties.csv. Stations should also be checked if they are included in the kernels. (earthstns_itrf93_201023.bsp). 
DSS-17 is added to the CSV and a seperate kernel was made for it. 
### Time limitations
Individual kernel coverage (union across all objects in the file): 
<pre>   
  de440.bsp                                [1549-12-30 23:59:18.815 -> 2650-01-24 23:58:50.815]    
  dss_17_prelim_itrf93_190814.bsp          [1949-12-31 23:59:18.816 -> 2049-12-31 23:58:50.816]     
  earth_200101_990827_predict.bpc          [2020-01-01 00:00:00.000 -> 2099-08-27 00:00:00.000]     
  earthstns_itrf93_201023.bsp              [1949-12-31 23:59:18.816 -> 2149-12-31 23:58:50.816]     
  moon_pa_de440_200625.bpc                 [1549-12-30 23:59:18.815 -> 2650-01-24 23:58:50.815]       

Time range covered by ALL kernels (safe for computation):     
  2020-01-01 00:00:00.000  ->  2049-12-31 23:58:50.816     
</pre>

## Setup
Python3.10  
Ensure that the following folder path exists:  
<pre> 
lunar-delay/  
├── resources/  
│   ├── dem_to_spice.py  
│   ├── LDEM_80S_20M.tif    
│   ├── LDEM_80S_20MPP_ADJ_ERR.tiff     
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
</pre>
Due to file size restrictions on GitHub, these files will need to be downloaded seperately:  
LDEM_80S_20m.tif (https://pgda.gsfc.nasa.gov/products/90)  
LDEM_80S_20MPP_ADJ_ERR.tiff (https://pgda.gsfc.nasa.gov/products/90)   
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
SPICE has a built in "LT" aberration correction, but this has certain limitations. The main limitation is that this works between two bodies that have their own SPK segments, such as the Moon's center and the Earth's center. Since a surface point is being queried, it needs to be reconstructed everytime since it is on a rotating and orbiting body. 
Functionally, this does the same thing as SPICE, but applied to the reconstructed point. 

## Uncertainty Rationale
### Ground stations
For all ground stations, except DSS-17, their positional uncertainty is pulled from Table 7 in this document (DSN No. 810-005, 301, Rev. M): https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301M.pdf. 
Since DSS-17 is not included in this table, another source had to be used. I struggled to find documented numbers for DSS-17, so I assumed a larger error bar of 50m in all directions.

### LDEM Values
Thankfully, there is an associated LDEM map that gives the associated uncertainties per pixel here, LDEM_80S_20MPP_ADJ_ERR.tiff (https://pgda.gsfc.nasa.gov/products/90). Since this is a discrete map, the same interpolation is used for this error map as the height map, bilinear interpolation. This provides a weighted sum between the 4 corners of the queried point, providing height from LDEM_80S_20m.tif and error from LDEM_80S_20MPP_ADJ_ERR.tiff. 
