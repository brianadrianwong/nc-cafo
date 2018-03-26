# North Carolinian CAFO Change Detection in Earth Engine
Used a combination of Sentinel-2 and Landsat-7's archive for change detections of North Carolinian concentrated animal feeding operations (CAFOs). Many thanks to Earth Engine for hosting the aforementioned archives. Analysis on such a large stack of imagery would not be possible otherwise.

Demo below was depoloyed for all of North Carolina. All in all, the workflow processed 1,497 Landsat-7 images in the final change-detection algorithm. The valid-mask creation from Sentinel-2 used 3,156 images from January 1, 2017 to December 31, 2017.

<br />

### Raw Scene Example:
Starting point exampled region of interest.
![alt](../master/images/nc_cafo_change_detection_demo_01.png?raw=true "Raw Scene")

<br />

### NDVI:
Used Sentinel-2 for NDVI due to higher resolution specs. Goal here was to isolate low-NDVI spaces that were most likely to be anthropogenic structures.
![alt](../master/images/nc_cafo_change_detection_demo_02.png?raw=true "Raw Scene")

<br />

### CAFO Buffer Overlay:
![alt](../master/images/nc_cafo_change_detection_demo_03.png?raw=true "Raw Scene")

<br />

### Clipped to Buffer Overlay:
![alt](../master/images/nc_cafo_change_detection_demo_04.png?raw=true "Raw Scene")

<br />

### CAFO Buffer Overlay:
Thresholded the imagery for final step in creating a valid-mask across the state. Cumuluative steps up to this point are essentially a simplified rural roof detection method.
![alt](../master/images/nc_cafo_change_detection_demo_05.png?raw=true "Raw Scene")


### Change Detection:
Temporal Analysis for the normalized difference of Tasseled-Cap Brightness & Tasseled-Cap Greenness. Have not seen this index before so currently referring to it as the Normalized TCB-TCG (NDTCBG) Index. Alternatively, also calling it the Rural Built-Up Index (RBUI). The primary observation of interest for using this normalized difference was because of confounding spectral signatures between rooftops and barren or non-vegetated fields, however, larger Tasseled Cap Greenness (TCG) responses of barren or non-vegetated fields relative to CAFO structures provided reason to exploit this relationship with a normalized difference approach.

![alt](../master/images/nc_cafo_change_detection_demo_06.png?raw=true "Raw Scene")
