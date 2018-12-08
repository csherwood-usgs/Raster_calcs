# Raster_calcs
This started out as examples of working with Geotiff files, but became the defacto location for processing the Sandwich mapping results.
Basically, there are three types of files here:
* Scripts to process Sandwich raster data...many of which lead to useful final-ish products.
* Scripts to process other Sandwich data...again, many of which are useful analyses.
* Tests or examples where I figure out / remind myself of how to do something.

Below, I briefly list filenames with brief notes to remind me of what they do, and where they stand. All are Jupyter notebooks (.ipynb) files. They require packages in my slightly enhanced version of the IOOS environment. GM == GlobalMapper v.19 or v.20.

```v1 maps``` - These are the maps collected in 2017_Karen_Sandwich_maps and do not reflect any reprocessing of the photogrammetry.

```v2 maps``` - These are reprocessed maps with associated spatial difference maps. These have not been included yet.

## Sandwich raster calcs
* ```Calc_Sandwich_volumes``` - Calculate volumes from 0.1-m maps in UTM space. Uses upper/lower beach polygons in ```Volume_polygons.geojson``` that were created in GM. These use delz = 0.08. Uses rasterio.

* ```Calc_Sandwich_variance``` - Calculate s.d. and deviations for "fixed" points from 0.1-m maps in UTM space. Includes March 2018, still v1 maps. *Final boxplot needs to be cleaned up. Also, needs to import renamed March 2018 map.* Uses rasterio. Includes my bi-linear interpolation routine.

* ```Calc_Sandwich_slopes_and_roughness``` - Includes functions for calculating slope and roughness. Demonstrates use of `scipy.ndimage.filters.generic_filter`.*Does not include March 2018; v1 maps.* Uses rasterio. Writes slope `.tif` files.

## Sandwich map data processing
*```Calc_Profile_retreat_rates``` - Loads profiles extracted from Global Mapper and calculates profile volumes. Updated to include Mrch 2018 using v1 maps, but not yet a complete analysis. Writes out `.mat` and `.csv` files with volumes, areas, and errors. *Need to update error estimates to 1.96xsqrt([2x0.08]^2)* 

## Sandwich oceanographic data

```CDIP221_to_nc``` - Processes a mix of wave and water-level data. Work in progress. Used to send TWL estimates to Patrick before 2018 Fall AGU.





## Sandwich processing in other repos
```Analysis of alongshore transport using groin images``` - https://github.com/csherwood-usgs/GE_groin_analysis/blob/master/analyze_GE_groin_asymmetry.ipynb

## Other stuff
