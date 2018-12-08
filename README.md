# Raster_calcs
This started out as examples of working with Geotiff files, but became the defacto location for processing the Sandwich mapping results.
Basically, there are three types of files here:
* Scripts to process Sandwich raster data...many of which lead to useful final-ish products.
* Scripts to process other Sandwich data...again, many of which are useful analyses.
* Tests or examples where I figure out / remind myself of how to do something.

Below, I briefly list filenames with brief notes to remind me of what they do, and where they stand. All are Jupyter notebooks (.ipynb) files. They require packages in my slightly enhanced version of the IOOS environment. GM == GlobalMapper v.19 or v.20.

```v1 maps``` - These are the maps collected in 2017_Karen_Sandwich_maps and do not reflect any reprocessing of the photogrammetry.

```v2 maps``` - These are reprocessed maps with associated spatial difference maps. These have not been included yet.

## Key Sandwich raster calcs
* ```Calc_Sandwich_volumes``` - Calculate volumes from 0.1-m maps in UTM space. Uses upper/lower beach polygons in ```Volume_polygons.geojson``` that were created in GM. These use delz = 0.08. Uses rasterio.

* ```Calc_Sandwich_variance``` - Calculate s.d. and deviations for "fixed" points from 0.1-m maps in UTM space. Includes March 2018, still v1 maps. *Final boxplot needs to be cleaned up. Also, needs to import renamed March 2018 map.* Uses rasterio. Includes my bi-linear interpolation routine.

* ```Sandwich_maps_local_coords``` - Load maps, rotate to local coords, and combine into 3d xarray Dataframe called `one_meter_test.nc`. Uses xarray. Uses v1 maps. *Need to update directory location and rename 2018 maps.* 

* ```Sandwich_analyze_ROI``` - Work on reduced, rotated region of interest after loading `one_meter_test.nc`. Uses v1 maps. Uses xarray. *Uses dx,dy = 1 m.* Demonstrates interpolation onto new grid in rotated (alonshore, cross-shore) coordinates. Masks data using `barrier_roi.geojson` after rotating. *Wanders off into transect analysis...should move this to another notebook.*

* ```Sandwich_compare_DEM_and_transect``` - Compares all v1 Sandwich DEMs with transect points, and writes the .txt files used in next notebook. *Needs to upgraded to reflect movement of 2017_Karen_Sandwich_maps directory.* Otherwise, final analysis.

* ```Sandwich_summarize_all_DEM_minus_transect``` - Reads the DEM minus transect .txt files (in `2017_Karen_Sandwich_maps\\dem_trans`; based on v1 maps, and makes a nice histogram and map of diffs on top of mean elevation map. *Needs to upgraded to reflect movement of 2017_Karen_Sandwich_maps directory.* Otherwise, final analysis.

## Other Sandwich raster calcs (not part of final analysis thread)
* ```Calc_Sandwich_slopes_and_roughness``` - Includes functions for calculating slope and roughness. Demonstrates use of `scipy.ndimage.filters.generic_filter`. *Does not include March 2018; v1 maps.* Uses rasterio. Writes slope `.tif` files.



## Sandwich map data processing
* ```Calc_Profile_retreat_rates``` - Loads profiles extracted from Global Mapper and calculates profile volumes. Updated to include Mrch 2018 using v1 maps, but not yet a complete analysis. Writes out `.mat` and `.csv` files with volumes, areas, and errors. *Need to update error estimates to 1.96xsqrt([2x0.08]^2).* 
* ```GCP_Experiment``` - Loads transect data from 2017-01-25 and several DEMs made with various arrays of GCPs to evaluate best GCP configuration.
* ```GCP_analysis``` - Experiments in analyzing spatial distribution of GCPs from 2017-04-28. Demos application of `scipy.spatial.ConvexHull` and `.Delauney`, as well as `mpl_toolkits.mplot3d.Axes3D`. Kind of a dead end.

## Sandwich oceanographic data
* ```CDIP221_to_nc``` - Processes a mix of wave and water-level data. Work in progress. Used to send TWL estimates to Patrick before 2018 Fall AGU.
* ```Sandwich_timeline``` - Load oceangraphic data and make little time-series plots. *Need to add Waverider Hs and second Groin DWave data. Maybe merge with CDIP_to_nc TWL calcs?*


## Sandwich processing in other repos
* ```Analysis of alongshore transport using groin images``` - https://github.com/csherwood-usgs/GE_groin_analysis/blob/master/analyze_GE_groin_asymmetry.ipynb

## Other stuff
* ```Compare_DEM_and_transect.ipynb``` - Example of loading .tif file and transect points and calculating a histogram of differences. Include bilinear interpolation routine. Uses rasterio. Demonstrates addding columns to Pandas dataframe and saving as .csv with values for NaNs and a fixed format.
* ```Compress_byte_array``` - Messing around with `zlib` compression. Completely off-topic. *Should move elsewhere.*
* ```Coordinate_conversion``` - Example of how to use `pyproj` and EPSG info to do coordinate conversions. 
