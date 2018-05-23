# paleotopography
pygplates code for reconstructing paleotopography

Requires the 'Paleogeography' repo https://github.com/siwill22/paleogeography for both data and code dependencies

Notebooks as follows:

#### tween-paleoshorelines-inplace: 
works on the paleogeography polygons and generates tweened multipoints. NOTE: this code produces multipoints where the shorelines are tweened to shift gradually within the area covered by the original paleogeography polygons - BUT - it does not deal with any gaps and overlaps at certain reconstruction times in between the original paleogeography snapshot times. This gets handled in the next notebook.....
#### tweened-paleotopography-clean: 
works on the multipoints generated from the above code and makes paleotopography grids. Includes a step to fill the gaps and overlaps that appear for time between the reconstruction time snapshot
#### make_present_day_environments
The youngest paleogeography in the sequence of Golonka/Cao++ represents 2 to 11 Ma (nominal age 6 Ma). To create a continuous time series to present day, we need to make a representation of present day that follows the same format (landmass, mountains, shallow marine and slope). How to do this is user-dependent - for example, if the user wants to use this to link to some paleotopography generated from the older paleogeographies, then the parameters used to generate these elevations (and at what spatial resolution) need to be considered to ensure some kind of self-consistency. This notebook allows present-day paleogeography polygons to be regenerated according to individual user-requirements. 
#### file_reformat
a short notebook to reformat between different formats, for example because gpmlz causes problems in parallelised processes
#

# dev folder
Some experimental notebooks that are even less documented:
- tweened-paleotopography-clean-dev.ipynb includes some testing
- 'watershed' and 'paleotopography-anisotropic-diffusion' contain some testing of image processing to smooth the 
paleotopography in different ways
- some very preliminary intermodel comparison, looking mainly at hypsometry 
