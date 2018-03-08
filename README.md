# paleotopography
pygplates code for reconstructing paleotopography

Requires the 'Paleogeography' repo, for both some of the python dependencies

Notebooks as follows:

- tween-paleoshorelines-inplace: works on the paleogeography polygons and generates tweened multipoints. NOTE: this code produces multipoints where the shorelines are tweened to shift gradually within the area covered by the original paleogeography polygons - BUT - it does not deal with any gaps and overlaps at certain reconstruction times in between the original paleogeography snapshot times. This gets handled in the next notebook.....
- tweened-paleotopography-clean: works on the multipoints generated from the above code and makes paleotopography grids. Includes a step to fill the gaps and overlaps that appear for time between the reconstruction time snapshot
#

- tweened-paleotopography-clean-dev.ipynb includes some testing
- 'watershed' and 'paleotopography-anisotropic-diffusion' contain some testing of image processing to smooth the 
paleotopography in different ways
