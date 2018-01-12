# paleotopography
pygplates code for reconstructing paleotopography

Requires the 'Paleogeography' repo, for both some of the python dependencies

Notebooks as follows:

--- tween-paleoshorelines-inplace: works on the paleogeography polygons and generates tweened multipoints
--- tweened-paleotopography-clean: works on the multipoints generated from the above code and makes paleotopography grids


--- tweened-paleotopography-clean-dev.ipynb includes some testing
--- 'watershed' and 'paleotopography-anisotropic-diffusion' contain some testing of image processing to smooth the 
paleotopography in different ways
