import pygplates
import numpy as np 
import pandas as pd
from sklearn.neighbors import KernelDensity
import paleogeography as pg
from pigplates import sphere_tools as pigsph


def reconstruct_dataframe(df,static_polygons,rotation_model,reconstruction_time,
                          longitude_field_name='lng',latitude_field_name='lat'):
    
    point_features = []
    # put the points into a feature collection, using Lat,Long coordinates from dataframe
    for index,row in df.iterrows():
        point = pygplates.PointOnSphere(float(row[latitude_field_name]),float(row[longitude_field_name]))
        point_feature = pygplates.Feature()
        point_feature.set_geometry(point)
        point_features.append(point_feature)

    # The partition points function can then be used as before
    partitioned_point_features = pygplates.partition_into_plates(static_polygons,
                                                                 rotation_model,
                                                                 point_features) 
    
    reconstructed_point_features = []
    pygplates.reconstruct(partitioned_point_features,rotation_model,
                          reconstructed_point_features,reconstruction_time)
    
    rlat = []
    rlon = []
    for point in reconstructed_point_features:
        rlat.append(point.get_reconstructed_geometry().to_lat_lon()[0])
        rlon.append(point.get_reconstructed_geometry().to_lat_lon()[1])
    
    return rlat,rlon
    

def kde_sklearn(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)


def indicator_z_distribution(df,static_polygons,rotation_model,comparison_time,
							 longitude_field_name='lng',latitude_field_name='lat'):

    rla,rlo = reconstruct_dataframe(df,static_polygons,rotation_model,comparison_time,
								   longitude_field_name,latitude_field_name)

    # interpolate heights at grid
    grdfile = '../paleotopography/paleotopo_grids/paleotopobathy_smooth_%0.2fMa.nc' % comparison_time

    topo_smoothX,topo_smoothY,topo_smoothZ = pg.load_netcdf(grdfile)

    topo_smoothXg, topo_smoothYg = np.meshgrid(topo_smoothX,topo_smoothY)
    d,l = pigsph.sampleOnSphere(topo_smoothXg.flatten(),
                                topo_smoothYg.flatten(), 
                                topo_smoothZ.flatten(),
                                np.array(rlo),
                                np.array(rla),
                                k=4)


    # based on http://earthpy.org/interpolation_between_grids_with_ckdtree.html
    w = 1./d**2
    pbdb_topo = np.sum(w * topo_smoothZ.flatten().ravel()[l],axis=1) / np.sum(w,axis=1)

    #age_range = np.array(df.early_age - df.late_age)

    #h,be = np.histogram(pbdb_topo,bins=40,range=(-1000,1000),weights=1./age_range)
    #bc = (be[:-1] + be[1:]) / 2
    bc=np.arange(-1000,3000,50)
    h = kde_sklearn(pbdb_topo,bc,bandwidth=100)

    # normalise so that histogram bin counts always sum to 1
    #h = h/np.sum(h)
        
    return bc,h