import pygplates
import numpy as np 
import pandas as pd
from sklearn.neighbors import KernelDensity
import paleogeography as pg
from pigplates import sphere_tools as pigsph

import sys
sys.path.append('/Users/Simon/GIT/GPlatesClassStruggle/')
from raster_reconstruction_classes import *


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


def indicator_z_distribution(df, static_polygons, rotation_model, grdfile, comparison_time,
							 longitude_field_name='lng', latitude_field_name='lat',
                             hist_type='kde', normalized=False):
    
    ptopo = GplatesRaster(grdfile)
    
    rla,rlo = reconstruct_dataframe(df, static_polygons, rotation_model, comparison_time,
								    longitude_field_name, latitude_field_name) 

    # sample grid at reconstructed points
    pbdb_topo = ptopo.sample(np.array(rlo), np.array(rla))
    
    if hist_type=='weighted':
        age_range = np.array(df.early_age - df.late_age)
        h,be = np.histogram(pbdb_topo,bins=40,range=(-1000,1000),weights=1./age_range)
        bc = (be[:-1] + be[1:]) / 2
    elif hist_type=='kde':
        bc=np.arange(-1000,3000,50)
        h = kde_sklearn(pbdb_topo,bc,bandwidth=100)

    if normalized:
        # normalise so that histogram bin counts always sum to 1
        h = h/np.sum(h)
        
    return bc,h


def plot_points_on_paleotopography(df, static_polygons, rotation_model, grdfile, comparison_time,
							       longitude_field_name='lng',latitude_field_name='lat',weighting=None):
    
    ptopo = GplatesRaster(grdfile)
    
    rla,rlo = reconstruct_dataframe(df, static_polygons, rotation_model, comparison_time,
								    longitude_field_name, latitude_field_name) 

    # sample grid at reconstructed points
    pbdb_topo = ptopo.sample(np.array(rlo), np.array(rla))
    
    # define some weights in range [0-1], based on the length of the valid age range for each sample
    age_range = np.array(df.early_age - df.late_age)
    weights = (1/age_range)/(1/age_range).max()
    #print weights
    
    ptopo.plot()
    if weighting=='alpha':
        # https://stackoverflow.com/questions/24767355/individual-alpha-values-in-scatter-plot-matplotlib
        rgba_colors = np.zeros((weights.shape[0],4))  # all zeros --> black
        # the fourth column needs to be your alphas
        rgba_colors[:, 3] = weights
        plt.scatter(rlo,rla,c=rgba_colors,edgecolors='',s=8)
    elif weighting=='size':
        plt.scatter(rlo,rla,c='k',edgecolors='',s=weights*20,alpha=0.5)
    elif weighting is None:
        plt.plot(rlo,rla,'k.')
    plt.show()
    
    
    
    
