import pygplates
import glob, re
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import xarray as xr

import polygon_processing as pp
import paleogeography as pg
import paleogeography_tweening as pgt

from proximity_query import *
from create_gpml import create_gpml_regular_long_lat_mesh
import points_in_polygons
from sphere_tools import sampleOnSphere
import points_spatial_tree



def write_xyz_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


# define a function that loads paleogeography multipoints at a specified time
# NOTE this time can be anything, not a time where the multipoints fit nicely together,
# hence the gaps and overlaps will be present
def add_reconstructed_points_to_xyz(points_file,rotation_model,reconstruction_time,zval,mask_file=None):
    
    reconstructed_points = []
    pygplates.reconstruct(points_file,rotation_model,reconstructed_points,reconstruction_time)
    
    mask_file = './OrogenMasks.gpml'
    if mask_file is not None:
        reconstructed_masks = []
        pygplates.reconstruct(mask_file,rotation_model,reconstructed_masks,reconstruction_time)

        for reconstructed_point in reconstructed_points:
            for reconstructed_mask in reconstructed_masks:
                if reconstructed_mask.get_reconstructed_geometry().is_point_in_polygon(reconstructed_point.get_reconstructed_geometry().get_centroid()):
                    reconstructed_point.get_feature().set_name(reconstructed_mask.get_feature().get_name())
        
    point_array = []
    zval_array = []
    for reconstructed_point in reconstructed_points:
        feature_coordinates = reconstructed_point.get_reconstructed_geometry().to_lat_lon_array()
        point_array.append(feature_coordinates)
        #if reconstructed_point.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('OrogenicBelt'):
        if reconstructed_point.get_feature().get_name() == 'OrogenicBelt':
            zval_array.append(np.ones((feature_coordinates.shape[0],1))*zval*3)
        elif reconstructed_point.get_feature().get_name() == 'Cordillera':
            zval_array.append(np.ones((feature_coordinates.shape[0],1))*zval*2)
        else:
            zval_array.append(np.ones((feature_coordinates.shape[0],1))*zval)
      
    xy_array = np.vstack(point_array)
    xyz_array = np.hstack((xy_array,np.vstack(zval_array)))
    
    return xyz_array


# function to facilitate the smoothing of topography
# at the edge of mountain range polygons
def get_distance_to_mountain_edge(point_array,reconstruction_basedir,rotation_model,time,area_threshold):
    
    distance_threshold_radians=None
    env_list = ['m']

    if time==0:
        pg_dir = './present_day_paleogeography.gmt'
        pg_features = pg.load_paleogeography(pg_dir,env_list,single_file=True,env_field='Layer')
    else:
        pg_dir = '%s/PresentDay_Paleogeog_Matthews2016_%dMa/' % (reconstruction_basedir,time)
        pg_features = pg.load_paleogeography(pg_dir,env_list)

    cf = pp.merge_polygons(pg_features,rotation_model,time=time,sampling=0.25)
    sieve_polygons_t1 = pp.polygon_area_threshold(cf,area_threshold)

    polygons_as_list = []
    for feature in sieve_polygons_t1:
        polygons_as_list.append(feature.get_geometry())
        
    res1 = find_closest_geometries_to_points([pygplates.PointOnSphere(point) for point in zip(point_array[:,0],point_array[:,1])],
                                             polygons_as_list,
                                             distance_threshold_radians = distance_threshold_radians)
    
    distance_to_polygon_boundary = np.degrees(np.array(zip(*res1)[0]))

    # Make a copy of list of distances.
    distance_to_polygon = list(distance_to_polygon_boundary)

    # Set distance to zero for any points inside a polygon (leave other points unchanged).
    res2 = points_in_polygons.find_polygons([pygplates.PointOnSphere(point) for point in zip(point_array[:,0],point_array[:,1])],
                                            polygons_as_list)

    for point_index, rpolygon in enumerate(res2):
        # If not inside any polygons then result will be None.
        if rpolygon is None:
            distance_to_polygon[point_index] = 0.0
            
    return distance_to_polygon


# This cell uses COB Terranes to make a masking polygon
# (which is called 'seive_polygons')
def get_merged_cob_terrane_polygons(COBterrane_file,rotation_model,reconstruction_time,
                                    sampling,area_threshold=None,return_raster=False):

    polygon_features = pygplates.FeatureCollection(COBterrane_file)

    cobter = pp.force_polygon_geometries(polygon_features)

    cf = pp.merge_polygons(cobter,rotation_model,time=reconstruction_time,sampling=sampling)
    
    if area_threshold is not None:
        sieve_polygons = pp.polygon_area_threshold(cf,area_threshold)
        return sieve_polygons

    else:
        return cf

# This cell uses COB Terranes to make a masking polygon
# (which is called 'seive_polygons')
def get_merged_cob_terrane_raster(COBterrane_file,rotation_model,reconstruction_time,
                                  sampling):

    polygon_features = pygplates.FeatureCollection(COBterrane_file)

    cobter = pp.force_polygon_geometries(polygon_features)

    mask = pp.merge_polygons(cobter,rotation_model,time=reconstruction_time,
                             sampling=sampling,return_raster=True)
    
    return mask


# use merged seive_polygons to get a regular lat-long multipoint that will contain points
# only within the COB Terranes (ie not within the 'deep ocean')
def get_land_sea_multipoints(sieve_polygons,sampling,depth_for_unknown_ocean,subdivision_depth=4):

    multipoints = create_gpml_regular_long_lat_mesh(sampling)
    grid_dims = (int(180/sampling)+1,int(360/sampling)+1)

    for multipoint in multipoints:
        for mp in multipoint.get_all_geometries():
            points = mp.to_lat_lon_point_list()

    #reconstructed_polygons = []
    #pygplates.reconstruct(cobter,rotation_model,reconstructed_polygons,reconstruction_time)

    rpolygons = []
    for polygon in sieve_polygons:
        if polygon.get_geometry():
            rpolygons.append(polygon.get_geometry())

    polygons_containing_points = points_in_polygons.find_polygons(points, rpolygons, subdivision_depth=subdivision_depth)

    lat = []
    lon = []
    zval = []

    lat_deep = []
    lon_deep = []
    zval_deep = []

    for pcp,point in zip(polygons_containing_points,points):
        if pcp is not None:
            lat.append(point.get_latitude())
            lon.append(point.get_longitude())
        else:
            lat_deep.append(point.get_latitude())
            lon_deep.append(point.get_longitude())
            zval_deep.append(depth_for_unknown_ocean)
            
    #plt.figure(figsize=(25,11))      
    #plt.plot(lon,lat,'.')
    
    #plt.figure(figsize=(25,11))      
    #plt.plot(lon_deep,lat_deep,'.')
    
    #plt.figure(figsize=(25,11))  
    #for polygon in rpolygons:
    #    plt.plot(polygon.to_lat_lon_array()[:,1],
    #             polygon.to_lat_lon_array()[:,0])
        
            
    return lat,lon,zval,lat_deep,lon_deep,zval_deep


