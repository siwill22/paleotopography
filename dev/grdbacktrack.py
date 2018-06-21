import pygplates
from raster_reconstruction_classes import *
import points_spatial_tree
import points_in_polygons
import tectonic_subsidence as ts


###################
rho_mantle = 3300
rho_water = 1030
phi = 0.56
c = 4.5
#q = 2500


def sample_seafloor_age_model(static_polygon_features, agegrid, healpix_resolution):

    equal_area_points = PointDistributionOnSphere(distribution_type='healpix',N=healpix_resolution)

    points = [point.to_lat_lon_point() for point in equal_area_points.multipoint.get_points()]


    spatial_tree_of_uniform_recon_points = points_spatial_tree.PointsSpatialTree(points)

    recon_static_polygons = []
    recon_static_polygon_plate_ids = []
    for static_polygon_feature in static_polygon_features:
        recon_plate_id = static_polygon_feature.get_reconstruction_plate_id()
        recon_conjugate_plate_id = static_polygon_feature.get_conjugate_plate_id()
        recon_polygon = static_polygon_feature.get_geometry()

        recon_static_polygon_plate_ids.append((recon_plate_id,recon_conjugate_plate_id))
        recon_static_polygons.append(recon_polygon)

    point_plate_pairs = points_in_polygons.find_polygons_using_points_spatial_tree(
        points, spatial_tree_of_uniform_recon_points, recon_static_polygons, recon_static_polygon_plate_ids)
    
    point_ages = agegrid.sample_using_gmt(equal_area_points.longitude,equal_area_points.latitude)
    
    return points, point_plate_pairs, point_ages


def return_conjugate_points(points, point_ages, point_plate_pairs, target_plate_pair, rotation_model):

    points_and_conjugate_points = []
    for point, point_age, plate_ids in zip(points, point_ages, point_plate_pairs):

        if plate_ids is None:
            continue
        elif plate_ids[0] not in target_plate_pair:
            continue
        elif plate_ids[1] not in target_plate_pair:
            continue
        else:
            
            if not np.isnan(point_age):

                finite_rotation = rotation_model.get_rotation(float(point_age), int(plate_ids[0]), 0, int(plate_ids[1]))

                reconstructed_point = finite_rotation * pygplates.PointOnSphere(point)

                points_and_conjugate_points.append([point.to_lat_lon()[1],
                                                    point.to_lat_lon()[0],
                                                    reconstructed_point.to_lat_lon()[1],
                                                    reconstructed_point.to_lat_lon()[0],
                                                    plate_ids[0],
                                                    point_age])
            
    return points_and_conjugate_points


def get_unloaded_bsmt_depth(lons,lats,topography,sediment_thickness):
    
    point_depth = topography.sample_using_gmt(lons, lats)    
    
    point_sed_thick = sediment_thickness.sample_using_gmt(lons, lats)  
    
    rhoSbar = ts.AverageSedimentDensity(point_sed_thick, phi=phi, c=c)
    
    q = ts.AverageDensityAboveBasement(rhoSbar,point_sed_thick,point_depth)
    
    delta_z = (rho_mantle*point_sed_thick - q) / (rho_mantle - rho_water)
    
    unloaded_bsmt_depth = point_depth - delta_z
    
    return unloaded_bsmt_depth

