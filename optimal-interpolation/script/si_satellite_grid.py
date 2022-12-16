'''
######################################################################
## SCRIPT:           si_satellite_grid.py
## LANGUAGE:         python
## DESCRIPTION:      1. Generate SI grid
#######################################################################
'''

import time
import parser
import argparse
from os import environ, path
from datetime import datetime
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from itertools import product
from calendar import monthrange

from core_functions import data_dl, create_ncfile, apply_interp_to_grid, rr_total_file_reader, haversine, correlation_model, domain_dict, map_type_dict

def main():
        
    # Start timing of code.

    start = time.time()

    # Command-Line input
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--date', dest='date_string', nargs='?', help='Date', default=None)
    args = parser.parse_args()

    # Setup parameters
    STATION_DIR = environ['station_dir']
    config_dir = environ['config_dir']
    grid_dir = environ['grid_dir']
    superobs_status = environ['SUPEROBS_STATUS']
    MAP_TYPE = environ['MAP_TYPE']
    EXTENT = environ['EXTENT']
    PERIOD = environ['ACCUM_MON']
    station_number_cutoff = int(environ['station_number_cutoff'])
    correlation_cutoff = float(environ['correlation_cutoff'])
    model_type = environ['model_type']
    resolution = float(environ['si_resolution'])
    runtime_dir=environ['runtime_dir']

    start_lat = domain_dict[EXTENT]['llcrnrlat']
    end_lat = domain_dict[EXTENT]['urcrnrlat']
    start_lon = domain_dict[EXTENT]['llcrnrlon']
    end_lon = domain_dict[EXTENT]['urcrnrlon']

    # Make the empty analysis grid.
    grid_lat = np.arange(start_lat,end_lat+resolution,resolution)
    grid_lon = np.arange(start_lon,end_lon+resolution,resolution)
    analysis_pcp = np.empty([len(grid_lat), len(grid_lon)])

    # Download the data
    date_dt = datetime.strptime(args.date_string, '%Y%m%d')
    dated_url = date_dt.strftime(map_type_dict[MAP_TYPE]['url'])
    print('Field filename: %s' % dated_url)
    satellite_data, _, _, _ = data_dl(dated_url, MAP_TYPE, runtime_dir)

    # Interpolate the background field.
    bg_precip, bg_lat, bg_lon = apply_interp_to_grid(satellite_data[2], satellite_data[1], satellite_data[0], start_lat, end_lat, start_lon, end_lon, resolution)

    # Read in station data.
    date_string = date_dt.strftime('%Y%m')
    year = date_string[:4]
    month = date_string[4:6]
    days_in_month = monthrange(int(year), int(month))[1]
    station_values, _, station_lat, station_lon = rr_total_file_reader(STATION_DIR, date_string, file_type = '', delimiter = ',')
    station_values = np.array(station_values)
    station_lat = np.array(station_lat)
    station_lon = np.array(station_lon)
    
    # Interpolate the background field to allow estimation of field at observation coordinates.
    bg_precip *= days_in_month
    bg_interp = RegularGridInterpolator((bg_lat, bg_lon), bg_precip, bounds_error = False, fill_value = None)   

    # Read in rz and L.
    with open('%s/si.conf' % config_dir, 'r') as in_file:
        content = in_file.readlines()
        for line_idx, line in enumerate(content):
            if line.split(',')[0] == 'Rz':
                rz_array = np.array([float(i) for i in line.split(',')[1:]])
            elif line.split(',')[0] == 'L':
                L_array = np.array([float(i) for i in line[1:].split(',')[1:]])
            elif line.split(',')[0] == 'Total variance':
                total_error_variance_array = np.array([float(i) for i in line.split(',')[1:]])

    month_idx = int(date_dt.strftime('%m')) - 1
    monthly_total_error_variance = total_error_variance_array[month_idx]
    rz = rz_array[month_idx]
    L = L_array[month_idx]

    # Calculate observation error and background error variances based on the total error variance. Daley's method.
    bg_error_variance = rz * monthly_total_error_variance
    obs_error_variance = monthly_total_error_variance - bg_error_variance   

    # Remove null observations.
    null_stations =  np.where(station_values < 0)[0]
    station_values = np.delete(station_values, null_stations)
    station_lat = np.delete(station_lat, null_stations)
    station_lon = np.delete(station_lon, null_stations)

    # Perform algorithm at each gridpoint.
    for lat_i, lon_j in product(range(len(grid_lat)), range(len(grid_lon))):

        # Initialise arrays needed.
        correlation_station_to_station_matrix = np.ones([station_number_cutoff, station_number_cutoff])
        correlation_point_to_station_vector = np.ones([station_number_cutoff,1])

        # Form the point to station distance matrix.
        distance_point_to_station_vector = haversine(np.array([station_lat, station_lon]), np.array([grid_lat[lat_i], grid_lon[lon_j]]))

        # Pick closest stations according to station number cutoff
        sorted_station_index = distance_point_to_station_vector.argsort()
        closest_station_index = sorted_station_index[0:station_number_cutoff]
        closest_stations = station_values[closest_station_index]
        closest_lats = station_lat[closest_station_index]
        closest_lons = station_lon[closest_station_index]

        reduced_stations = np.copy(closest_stations)
        reduced_lats = np.copy(closest_lats)
        reduced_lons = np.copy(closest_lons)
       
        # Form the background error covariance matrix. Based on bg_error_variance and correlation of station to station matrix.
        # To form correlation of station to station matrix, use only closest stations, going through row by row for each station.
        for station_i in range(len(closest_stations)):
            distance_station_to_station_vector = haversine(np.array([closest_lats, closest_lons]), np.array([closest_lats[station_i], closest_lons[station_i]]))
            correlation_station_to_station_matrix[station_i,:] = correlation_model(distance_station_to_station_vector, L, model_type=model_type)

        if superobs_status==True:
            # Form superobs
            while np.any(np.logical_and(correlation_station_to_station_matrix >= correlation_cutoff, correlation_station_to_station_matrix < 1)):
                for station_i in range(len(reduced_stations)):
                    supobs_index = np.where(np.logical_and(correlation_station_to_station_matrix[station_i,:] >= correlation_cutoff, correlation_station_to_station_matrix[station_i,:] < 1))
                    if len(supobs_index[0]) != 0:
                        supobs_filter = np.logical_and(correlation_station_to_station_matrix[station_i,:] >= correlation_cutoff, correlation_station_to_station_matrix[station_i,:] <= 1)

                        supob_lat = np.mean(reduced_lats[supobs_filter])
                        supob_lon = np.mean(reduced_lons[supobs_filter])
                        supob_station = np.mean(reduced_stations[supobs_filter])

                        reduced_lats = np.delete(reduced_lats, supobs_filter)
                        reduced_lons = np.delete(reduced_lons, supobs_filter)
                        reduced_stations = np.delete(reduced_stations, supobs_filter)

                        reduced_lats = np.insert(reduced_lats, 0, supob_lat)
                        reduced_lons = np.insert(reduced_lons, 0, supob_lon)
                        reduced_stations = np.insert(reduced_stations, 0, supob_station)

                #         Form the station to station distance matrix.
                        correlation_station_to_station_matrix = np.ones([len(reduced_stations), len(reduced_stations)])
                        for station_i in range(len(reduced_stations)):
                            distance_station_to_station_vector = haversine(np.array([reduced_lats, reduced_lons]), np.array([reduced_lats[station_i], reduced_lons[station_i]]))
                            correlation_station_to_station_matrix[station_i,:] = correlation_model(distance_station_to_station_vector, L, model_type=model_type)

                        break

        background_error_covarince_matrix = bg_error_variance * correlation_station_to_station_matrix

        # Form the point to station covariance vector                        
        point_to_station_vector = haversine(np.array([reduced_lats, reduced_lons]), np.array([grid_lat[lat_i], grid_lon[lon_j]]))
        correlation_point_to_station_vector = correlation_model(point_to_station_vector, L, model_type=model_type)
        covariance_point_to_station_vector = bg_error_variance * correlation_point_to_station_vector

        # Form observation error covariance matrix. Based on obs_error variance and identity matrix.
        obs_error_covariance_matrix = obs_error_variance * np.identity(len(reduced_stations))

        # Calculate the innovation (difference between stations and background field at station locations)
        bg_at_station_vector = bg_interp((np.array([reduced_lats, reduced_lons]).T))
        innovation = (reduced_stations - bg_at_station_vector ).reshape(len(reduced_stations),1)

       # Weight matrix based on distance from station correlations.
        weight = np.matmul(covariance_point_to_station_vector, np.linalg.inv(background_error_covarince_matrix + obs_error_covariance_matrix))

        # Analysis point is equal to background field plus weights * innovation.
        analysis_pcp[lat_i, lon_j] = bg_precip[lat_i, lon_j] + np.matmul(weight, innovation)

    # Zero very small values (and negative values) from final analysis.
    analysis_pcp[analysis_pcp < 0.1] = 0
       
    # Create the NetCDF file.
    grid_output_path = '%s/si' % grid_dir
    create_ncfile(grid_output_path, '%s.si.%s.month%s.%s.nc' % (MAP_TYPE, EXTENT, str(PERIOD), date_dt.strftime('%Y%m')), analysis_pcp, 'precip', 'Statistical interpolation analysis, based on totals', 'mm', grid_lat, grid_lon, date_dt)

# To allow running on Unix.

if __name__ == "__main__":
   main()