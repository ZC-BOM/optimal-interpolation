######################################################################
## SCRIPT:       core_functions.py
## LANGUAGE:     python
## PURPOSE:      Core functions and dictionary
######################################################################

domain_dict = {
                'aus' : {'llcrnrlat' : -43.9, 'llcrnrlon' : 111.5, 'urcrnrlat' : -8, 'urcrnrlon' : 155,},
}

map_type_dict = {
            'gsmap': {'url': 'ftp://swcem:SEMaP+2004@hokusai.eorc.jaxa.jp/EAWP/GSMaP_GNRT/DATA/%Y/%Y%m/SEMDP_GSMaP_GNRT6_0.10deg-MON_%Y%m.nc', 'nc_name': 'gsmap', 'lat_name': 'lat', 'lon_name': 'lon',
            'map_path': 'ftp://swcem:SEMaP+2004@hokusai.eorc.jaxa.jp/EAWP/GSMaP_GNRT/DATA/%Y/%Y%m/SEMDP_GSMaP_GNRT6_0.10deg-', 'ftp': 'ftp://swcem:SEMaP+2004@hokusai.eorc.jaxa.jp/', 'res': 0.1, 'climo_nc_name' : 'precip'},
            }

def check_directory(directory, create=False):

    '''
    Function to check if output folder exists
    '''
    import os
    import logging

    log = logging.getLogger()

    if not os.path.exists(directory):
        log.info("Directory %s was not found" % directory)

        if create:
            os.makedirs(directory)
            log.info("Directory %s has been created" % directory)

def url_downloader(url, local_file):

    try:  
        from urllib2 import urlopen
    except ImportError:  
        from urllib.request import urlopen

    downloaded_file = urlopen(url)
    
    with open(local_file, 'wb') as f:
        f.write(downloaded_file.read())

def data_dl(url, data_variable_name, runtime_dir, *args, **kwargs):
    
    from os import remove, environ
    import netCDF4
    from shutil import copyfile
    from datetime import datetime
    
    download_status = kwargs.get('download_status', True)
    
    anom_data = None
    pct_data = None
    
    # Download file

    destination = '%s/data_dl_temp_%s.nc' % (runtime_dir, datetime.strftime(datetime.now(),'%H%M%S%d%m%Y'))
        
    if download_status:
        url_downloader(url, destination)
    else:
        copyfile(url, destination)
    
    nc = netCDF4.Dataset(destination, 'r')
    
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    data = nc.variables[data_variable_name][0, :, :]
    time = nc.variables['time'][0]

    time_units = nc.variables['time'].getncattr('units')
    
    try:
        time_cal = nc.variables['time'].getncattr('calendar')
    except:
        time_cal = 'standard'
        
    time_dt = netCDF4.num2date(time, units = time_units, calendar = time_cal)
    
    plot_data = (lon, lat, data)
        
    nc.close()
    remove(destination)
    
    return plot_data, time_dt, anom_data, pct_data

def create_ncfile(file_path, file_name, data, data_name, data_description, data_units, lat_in, lon_in, time_dt, *args, **kwargs):
    # Makes a ncfile with data, lat, lon, time.
    
    import netCDF4
    import datetime
    import logging

    log = logging.getLogger()

    missing_value = kwargs.get('missing_value', None)

    check_directory(file_path, create = True)

    dataset = netCDF4.Dataset('%s/%s' % (file_path, file_name), 'w', format='NETCDF4_CLASSIC')
    
    # Make dimensions
    
    time = dataset.createDimension('time', 1)
    lat = dataset.createDimension('lat',len(lat_in))
    lon = dataset.createDimension('lon',len(lon_in))
    
    # Make variables. Set up based on time and station numbers as dimensions.
    
    times = dataset.createVariable('time','f4',('time',))
    latitudes = dataset.createVariable('lat','f4',('lat',))
    longitudes = dataset.createVariable('lon','f4',('lon',))
    netcdf_data = dataset.createVariable(data_name,'f4', ('time','lat','lon',))

    # Create time attributes.
    
    times.calendar = 'standard'
    times.units = 'days since 1995-01-01 00:00'
    
    # Fill up the variables.
    
    latitudes[:] = lat_in
    longitudes[:] = lon_in
    netcdf_data[0,:,:] = data
    
    # Convert the datetime to netCDF time.

    netcdf_time = netCDF4.date2num(time_dt,units=times.units,calendar=times.calendar) 
    times[0] = netcdf_time
    
    # Create further attributes.
    
    dataset.description = data_description
    dataset.history = 'Created ' + datetime.datetime.strftime(datetime.datetime.now(),'%H:%M:%S %d/%m/%Y')
    
    latitudes.units = 'degrees_north'
    latitudes.axis = 'y'
    latitudes.long_name = 'latitude of center of grid box'
    longitudes.units = 'degrees_east'
    longitudes.axis = 'x'
    longitudes.long_name = 'longitude of center of grid box'
    netcdf_data.units = data_units

    if missing_value:
        netcdf_data.missing_value = missing_value
    
    dataset.close()

    log.info('%s/%s created' % (file_path, file_name))

def apply_interp_to_grid(array, lat, lon, lat_min, lat_max, lon_min, lon_max, res):
    
    from scipy.interpolate import RegularGridInterpolator
    import numpy as np

    interp_lat = np.arange(lat_min,lat_max+res,res)
    interp_lon = np.arange(lon_min,lon_max+res,res)

    obs_regrid = RegularGridInterpolator((lat, lon), array, bounds_error = False, fill_value = None)
    xx,yy = np.meshgrid(interp_lat, interp_lon, indexing='ij')
    out_mesh=(xx,yy)
    interp_array = obs_regrid(out_mesh, method='linear')

    return interp_array, interp_lat, interp_lon

def rr_total_file_reader(input_path, date, elevation_status=False, file_type = 'total.', delimiter = ' '):

    from calendar import monthrange
    from datetime import datetime
    
    append_status = False
    station_names = []
    station_lat = []
    station_lon = []
    station_values = []
    station_elevation = []
    
    date_dt = datetime.strptime(date, '%Y%m')
    year = date[:4]
    month = date[4:6]
    days_in_mon = monthrange(int(year), int(month))[1]
    
    input_filename = date_dt.strftime('rr.' + file_type  + '%Y%m01%Y%m' + str(days_in_mon))
    print('Working on: %s/%s' % (input_path, input_filename))
    with open('%s/%s' % (input_path, input_filename), 'r') as in_file:
        content = in_file.readlines()
        for line_idx, line in enumerate(content):
            if file_type in ['ratio.', '']:
                start_mark = '[$]'
            else:
                start_mark = '[$]\n'
            if line.strip() == start_mark:
                if append_status == True:
                    append_status = False
                else:
                    append_status = True
                    continue
            
            if append_status == True:
                if file_type == 'ratio.':
                    station_values.append(float(line.split()[4].strip()))
                    station_names.append(line.split()[9].strip())
                    station_lat.append(float(line.split()[2].strip()))
                    station_lon.append(float(line.split()[3].strip()))
                else:
                    if float(line.split(delimiter)[4].strip()) >= 0:
                        station_values.append(float(line.split(delimiter)[4].strip()))
                        station_names.append(line.split(delimiter)[9].strip())
                        station_lat.append(float(line.split(delimiter)[2].strip()))
                        station_lon.append(float(line.split(delimiter)[3].strip()))              
                if elevation_status:
                    station_elevation.append(float(line.split()[8]))
                    
    if elevation_status:
        return station_values, station_names, station_lat, station_lon, station_elevation
    else:
        return station_values, station_names, station_lat, station_lon

def haversine(coords1, coords2):
    import numpy as np
    R = 6371
    phi1 = coords1[0] * np.pi/180
    phi2 = coords2[0] * np.pi/180
    diff_phi = (coords2[0]-coords1[0]) * np.pi/180
    diff_lamda = (coords2[1]-coords1[1]) * np.pi/180
    
    a = np.sin(diff_phi/2)*np.sin(diff_phi/2) + np.cos(phi1)*np.cos(phi2)*np.sin(diff_lamda/2)*np.sin(diff_lamda/2)
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))

    d = R*c
    
    return d

def correlation_model(r, L, model_type = 'thiebaux'):
    import numpy as np
    # r is distance, L is length scale
    if model_type == 'gaussian':
        correlation = np.exp(-(r**2)/(2*L**2))
    elif model_type == 'thiebaux':
        correlation = (1+np.divide(r,L))*np.exp(-(np.divide(r,L)))
    return correlation