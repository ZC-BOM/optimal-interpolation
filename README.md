# optimal-interpolation
Python scripts to perform Optimal interpolation algorithm

This repository contains code to perform Optimal Interpolation, a geospatial data assimilation algorithm.

**Instructions to run:**

1. Ensure all required Python modules are installed.

2. Adjust si.conf files if desired.

si.conf contains three lines that correspond to the a priori parameters that specify the error and correlation information.
First line - Rz - Ratio of background error variance to total error variance. Greater values will result in the background field being weighed less.
Second line - L - Length scale for the decorrelation model. Greater values will result in station information being transferred across a greater area.
Third line - Total variance - Total error variance of background and insitu values.

3. Run either by:

a) Through Jupyter Notebooks

  i) Copying all the individual files (except run_si.bash) into individual cells in a Python notebook. env.sh has to be in the first cell, si_satellite_grid.py has to be in the last cell. In the env.sh cell, remove references to exporting and define paths as Python variables instead. 
  
  ii) Ensure all paths in env.sh exist.
  
  iii) As you are not going to parse in a date through a command line, remove the parser in si_satellite_grid.py. Replace 'args.date_string' with 'date_string' and set        this to a string of format 'YYYYMMDD'.
  
  iv) Run the si_satellite_grid.py cell.

b) Through Unix:

  i) In run_si.bash, set the paths to env.sh and config.conf. Ensure all paths in env.sh exist.
  
  ii) Ensure input station file exist and is in 'grid_dir'. The format used in the study is given. If your station data is in a different format, it might be better to         develop your own code for reading the data file and return the station rainfall values, station latitude values and station longitude values as the arrays               'station_values', 'station_lat' and 'station_lon'. The line 'station_values, _, station_lat, station_lon = rr_total_file_reader(STATION_DIR, date_string, file_type         = '', delimiter = ',')' in si_satellite_grid.py should be replaced in this case.
  
  ii) Ensure run_si.bash is executable.
  
  iii) Use the command 'run_si.bash -d YYYYDDMM si_satellite_grid'.

