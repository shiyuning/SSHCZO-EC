SSHCZO-EC
=========

Eddy covariance flux code for Shale Hills Critical Zone Observatory
Contact: Yuning Shi (yshi@psu.edu)

Type "make all" to compile all executables.

INSTRUCTION:
1. Read raw data of a single month into one file:
  When the calculate_flux program processes the flux data of a certain month, it will look into the "Data/YYYY-MM/" directory for the "YYYY-MM.dat" file.
  To produce the "YYYY-MM.dat" file, we must run the read_data program first.
  The syntax is:

  "./read_data YYYY-MM filename1 [filename2 ...]".

  EXAMPLE:

  If we want to process the data of May 2014, and the raw data come in from three files, ./incoming/EC_ts_data1.dat, ./incoming/EC_ts_data2.dat, and ./incoming/EC_ts_data3.dat, we run:

  "./read_data 2014-05 ./incoming/EC_ts_data1.dat ./incoiming/EC_ts_data2.dat ./incoming/EC_ts_data3.dat"
 
To split one month data into single days, run "./split_data YYYY-MM".
To process flux data of a single month, run "./calculate_flux YYYY-MM [-WPL pressure_filename]".
