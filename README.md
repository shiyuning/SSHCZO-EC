SSHCZO-EC
=========

Eddy covariance flux code for Shale Hills Critical Zone Observatory
Contact: Yuning Shi (yshi@psu.edu)

Type "make all" to compile all executables.

INSTRUCTION:
------------

1. Read raw data of a single month into one file:  
   When the calculate_flux program processes the flux data of a certain month, it will look into the "Data/YYYY-MM/" directory for the "YYYY-MM.dat" file.To produce the "YYYY-MM.dat" file, we must run the read_data program first. The syntax is:
~~~
./read_data YYYY-MM filename1 [filename2 ...]
~~~
   EXAMPLE:
   If we want to process the data of May 2014, and the raw data come in from three files, ./incoming/EC_ts_data1.dat, ./incoming/EC_ts_data2.dat, and ./incoming/EC_ts_data3.dat, we run:
~~~
./read_data 2014-05 ./incoming/EC_ts_data1.dat ./incoiming/EC_ts_data2.dat ./incoming/EC_ts_data3.dat
~~~
2. Split one month data into single days:
   This step is optional. This should be executed after the first step, if needed. The syntax is:
~~~
./split_data YYYY-MM
~~~
3. Process flux data of a single month:
   This is the core of the EC process. When the program begins processing, it will look into the "Data" directory for the "YYYY-MM.dat" file. So make sure the first step has been done. The syntax for flux processing program is:
~~~
./calculate_flux YYYY-MM [-WPL pressure_filename]
~~~
   -WPL is an optional parameter. When we use the -WPL parameter, the Webb-Pearman-Leuning correction will be applied. The WPL correction requires the surface pressure data, so the location of the surface pressure records (in the ten-minute tower top file) must be specified.  
   EXAMPLE:
   If we want to process the data of May 2014 with WPL correction, and the ten-minute tower top file is at "incoming/EC_ten_min_data.dat", we run
~~~
./calculate_flux 2014-05 -WPL ./incoming/EC_ten_min_data.dat
~~~
