SSHCZO-EC
=========

Eddy covariance flux code for Shale Hills Critical Zone Observatory

Contact: Yuning Shi [(Send Email)](mailto:yshi@psu.edu)


INSTRUCTION:
------------

1. Type `make all` to compile all executables.

1. Read raw data of a single month into one file:  
   When the calculate_flux program processes the flux data of a certain month, it will look into the "Data/YYYY-MM/" directory for the "YYYY-MM.dat" file.
   To produce the "YYYY-MM.dat" file, we must run the read_data program first.
   The syntax is:

   ~~~
   ./read_data YYYY-MM filename1 [filename2 ...]
   ~~~

   **EXAMPLE:**

   If we want to process the data of May 2014, and the raw data come in from three files, ./incoming/EC_ts_data1.dat, ./incoming/EC_ts_data2.dat, and ./incoming/EC_ts_data3.dat, we run:

   ~~~
   ./read_data 2014-05 ./incoming/EC_ts_data1.dat ./incoiming/EC_ts_data2.dat ./incoming/EC_ts_data3.dat
   ~~~

2. Split one month data into single days:
   This step is optional. This should be executed after the first step, if needed.
   The syntax is:

   ~~~
   ./split_data YYYY-MM
   ~~~

3. Process flux data of a single month:
   This is the core of the EC process.
   When the program begins processing, it will look into the "Data" directory for the "YYYY-MM.dat" file.
   So make sure the first step has been done. The syntax for flux processing program is:

   ~~~
   ./calculate_flux YYYY-MM [-WPL pressure_filename]
   ~~~

   -WPL is an optional parameter.
   When we use the -WPL parameter, the Webb-Pearman-Leuning correction will be applied.
   The WPL correction requires the surface pressure data, so the location of the surface pressure records (in the ten-minute tower top file) must be specified.  

   **EXAMPLE:**

   If we want to process the data of May 2014 with WPL correction, and the ten-minute tower top file is at "incoming/EC_ten_min_data.dat", we run
   ~~~
   ./calculate_flux 2014-05 -WPL ./incoming/EC_ten_min_data.dat
   ~~~

CSAT3 AND IRGA DIAGNOSTIC INFORMATION
-------------------------------------

The diagnostic value in the 10-Hz data is a 12-bit integer.

bit 11&nbsp;&nbsp;&nbsp;bit 10&nbsp;&nbsp;&nbsp;bit 9&nbsp;&nbsp;&nbsp;bit 8|bit 7&nbsp;&nbsp;&nbsp;bit 6&nbsp;&nbsp;&nbsp;bit 5&nbsp;&nbsp;&nbsp;bit 4|bit 3&nbsp;&nbsp;&nbsp;bit 2&nbsp;&nbsp;&nbsp;bit 1&nbsp;&nbsp;&nbsp;bit 0
:--------------------------:|:------------------------:|:------------------------:
CSAT3 flags|IRGA flags|AGC/6.25


### CSAT3 flags:

* 1000: Difference in the speed of sound between the three non-orthogonal axes is greater than 2.360 ms-1
* 0100: Poor signal lock
* 0010: Sonic signal amplitude too high
* 0001: Sonic signal amplitude too low
* 1001: Lost trigger special case
* 1010: No data special case
* 1011: Wrong CSAT3 embedded code special case
* 1100: SDM error special case
* 1101: NaN special case

### IRGA flags:

* 1000: chopper
* 0100: detector
* 0010: pll
* 0001: sync

### Automatic Gain Control (AGC) value
The automatic gain control (AGC) valude indicates how dirty the sensor head windows are.
Typical AGC values are 50-60%.
As dirt accumulates on the sensor head windows the value of AGC will increase.
Droplets on the window can also increase AGC value.
The AGC value should be monitored and the LI-7500 optical windows should be cleaned when necessary (when the AGC value approaches 100%).
