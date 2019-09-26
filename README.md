# SSHCZO-EC

Eddy covariance flux code for Shale Hills Critical Zone Observatory

Contact: Yuning Shi [(Send Email)](mailto:yshi@psu.edu)

## INSTRUCTION

1. Download the source files and the Makefile to your work directory. In your work directory, type `make all` to compile all executables. Three executables will be generated: `read_data`, `split_data`, and `calculate_flux`.

2. Read raw data of a single month into one file:
   When the `calculate_flux` program processes the flux data of a certain month, it searches the `Data/YYYY-MM/` directory for the `YYYY-MM.dat` file.
   To produce the `YYYY-MM.dat` file, we must run the `read_data` executable first.
   The syntax is:

   ```shell
   ./read_data YYYY-MM filename1 [filename2 ...]
   ```

   **EXAMPLE:**

   If we want to process the data from May 2014, and the raw data come in from three files, `./incoming/EC_ts_data1.dat`, `./incoming/EC_ts_data2.dat`, and `./incoming/EC_ts_data3.dat`, we run:

   ```shell
   ./read_data 2014-05 ./incoming/EC_ts_data1.dat ./incoming/EC_ts_data2.dat ./incoming/EC_ts_data3.dat
   ```

3. Split one month data into single days:
   This step is optional. This should be executed after the previous step, if needed.
   The syntax is:

   ```shell
   ./split_data YYYY-MM
   ```

4. Process flux data of a single month:
   This is the core of the EC process.
   When the code begins processing, it searches the `Data` directory for the `YYYY-MM.dat` file.
   Make sure `read_data` has been executed. The syntax for flux processing program is:

   ```shell
   ./calculate_flux YYYY-MM [-WPL pressure_filename]
   ```

   `-WPL` is an optional parameter.
   When the `-WPL` parameter is used, the Webb-Pearman-Leuning correction will be applied.
   The WPL correction, however, requires the surface pressure data.
   Therefore the path to the surface pressure records (i.e., the ten-minute tower top file) must be specified.

   **EXAMPLE:**

   If we want to process the data from May 2014 with WPL correction, and the ten-minute tower top file is located at `incoming/EC_ten_min_data.dat`, we run

   ```shell
   ./calculate_flux 2014-05 -WPL ./incoming/EC_ten_min_data.dat
   ```

   **CAUTION:**

   If you have problems reading pressure data, please check if the code and the data format match.

   Current format of Pressure file (as of June 01, 2014) is below:

   ```
   "TIMESTAMP","RECORD","pressure_irga_mean","T_hmp_mean","T_hmp_current","RH_hmp_current","h2o_irga_mean","h2o_hmp_mean"
   ```

   Code in `calculate_flux.f90` reading these data:

   ```Fortran
   READ(600,*,IOSTAT = error) buffer, RECORD, P, T1, T2, RH, H2O1, H2O2
   ```

## CSAT3 AND IRGA DIAGNOSTIC INFORMATION

The diagnostic value in the 10-Hz data is a 12-bit integer.

bit 11&nbsp;&nbsp;&nbsp;bit 10&nbsp;&nbsp;&nbsp;bit 9&nbsp;&nbsp;&nbsp;bit 8|bit 7&nbsp;&nbsp;&nbsp;bit 6&nbsp;&nbsp;&nbsp;bit 5&nbsp;&nbsp;&nbsp;bit 4|bit 3&nbsp;&nbsp;&nbsp;bit 2&nbsp;&nbsp;&nbsp;bit 1&nbsp;&nbsp;&nbsp;bit 0
:--------------------------:|:------------------------:|:------------------------:
CSAT3 flags|IRGA flags|AGC/6.25


### CSAT3 flags

* 1000: Difference in the speed of sound between the three non-orthogonal axes is greater than 2.360&nbsp;m&nbsp;s<sup>-1</sup>
* 0100: Poor signal lock
* 0010: Sonic signal amplitude too high
* 0001: Sonic signal amplitude too low
* 1001: Lost trigger special case
* 1010: No data special case
* 1011: Wrong CSAT3 embedded code special case
* 1100: SDM error special case
* 1101: NaN special case

### IRGA flags

* 1000: chopper
* 0100: detector
* 0010: pll
* 0001: sync

### Automatic Gain Control (AGC) value

The automatic gain control (AGC) value indicates how dirty the sensor head windows are.
Typical AGC values are 50-60%.
As dirt accumulates on the sensor head windows the value of AGC will increase.
Droplets on the window can also increase AGC value.
The AGC value should be monitored and the LI-7500 optical windows should be cleaned when necessary (when the AGC value approaches 100%).
