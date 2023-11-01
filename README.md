# SSHCZO-EC

Eddy covariance flux code for Shale Hills Critical Zone Observatory

Contact: Yuning Shi [(Send Email)](mailto:yshi@psu.edu)

## INSTRUCTION

The code should be run using Python version 3, and requires the `pandas` and `numpy` packages.

1. Define site-specific parameters:

   Site-specific parameters should be defined in the `bin/site_parameters.py` file.
   This includes site name, CSAT3 y-axis direction, instrument sampling frequency, data column names, quality control threshold values, etc.

2. Process flux data of a single month:

   The syntax for flux processing program is:
   ```shell
   python3 ./bin/ec_flux.py --month YYYY-MM --files filename1 [filename2 ...] [--WPL pressure_filename]
   ```

   **EXAMPLE:**

   If we want to process the data from May 2014, and the raw data come in from three files, `./incoming/EC_ts_data1.dat`, `./incoming/EC_ts_data2.dat`, and `./incoming/EC_ts_data3.dat`, we run:
   ```shell
   python3 ./bin/ec_flux.py --month 2014-05 --files ./incoming/EC_ts_data1.dat ./incoming/EC_ts_data2.dat ./incoming/EC_ts_data3.dat
   ```
   `-WPL` is an optional parameter.
   When the `-WPL` parameter is used, the Webb-Pearman-Leuning correction will be applied.
   The WPL correction, however, requires the surface pressure data.
   Therefore the path to the surface pressure records (i.e., the ten-minute tower top file) must be specified.

   For each month, three csv output files are generated.
   The main flux file is named `<site>_<averaging period>_<start timestamp>_<end timestamp>.csv`.
   For example, Shale Hills (Ameriflux site US-SSH) May 2014 30-min flux data file will be named `US-SSH_HH_201405010000_201406010000.csv`
   Two diagnostic files will be generated along with the flux file.
   A `*_diag.csv`, which contains all the diagnostic parameters, and a `*_flag.csv`, which contains the quality control flags.

3. Reprocess using new quality control threshold values:

   If different quality control threshold values need to be applied to change the filter of bad flux data, we don't have to reprocess everything.
   Because the `*_diag.csv` files contain all the diagnostic parameters needed to flag the flux data, we only need to read the `*_diag.csv` file, apply new quality control threshold values, and generate new flux data files.
   To do this, simply run:
   ```shell
   python3 ./bin/write_flux_csv.py --month YYYY-MM
   ```

   New `*_flag.csv` and main flux files will be generated for the specified month.

### Site-specific parameters

The following site-specific parameters or functions can be defined in `site_parameters.py`.

#### `SITE`
This is the site name, which is used to name output csv files.

#### `CSAT3_AZIMUTH`
Direction between CSAT3 (or other anemometer) y-axis and true north direction. This parameter is being used for wind direction calculation.

#### `AVERAGING_PERIOD_MINUTES`
Eddy covariance flux averaging time period in minutes.

#### `FREQUENCY_HZ`
Eddy covariance data sampling frequency in Hz.

#### Eddy covariance data file column names
`TIME` is the name of the time column; `U` is the wind speed in x direction, `V` is the wind speed in y direction, `W` is the wind speed in z direction; `T_SONIC` is the sonic temperature; `H2O` is water vapor concentration; and `CO2` is the carbon dioxide concentration.

#### Other parameters for eddy covariance data files
To help read the eddy covariance data files, the `SKIP_ROW` and `COMMENT` parameters can be used. Both parameters are correspondent to the same parameters for the Python `pandas` [`read_csv`](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html) function.

#### Unit conversions
To correctly calculate fluxes, we require wind speed in m s<sup>-1</sup>, sonic temperature in Celsius, water vapor in mg m<sup>-3</sup>, and carbon dioxide in g m<sup>-3</sup>. If the input is different from those units, the `wind_speed_m_per_s`, `tsonic_celcius`, `h2o_mg_per_m3`, and `co2_g_per_m3` functions can be used for unit conversion.
For example, if the sonic temperature in the eddy covariance data files is in Kelvin, then we can define

```Python
def tsonic_celsiu(t): return t - 273.15
```

Bad data filters: Callable `ANEMOMETER_FILTER` and `IRGA_FILTER` can be defined to flag bad data from anemomters and gas analyzers.
For example, at Shale Hills CZO, The diagnostic value in the 10-Hz data is a 12-bit integer.

bit 11&nbsp;&nbsp;&nbsp;bit 10&nbsp;&nbsp;&nbsp;bit 9&nbsp;&nbsp;&nbsp;bit 8|bit 7&nbsp;&nbsp;&nbsp;bit 6&nbsp;&nbsp;&nbsp;bit 5&nbsp;&nbsp;&nbsp;bit 4|bit 3&nbsp;&nbsp;&nbsp;bit 2&nbsp;&nbsp;&nbsp;bit 1&nbsp;&nbsp;&nbsp;bit 0
:--------------------------:|:------------------------:|:------------------------:
CSAT3 flags|IRGA flags|AGC/6.25

**CSAT3 flags**

* 1000: Difference in the speed of sound between the three non-orthogonal axes is greater than 2.360&nbsp;m&nbsp;s<sup>-1</sup>
* 0100: Poor signal lock
* 0010: Sonic signal amplitude too high
* 0001: Sonic signal amplitude too low
* 1001: Lost trigger special case
* 1010: No data special case
* 1011: Wrong CSAT3 embedded code special case
* 1100: SDM error special case
* 1101: NaN special case

**IRGA flags**

* 1000: chopper
* 0100: detector
* 0010: pll
* 0001: sync

**Automatic Gain Control (AGC) value**

The automatic gain control (AGC) value indicates how dirty the sensor head windows are.
Typical AGC values are 50-60%.
As dirt accumulates on the sensor head windows the value of AGC will increase.
Droplets on the window can also increase AGC value.
The AGC value should be monitored and the LI-7500 optical windows should be cleaned when necessary (when the AGC value approaches 100%).

Therefore, to filter out bad data from CSAT-3, we define
```Python
ANEMOMETER_FILTER = lambda x: (x['diag'] & 3840) == 0   # 3840 = 1111 0000 0000
IRGA_FILTER = lambda x: (x['diag'] & 240) == 0          # 240 = 0000 1111 0000
```

#### Pressure data file column names
`PRESSURE_TIME` is the name of the time column; `PRESSURE` is the air pressure, `T_AIR` is the air temperature.

#### Other parameters for pressure data files
Similar to eddy covariance data files, `PRESSURE_SKIP_ROWS` and `PRESSURE_COMMENT` parameters can be used.

#### Unit conversions
To correctly calculate fluxes, we require air pressure in Pascal, and air temperature in Celsius.
If the input is different from those units, the `pressure_pa` and `tair_celsius` functions can be used for conversion.

#### `QC_THRESHOLDS`
These are the quality control diagnostic thresholds as defined in [Vickers and Mahrt (1997)](https://doi.org/10.1175/1520-0426(1997)014<0512:QCAFSP>2.0.CO;2).
Please refer to the original paper for the definition of those thresholds.

## References:
Burba, G. and Anderson, D., 2010. *A brief practical guide to eddy covariance flux measurements: Principles and workflow examples for scientific and industrial applications*. Li-Cor Biosciences.

Lee, X., Massman, W. and Law, B. eds., 2004. *Handbook of micrometeorology: A guide for surface flux measurement and analysis*. Springer Science & Business Media.

Vickers, D., and L. Mahrt, 1997: Quality control and flux sampling problems for tower and aircraft data. *Journal of Atmospheric and Oceanic Technology*, **14**, 512–526, [doi:10.1175/1520-0426(1997)014<0512:QCAFSP>2.0.CO;2](https://doi.org/10.1175/1520-0426(1997)014<0512:QCAFSP>2.0.CO;2).
