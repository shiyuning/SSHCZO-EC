# SSHCZO-EC

Eddy covariance flux code for Shale Hills Critical Zone Observatory

Contact: Yuning Shi [(Send Email)](mailto:yshi@psu.edu)

## INSTRUCTION

Ensure you have Python version 3 installed along with the required packages, `pandas` and `numpy`.

1. Define site-specific parameters:

   Site-specific parameters should be configured in the `bin/site_parameters.py` file.
   This includes information such as the site name, CSAT3 y-axis direction, instrument sampling frequency, data column names, quality control threshold values, and more.

2. Process flux data of a single month:

   To process flux data for a specific month, use the following command:

   ```shell
   python3 ./bin/ec_flux.py --month YYYY-MM --files filename1 [filename2 ...] [--WPL pressure_filename]
   ```

   **EXAMPLE:**

   Suppose you want to process data from May 2014, and you have three raw data files: `./incoming/EC_ts_data1.dat`, `./incoming/EC_ts_data2.dat`, and `./incoming/EC_ts_data3.dat`.
   You can run the following command:

   ```shell
   python3 ./bin/ec_flux.py --month 2014-05 --files ./incoming/EC_ts_data1.dat ./incoming/EC_ts_data2.dat ./incoming/EC_ts_data3.dat --WPL ./incoming/ten_min_data.dat
   ```

   The `-WPL` parameter is optional and is used to apply the Webb-Pearman-Leuning correction.
   This correction requires surface pressure and air temperature data, so you must specify the path to the surface pressure and air temperature records when using this parameter.

   For each month, three CSV output files will be generated.
   The main flux file will be named `<site>_<averaging period>_<start timestamp>_<end timestamp>.csv`.
   For example, the 30-minute flux data file for Shale Hills (Ameriflux site US-SSH) in May 2014 will be named `US-SSH_HH_201405010000_201406010000.csv`.
   Additionally, two diagnostic files will be generated: `*_diag.csv` containing diagnostic parameters and `*_flag.csv` containing quality control flags.

3. Reprocess using new quality control threshold values:

   If you need to apply different quality control threshold values to filter out bad flux data without reprocessing everything, you can follow these steps:

   ```shell
   python3 ./bin/write_flux_csv.py --month YYYY-MM
   ```

   New `*_flag.csv` and main flux files will be generated for the specified month with the updated quality control thresholds.

### Site-specific parameters

In the `site_parameters.py` file, you can define the following site-specific parameters:

#### `SITE`
The site name used for naming output CSV files.

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

## Output variables

| Variable        | Description                                                  | Units                                            |
| --------------- | ------------------------------------------------------------ | ------------------------------------------------ |
| TIMESTAMP_START | ISO timestamp start of averaging period                      | YYYYMMDDHHMM                                     |
| TIMESTAMP_END   | ISO timestamp end of averaging period                        | YYYYMMDDHHMM                                     |
| USTAR           | Friction velocity                                            | m s<sup>-1</sup>                                 |
| WD              | Wind direction                                               | Decimal degrees                                  |
| WS              | Wind speed                                                   | m s<sup>-1</sup>                                 |
| FC              | Carbon Dioxide (CO2) turbulent flux (no storage correction)  | µmolCO<sub>2</sub> m<sup>-2</sup> s<sup>-1</sup> |
| H               | Sensible heat turbulent flux (no storage correction)         | W m<sup>-2</sup>                                 |
| LE              | Latent heat turbulent flux (no storage correction)           | W m<sup>-2</sup>                                 |
| CO2             | Carbon Dioxide (CO2) mole fraction in wet air                | mg CO<sub>2</sub> m<sup>-3</sup>                 |
| H2O             | Water (H2O) vapor in mole fraction of wet air                | g H<sub>2</sub>O m<sup>-3</sup>                  |
| PA              | Atmospheric pressure                                         | kPa                                              |
| T_SONIC         | Sonic temperature                                            | deg C                                            |
| TA              | Air temperature                                              | deg C                                            |


## QC flags

QC flags reported in the `_flag.csv` file are interpreted as binary numbers, where `1` in each bit represents a failed QC check.

- bit 0: missing data
- bit 1: spikes
- bit 2: amplitude resolution
- bit 3: dropouts
- bit 4: absolute limits
- bit 5: higher moment statistics
- bit 6: discontinuities
- bit 7: nonstationarity

## References:
Burba, G. and Anderson, D., 2010. *A brief practical guide to eddy covariance flux measurements: Principles and workflow examples for scientific and industrial applications*. Li-Cor Biosciences.

Lee, X., Massman, W. and Law, B. eds., 2004. *Handbook of micrometeorology: A guide for surface flux measurement and analysis*. Springer Science & Business Media.

Vickers, D., and L. Mahrt, 1997: Quality control and flux sampling problems for tower and aircraft data. *Journal of Atmospheric and Oceanic Technology*, **14**, 512–526, [doi:10.1175/1520-0426(1997)014<0512:QCAFSP>2.0.CO;2](https://doi.org/10.1175/1520-0426(1997)014<0512:QCAFSP>2.0.CO;2).
