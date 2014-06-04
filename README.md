SSHCZO-EC
=========

Eddy covariance flux code for Shale Hills Critical Zone Observatory
Contact: Yuning Shi (yshi@psu.edu)

Type "make all" to compile all executables.

To write raw data of a single month into one file, run "./read_data YYYY-MM filename1 [filename2 ...]".
To split one month data into single days, run "./split_data YYYY-MM".
To process flux data of a single month, run "./calculate_flux YYYY-MM [-WPL pressure_filename]".
