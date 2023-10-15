#!/usr/bin/env python3

def initialize_csv_files(flux_file, diag_file):
    """Initialize the output CSV files
    """
    with open(diag_file, 'w') as f:
        f.write('YEAR,MONTH,MDAY,DOY,HRMIN,DTIME,UST,TA,WD,WS,VWS,FC,H,LE,CO2,H2O,U_FLAG,W_FLAG,T_FLAG,H2O_FLAG,CO2_FLAG\n')
        f.write('-,-,-,-,-,-,m/s,degC,deg,m/s,m/s,umol/m2/s,w/m2,w/m2,mg/m3,g/m3,-,-,-,-,-\n')

    with open(flux_file, 'w') as f:
        f.write('YEAR,MONTH,MDAY,DOY,HRMIN,DTIME,UST,TA,WD,WS,VWS,FC,H,LE,CO2,H2O\n')
        f.write('-,-,-,-,-,-,m/s,degC,deg,m/s,m/s,umol/m2/s,w/m2,w/m2,mg/m3,g/m3\n')