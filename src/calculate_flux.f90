!-------------------------------------------------------------------------------
! Program for processing Shale Hills flux tower eddy covariance data
! Yuning Shi, Penn State Meteorology
! Email: yshi@psu.edu
!-------------------------------------------------------------------------------


PROGRAM SHEC

  USE ECtype
  USE fluxtype
  USE time

  IMPLICIT NONE
 
  CHARACTER(100) :: buffer, inputname, output_flag, output_no_flag
  CHARACTER(200) :: pressure_filename
  INTEGER :: error

  INTEGER :: datayear, datamonth
  INTEGER :: ind, diag
  INTEGER :: time_temp
  REAL :: second, Ux, Uy, Uz, Ts, co2, h2o, h2o_hmp, T_hmp
  INTEGER :: recordtime, time_pointer
  INTEGER :: i,jout
  INTEGER :: infid = 100, outfid_flag = 200, outfid_no_flag = 300
  REAL :: unit_k(3), b0

  LOGICAL :: WPL

  TYPE(ECdata) :: record
  TYPE(fluxdata) :: flux_record
  TYPE(tm) :: timeinfo

  CALL get_command_argument(1, buffer)      ! Read year and month from command line

  READ(buffer(1:4),*) datayear
  READ(buffer(6:7),*) datamonth

  WRITE(*,"('\n  Processing ',I4.4,'-',I2.2,' data')") datayear, datamonth

  CALL get_command_argument(2, buffer)

  IF (TRIM(buffer) == "-WPL") THEN
    WPL = 1
    WRITE(*,"('\n  WPL correction will be applied.')")
  ELSE
    WPL = 0
    WRITE(*,"('\n  WPL correction will not be applied.')")
  END IF

  IF (WPL == 1) CALL get_command_argument(3, pressure_filename)

  IF (LEN_TRIM(pressure_filename) == 0) THEN
    WRITE(*,"('\n  ERROR! Must specify the file that contains pressure &
measurements when WPL correction is applied.')")
    STOP
  END IF

  !-------------------------------------------------------------------------------
  ! Define the start point (YYYY-MM-01 00:30:00)
  !-------------------------------------------------------------------------------

  timeinfo = tm(datayear, datamonth, 1, 0, 30, 0)

  CALL timegm(timeinfo, time_pointer)

  !-------------------------------------------------------------------------------
  ! Open eddy-covariance data
  !-------------------------------------------------------------------------------

  WRITE(inputname,"('Data/',I4.4,'-',I2.2,'/',I4.4,'-',I2.2,'.dat')") datayear, datamonth, datayear, datamonth
!  WRITE(inputname,"(I4.4,'-',I2.2,'/','2009-12-31.dat')") datayear, datamonth
  OPEN(UNIT=infid, FILE=inputname, ACTION='read', IOSTAT=error)
  IF (error/=0) GOTO 1000

  WRITE(*,*)"\n  File opened ...\n"

  !-------------------------------------------------------------------------------
  ! Open output files:
  ! *-flux-flag.dat contain diagnostic flags
  ! *-flux.dat is to be posted online
  !-------------------------------------------------------------------------------

  WRITE(output_flag,"('Data/',I4.4,'-',I2.2,'/',I4.4,'-',I2.2,'-flux-flag.dat')") datayear, datamonth, datayear, datamonth
  OPEN(UNIT=outfid_flag, FILE=output_flag, ACTION='WRITE', IOSTAT=error)

  WRITE(output_no_flag,"('Data/',I4.4,'-',I2.2,'/',I4.4,'-',I2.2,'-flux.dat')") datayear, datamonth, datayear, datamonth
  OPEN(UNIT=outfid_no_flag, FILE=output_no_flag, ACTION='WRITE', IOSTAT=error)

  !-------------------------------------------------------------------------------
  ! Write file headers
  !-------------------------------------------------------------------------------
  
  WRITE(outfid_flag,"(A)") "Shale Hills CZO flux tower data"
  WRITE(outfid_flag,"(A)") "Contact: Yuning Shi (yshi@psu.edu)"
  WRITE(outfid_flag,"(A)") "YEAR,MONTH,MDAY,DOY,HRMIN,DTIME,UST,TA,WD,WS,VWS,FC,H,LE,CO2,H2O,U_FLAG,W_FLAG,T_FLAG,H2O_FLAG,CO2_FLAG"
  WRITE(outfid_flag,"(A)") "-,-,-,-,-,-,m/s,degC,deg,m/s,m/s,umol/m2/s,w/m2,w/m2,mg/m3,g/m3,-,-,-,-,-"

  WRITE(outfid_no_flag,"(A)") "Shale Hills CZO flux tower data"
  WRITE(outfid_no_flag,"(A)") "Contact: Yuning Shi (yshi@psu.edu)"
  WRITE(outfid_no_flag,"(A)") "YEAR,MONTH,MDAY,DOY,HRMIN,DTIME,UST,TA,WD,WS,VWS,FC,H,LE,CO2,H2O"
  WRITE(outfid_no_flag,"(A)") "-,-,-,-,-,-,m/s,degC,deg,m/s,m/s,umol/m2/s,w/m2,w/m2,mg/m3,g/m3"


  !-------------------------------------------------------------------------------
  ! Determine unit vector k of planar fit coordinate
  ! (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 4)
  ! unit_k is unit vector parallel to new coordinate z axis
  ! bo is instrument offset in w1
  !-------------------------------------------------------------------------------

  WRITE(*,*)"\n  Calculate unit vector k of planar fit coordinate ...\n"

  unit_k(1) = 0
  unit_k(2) = 0
  unit_k(3) = 1
  b0 = 0

  ! Comment the following statement out if do not want to do coordinate correction

  CALL unit_vector_k(infid, unit_k, b0)

  WRITE(*,"('  k_vector = [',F,',',F,',',F,']')") unit_k(1),unit_k(2),unit_k(3)
  WRITE(*,"('  b0 = ',F)") b0

  !-------------------------------------------------------------------------------
  ! Read input data
  !-------------------------------------------------------------------------------

  WRITE(*,*)"\n  Start reading in file ... \n"

  REWIND(infid,iostat = error)

  READ(infid,*,IOSTAT = error) buffer      ! Skip head-lines
  READ(infid,*,IOSTAT = error) buffer
  READ(infid,*,IOSTAT = error) buffer
  READ(infid,*,IOSTAT = error) buffer

  record%record_no = 0          ! Initiate ECdata

  DO WHILE (.true.)
    READ(infid,*,IOSTAT = error) buffer, ind, Ux, Uy, Uz, Ts, co2, h2o, diag, h2o_hmp, T_hmp
    IF (error < 0 ) THEN        ! If reaches end of file, calculate flux using existing record
      CALL flux(time_pointer, record, unit_k, b0, flux_record, outfid_flag, outfid_no_flag, WPL, pressure_filename)
      GOTO 1000
    ELSE IF (error > 0) THEN        ! If error occurs reading in a record, skip
      CYCLE
    END IF

    !-------------------------------------------------------------------------------
    ! Source file Diag configuration
    ! ------------------------------------------------------------------------------
    ! 11 | 10 |  9 |  8 |  7 |  6 |  5 |  4 |  3 |  2 |  1 |  0 |
    !    CSAT3 flags    |    IRGA flags     |    AGC/6.25       |
    ! CSAT3:
    ! 9: lost trigger special case
    ! 10: no data special case
    ! 11: wrong CSAT3 embedded code special case
    ! 12: SDM error special case
    ! 13: NaN special case
    !
    ! IRGA:
    ! 1000: chopper
    ! 0100: detector
    ! 0010: pll
    ! 0001: sync
    !-------------------------------------------------------------------------------
    ! Determine time for current record
    !-------------------------------------------------------------------------------

    READ(buffer(:4),*) timeinfo%tm_year
    READ(buffer(6:7),*) timeinfo%tm_mon
    READ(buffer(9:10),*) timeinfo%tm_mday
    READ(buffer(12:13),*) timeinfo%tm_hour
    READ(buffer(15:16),*) timeinfo%tm_min
    READ(buffer(18:),*) second

    timeinfo%tm_sec = 0

    CALL timegm(timeinfo, recordtime)      ! Convert timestamp to integer value

    !-------------------------------------------------------------------------------
    ! Copy data and calculate fluxes
    !-------------------------------------------------------------------------------

    DO WHILE (recordtime+CEILING(second)>time_pointer)  ! If 30-min records are collected, calculate fluxes
      CALL flux(time_pointer, record, unit_k, b0, flux_record, outfid_flag, outfid_no_flag, WPL, pressure_filename)
      record%record_no = 0        ! Re-initiate ECdata
      time_pointer = time_pointer + 30*60     ! Time pointer points to the next 30-min
    END DO

    IF (recordtime+CEILING(second)>time_pointer - 30*60 .AND. diag<16) THEN  ! Copy data to ECdata
      record%record_no = record%record_no +1
      record%time(record%record_no) = REAL(recordtime - time_pointer + 30*60) + second
      record%Ux(record%record_no) = Ux
      record%Uy(record%record_no) = Uy
      record%Uz(record%record_no) = Uz
      record%Ts(record%record_no) = Ts
      record%co2(record%record_no) = co2
      record%h2o(record%record_no) = h2o
      record%diag(record%record_no) = diag
      record%h2o_hmp(record%record_no) = h2o_hmp
      record%T_hmp(record%record_no) = T_hmp
    END IF

  END DO

1000 CONTINUE

  CLOSE(infid)
  CLOSE(outfid_flag)
  CLOSE(outfid_no_flag)

  WRITE(*,*)"\n  done.\n"


END PROGRAM SHEC

SUBROUTINE unit_vector_k(fid, unit_k, b0)

  !-------------------------------------------------------------------------------
  ! Determine unit vector k of planar fit coordinate
  ! (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 4)
  !-------------------------------------------------------------------------------

  USE, INTRINSIC :: IEEE_ARITHMETIC
  CHARACTER(100) :: buffer
  INTEGER :: diag, ind
  REAL :: Ux, Uy, Uz, Ts, co2, h2o, h2o_hmp, T_hmp
  REAL :: su, sv, sw, suv, suw, svw, su2, sv2
  REAL :: flen
  INTEGER :: i, fid
  INTEGER :: errorflag, error
  REAL :: H1(3,3), invH1(3,3),g(3),x(3)
  REAL :: b0, b1, b2
  REAL :: unit_k(3)

  !-------------------------------------------------------------------------------
  ! Initiate variables
  !-------------------------------------------------------------------------------

  su = 0
  sv = 0
  sw = 0
  suv = 0
  suw = 0
  svw = 0
  su2 = 0
  sv2 = 0
  flen = 0

  !-------------------------------------------------------------------------------
  ! Skip headlines
  !-------------------------------------------------------------------------------

  READ(fid,*,IOSTAT = error) buffer
  READ(fid,*,IOSTAT = error) buffer
  READ(fid,*,IOSTAT = error) buffer
  READ(fid,*,IOSTAT = error) buffer

  !-------------------------------------------------------------------------------
  ! Read data
  !-------------------------------------------------------------------------------

  DO WHILE (.TRUE.)
    READ(fid,*,iostat = error) buffer, ind, Ux, Uy, Uz, Ts, co2, h2o, diag, h2o_hmp, T_hmp
    IF (error < 0 ) EXIT
    IF (error > 0 .OR. diag>255) CYCLE
    su = su + Ux
    sv = sv + Uy
    sw = sw + Uz
    suv = suv + Ux*Uy
    suw = suw + Ux*Uz
    svw = svw + Uy*Uz
    su2 = su2 + Ux*Ux
    sv2 = sv2 + Uy*Uy
    flen = flen +1
!    IF(IEEE_IS_NAN(su)) THEN
!      WRITE(*,*) buffer, ind, Ux, Uy, Uz, error
!      WRITE(*,*) IBITS(diag,4,4)
!      EXIT
!    END IF
  END DO


  !-------------------------------------------------------------------------------
  ! Calculate inverse matrix
  !-------------------------------------------------------------------------------

  H1(1,:) = (/flen,su,sv/)
  H1(2,:) = (/su, su2, suv/)
  H1(3,:) = (/sv, suv, sv2/)

  g = (/sw, suw, svw/)

  CALL FindInv(H1, invH1, 3, ErrorFlag)

  !-------------------------------------------------------------------------------
  ! Calculate unit k vector
  !-------------------------------------------------------------------------------

  DO i = 1,3
    x(i) = SUM(invH1(i,:)*g(:))
  END DO

  b0 = x(1)
  b1 = x(2)
  b2 = x(3)

  unit_k(3) = 1/(1+b1*b1+b2*b2)
  unit_k(1) = -b1*unit_k(3)
  unit_k(2) = -b2*unit_k(3)

END SUBROUTINE unit_vector_k

SUBROUTINE flux(rawtime, record, unit_k, b0, flux_record, outfid_flag, outfid_no_flag, WPL, pressure_file)

!-------------------------------------------------------------------------------
! Calculate surface fluxes
!-------------------------------------------------------------------------------

  USE ECtype
  USE fluxtype
  USE QC
  USE time
  USE Pressure

  IMPLICIT NONE

  TYPE(fluxdata) :: flux_record
  TYPE(ECdata) :: record
  TYPE(ECdata) :: perturb
  TYPE(tm) :: timeinfo
  TYPE(Pressuredata) :: Precord

  INTEGER :: rawtime

  LOGICAL :: WPL

  REAL, PARAMETER :: Rd = 287.05
  REAL, PARAMETER :: Lv = 2503000.0
  REAL, PARAMETER :: c_air = 1004.0
  REAL, PARAMETER :: Pi = 3.1415927
  REAL, PARAMETER :: CSAT3_AZIMUTH = 199.0  ! Direction between CSAT3 y axis and true north direction
  INTEGER :: outfid_flag, outfid_no_flag
  REAL :: u_bar, v_bar, w_bar, h2o_bar, T_bar, eta, co2_bar
  REAL :: wT_bar, wh2o_bar, uT_bar, uh2o_bar, vT_bar, vh2o_bar, uco2_bar, vco2_bar, wco2_bar
  REAL :: uw_bar, vw_bar
  REAL :: uw(record%record_no), vw(record%record_no)
  REAL :: U_trend(record%record_no), V_trend(record%record_no), W_trend(record%record_no), CO2_trend(record%record_no), H2O_trend(record%record_no), T_trend(record%record_no)
  REAL :: U1(3)
  REAL :: b0
  REAL :: CE, SE
  REAL :: P, rho, rho_d, rho_v, rho_c, vp, q, E0, F0, Ta, Fc, E, H
  INTEGER :: i,j
  INTEGER :: INSTRUMENT_FLAG, NONSTATIONARITY_FLAG
  INTEGER :: SPIKES_FLAG(5), AMPLITUDE_RESOLUTION_FLAG(5), DROPOUTS_FLAG(5), ABSOLUTE_LIMITS_FLAG(5), HIGHER_MOMENT_FLAG(5), DISCONTINUITIES_FLAG(5)
  INTEGER :: ErrorFlag
  REAL :: unit_k(3), unit_j(3), unit_i(3)
  CHARACTER(200) :: pressure_file

  !-------------------------------------------------------------------------------
  ! Read pressure file for WPL correction
  !-------------------------------------------------------------------------------

  IF (WPL == 1) CALL read_pressure_file(Precord, pressure_file)

  !-------------------------------------------------------------------------------
  ! Determine timestamp for surface fluxes
  !-------------------------------------------------------------------------------

  CALL gmtime(rawtime, timeinfo)
  flux_record%YEAR = timeinfo%tm_year
  flux_record%MONTH = timeinfo%tm_mon
  flux_record%MDAY = timeinfo%tm_mday
  CALL doy(timeinfo%tm_year, timeinfo%tm_mon, timeinfo%tm_mday, flux_record%DOY)
  flux_record%HRMIN = 100*timeinfo%tm_hour + timeinfo%tm_min
  flux_record%DTIME = (REAL(timeinfo%tm_hour) + REAL(timeinfo%tm_min)/60)/24

  perturb%record_no = record%record_no
  perturb%time = record%time

  !-------------------------------------------------------------------------------
  ! Initialize diagnostic flags
  !-------------------------------------------------------------------------------

  SPIKES_FLAG = 0
  ABSOLUTE_LIMITS_FLAG = 0
  NONSTATIONARITY_FLAG = 0
  AMPLITUDE_RESOLUTION_FLAG = 0
  DROPOUTS_FLAG = 0
  HIGHER_MOMENT_FLAG = 0
  DISCONTINUITIES_FLAG = 0
  INSTRUMENT_FLAG = 0

  WRITE(*,"('')")
  WRITE(*,"(I4,'-',I2.2,'-',I2.2,' ',I4.4, I)") flux_record%YEAR, flux_record%MONTH, flux_record%MDAY, flux_record%HRMIN, record%record_no

  !-------------------------------------------------------------------------------
  ! Quality Control ... see QC.f90
  !-------------------------------------------------------------------------------

  IF(record%record_no > 0) THEN

    CALL instrument(record, INSTRUMENT_FLAG, flux_record)

    CALL spikes(record, SPIKES_FLAG, flux_record)

    WRITE(*,"('  Ux range = [',F7.2,',',F7.2,']')") MINVAL(ABS(record%Ux(1:record%record_no))), MAXVAL(ABS(record%Ux(1:record%record_no)))
    WRITE(*,"('  Uy range = [',F7.2,',',F7.2,']')") MINVAL(ABS(record%Uy(1:record%record_no))), MAXVAL(ABS(record%Uy(1:record%record_no)))
    WRITE(*,"('  Uz range = [',F7.2,',',F7.2,']')") MINVAL(ABS(record%Uz(1:record%record_no))), MAXVAL(ABS(record%Uz(1:record%record_no)))
    WRITE(*,"('  Ts range = [',F7.2,',',F7.2,']')") MINVAL(record%Ts(1:record%record_no)), MAXVAL(record%Ts(1:record%record_no))
    WRITE(*,"('  CO2 range = [',F7.2,',',F7.2,']')") MINVAL(record%co2(1:record%record_no)), MAXVAL(record%co2(1:record%record_no))
    WRITE(*,"('  H2O range = [',F7.2,',',F7.2,']')") MINVAL(record%h2o(1:record%record_no)), MAXVAL(record%h2o(1:record%record_no))

    CALL amplitude_resolution(record, AMPLITUDE_RESOLUTION_FLAG, flux_record)

    CALL dropouts(record, DROPOUTS_FLAG, flux_record)

    CALL absolute_limits(record, ABSOLUTE_LIMITS_FLAG)

    !-------------------------------------------------------------------------------
    ! Transformation to the natural wind coordinate system
    ! (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 4)
    !-------------------------------------------------------------------------------

    record%Uz(1:record%record_no) = record%Uz(1:record%record_no) - b0
    u_bar = SUM(record%Ux(1:record%record_no))/REAL(record%record_no)
    v_bar = SUM(record%Uy(1:record%record_no))/REAL(record%record_no)
    w_bar = SUM(record%Uz(1:record%record_no))/REAL(record%record_no)
    h2o_bar = SUM(record%h2o(1:record%record_no))/REAL(record%record_no)
    co2_bar = SUM(record%co2(1:record%record_no))/REAL(record%record_no)
    T_bar = SUM(record%Ts(1:record%record_no))/REAL(record%record_no)

    !-------------------------------------------------------------------------------
    ! Calculate wind direction and speed
    !-------------------------------------------------------------------------------

    CE = u_bar/SQRT(u_bar*u_bar+v_bar*v_bar)
    SE = v_bar/SQRT(u_bar*u_bar+v_bar*v_bar)

    eta = ACOS(CE)/Pi*180
    IF (SE<0) eta = 360 - eta
    eta = MODULO(360.0 - eta + CSAT3_AZIMUTH, 360.0)

    U1 = (/u_bar,v_bar,w_bar/)

    !-------------------------------------------------------------------------------
    ! Calculate unit vector i and j
    !-------------------------------------------------------------------------------

    CALL cross(unit_k,U1,unit_j)
    unit_j = unit_j/SQRT(SUM(unit_j*unit_j))
    CALL cross(unit_j,unit_k,unit_i)

    !-------------------------------------------------------------------------------
    ! Linear detrend
    !-------------------------------------------------------------------------------

    CALL detrend(record%time,record%Ux, U_trend, perturb%Ux, record%record_no)
    CALL detrend(record%time,record%Uy, V_trend, perturb%Uy, record%record_no)
    CALL detrend(record%time,record%Uz, W_trend, perturb%Uz, record%record_no)
    CALL detrend(record%time,record%co2, CO2_trend, perturb%co2, record%record_no)
    CALL detrend(record%time,record%h2o, H2O_trend, perturb%h2o, record%record_no)
    CALL detrend(record%time,record%Ts, T_trend, perturb%Ts, record%record_no)

    CALL higher_moment(perturb,HIGHER_MOMENT_FLAG, flux_record)

    CALL discontinuities(record, DISCONTINUITIES_FLAG, flux_record)

    CALL nonstationary(record%Ux, record%Uy, record%record_no, NONSTATIONARITY_FLAG, flux_record)

    !-------------------------------------------------------------------------------
    ! Calculate fluxes
    !-------------------------------------------------------------------------------

    uw = unit_i(1)*unit_k(1)*(perturb%Ux(1:record%record_no)*perturb%Ux(1:record%record_no))
    uw = uw + unit_i(2)*unit_k(2)*(perturb%Uy(1:record%record_no)*perturb%Uy(1:record%record_no))
    uw = uw + unit_i(3)*unit_k(3)*(perturb%Uz(1:record%record_no)*perturb%Uz(1:record%record_no))
    uw = uw + (unit_i(1)*unit_k(2)+unit_i(2)*unit_k(1))*(perturb%Ux(1:record%record_no)*perturb%Uy(1:record%record_no))
    uw = uw + (unit_i(1)*unit_k(3)+unit_i(3)*unit_k(1))*(perturb%Ux(1:record%record_no)*perturb%Uz(1:record%record_no))
    uw = uw + (unit_i(2)*unit_k(3)+unit_i(3)*unit_k(2))*(perturb%Uy(1:record%record_no)*perturb%Uz(1:record%record_no))

    uw_bar = SUM(uw)/REAL(record%record_no)

    vw = unit_j(1)*unit_k(1)*(perturb%Ux(1:record%record_no)*perturb%Ux(1:record%record_no))
    vw = vw + unit_j(2)*unit_k(2)*(perturb%Uy(1:record%record_no)*perturb%Uy(1:record%record_no))
    vw = vw + unit_j(3)*unit_k(3)*(perturb%Uz(1:record%record_no)*perturb%Uz(1:record%record_no))
    vw = vw + (unit_j(1)*unit_k(2)+unit_j(2)*unit_k(1))*(perturb%Ux(1:record%record_no)*perturb%Uy(1:record%record_no))
    vw = vw + (unit_j(1)*unit_k(3)+unit_j(3)*unit_k(1))*(perturb%Ux(1:record%record_no)*perturb%Uz(1:record%record_no))
    vw = vw + (unit_j(2)*unit_k(3)+unit_j(3)*unit_k(2))*(perturb%Uy(1:record%record_no)*perturb%Uz(1:record%record_no))

    vw_bar = SUM(vw)/REAL(record%record_no)

    flux_record%UST = SQRT(SQRT(uw_bar*uw_bar+vw_bar*vw_bar))

    uT_bar = SUM(perturb%Ux(1:record%record_no)*perturb%Ts(1:record%record_no))/REAL(record%record_no)
    vT_bar = SUM(perturb%Uy(1:record%record_no)*perturb%Ts(1:record%record_no))/REAL(record%record_no)
    wT_bar = SUM(perturb%Uz(1:record%record_no)*perturb%Ts(1:record%record_no))/REAL(record%record_no)
    uh2o_bar = SUM(perturb%Ux(1:record%record_no)*perturb%h2o(1:record%record_no))/REAL(record%record_no)
    vh2o_bar = SUM(perturb%Uy(1:record%record_no)*perturb%h2o(1:record%record_no))/REAL(record%record_no)
    wh2o_bar = SUM(perturb%Uz(1:record%record_no)*perturb%h2o(1:record%record_no))/REAL(record%record_no)
    uco2_bar = SUM(perturb%Ux(1:record%record_no)*perturb%co2(1:record%record_no))/REAL(record%record_no)
    vco2_bar = SUM(perturb%Uy(1:record%record_no)*perturb%co2(1:record%record_no))/REAL(record%record_no)
    wco2_bar = SUM(perturb%Uz(1:record%record_no)*perturb%co2(1:record%record_no))/REAL(record%record_no)

    F0 = (uco2_bar*unit_k(1)+vco2_bar*unit_k(2)+wco2_bar*unit_k(3))
    E0 = (uh2o_bar*unit_k(1)+vh2o_bar*unit_k(2)+wh2o_bar*unit_k(3))/1000

    !-------------------------------------------------------------------------------
    ! WPL correction (Webb et al. 1980)
    !-------------------------------------------------------------------------------
    IF (WPL == 1) THEN
      CALL read_pressure(Precord, rawtime, P)
      P = P*1000
      Ta = T_bar + 273.15
      rho_v = h2o_bar/1000
      rho_c = co2_bar

      vp = rho_v*Rd*Ta/0.622
      q = 0.622*vp/(P-0.378*vp)
      rho = P/(Rd*(1+0.608*q)*Ta)
      rho_d = rho - rho_v

      H = rho*c_air*(uT_bar*unit_k(1)+vT_bar*unit_k(2)+wT_bar*unit_k(3))

      ! Sonic correction is commented out because it is already done in LI-7500

!      H = H+rho*c_air*(-0.51*Ta*wh2o_bar/1000)/rho  ! Sonic correction

      E = (1+ 1.6077*rho_v/rho_d)*(E0+H/rho/c_air*rho_v/Ta)
      Fc = F0 + 1.6077*E/rho_d*rho_c/(1+1.6077*(rho_v/rho_d)) + H/rho/c_air*rho_c/Ta    

      WRITE(*,"('  Air pressure = ', F9.2,' Pa')") P
      WRITE(*,"('  Air density = ', F9.2, ' kg m-3')") rho
      WRITE(*,"('  Sensible heat flux = ', F9.2, ' W m-2')") H
      WRITE(*,"('  H2O flux before WPL correction = ', F9.2, ' kg m-2')") E0*Lv
      WRITE(*,"('  H2O flux after WPL correction = ', F9.2, ' kg m-2')") E*Lv
      WRITE(*,"('  CO2 flux before WPL correction = ', F9.2, ' ug m-2 s-1')") F0*1000/44
      WRITE(*,"('  CO2 flux after WPL correction = ', F9.2, ' ug m-2 s-1')") Fc*1000/44
    ELSE
      rho = 1.20

      F0 = (uco2_bar*unit_k(1)+vco2_bar*unit_k(2)+wco2_bar*unit_k(3))
      E0 = (uh2o_bar*unit_k(1)+vh2o_bar*unit_k(2)+wh2o_bar*unit_k(3))/1000

      H = rho*c_air*(uT_bar*unit_k(1)+vT_bar*unit_k(2)+wT_bar*unit_k(3))

      ! Sonic correction is commented out because it is already done in LI-7500

!      H = H+rho*c_air*(-0.51*Ta*wh2o_bar/1000)/rho  ! Sonic correction

      E = E0
      Fc = F0

      WRITE(*,"('  Sensible heat flux = ', F9.2, ' W m-2')") H
      WRITE(*,"('  H2O flux = ', F9.2, ' kg m-2')") E*Lv
      WRITE(*,"('  CO2 flux = ', F9.2, ' ug m-2 s-1')") Fc*1000/44
    END IF      

    flux_record%TA = T_bar
    flux_record%WD = eta
    flux_record%WS = SQRT(SUM(unit_i*U1)*SUM(unit_i*U1)+SUM(unit_j*U1)*SUM(unit_j*U1))
    flux_record%VWS = w_bar
    flux_record%FC = 1000/44*Fc
    flux_record%H = H
    flux_record%LE = Lv*E
    flux_record%CO2 = co2_bar
    flux_record%H2O = h2o_bar

  ELSE
    INSTRUMENT_FLAG = 999
  END IF

!  CALL SEB(rawtime, flux_record, RN_FLAG)

  !-------------------------------------------------------------------------------
  ! Flags:
  ! Bit 1:
  ! Bit 2: Spikes (2)
  ! Bit 3: Amplitude resolution (4)
  ! Bit 4: Dropouts (8)
  ! Bit 5: Absolute limits (16)
  ! Bit 6: Higher-moment statistics (32)
  ! Bit 7: Discontinuities (64)
  ! Bit 8: Nonstationarity of the horizontal wind (128)
  ! Bit 9: Surface energy balance (256)
  !-------------------------------------------------------------------------------

  flux_record%U_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(1)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(1)*(2**2)+DROPOUTS_FLAG(1)*(2**3)+ABSOLUTE_LIMITS_FLAG(1)*(2**4)+HIGHER_MOMENT_FLAG(1)*(2**5)+DISCONTINUITIES_FLAG(1)*(2**6)+ NONSTATIONARITY_FLAG*(2**7) !+ RN_FLAG*(2**8)

  flux_record%W_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(2)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(2)*(2**2)+DROPOUTS_FLAG(2)*(2**3)+ABSOLUTE_LIMITS_FLAG(2)*(2**4)+HIGHER_MOMENT_FLAG(2)*(2**5)+DISCONTINUITIES_FLAG(2)*(2**6)+ NONSTATIONARITY_FLAG*(2**7) !+ RN_FLAG*(2**8)

  flux_record%T_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(3)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(3)*(2**2)+DROPOUTS_FLAG(3)*(2**3)+ABSOLUTE_LIMITS_FLAG(3)*(2**4)+HIGHER_MOMENT_FLAG(3)*(2**5)+DISCONTINUITIES_FLAG(3)*(2**6)+ NONSTATIONARITY_FLAG*(2**7) !+ RN_FLAG*(2**8)

  flux_record%H2O_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(4)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(4)*(2**2)+DROPOUTS_FLAG(4)*(2**3)+ABSOLUTE_LIMITS_FLAG(4)*(2**4)+HIGHER_MOMENT_FLAG(4)*(2**5)+DISCONTINUITIES_FLAG(4)*(2**6)+ NONSTATIONARITY_FLAG*(2**7) !+ RN_FLAG*(2**8)

  flux_record%CO2_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(5)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(5)*(2**2)+DROPOUTS_FLAG(5)*(2**3)+ABSOLUTE_LIMITS_FLAG(5)*(2**4)+HIGHER_MOMENT_FLAG(5)*(2**5)+DISCONTINUITIES_FLAG(5)*(2**6)+ NONSTATIONARITY_FLAG*(2**7) !+ RN_FLAG*(2**8)

  CALL write_to_file_flag(outfid_flag, flux_record)
  CALL write_to_file_no_flag(outfid_no_flag, flux_record)

END SUBROUTINE flux


SUBROUTINE cross(A, B, C)

  IMPLICIT NONE

  REAL :: A(3),B(3),C(3)

  C(1) = A(2)*B(3)-A(3)*B(2)
  C(2) = A(3)*B(1)-A(1)*B(3)
  C(3) = A(1)*B(2)-A(2)*B(1)

END SUBROUTINE cross


SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
  IMPLICIT NONE

  !-------------------------------------------------------------------------------
  !Declarations
  !-------------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: errorflag      !Return error status. -1 for error, 0 for normal
  REAL, INTENT(IN), DIMENSION(n,n) :: matrix    !Input matrix
  REAL, INTENT(OUT), DIMENSION(n,n) :: inverse    !Inverted matrix
  
  LOGICAL :: FLAG = .TRUE.
  INTEGER :: i, j, k, l
  REAL :: m
  REAL, DIMENSION(n,2*n) :: augmatrix      !augmented matrix

  !-------------------------------------------------------------------------------  
  !Augment input matrix with an identity matrix
  !-------------------------------------------------------------------------------

  DO i = 1, n
    DO j = 1, 2*n
      IF (j <= n ) THEN
        augmatrix(i,j) = matrix(i,j)
      ELSE IF ((i+n) == j) THEN
        augmatrix(i,j) = 1
      Else
        augmatrix(i,j) = 0
      ENDIF
    END DO
  END DO

  !-------------------------------------------------------------------------------  
  !Reduce augmented matrix to upper traingular form
  !-------------------------------------------------------------------------------

  DO k =1, n-1
    IF (augmatrix(k,k) == 0) THEN
      FLAG = .FALSE.
      DO i = k+1, n
        IF (augmatrix(i,k) /= 0) THEN
          DO j = 1,2*n
            augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
          END DO
          FLAG = .TRUE.
          EXIT
        ENDIF
        IF (FLAG .EQV. .FALSE.) THEN
          PRINT*, "Matrix is non - invertible"
          inverse = 0
          errorflag = -1
          return
        ENDIF
      END DO
    ENDIF
    DO j = k+1, n      
      m = augmatrix(j,k)/augmatrix(k,k)
      DO i = k, 2*n
        augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
      END DO
    END DO
  END DO

  !-------------------------------------------------------------------------------  
  !Test for invertibility
  !-------------------------------------------------------------------------------

  DO i = 1, n
    IF (augmatrix(i,i) == 0) THEN
      PRINT*, "Matrix is non - invertible"
      inverse = 0
      errorflag = -1
      return
    ENDIF
  END DO
  !-------------------------------------------------------------------------------  
  !Make diagonal elements as 1
  !-------------------------------------------------------------------------------

  DO i = 1 , n
    m = augmatrix(i,i)
    DO j = i , (2 * n)        
         augmatrix(i,j) = (augmatrix(i,j) / m)
    END DO
  END DO
  !-------------------------------------------------------------------------------  
  !Reduced right side half of augmented matrix to identity matrix
  !-------------------------------------------------------------------------------

  DO k = n-1, 1, -1
    DO i =1, k
    m = augmatrix(i,k+1)
      DO j = k, (2*n)
        augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
      END DO
    END DO
  END DO        

  !-------------------------------------------------------------------------------  
  !store answer
  !-------------------------------------------------------------------------------

  DO i =1, n
    DO j = 1, n
      inverse(i,j) = augmatrix(i,j+n)
    END DO
  END DO
  errorflag = 0
END SUBROUTINE FINDinv

SUBROUTINE detrend(x, y, trend, perturb, N)

  !-------------------------------------------------------------------------------
  ! Linear detrend
  !-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: N
  REAL :: x(N), y(N), trend(N), perturb(N)
  REAL :: s,t
  REAL :: A,B,C,D,E,F
  INTEGER :: i

  A = 0; B = 0; E = 0; F = 0
  DO i=1,N
    A = A + x(i)*x(i) 
    B = B + x(i)              
    E = E + x(i)*y(i) 
    F = F + y(i)              
  END DO

  C = B
  D = N

  s=(B*F-E*D)/(B*C-A*D)
  t=(E*C-A*F)/(B*C-A*D)

  trend = s*x + t
  perturb = y-trend

END SUBROUTINE detrend

SUBROUTINE write_to_file_flag(fid, flux_record)

  !-------------------------------------------------------------------------------
  ! Write to output file with flag information
  !-------------------------------------------------------------------------------

  USE fluxtype

  IMPLICIT NONE
  INTEGER :: fid
  CHARACTER(15) :: year_char, month_char, mday_char, doy_char, hrmin_char, dtime_char
  CHARACTER(15) :: ust_char, t_bar_char, eta_char, WS_char, FC_char, H_char, LE_char, co2_bar_char, h2o_bar_char, vws_char
  CHARACTER(15) :: u_flag_char, w_flag_char, t_flag_char, h2o_flag_char, co2_flag_char

  TYPE(fluxdata) :: flux_record

  WRITE(year_char, "(I)") flux_record%YEAR
  WRITE(month_char, "(I)") flux_record%MONTH
  WRITE(mday_char, "(I)") flux_record%MDAY
  WRITE(doy_char, "(I)") flux_record%DOY
  WRITE(hrmin_char, "(I4.4)") flux_record%HRMIN
  WRITE(dtime_char, "(F10.5)") flux_record%DTIME
  WRITE(ust_char, "(F10.3)") flux_record%UST
  WRITE(t_bar_char, "(F10.2)") flux_record%TA
  WRITE(eta_char, "(F10.3)") flux_record%WD
  WRITE(WS_char, "(F10.3)") flux_record%WS
  WRITE(FC_char, "(F10.3)") flux_record%FC
  WRITE(H_char, "(F10.3)") flux_record%H
  WRITE(LE_char, "(F10.3)") flux_record%LE
  WRITE(co2_bar_char, "(F10.3)") flux_record%CO2
  WRITE(h2o_bar_char, "(F10.3)") flux_record%H2O
  WRITE(vws_char, "(F10.3)") flux_record%VWS
  WRITE(u_flag_char, "(I)") flux_record%U_FLAG
  WRITE(w_flag_char, "(I)") flux_record%W_FLAG
  WRITE(t_flag_char, "(I)") flux_record%T_FLAG
  WRITE(h2o_flag_char, "(I)") flux_record%H2O_FLAG
  WRITE(co2_flag_char, "(I)") flux_record%CO2_FLAG

  WRITE(fid,"(A)") trim(adjustl(year_char))//","//trim(adjustl(month_char))//","//trim(adjustl(mday_char))//","//trim(adjustl(doy_char))//","//&
    trim(adjustl(hrmin_char))//","//trim(adjustl(dtime_char))//","//trim(adjustl(ust_char))//","//trim(adjustl(t_bar_char))//","//&
    trim(adjustl(eta_char))//","//trim(adjustl(WS_char))//","//trim(adjustl(vws_char))//","//trim(adjustl(FC_char))//","//trim(adjustl(H_char))//","//&
    trim(adjustl(LE_char))//","//trim(adjustl(co2_bar_char))//","//trim(adjustl(h2o_bar_char))//","//trim(adjustl(u_flag_char))//","//&
    trim(adjustl(w_flag_char))//","//trim(adjustl(t_flag_char))//","//trim(adjustl(h2o_flag_char))//","//trim(adjustl(co2_flag_char))

  WRITE(*,"(A)") "  "//trim(adjustl(year_char))//","//trim(adjustl(month_char))//","//trim(adjustl(mday_char))//","//trim(adjustl(doy_char))//","//&
    trim(adjustl(hrmin_char))//","//trim(adjustl(dtime_char))//","//trim(adjustl(ust_char))//","//trim(adjustl(t_bar_char))//","//&
    trim(adjustl(eta_char))//","//trim(adjustl(WS_char))//","//trim(adjustl(vws_char))//","//trim(adjustl(FC_char))//","//trim(adjustl(H_char))//","//&
    trim(adjustl(LE_char))//","//trim(adjustl(co2_bar_char))//","//trim(adjustl(h2o_bar_char))//","//trim(adjustl(u_flag_char))//","//&
    trim(adjustl(w_flag_char))//","//trim(adjustl(t_flag_char))//","//trim(adjustl(h2o_flag_char))//","//trim(adjustl(co2_flag_char))
!  WRITE(*,"('  Instrument = ',F,', Nonstationarity = ',F, ', SEB = ',F)") flux_record%instrument_par, flux_record%nonstationarity_par, flux_record%rn_par
  WRITE(*,"('  Instrument = ',F,', Nonstationarity = ',F)") flux_record%instrument_par, flux_record%nonstationarity_par
  WRITE(*,"('  U  :',F,F,F,F,F,F,F,F)") flux_record%U_spikes_par, flux_record%U_amplitude_resolution_par, flux_record%U_dropouts_par, flux_record%U_extreme_dropouts_par, &
    flux_record%U_skewness_par, flux_record%U_kurtosis_par, flux_record%U_haar_mean_par, flux_record%U_haar_variance_par
  WRITE(*,"('  W  :',F,F,F,F,F,F,F,F)") flux_record%W_spikes_par, flux_record%W_amplitude_resolution_par, flux_record%W_dropouts_par, flux_record%W_extreme_dropouts_par, &
    flux_record%W_skewness_par, flux_record%W_kurtosis_par, flux_record%W_haar_mean_par, flux_record%W_haar_variance_par
  WRITE(*,"('  T  :',F,F,F,F,F,F,F,F)") flux_record%T_spikes_par, flux_record%T_amplitude_resolution_par, flux_record%T_dropouts_par, flux_record%T_extreme_dropouts_par, &
    flux_record%T_skewness_par, flux_record%T_kurtosis_par, flux_record%T_haar_mean_par, flux_record%T_haar_variance_par
  WRITE(*,"('  H2O:',F,F,F,F,F,F,F,F)") flux_record%H2O_spikes_par, flux_record%H2O_amplitude_resolution_par, flux_record%H2O_dropouts_par, flux_record%H2O_extreme_dropouts_par, &
    flux_record%H2O_skewness_par, flux_record%H2O_kurtosis_par, flux_record%H2O_haar_mean_par, flux_record%H2O_haar_variance_par
  WRITE(*,"('  CO2:',F,F,F,F,F,F,F,F)") flux_record%CO2_spikes_par, flux_record%CO2_amplitude_resolution_par, flux_record%CO2_dropouts_par, flux_record%CO2_extreme_dropouts_par, &
    flux_record%CO2_skewness_par, flux_record%CO2_kurtosis_par, flux_record%CO2_haar_mean_par, flux_record%CO2_haar_variance_par

END SUBROUTINE write_to_file_flag

SUBROUTINE write_to_file_no_flag(fid, flux_record)

  !-------------------------------------------------------------------------------
  ! Write to output file WITHOUT flag information
  !-------------------------------------------------------------------------------

  USE fluxtype

  IMPLICIT NONE
  INTEGER :: fid
  CHARACTER(15) :: year_char, month_char, mday_char, doy_char, hrmin_char, dtime_char
  CHARACTER(15) :: ust_char, t_bar_char, eta_char, WS_char, FC_char, H_char, LE_char, co2_bar_char, h2o_bar_char, vws_char

  TYPE(fluxdata) :: flux_record

  WRITE(year_char, "(I)") flux_record%YEAR
  WRITE(month_char, "(I)") flux_record%MONTH
  WRITE(mday_char, "(I)") flux_record%MDAY
  WRITE(doy_char, "(I)") flux_record%DOY
  WRITE(hrmin_char, "(I4.4)") flux_record%HRMIN
  WRITE(dtime_char, "(F10.5)") flux_record%DTIME


  IF (flux_record%U_FLAG + flux_record%W_FLAG > 0) flux_record%UST = -999
  IF (flux_record%T_FLAG > 0) flux_record%TA = -999
  IF (flux_record%U_FLAG > 0) THEN
    flux_record%WD = -999
    flux_record%WS = -999
  END IF
  IF (flux_record%W_FLAG + flux_record%CO2_FLAG + flux_record%H2O_FLAG >0) THEN
    flux_record%FC = -999
    flux_record%LE = -999
  END IF
  IF (flux_record%W_FLAG + flux_record%T_FLAG > 0) flux_record%H = -999
  IF (flux_record%CO2_FLAG + flux_record%H2O_FLAG >0) THEN
    flux_record%CO2 = -999
    flux_record%H2O = -999
  END IF
  IF (flux_record%W_FLAG > 0) flux_record%VWS = -999


  IF (flux_record%UST == -999) THEN
    WRITE(ust_char, "(I)") INT(flux_record%UST)
  ELSE  
    WRITE(ust_char, "(F10.3)") flux_record%UST
  END IF

  IF (flux_record%TA == -999) THEN
    WRITE(t_bar_char, "(I)") INT(flux_record%TA)
  ELSE
    WRITE(t_bar_char, "(F10.2)") flux_record%TA
  END IF

  IF (flux_record%WD == -999) THEN
    WRITE(eta_char, "(I)") INT(flux_record%WD)
  ELSE
    WRITE(eta_char, "(F10.3)") flux_record%WD
  END IF

  IF (flux_record%WS == -999) THEN
    WRITE(WS_char, "(I)") INT(flux_record%WS)
  ELSE
    WRITE(WS_char, "(F10.3)") flux_record%WS
  END IF

  IF (flux_record%FC == -999) THEN
    WRITE(FC_char, "(I)") INT(flux_record%FC)
  ELSE
    WRITE(FC_char, "(F10.3)") flux_record%FC
  END IF

  IF (flux_record%H == -999) THEN
    WRITE(H_char, "(I)") INT(flux_record%H)
  ELSE
    WRITE(H_char, "(F10.3)") flux_record%H
  END IF

  IF (flux_record%LE == -999) THEN
    WRITE(LE_char, "(I)") INT(flux_record%LE)
  ELSE
    WRITE(LE_char, "(F10.3)") flux_record%LE
  END IF

  IF (flux_record%CO2 == -999) THEN
    WRITE(co2_bar_char, "(I)") INT(flux_record%CO2)
  ELSE
    WRITE(co2_bar_char, "(F10.3)") flux_record%CO2
  END IF

  IF (flux_record%H2O == -999) THEN
    WRITE(h2o_bar_char, "(I)") INT(flux_record%H2O)
  ELSE
    WRITE(h2o_bar_char, "(F10.3)") flux_record%H2O
  END IF

  IF (flux_record%VWS == -999) THEN
    WRITE(vws_char, "(I)") INT(flux_record%VWS)
  ELSE
    WRITE(vws_char, "(F10.3)") flux_record%VWS
  END IF


  WRITE(fid,"(A)") trim(adjustl(year_char))//","//trim(adjustl(month_char))//","//trim(adjustl(mday_char))//","//trim(adjustl(doy_char))//","//&
    trim(adjustl(hrmin_char))//","//trim(adjustl(dtime_char))//","//trim(adjustl(ust_char))//","//trim(adjustl(t_bar_char))//","//&
    trim(adjustl(eta_char))//","//trim(adjustl(WS_char))//","//trim(adjustl(vws_char))//","//trim(adjustl(FC_char))//","//trim(adjustl(H_char))//","//&
    trim(adjustl(LE_char))//","//trim(adjustl(co2_bar_char))//","//trim(adjustl(h2o_bar_char))

END SUBROUTINE write_to_file_no_flag

SUBROUTINE read_pressure_file(Precord, pressure_file)
  USE time
  USE Pressure

  IMPLICIT NONE

  CHARACTER(100) :: buffer
  CHARACTER(200) :: pressure_file
  INTEGER :: RECORD, error
  REAL :: P, H2O1, H2O2, T1, T2, RH, Rn !PAR, LWS
  INTEGER :: ind
  TYPE(tm) :: timeinfo
  TYPE(Pressuredata) :: Precord

  ind = 1
  error = 0

  OPEN(UNIT=600, FILE=TRIM(ADJUSTL(pressure_file)), ACTION='read',IOSTAT=error)

  IF (error/=0) GOTO 1000

  DO WHILE(.TRUE.)
!    READ(600,*,IOSTAT = error) buffer, RECORD, P, H2O1, H2O2, LWS, T1, T2, RH, Rn, PAR
    READ(600,*,IOSTAT = error) buffer, RECORD, P, T1, T2, RH, H2O1, H2O2
    IF (error < 0 ) THEN
      EXIT
    ELSE IF (error > 0) THEN
      CYCLE
    END IF

    Precord%pressure(ind) = P

    READ(buffer(:4),*) timeinfo%tm_year
    READ(buffer(6:7),*) timeinfo%tm_mon
    READ(buffer(9:10),*) timeinfo%tm_mday
    READ(buffer(12:13),*) timeinfo%tm_hour
    READ(buffer(15:16),*) timeinfo%tm_min
    READ(buffer(18:),*) timeinfo%tm_sec

    CALL timegm(timeinfo, Precord%time(ind))
    ind = ind + 1
  END DO
  CLOSE(600)
1000 CONTINUE
END SUBROUTINE read_pressure_file

SUBROUTINE read_pressure(Precord, time, P)

  USE pressure

  IMPLICIT NONE

  INTEGER :: time
  TYPE(Pressuredata) :: Precord
  REAL :: P
  INTEGER :: counter, i

  counter = 0
  P = 0
  DO i = 1, 100000
    IF (Precord%time(i)>time-60*30 .AND. Precord%time(i)<=time) THEN
      P = P + Precord%pressure(i)
      counter = counter + 1
    ELSE
      IF (Precord%time(i) > time) THEN
        EXIT
      END IF
    END IF
  END DO


  P = P/REAL(counter)

END SUBROUTINE read_pressure
