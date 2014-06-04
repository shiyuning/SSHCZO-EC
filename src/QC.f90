!-------------------------------------------------------------------------------
! Module for quality control of eddy-covariance data
! Yuning Shi, Penn State Meteorology
! Email: yshi@psu.edu
!-------------------------------------------------------------------------------

MODULE QC

  USE ECtype
  USE fluxtype
  USE time
  IMPLICIT NONE

  !-------------------------------------------------------------------------------
  ! Thresholds for quality control
  !-------------------------------------------------------------------------------

  REAL, PARAMETER :: INSTRUMENT_THRESHOLD = 0.95
  REAL, PARAMETER :: SPIKE_THRESHOLD = 0.03
  REAL, PARAMETER :: EMPTY_BINS_THRESHOLD = 0.7
  REAL, PARAMETER :: DROPOUTS_THRESHOLD = 0.1, EXTREME_DROPOUTS_THRESHOLD = 0.06
  REAL, PARAMETER :: SKEWNESS_THRESHOLD = 3.0, KURTOSIS_THRESHOLD = 5.0
  REAL, PARAMETER :: DISCONTINUITIES_THRESHOLD = 3.0
  REAL, PARAMETER :: NONSTATIONARITY_THRESHOLD = 0.3
  

CONTAINS

  SUBROUTINE instrument(record, INSTRUMENT_FLAG, flux_record)

  !-------------------------------------------------------------------------------
  ! If available records for current time period is less than a threshold, or is 
  ! more than 18000, a flag is placed
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: no_of_diag, no_of_co2, no_of_h2o, no_of_instrument
    INTEGER :: INSTRUMENT_FLAG
    INTEGER :: i, j, k
    TYPE(ECdata) :: record
    TYPE(fluxdata) :: flux_record

    IF (REAL(record%record_no)/(30*60*10.0) < INSTRUMENT_THRESHOLD) INSTRUMENT_FLAG = 1
    IF (REAL(record%record_no)/(30*60*10.0) > 1.0) THEN
      INSTRUMENT_FLAG = 1
      WRITE(*,*) "\nWARNING: Timestamp error!"
    END IF
    flux_record%instrument_par = REAL(record%record_no)/(30*60*10.0)

  END SUBROUTINE instrument

  SUBROUTINE spikes(record, SPIKES_FLAG, flux_record)

    !-------------------------------------------------------------------------------
    ! Quality controal
    ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
    ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
    !-------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(ECdata) :: record
    TYPE(fluxdata) :: flux_record
    INTEGER :: SPIKES_FLAG(5)
    REAL :: temp

    flux_record%U_spikes_par = 0

    CALL despike(record%time, record%Ux, record%record_no, SPIKES_FLAG(1), temp)
    flux_record%U_spikes_par = MAX(temp, flux_record%U_spikes_par)
    CALL despike(record%time, record%Uy, record%record_no, SPIKES_FLAG(1), temp)
    flux_record%U_spikes_par = MAX(temp, flux_record%U_spikes_par)
    CALL despike(record%time, record%Uz, record%record_no, SPIKES_FLAG(2), flux_record%W_spikes_par)
    CALL despike(record%time, record%co2, record%record_no, SPIKES_FLAG(5), flux_record%CO2_spikes_par)
    CALL despike(record%time, record%h2o, record%record_no, SPIKES_FLAG(4), flux_record%H2O_spikes_par)
    CALL despike(record%time, record%Ts, record%record_no, SPIKES_FLAG(3), flux_record%T_spikes_par)

  END SUBROUTINE spikes

  SUBROUTINE despike(time, rawdata, record_no, SPIKES_FLAG, spikes_par)
  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: record_no
    REAL :: rawdata(record_no), time(record_no)
    INTEGER :: L1, window_starter 
    REAL,ALLOCATABLE :: window(:)
    REAL :: window_mean, window_sd
    REAL :: SD_THRESHOLD
    INTEGER, PARAMETER :: CONSECUTIVE_SPIKES = 3
    INTEGER :: i,j,k
    INTEGER :: no_of_spikes
    INTEGER :: SPIKES_FLAG
    REAL :: spikes_par

    L1 = 3000
    window_starter = 1
    SD_THRESHOLD = 4.5
    spikes_par = 0
    no_of_spikes = 0

    L1 = MIN(L1, record_no)
    ALLOCATE(window(L1))

    DO WHILE (.TRUE.)
      DO WHILE (window_starter+L1-1<=record_no)
        window = rawdata(window_starter:window_starter+L1-1)
        window_mean = SUM(window)/REAL(L1)
        window_sd = SQRT(SUM((window-window_mean)*(window-window_mean))/REAL(L1))

        i = 2

        DO WHILE (i < L1)
          IF (ABS(window(i)-window_mean) > SD_THRESHOLD*window_sd) THEN
            j = i + 1
            DO WHILE (j< L1)
              IF (ABS(window(j)-window_mean) < SD_THRESHOLD*window_sd) EXIT
              j = j + 1
            END DO
            IF (j>i+CONSECUTIVE_SPIKES) THEN
              i = j
            ELSE
              DO k = i, j-1
                no_of_spikes = no_of_spikes + 1
                window(k) = window(i-1)+(time(k)-time(i-1))*(window(j)-window(i-1))/(time(j)-time(i-1))
                rawdata(window_starter+k-1) = window(k)
              END DO
              i = j
            END IF
          ELSE
            i = i+1
          END IF
        END DO
        window_starter = window_starter + 1
      END DO

      IF (no_of_spikes == 0) THEN
        EXIT
      ELSE
        IF (REAL(no_of_spikes)/REAL(record_no)>SPIKE_THRESHOLD) SPIKES_FLAG = 1
        spikes_par = MAX(spikes_par, REAL(no_of_spikes)/REAL(record_no))
        window_starter = 1
        no_of_spikes = 0
        SD_THRESHOLD = SD_THRESHOLD + 0.1
      END IF
    END DO
  END SUBROUTINE despike

  SUBROUTINE amplitude_resolution(record, AMPLITUDE_RESOLUTION_FLAG, flux_record)

    IMPLICIT NONE

    TYPE(ECdata) :: record
    TYPE(fluxdata) :: flux_record
    INTEGER :: AMPLITUDE_RESOLUTION_FLAG(5)
    REAL :: temp

    flux_record%U_amplitude_resolution_par = 0

    CALL check_amplitude_resolution(record%Ux, record%record_no, AMPLITUDE_RESOLUTION_FLAG(1), temp)
    flux_record%U_amplitude_resolution_par = MAX(temp, flux_record%U_amplitude_resolution_par)
    CALL check_amplitude_resolution(record%Uy, record%record_no, AMPLITUDE_RESOLUTION_FLAG(1), temp)
    flux_record%U_amplitude_resolution_par = MAX(temp, flux_record%U_amplitude_resolution_par)
    CALL check_amplitude_resolution(record%Uz, record%record_no, AMPLITUDE_RESOLUTION_FLAG(2), flux_record%W_amplitude_resolution_par)
    CALL check_amplitude_resolution(record%co2, record%record_no, AMPLITUDE_RESOLUTION_FLAG(5), flux_record%CO2_amplitude_resolution_par)
    CALL check_amplitude_resolution(record%h2o, record%record_no, AMPLITUDE_RESOLUTION_FLAG(4), flux_record%H2O_amplitude_resolution_par)
    CALL check_amplitude_resolution(record%Ts, record%record_no, AMPLITUDE_RESOLUTION_FLAG(3), flux_record%T_amplitude_resolution_par)

  END SUBROUTINE amplitude_resolution

  SUBROUTINE check_amplitude_resolution(rawdata, record_no, AMPLITUDE_RESOLUTION_FLAG, amplitude_resolution_par)

  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: record_no
    REAL :: rawdata(record_no)
    INTEGER :: L1
    INTEGER :: window_starter
    REAL,ALLOCATABLE :: window(:)
    REAL :: window_mean, window_sd
    INTEGER :: i,j,k
    INTEGER :: AMPLITUDE_RESOLUTION_FLAG
    REAL :: empty_bins
    INTEGER,ALLOCATABLE :: bin(:)
    REAL :: bin_range, bin_interval
    REAL :: amplitude_resolution_par

    L1 = 1000
    amplitude_resolution_par = 0
    L1 = MIN(L1, record_no)
    ALLOCATE(window(L1))
    ALLOCATE(bin(L1))
    window_starter = 1
    bin = 0

    DO WHILE (window_starter+L1-1<=record_no)
      window = rawdata(window_starter:window_starter+L1-1)
      window_mean = SUM(window)/REAL(L1)
      window_sd = SQRT(SUM((window-window_mean)*(window-window_mean))/REAL(L1))
      bin_range = MIN(7*window_sd, MAXVAL(window)-MINVAL(window))
      bin_interval = bin_range/100
      bin = CEILING((window-((MAXVAL(window)+MINVAL(window))/2-bin_range/2))/bin_interval)

      empty_bins = 0

      DO i = 1,100
        IF (ANY(bin == i)) THEN
          CYCLE
        ELSE
          empty_bins = empty_bins + 1
        END IF
      END DO
      amplitude_resolution_par = MAX(amplitude_resolution_par, empty_bins/100)
      IF (empty_bins > EMPTY_BINS_THRESHOLD*100) AMPLITUDE_RESOLUTION_FLAG = 1
      window_starter = window_starter + 250
    END DO

  END SUBROUTINE check_amplitude_resolution

  SUBROUTINE dropouts(record, DROPOUTS_FLAG, flux_record)

    IMPLICIT NONE

    TYPE(ECdata) :: record
    TYPE(fluxdata) :: flux_record
    INTEGER :: DROPOUTS_FLAG(5)
    REAL :: temp1, temp2


    CALL check_dropouts(record%Ux, record%record_no, DROPOUTS_FLAG(1), temp1, temp2)
    CALL check_dropouts(record%Uy, record%record_no, DROPOUTS_FLAG(1), flux_record%U_dropouts_par, flux_record%U_extreme_dropouts_par)
    flux_record%U_dropouts_par = MAX(temp1, flux_record%U_dropouts_par)
    flux_record%U_extreme_dropouts_par = MAX(temp2, flux_record%U_extreme_dropouts_par)
    CALL check_dropouts(record%Uz, record%record_no, DROPOUTS_FLAG(2), flux_record%W_dropouts_par, flux_record%W_extreme_dropouts_par)
    CALL check_dropouts(record%h2o, record%record_no, DROPOUTS_FLAG(4), flux_record%H2O_dropouts_par, flux_record%H2O_extreme_dropouts_par)
    CALL check_dropouts(record%co2, record%record_no, DROPOUTS_FLAG(5), flux_record%CO2_dropouts_par, flux_record%CO2_extreme_dropouts_par)
    CALL check_dropouts(record%Ts, record%record_no, DROPOUTS_FLAG(3), flux_record%T_dropouts_par, flux_record%T_extreme_dropouts_par)

  END SUBROUTINE dropouts

  SUBROUTINE check_dropouts(rawdata, record_no, DROPOUTS_FLAG, dropouts_par, extreme_dropouts_par)
  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: record_no
    REAL :: rawdata(record_no)
    INTEGER :: L1, window_starter
    REAL,ALLOCATABLE :: window(:)
    INTEGER,ALLOCATABLE :: bin(:)
    REAL :: window_mean, window_sd
    INTEGER :: i,j,k
    INTEGER :: DROPOUTS_FLAG
    REAL :: no_of_dropouts, no_of_extreme_dropouts
    REAL :: bin_range, bin_interval
    REAL :: dropouts_par, extreme_dropouts_par

    dropouts_par = 0
    extreme_dropouts_par = 0

    L1 = 1000
    L1 = MIN(L1, record_no)
    ALLOCATE(window(L1))
    ALLOCATE(bin(L1))

    bin = 0
    window_starter = 1

    DO WHILE (window_starter+L1-1<=record_no)
      window = rawdata(window_starter:window_starter+L1-1)
      window_mean = SUM(window)/REAL(L1)
      window_sd = SQRT(SUM((window-window_mean)*(window-window_mean))/REAL(L1))
      bin_range = MIN(7*window_sd, MAXVAL(window)-MINVAL(window))
      bin_interval = bin_range/100
      bin = CEILING((window-((MAXVAL(window)+MINVAL(window))/2-bin_range/2))/bin_interval)

      no_of_dropouts = 0
      no_of_extreme_dropouts = 0

      i = 2

      DO WHILE (i < L1)
        IF (bin(i) == bin(i-1)) THEN
          j = i + 1
          DO WHILE (j <= L1)
            IF (bin(j)/=bin(j-1)) EXIT
            j = j+1
          END DO
          no_of_dropouts = MAX(no_of_dropouts, REAL(j-i+1))
          IF (bin(i)<10 .OR. bin(i)>90) no_of_extreme_dropouts = no_of_dropouts
          i = j + 1
        ELSE
          i = i + 1
        END IF
      END DO

!      write(*,"('Number of extreme dropouts',F)") no_of_extreme_dropouts
      dropouts_par = MAX(dropouts_par, no_of_dropouts/REAL(L1))
      extreme_dropouts_par = MAX(extreme_dropouts_par, no_of_extreme_dropouts/REAL(L1))

      IF (no_of_dropouts>DROPOUTS_THRESHOLD*REAL(L1) .OR. no_of_extreme_dropouts>EXTREME_DROPOUTS_THRESHOLD*REAL(L1)) DROPOUTS_FLAG = 1
      window_starter = window_starter + 250
    END DO


  END SUBROUTINE check_dropouts


  SUBROUTINE absolute_limits(record, ABSOLUTE_LIMITS_FLAG)

  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: ABSOLUTE_LIMITS_FLAG(5)

    TYPE(ECdata) :: record

    IF (ANY(ABS(record%Ux(1:record%record_no)) > 30)) THEN
      ABSOLUTE_LIMITS_FLAG(1) = 1
    END IF
    IF (ANY(ABS(record%Uy(1:record%record_no)) > 30)) THEN
      ABSOLUTE_LIMITS_FLAG(1) = 1
    END IF
    IF (ANY(ABS(record%Uz(1:record%record_no)) > 10)) THEN
      ABSOLUTE_LIMITS_FLAG(2) = 1
    END IF
    IF (ANY(record%Ts(1:record%record_no) > 60) .OR. ANY(record%Ts(1:record%record_no)<-50) .OR. MAXVAL(record%Ts(1:record%record_no))-MINVAL(record%Ts(1:record%record_no))>10) THEN
      ABSOLUTE_LIMITS_FLAG(3) = 1
    END IF
    IF (ANY(record%h2o(1:record%record_no) > 35) .OR. ANY(record%h2o(1:record%record_no)<2.5) .OR. MAXVAL(record%h2o(1:record%record_no))-MINVAL(record%h2o(1:record%record_no))>8) THEN
      ABSOLUTE_LIMITS_FLAG(4) = 1
!      WRITE(*,*) ABSOLUTE_LIMITS_FLAG(4)
    END IF
    IF (ANY(record%co2(1:record%record_no) > 950) .OR. ANY(record%co2(1:record%record_no)<550) .OR. MAXVAL(record%co2(1:record%record_no))-MINVAL(record%co2(1:record%record_no))>120) THEN
      ABSOLUTE_LIMITS_FLAG(5) = 1
    END IF

  END SUBROUTINE absolute_limits


  SUBROUTINE higher_moment(record, HIGHER_MOMENT_FLAG, flux_record)

  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(ECdata) :: record
    TYPE(fluxdata) :: flux_record
    INTEGER :: HIGHER_MOMENT_FLAG(5)
    REAL :: temp1, temp2

    CALL check_higher_moment(record%Ux, record%record_no, HIGHER_MOMENT_FLAG(1), temp1, temp2)
    CALL check_higher_moment(record%Uy, record%record_no, HIGHER_MOMENT_FLAG(1), flux_record%U_skewness_par, flux_record%U_kurtosis_par)
    flux_record%U_skewness_par = MAX(temp1, flux_record%U_skewness_par)
    flux_record%U_kurtosis_par = MAX(temp2, flux_record%U_kurtosis_par)
    CALL check_higher_moment(record%Uz, record%record_no, HIGHER_MOMENT_FLAG(2), flux_record%W_skewness_par, flux_record%W_kurtosis_par)
    CALL check_higher_moment(record%co2, record%record_no, HIGHER_MOMENT_FLAG(5), flux_record%CO2_skewness_par, flux_record%CO2_kurtosis_par)
    CALL check_higher_moment(record%h2o, record%record_no, HIGHER_MOMENT_FLAG(4), flux_record%H2O_skewness_par, flux_record%H2O_kurtosis_par)
    CALL check_higher_moment(record%Ts, record%record_no, HIGHER_MOMENT_FLAG(3), flux_record%T_skewness_par, flux_record%T_kurtosis_par)

  END SUBROUTINE higher_moment

  SUBROUTINE check_higher_moment(rawdata, record_no, HIGHER_MOMENT_FLAG, skewness_par, kurtosis_par)

  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: record_no
    REAL :: rawdata(record_no)
    INTEGER :: i,j,k
    INTEGER :: HIGHER_MOMENT_FLAG
    REAL :: miu, sigma
    REAL :: skewness, kurtosis
    REAL :: skewness_par, kurtosis_par

    skewness_par = 0
    kurtosis_par = 0

    miu = SUM(rawdata)/REAL(record_no)
    sigma = SQRT(SUM((rawdata-miu)*(rawdata-miu))/REAL(record_no))

    skewness = SUM((rawdata - miu)**3)/REAL(record_no)/(sigma**3)
    kurtosis = SUM((rawdata - miu)**4)/REAL(record_no)/(sigma**4)

    skewness_par = ABS(skewness)
    kurtosis_par = ABS(kurtosis - 5)
    IF (ABS(skewness) > SKEWNESS_THRESHOLD) HIGHER_MOMENT_FLAG = 1
    IF (ABS(kurtosis -5) > KURTOSIS_THRESHOLD) HIGHER_MOMENT_FLAG = 1

  END SUBROUTINE check_higher_moment

  SUBROUTINE discontinuities(record, DISCONTINUITIES_FLAG, flux_record)

  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(ECdata) :: record
    TYPE(fluxdata) :: flux_record
    INTEGER :: DISCONTINUITIES_FLAG(5)
    REAL :: temp1, temp2
!    WRITE(*,*) "Ux"
    CALL check_discontinuities(record%Ux, record%record_no, DISCONTINUITIES_FLAG(1), flux_record%U_haar_mean_par, flux_record%U_haar_variance_par)
!    WRITE(*,*) "Uy"
    CALL check_discontinuities(record%Uy, record%record_no, DISCONTINUITIES_FLAG(1), temp1, temp2)
    flux_record%U_haar_mean_par = MAX(temp1,  flux_record%U_haar_mean_par)
    flux_record%U_haar_variance_par = MAX(temp2,  flux_record%U_haar_variance_par)
!    WRITE(*,*) "Uz"
    CALL check_discontinuities(record%Uz, record%record_no, DISCONTINUITIES_FLAG(2), flux_record%W_haar_mean_par, flux_record%W_haar_variance_par)
!    WRITE(*,*) "H2O"
    CALL check_discontinuities(record%h2o, record%record_no, DISCONTINUITIES_FLAG(4), flux_record%H2O_haar_mean_par, flux_record%H2O_haar_variance_par)
!    WRITE(*,*) "CO2"
    CALL check_discontinuities(record%co2, record%record_no, DISCONTINUITIES_FLAG(5), flux_record%CO2_haar_mean_par, flux_record%CO2_haar_variance_par)
!    WRITE(*,*) "Ts"
    CALL check_discontinuities(record%Ts, record%record_no, DISCONTINUITIES_FLAG(3), flux_record%T_haar_mean_par, flux_record%T_haar_variance_par)

  END SUBROUTINE discontinuities

  SUBROUTINE check_discontinuities(rawdata, record_no, DISCONTINUITIES_FLAG, haar_mean_par, haar_variance_par)
  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: record_no
    REAL :: rawdata(record_no)
    INTEGER :: L1, window_starter
    REAL,ALLOCATABLE :: window(:)
    REAL :: record_mean, record_sd, record_range
    REAL :: half1_mean, half2_mean, half1_variance, half2_variance, A
    INTEGER :: i,j,k
    INTEGER :: DISCONTINUITIES_FLAG
    REAL :: haar_mean_par, haar_variance_par


    L1 = 5*60*10
    L1 = MIN(L1, record_no)
    ALLOCATE(window(L1))

    window_starter = 1

    record_mean = SUM(rawdata)/REAL(record_no)
    record_sd = SQRT(SUM((rawdata-record_mean)*(rawdata-record_mean))/REAL(record_no))
    record_range = MAXVAL(rawdata)-MINVAL(rawdata)
    A = MIN(record_sd, record_range/4)

!    WRITE(*,*) A

    haar_mean_par = 0
    haar_variance_par = 0

    DO WHILE (window_starter+L1-1<=record_no)
      window = rawdata(window_starter:window_starter+L1-1)

      half1_mean = SUM(window(1:L1/2))/REAL(L1/2)
      half2_mean = SUM(window(L1/2+1:L1))/REAL(L1/2)
!      WRITE(*,*) ABS(half2_mean-half1_mean)
      half1_variance = SUM((window(1:L1/2)-half1_mean)*(window(1:L1/2)-half1_mean))/REAL(L1/2)
      half2_variance = SUM((window(L1/2+1:L1)-half2_mean)*(window(L1/2+1:L1)-half2_mean))/REAL(L1/2)
      haar_mean_par = MAX(haar_mean_par,ABS((half2_mean-half1_mean)/A))
      haar_variance_par = MAX(haar_variance_par,ABS((half2_variance-half1_variance)/record_sd/record_sd))
      IF (ABS((half2_mean-half1_mean)/A)>DISCONTINUITIES_THRESHOLD) DISCONTINUITIES_FLAG = 1
      IF (ABS((half2_variance-half1_variance)/record_sd/record_sd)>DISCONTINUITIES_THRESHOLD) DISCONTINUITIES_FLAG = 1
      window_starter = window_starter + 1
    END DO
  END SUBROUTINE check_discontinuities


  SUBROUTINE nonstationary(u, v, record_no, NONSTATIONARITY_FLAG, flux_record)

  !-------------------------------------------------------------------------------
  ! Quality controal
  ! (Vickers, D., and L. Mahrt, 1997: Quality control and flx sampling problems
  ! for tower and aircraft data. J. Atmos. Oceanic tech., 14, 512-526)
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(fluxdata) :: flux_record
    INTEGER :: record_no
    REAL :: u(record_no), v(record_no)
    REAL :: u_bar, v_bar
    INTEGER :: NONSTATIONARITY_FLAG
    REAL :: vector_avg, instant_avg

    !-------------------------------------------------------------------------------
    ! wind speed reduction
    !-------------------------------------------------------------------------------

    u_bar = SUM(u)/REAL(record_no)
    v_bar = SUM(v)/REAL(record_no)

    vector_avg = SQRT(u_bar*u_bar+v_bar*v_bar)
    instant_avg = SUM(SQRT(u*u+v*v))/REAL(record_no)

    IF (vector_avg < instant_avg*NONSTATIONARITY_THRESHOLD) NONSTATIONARITY_FLAG = 1
    flux_record%nonstationarity_par = vector_avg/instant_avg

  END SUBROUTINE nonstationary

  SUBROUTINE SEB(rawtime, flux_record, RN_FLAG)

  !-------------------------------------------------------------------------------
  ! Quality controal
  ! Constrain surface heat fluxes using net radiation measurements
  !-------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(fluxdata) :: flux_record
    TYPE(tm) :: timeinfo
    INTEGER :: record_no, rntime, rawtime

    CHARACTER(100) :: buffer, buffer1
    INTEGER :: error, RN_FLAG
    REAL :: Rn, threshold

    error = 0
    OPEN(UNIT=500, FILE='rn.dat', ACTION='read', IOSTAT=error)
    IF (error/=0) THEN
      WRITE(*,*) "\nWARNING: Net radiation data cannot be found!"
      GOTO 1000
    END IF
!    WRITE(*,*)"\n  File opened ...\n"
    DO WHILE (.TRUE.)
      READ(500,*,IOSTAT = error) buffer, buffer1, Rn
      IF (error/=0) THEN
        WRITE(*,*) "\nWARNING: Net radiation read error!"
        GOTO 1000
      END IF
      READ(buffer(:4),*) timeinfo%tm_year
      READ(buffer(6:7),*) timeinfo%tm_mon
      READ(buffer(9:10),*) timeinfo%tm_mday
      READ(buffer1(1:2),*) timeinfo%tm_hour
      READ(buffer1(4:5),*) timeinfo%tm_min
      timeinfo%tm_sec = 0
      CALL timegm(timeinfo, rntime)
      IF (rntime == rawtime) THEN
        WRITE(*,"('  H + LE = ',F,', Rn = ',F)") flux_record%H + flux_record%LE, Rn
        threshold = ABS(Rn)*0.4
        threshold = MAX(threshold, 100.0)
        flux_record%rn_par = flux_record%H + flux_record%LE - Rn
        IF (ABS(flux_record%rn_par)>threshold) THEN
          RN_FLAG = 1
        END IF
        EXIT
      END IF
    END DO
1000 CONTINUE
    CLOSE (500)
  END SUBROUTINE SEB

END MODULE QC
