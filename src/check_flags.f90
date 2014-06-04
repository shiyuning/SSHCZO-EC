PROGRAM check_flags

  IMPLICIT NONE
 
  CHARACTER(100) :: buffer, inputname
  INTEGER :: error

  INTEGER :: datayear, datamonth
  INTEGER :: year, month, day, hour, minute, ind, diag
  REAL :: second, time, Ux, Uy, Uz, Ts, co2, h2o, h2o_hmp, T_hmp
  INTEGER :: i,j
  INTEGER :: totalrecord = 0
  INTEGER :: infid = 100
  INTEGER :: flags(13) = 0
  REAL :: sync(50000)

  CALL get_command_argument(1, buffer)			! Read year and month from command line

  READ(buffer(1:4),*) datayear
  READ(buffer(6:7),*) datamonth

  WRITE(*,"('\n  Checking ',I4.4,'-',I2.2,' data')") datayear, datamonth
  WRITE(inputname,"(I4.4,'-',I2.2,'/',I4.4,'-',I2.2,'.dat')") datayear, datamonth, datayear, datamonth
  OPEN(UNIT=infid, FILE=inputname, ACTION='read', IOSTAT=error)
  IF (error/=0) GOTO 1000

  WRITE(*,*)"\n  File opened ...\n"

  DO WHILE (.TRUE.)
    READ(infid,*,iostat = error) buffer, ind, Ux, Uy, Uz, Ts, co2, h2o, diag, h2o_hmp, T_hmp
    IF (error < 0 ) THEN				! If reaches end of file, calculate flux using existing record
      EXIT
    ELSE IF (error > 0) THEN				! If error occurs reading in a record, skip
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

    totalrecord = totalrecord + 1
    DO i = 0,3
      IF(BTEST(diag,i+4)) flags(i) = flags(i) + 1
    END DO
    IF(diag/256 == 9) THEN      ! lost trigger
      flags(9) = flags(9) + 1
    ELSE IF (diag/256 == 10) THEN       ! no data
      flags(10) = flags(10) + 1
    ELSE IF (diag/256 == 11) THEN       ! wrong CSAT3 embedded code
      flags(11) = flags(11) + 1
    ELSE IF (diag/256 == 12) THEN        ! SDM error
      flags(12) = flags(12) + 1
    ELSE IF (diag/256 == 13) THEN        ! NaN
      flags(13) = flags(13) + 1
    ELSE
      DO i = 5,8
        IF (BTEST(diag,i+4)) flags(i) = flags(i) + 1
      END DO
    END IF
  END DO

  WRITE(*,"('Total records:',I)") totalrecord
  WRITE(*,"('Delta temperature warning flags:',I,',',F6.3,'%')") flags(1), REAL(flags(1))/REAL(totalrecord)*100
  WRITE(*,"('Poor signal lock warning flags:',I,',',F6.3,'%')") flags(2),REAL(flags(2))/REAL(totalrecord)*100
  WRITE(*,"('Amplitude high warning flags:',I,',',F6.3,'%')") flags(3),REAL(flags(3))/REAL(totalrecord)*100
  WRITE(*,"('Amplitude low warning flags:',I,',',F6.3,'%')") flags(4),REAL(flags(4))/REAL(totalrecord)*100
  WRITE(*,"('Chopper warning flags:',I,',',F6.3,'%')") flags(5),REAL(flags(5))/REAL(totalrecord)*100
  WRITE(*,"('Detector warning flags:',I,',',F6.3,'%')") flags(6),REAL(flags(6))/REAL(totalrecord)*100
  WRITE(*,"('PLL warning flags:',I,',',F6.3,'%')") flags(7),REAL(flags(7))/REAL(totalrecord)*100
  WRITE(*,"('Synchronization flags:',I,',',F6.3,'%')") flags(8),REAL(flags(8))/REAL(totalrecord)*100
  WRITE(*,"('Lost trigger flags:',I,',',F6.3,'%')") flags(9),REAL(flags(9))/REAL(totalrecord)*100
  WRITE(*,"('No data flags:',I,',',F6.3,'%')") flags(10),REAL(flags(10))/REAL(totalrecord)*100
  WRITE(*,"('Wrong CSAT3 embedded code flags:',I,',',F6.3,'%')") flags(11),REAL(flags(11))/REAL(totalrecord)*100
  WRITE(*,"('SDM error flags:',I,',',F6.3,'%')") flags(12),REAL(flags(12))/REAL(totalrecord)*100
  WRITE(*,"('NaN flags:',I,',',F6.3,'%')") flags(13),REAL(flags(13))/REAL(totalrecord)*100
 
1000 CONTINUE

END PROGRAM
