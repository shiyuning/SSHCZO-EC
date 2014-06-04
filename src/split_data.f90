PROGRAM split_data

  IMPLICIT NONE
 
  CHARACTER(200) :: filename, buffer
  INTEGER :: error
  INTEGER :: i,j
  INTEGER:: year, month
  INTEGER :: day
  INTEGER :: eom(12)
  INTEGER :: fid

  TYPE :: ECdata
    CHARACTER(200) :: buffer
    INTEGER :: year, month, day, hour, minute
    CHARACTER(200) :: Ux, Uy, Uz, Ts, co2, h2o, ind
    REAL :: second
    CHARACTER(200) :: diag
    CHARACTER(200) :: h2o_hmp, T_hmp
  END TYPE ECdata

  INTEGER :: daytemp = 0

  TYPE(ECdata) :: record

  CALL get_command_argument(1, buffer)    ! Read year and month from command line

  READ(buffer(1:4),*) year
  READ(buffer(6:7),*) month

  WRITE(*,*)"\n  Start reading in file ... \n"

  WRITE(filename,"('Data/',I4.4,'-',I2.2,'/',I4.4,'-',I2.2,'.dat')") year, month, year, month
  OPEN(UNIT=100,FILE=filename,ACTION='read',IOSTAT=error)
  IF (error/=0) GOTO 1000

  WRITE(*,*)"\n  File opened ...\n"

  IF ((mod(year,4)==0 .AND. mod(year,100)/=0) .or. mod(year,400)==0) THEN
    eom=(/31,29,31,30,31,30,31,31,30,31,30,31/)
  ELSE
    eom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  END IF

  READ(100,*,IOSTAT = error) buffer
  READ(100,*,IOSTAT = error) buffer
  READ(100,*,IOSTAT = error) buffer
  READ(100,*,IOSTAT = error) buffer

  DO WHILE(.TRUE.)
    READ(100,*,IOSTAT = error) record%buffer, record%ind, record%Ux, record%Uy, record%Uz, record%Ts, record%co2, record%h2o, record%diag, record%h2o_hmp, record%T_hmp
    IF (error /= 0 ) EXIT
    READ(record%buffer(:4),*) record%year
    READ(record%buffer(6:7),*) record%month
    READ(record%buffer(9:10),*) record%day
    READ(record%buffer(12:13),*) record%hour
    READ(record%buffer(15:16),*) record%minute
    READ(record%buffer(18:),*) record%second

!    write(*,*) record%day, record%second

    IF (record%day /= daytemp) THEN
      WRITE(*,"('\n  Now copying ',I2.2,'-',I2.2,' data...')") record%month,record%day
      IF (record%day>1) CLOSE(200)
      WRITE(filename,"('Data/',I4.4,'-',I2.2,'/',I4.4,'-',I2.2,'-',I2.2,'.dat')") year, month, year, month, record%day
      OPEN(UNIT=200,FILE=filename,ACTION='write',IOSTAT=error)
      IF (error/=0) GOTO 1000
      WRITE(200,"(A)") """TOA5"",""CZO_EC1"",""CR1000"",""4745"",""CR1000.Std.16"",""CPU:CZO_EC_d2.cr1"",""61052"",""ts_data"""
      WRITE(200,"(A)") """TIMESTAMP"",""RECORD"",""Ux"",""Uy"",""Uz"",""Ts"",""co2"",""h2o"",""diag"",""h2o_hmp"",""T_hmp"""
      WRITE(200,"(A)") """TS"",""RN"",""m/s"",""m/s"",""m/s"",""C"",""mg/(m^3)"",""g/(m^3)"",""unitless"",""g/(m^3)"",""C"""
      WRITE(200,"(A)") """"","""",""Smp"",""Smp"",""Smp"",""Smp"",""Smp"",""Smp"",""Smp"",""Smp"",""Smp"""
      daytemp = record%day
    END IF
    CALL write_to_file(200, record)
  END DO

  CLOSE(100)
  CLOSE(200)

  WRITE(*,*)"\n  done.\n"

1000 CONTINUE

END

SUBROUTINE write_to_file(fid, record)

  IMPLICIT NONE
  INTEGER :: fid
  TYPE :: ECdata
    CHARACTER(200) :: buffer
    INTEGER :: year, month, day, hour, minute
    CHARACTER(200) :: Ux, Uy, Uz, Ts, co2, h2o, ind
    REAL :: second
    CHARACTER(200) :: diag
    CHARACTER(200) :: h2o_hmp, T_hmp
  END TYPE ECdata

  TYPE(ECdata) :: record

  WRITE(fid,"(A)") """"//trim(adjustl(record%buffer))//""""//","//trim(adjustl(record%ind))//","//trim(adjustl(record%Ux))//","//trim(adjustl(record%Uy))//","//&
    trim(adjustl(record%Uz))//","//trim(adjustl(record%Ts))//","//trim(adjustl(record%co2))//","//trim(adjustl(record%h2o))//","//&
    trim(adjustl(record%diag))//","//trim(adjustl(record%h2o_hmp))//","//trim(adjustl(record%T_hmp))

END SUBROUTINE write_to_file
