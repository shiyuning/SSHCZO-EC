PROGRAM fill_gap

	USE time

	IMPLICIT NONE

	INTEGER :: inputfluxid, inputrnid, outputid

	REAL:: dtime, ust, ta, wd, ws, vws, fc, co2, h2o, h, le, Rn, H_sum, LE_sum, Rn_sum
	INTEGER :: year, month, mday, jday, hrmin, hour, minute, rntime
	INTEGER datasize, maxind, icount
	REAL, ALLOCATABLE :: avg_H(:), avg_LE(:)
	INTEGER, ALLOCATABLE:: avg_time(:)

	TYPE :: ecdata
		REAL :: H, LE, Rn
		INTEGER :: time, H_FLAG, LE_FLAG
	END TYPE ecdata

	TYPE(ecdata) :: record(20000)

	INTEGER :: rawtime, time_start, time_end, interval
	TYPE(tm) :: timeinfo_start, timeinfo_end, timeinfo
	INTEGER :: i,j,k, temp
	CHARACTER(60) :: buffer, buffer1
	INTEGER :: error

	inputfluxid = 100
	inputrnid = 200
	outputid = 300

	OPEN(UNIT=inputfluxid,FILE='2009-flux.dat',ACTION='read',IOSTAT=error)

	IF (error/=0) goto 1000

	READ(inputfluxid,*) buffer
	READ(inputfluxid,*) buffer
	READ(inputfluxid,*) buffer
	READ(inputfluxid,*) buffer

!	OPEN(UNIT=outputid,FILE='wind.dat',ACTION='WRITE',IOSTAT=error)

	i = 1

	DO WHILE(.TRUE.)
		READ(inputfluxid,*,IOSTAT = error) year, month, mday, jday, hrmin, dtime, ust, ta, wd, ws, vws, fc, record(i)%H, record(i)%LE, co2, h2o
		IF (error/=0) EXIT
!		CALL doy2date(year, jday, month, day)
		hour = INT(hrmin/100)
		minute = hrmin - hour*100
		timeinfo = tm(year,month,mday,hour,minute,0)
		CALL timegm(timeinfo, rawtime)
		record(i)%time = rawtime
		i = i+1
	END DO

	CLOSE(inputfluxid)

	maxind = i - 1
	i = 1

	OPEN(UNIT=inputrnid, FILE='rn.dat', ACTION='read', IOSTAT=error)
	IF (error/=0) THEN
		WRITE(*,*) "\nWARNING: Net radiation data cannot be found!"
		GOTO 1000
	END IF
	WRITE(*,*)"\n  File opened ...\n"

	DO WHILE (.TRUE.)
		READ(inputrnid,*,IOSTAT = error) buffer, buffer1, Rn
		IF (error/=0) EXIT
		READ(buffer(:4),*) timeinfo%tm_year
		READ(buffer(6:7),*) timeinfo%tm_mon
		READ(buffer(9:10),*) timeinfo%tm_mday
		READ(buffer1(1:2),*) timeinfo%tm_hour
		READ(buffer1(4:5),*) timeinfo%tm_min
		timeinfo%tm_sec = 0
		CALL timegm(timeinfo, rntime)

		temp = 1
		IF (rntime <= record(maxind)%time) THEN
			DO i = temp, maxind
				IF (rntime == record(i)%time) THEN
					record(i)%Rn = Rn
					WRITE(*,*) record(i)%H, record(i)%Rn
					temp = i
					EXIT
				END IF
			END DO
		END IF
	END DO

	CLOSE(inputrnid)

	DO i = 8, maxind-1
		CALL gmtime(record(i)%time, timeinfo)
		H_sum = 0
		Rn_sum = 0
		icount = 0
		IF (record(i)%H == -999) THEN
			IF (record(i-1)%H>-999 .AND. record(i+1)%H>-999) THEN
				IF (timeinfo%tm_hour >=8 .AND. timeinfo%tm_hour <=17) THEN
					record(i)%H = (record(i-1)%H+record(i+1)%H)/(record(i-1)%Rn+record(i+1)%Rn)*record(i)%Rn
				ELSE
					record(i)%H = (record(i-1)%H+record(i+1)%H)/2
				END IF
!			WRITE(*,*) record(i)%H, record(i-1)%H, record(i+1)%H, record(i-1)%Rn, record(i+1)%Rn
				record(i)%H_FLAG = 1
			ELSE
				DO j = -7, 7
					IF (i-j*48>0 .AND. i-j*48<=maxind) THEN
						IF (record(i-j*48)%H>-999) THEN
							icount = icount + 1
							H_sum = H_sum + record(i-j*48)%H
							Rn_sum = Rn_sum + record(i-j*48)%Rn
						END IF
					END IF
				END DO

				IF (icount >0) THEN
					IF (timeinfo%tm_hour >=8 .AND. timeinfo%tm_hour <=17) THEN
						record(i)%H = H_sum/Rn_sum*record(i)%Rn
					ELSE
						record(i)%H = H_sum/icount
					END IF
!					WRITE(*,*) record(i)%H
					record(i)%H_FLAG = 2
				END IF
			END IF
		END IF

		LE_sum = 0
		Rn_sum = 0
		icount = 0
		IF (record(i)%LE == -999) THEN
			IF (record(i-1)%LE>-999 .AND. record(i+1)%LE>-999) THEN
				IF (timeinfo%tm_hour >=8 .AND. timeinfo%tm_hour <=17) THEN
					record(i)%LE = (record(i-1)%LE+record(i+1)%LE)/(record(i-1)%Rn+record(i+1)%Rn)*record(i)%Rn
				ELSE
					record(i)%LE = (record(i-1)%LE+record(i+1)%LE)/2
				END IF
				record(i)%LE_FLAG = 1
			ELSE
				DO j = -7, 7
					IF (i-j*48>0 .AND. i-j*48<=maxind) THEN
						IF (record(i-j*48)%LE>-999) THEN
							icount = icount + 1
							LE_sum = LE_sum + record(i-j*48)%LE
							Rn_sum = Rn_sum + record(i-j*48)%Rn
						END IF
					END IF
				END DO

				IF (icount >0) THEN
					IF (timeinfo%tm_hour >=8 .AND. timeinfo%tm_hour <=17) THEN
						record(i)%LE = LE_sum/Rn_sum*record(i)%Rn
					ELSE
						record(i)%LE = LE_sum/icount
					END IF
					record(i)%LE_FLAG = 2
					WRITE(*,*) timeinfo, record(i)%LE
				END IF
			END IF
		END IF
	END DO

	timeinfo_start = tm(2009,4,1,0,0,0)
	timeinfo_end = tm(2010,1,1,0,0,0)

	CALL timegm(timeinfo_start, time_start)
	CALL timegm(timeinfo_end, time_end)

	interval = 60*60

	datasize = (time_end-time_start)/interval + 1

	ALLOCATE(avg_time(datasize+1))
	ALLOCATE(avg_H(datasize))
	ALLOCATE(avg_LE(datasize))

	avg_time(0:datasize) = 0
	avg_H(1:datasize) = -999
	avg_LE(1:datasize) = -999

	DO i = 0, datasize - 1
		avg_time(i) = time_start + i*interval
	END DO

	CALL average_data(record%time, record%H, 1, avg_time, avg_H, 0)
	CALL average_data(record%time, record%LE, 1, avg_time, avg_LE, 0)

	OPEN(UNIT=outputid,FILE='flux_gap_filled.dat',ACTION='write')

	DO i = 1, datasize-1
		CALL gmtime(avg_time(i), timeinfo)
		WRITE(outputid,"(I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':00',F10.2, F10.2") timeinfo%tm_year,timeinfo%tm_mon,timeinfo%tm_mday,timeinfo%tm_hour,timeinfo%tm_min,avg_H(i), avg_LE(i)
	END DO

	CLOSE(outputid)


	WRITE(*,*)"\n  done.\n"

1000 CONTINUE

END PROGRAM fill_gap
