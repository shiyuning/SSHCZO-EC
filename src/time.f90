MODULE time

	IMPLICIT NONE

	TYPE :: tm
		INTEGER :: tm_year, tm_mon, tm_mday, tm_hour, tm_min, tm_sec
	END TYPE tm

	INTEGER, PARAMETER :: year_org = 70
	INTEGER, PARAMETER :: mon_org = 1
	INTEGER, PARAMETER :: mday_org = 1
	INTEGER, PARAMETER :: hour_org = 0
	INTEGER, PARAMETER :: min_org = 0
	INTEGER, PARAMETER :: sec_org = 0

CONTAINS

	SUBROUTINE timegm(timeinfo, rawtime)

		IMPLICIT NONE

		TYPE(tm) :: timeinfo

		INTEGER :: rawtime
		INTEGER :: eom(12)
		INTEGER :: year, month, mday

		rawtime = 0

		timeinfo%tm_year = timeinfo%tm_year - 1900
		timeinfo%tm_mon = timeinfo%tm_mon

		IF ((mod(timeinfo%tm_year + 1900,4)==0 .AND. mod(timeinfo%tm_year + 1900,100)/=0) .OR. mod(timeinfo%tm_year + 1900,400)==0) THEN
			eom=(/31,29,31,30,31,30,31,31,30,31,30,31/)
		ELSE
			eom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
		END IF

		!----------------------------------------------------!
		! Check input format                                 !
		!----------------------------------------------------!

		IF (timeinfo%tm_mon>12 .OR. timeinfo%tm_mon<1) THEN
			WRITE(*,*) "\nERROR: Month out of range!!!"
			STOP
		END IF

		IF (timeinfo%tm_mday>eom(timeinfo%tm_mon) .OR. timeinfo%tm_mday<1) THEN
			WRITE(*,*) "\nERROR: Day out of range!!!"
			STOP
		END IF

		IF (timeinfo%tm_hour>23 .OR. timeinfo%tm_hour<0) THEN
			WRITE(*,*) "\nERROR: Hour out of range!!!"
			STOP
		END IF

		IF (timeinfo%tm_min>59 .OR. timeinfo%tm_min<0) THEN
			WRITE(*,*) "\nERROR: Minute out of range!!!"
			STOP
		END IF

		IF (timeinfo%tm_sec>59 .OR. timeinfo%tm_sec<0) THEN
			WRITE(*,*) "\nERROR: Second out of range!!!"
			STOP
		END IF

		IF (timeinfo%tm_year >= year_org) THEN
			IF (timeinfo%tm_year > year_org) THEN
				DO year = year_org, timeinfo%tm_year - 1
					IF ((mod(year + 1900,4)==0 .AND. mod(year + 1900,100)/=0) .OR. mod(year + 1900,400)==0) THEN
						rawtime = rawtime + 366*24*60*60
					ELSE
						rawtime = rawtime + 365*24*60*60
					END IF
				END DO
			END IF

			IF (timeinfo%tm_mon > mon_org) THEN
				DO month = 1, timeinfo%tm_mon - 1
					rawtime = rawtime + eom(month)*24*60*60
				END DO
			END IF

			rawtime = rawtime + (timeinfo%tm_mday - mday_org)*24*60*60 + timeinfo%tm_hour*60*60 + timeinfo%tm_min*60 +timeinfo%tm_sec
		END IF

		IF (timeinfo%tm_year < year_org) THEN
			DO year = timeinfo%tm_year, year_org - 1
				IF ((mod(year + 1900,4)==0 .AND. mod(year + 1900,100)/=0) .OR. mod(year + 1900,400)==0) THEN
					rawtime = rawtime - 366*24*60*60
				ELSE
					rawtime = rawtime - 365*24*60*60
				END IF
			END DO

			IF (timeinfo%tm_mon > mon_org) THEN
				DO month = 1, timeinfo%tm_mon - 1
					rawtime = rawtime + eom(month)*24*60*60
				END DO
			END IF
			rawtime = rawtime + (timeinfo%tm_mday - mday_org)*24*60*60 + timeinfo%tm_hour*60*60 + timeinfo%tm_min*60 +timeinfo%tm_sec
		END IF

		timeinfo%tm_year = timeinfo%tm_year + 1900
	
	END SUBROUTINE timegm

	SUBROUTINE gmtime(rawtime, timeinfo)
	
		IMPLICIT NONE

		TYPE(tm) :: timeinfo

		INTEGER :: rawtime, time, time_temp
		INTEGER :: eom(12)
		INTEGER :: year, month, mday

		time_temp = 0
		time = 0

		year = year_org
		month = 1
		mday = 1

		IF (rawtime>=0) THEN
			DO WHILE(.TRUE.)
				IF ((mod(year + 1900,4)==0 .AND. mod(year + 1900,100)/=0) .OR. mod(year + 1900,400)==0) THEN
					time_temp = time + 366*24*60*60
				ELSE
					time_temp = time + 365*24*60*60
				END IF

				IF(time_temp>rawtime) EXIT

				year = year + 1
				time = time_temp
			END DO

			timeinfo%tm_year = year


			IF ((mod(timeinfo%tm_year + 1900,4)==0 .AND. mod(timeinfo%tm_year + 1900,100)/=0) .OR. mod(timeinfo%tm_year + 1900,400)==0) THEN
				eom=(/31,29,31,30,31,30,31,31,30,31,30,31/)
			ELSE
				eom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
			END IF

			DO WHILE(.TRUE.)
		
				time_temp = time + eom(month)*24*60*60
				IF(time_temp>rawtime) EXIT
				month = month + 1
				time = time_temp
			END DO

			timeinfo%tm_mon = month
			timeinfo%tm_mday = FLOOR(REAL(rawtime - time)/24/60/60)+1
			time = time + (timeinfo%tm_mday-1)*24*60*60
			timeinfo%tm_hour = INT((rawtime - time)/60/60)
			time = time + timeinfo%tm_hour*60*60
			timeinfo%tm_min = INT((rawtime - time)/60)
			time = time + timeinfo%tm_min*60
			timeinfo%tm_sec = rawtime - time
			timeinfo%tm_year = timeinfo%tm_year + 1900
		END IF

		IF (rawtime<0) THEN
			year = year - 1
			DO WHILE(.TRUE.)
				IF ((mod(year + 1900,4)==0 .AND. mod(year + 1900,100)/=0) .OR. mod(year + 1900,400)==0) THEN
					time = time - 366*24*60*60
				ELSE
					time = time - 365*24*60*60
				END IF

				IF(time<=rawtime) EXIT

				year = year - 1
			END DO

			timeinfo%tm_year = year


			IF ((mod(timeinfo%tm_year + 1900,4)==0 .AND. mod(timeinfo%tm_year + 1900,100)/=0) .OR. mod(timeinfo%tm_year + 1900,400)==0) THEN
				eom=(/31,29,31,30,31,30,31,31,30,31,30,31/)
			ELSE
				eom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
			END IF

			DO WHILE(.TRUE.)
		
				time_temp = time + eom(month)*24*60*60
				IF(time_temp>rawtime) EXIT
				month = month + 1
				time = time_temp
			END DO

			timeinfo%tm_mon = month
			timeinfo%tm_mday = FLOOR(REAL(rawtime - time)/24/60/60)+1
			time = time + (timeinfo%tm_mday-1)*24*60*60
			timeinfo%tm_hour = INT((rawtime - time)/60/60)
			time = time + timeinfo%tm_hour*60*60
			timeinfo%tm_min = INT((rawtime - time)/60)
			time = time + timeinfo%tm_min*60
			timeinfo%tm_sec = rawtime - time
			timeinfo%tm_year = timeinfo%tm_year + 1900
		END IF
	END SUBROUTINE gmtime

	SUBROUTINE doy(year,month,day,jday)

		IMPLICIT NONE
		INTEGER :: year,month,day,jday
		INTEGER :: eom(12)
		INTEGER :: daysum,i

		daysum = 0

		IF ((mod(year,4)==0 .AND. mod(year,100)/=0) .OR. mod(year,400)==0) THEN
			eom=(/31,29,31,30,31,30,31,31,30,31,30,31/)
		ELSE
			eom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
		END IF

		IF (month==1) THEN
			jday=day
		ELSE
			daysum=eom(1)
		END IF

		DO i=2,month-1
			daysum=daysum+eom(i)
		END DO
		jday=daysum+day
	END SUBROUTINE doy

	SUBROUTINE doy2date(year,jday,month,day)

		IMPLICIT NONE
		INTEGER :: daysum, year, month, day, jday
		INTEGER :: i
		INTEGER :: eom(12)

		daysum = 0
		i = 1

		IF ((mod(year,4)==0 .AND. mod(year,100)/=0) .OR. mod(year,400)==0) THEN
			eom=(/31,29,31,30,31,30,31,31,30,31,30,31/)
		ELSE
			eom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
		END IF

		DO WHILE (i<=12)
			IF (daysum+eom(i)>=jday) THEN
				month = i
			        day=jday-daysum
				EXIT
			ELSE
				daysum=daysum+eom(i)
				i = i+1
			END IF
		END DO

	END SUBROUTINE doy2date

END MODULE time
