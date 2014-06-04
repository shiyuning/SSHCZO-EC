!-------------------------------------------------------------------------------
! Program for processing Shale Hills flux tower eddy covariance data
! Yuning Shi, Penn State Meteorology
! Email: yshi@psu.edu
!-------------------------------------------------------------------------------

MODULE cospectratype

	REAL, PARAMETER :: bin_min = 0.001
	REAL, PARAMETER :: bin_max = 5.0
	INTEGER, PARAMETER :: bin_no = 23

	TYPE :: cospectradata
	 	REAL :: Co_H(bin_no), Co_LE(bin_no), Co_FC(bin_no)
		INTEGER :: members
	END TYPE cospectradata
	
END MODULE cospectratype

MODULE ECtype

	TYPE :: ECdata
		INTEGER :: record_no
		REAL :: Ux(20000), Uy(20000), Uz(20000), Ts(20000), co2(20000), h2o(20000), time(20000)
		INTEGER :: diag(20000)
		REAL :: h2o_hmp(20000), T_hmp(20000)
	END TYPE ECdata

END MODULE ECtype

MODULE fluxtype

	TYPE :: fluxdata
		INTEGER :: YEAR, DOY, HRMIN
		REAL :: DTIME, UST, TA, WD, WS, VWS, FC, H, LE, CO2, H2O
		INTEGER :: U_FLAG, W_FLAG, T_FLAG, H2O_FLAG, CO2_FLAG
		REAL :: instrument_par, spikes_par, amplitude_resolution_par, dropouts_par, extreme_dropouts_par, nonstationarity_par
		REAL :: higher_moment_par1, higher_moment_par2, discontinuities_par1, discontinuities_par2
	END TYPE fluxdata
END MODULE fluxtype

PROGRAM cospectrum

	USE ECtype
	USE cospectratype

	IMPLICIT NONE
 
	CHARACTER(100) :: buffer
	INTEGER :: error

	INTEGER :: year, month, day, hour, minute, ind, diag
	REAL :: second, record_time, Ux, Uy, Uz, Ts, co2, h2o, h2o_hmp, T_hmp
	INTEGER :: i,j
	INTEGER :: infid = 100, outfid = 200
	REAL :: unit_k(3) =  (/-0.1395917, 0.0250885, 0.9794629/), b0 =  -0.2013965
	
	TYPE(ECdata) :: record
	TYPE(cospectradata) :: avg_cospectra

	avg_cospectra%members = 0
	avg_cospectra%Co_H = 0
	avg_cospectra%Co_LE = 0
	avg_cospectra%Co_FC = 0

	record%record_no = 0

	WRITE(*,*)"\n  Start reading in file ... \n"
	OPEN(UNIT=infid, FILE='2009-05/2009-05.dat', ACTION='read', IOSTAT=error)
	IF (error/=0) GOTO 1000

	OPEN(UNIT=outfid, FILE='2009-05/cospectra.dat', ACTION='WRITE', IOSTAT=error)

	WRITE(*,*)"\n  File opened ...\n"

	READ(100,*,IOSTAT = error) buffer
	READ(100,*,IOSTAT = error) buffer
	READ(100,*,IOSTAT = error) buffer
	READ(100,*,IOSTAT = error) buffer

	DO WHILE (.true.)
		READ(100,*,IOSTAT = error) buffer, ind, Ux, Uy, Uz, Ts, co2, h2o, diag, h2o_hmp, T_hmp
		IF (error /= 0 ) THEN
			EXIT
		END IF

		READ(buffer(:4),*) year
		READ(buffer(6:7),*) month
		READ(buffer(9:10),*) day
		READ(buffer(12:13),*) hour
		READ(buffer(15:16),*) minute
		READ(buffer(18:),*) second

		record_time = hour + REAL(minute)/60 + second/60/60

		IF (record_time >=time_range(1) .AND. record_time<time_range(2)) THEN
			WRITE(*,*) "\n  Begin reading in data...\n"
			DO WHILE (.TRUE.)
				READ(100,*,IOSTAT = error) buffer, ind, Ux, Uy, Uz, Ts, co2, h2o, diag, h2o_hmp, T_hmp
				READ(buffer(:4),*) year
				READ(buffer(6:7),*) month
				READ(buffer(9:10),*) day
				READ(buffer(12:13),*) hour
				READ(buffer(15:16),*) minute
				READ(buffer(18:),*) second
				record_time = hour + REAL(minute)/60 + second/60/60
				IF (record_time>time_range(2)) THEN
					CALL cospectra(record, unit_k, b0, avg_cospectra)
					record%record_no = 0
					EXIT
				END IF

				IF ( h2o>0 .AND. co2>0 .AND. IBITS(diag,4,4)==0) THEN
					record%record_no = record%record_no +1
					record%time(record%record_no) = ((REAL(day)-1)*24.0 + REAL(hour) - (time_end - 0.5))*60.0*60.0 + REAL(minute)*60.0 + second 
!					record%ind(record%record_no) = ind
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
		END IF
	END DO

	CLOSE(infid)
	CLOSE(outfid)

	WRITE(*,*)"\n  done.\n"

1000 CONTINUE

END PROGRAM cospectrum

SUBROUTINE unit_vector_k(fid, unit_k, b0)

	!-------------------------------------------------------------------------------
	! Determine unit vector k of planar fit coordinate
	! (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 4)
	!-------------------------------------------------------------------------------

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

	su = 0
	sv = 0
	sw = 0
	suv = 0
	suw = 0
	svw = 0
	su2 = 0
	sv2 = 0
	flen = 0

	READ(fid,*,IOSTAT = error) buffer
	READ(fid,*,IOSTAT = error) buffer
	READ(fid,*,IOSTAT = error) buffer
	READ(fid,*,IOSTAT = error) buffer

	DO WHILE (.TRUE.)
!	DO i = 1,100
		READ(fid,*,iostat = error) buffer, ind, Ux, Uy, Uz, Ts, co2, h2o, diag, h2o_hmp, T_hmp
		IF (error /= 0 ) EXIT
		su = su + Ux
		sv = sv + Uy
		sw = sw + Uz
		suv = suv + Ux*Uy
		suw = suw + Ux*Uz
		svw = svw + Uy*Uz
		su2 = su2 + Ux*Ux
		sv2 = sv2 + Uy*Uy
		flen = flen +1
	END DO

	H1(1,:) = (/flen,su,sv/)
	H1(2,:) = (/su, su2, suv/)
	H1(3,:) = (/sv, suv, sv2/)

	g = (/sw, suw, svw/)

	CALL FindInv(H1, invH1, 3, ErrorFlag)

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

SUBROUTINE flux(record, unit_k, b0, flux_record, outfid)

	USE ECtype
	USE fluxtype
	USE QC

	IMPLICIT NONE

	TYPE(fluxdata) :: flux_record

	TYPE(ECdata) :: record
	TYPE(ECdata) :: perturb

	REAL, PARAMETER :: rho = 1.2
	REAL, PARAMETER :: Lv = 2503000.0
	REAL, PARAMETER :: c_air = 1004.0
	REAL, PARAMETER :: Pi = 3.1415927
	REAL, PARAMETER :: CSAT3_AZIMUTH = 199.0
	INTEGER :: outfid
	REAL :: u_bar, v_bar, w_bar, h2o_bar, T_bar, eta, co2_bar
	REAL :: wT_bar, wh2o_bar, uT_bar, uh2o_bar, vT_bar, vh2o_bar, uco2_bar, vco2_bar, wco2_bar
	REAL :: uv_bar, uw_bar, vw_bar
	REAL :: uw(record%record_no), vw(record%record_no)
!	REAL :: u2(record%record_no), v2(record%record_no), w2(record%record_no)
	REAL :: U_trend(record%record_no), V_trend(record%record_no), W_trend(record%record_no), CO2_trend(record%record_no), H2O_trend(record%record_no), T_trend(record%record_no)
	REAL :: U1(3)
	REAL :: b0
	REAL :: CE, SE
	INTEGER :: i,j
	INTEGER :: INSTRUMENT_FLAG, NONSTATIONARITY_FLAG
	INTEGER :: SPIKES_FLAG(5), AMPLITUDE_RESOLUTION_FLAG(5), DROPOUTS_FLAG(5), ABSOLUTE_LIMITS_FLAG(5), HIGHER_MOMENT_FLAG(5), DISCONTINUITIES_FLAG(5)
	INTEGER :: ErrorFlag
	REAL :: unit_k(3), unit_j(3), unit_i(3)

	perturb%record_no = record%record_no
	perturb%time = record%time

	SPIKES_FLAG = 0
	ABSOLUTE_LIMITS_FLAG = 0
	NONSTATIONARITY_FLAG = 0
	AMPLITUDE_RESOLUTION_FLAG = 0
	DROPOUTS_FLAG = 0
	HIGHER_MOMENT_FLAG = 0
	DISCONTINUITIES_FLAG = 0
	INSTRUMENT_FLAG = 0

	WRITE(*,"('  ',I4,'-',I3.3,'-',I4.4, I)") flux_record%YEAR, flux_record%DOY, flux_record%HRMIN, record%record_no

!	flux_record%YEAR = record(record_no)%year
!	CALL doy(record(record_no)%year, record(record_no)%month, record(record_no)%day, flux_record%DOY)
!	flux_record%HRMIN = record(record_no)%hour*100+record(record_no)%minute
!	flux_record%DTIME = (REAL(record(record_no)%hour) + REAL(record(record_no)%minute)/60)/24

	CALL instrument(record, INSTRUMENT_FLAG, flux_record%instrument_par)

	CALL spikes(record, SPIKES_FLAG, flux_record%spikes_par)

	WRITE(*,"('  Ux range = [',F7.2,',',F7.2,']')") MINVAL(ABS(record%Ux(1:record%record_no))), MAXVAL(ABS(record%Ux(1:record%record_no)))
	WRITE(*,"('  Uy range = [',F7.2,',',F7.2,']')") MINVAL(ABS(record%Uy(1:record%record_no))), MAXVAL(ABS(record%Uy(1:record%record_no)))
	WRITE(*,"('  Uz range = [',F7.2,',',F7.2,']')") MINVAL(ABS(record%Uz(1:record%record_no))), MAXVAL(ABS(record%Uz(1:record%record_no)))
	WRITE(*,"('  Ts range = [',F7.2,',',F7.2,']')") MINVAL(record%Ts(1:record%record_no)), MAXVAL(record%Ts(1:record%record_no))
	WRITE(*,"('  CO2 range = [',F7.2,',',F7.2,']')") MINVAL(record%co2(1:record%record_no)), MAXVAL(record%co2(1:record%record_no))
	WRITE(*,"('  H2O range = [',F7.2,',',F7.2,']')") MINVAL(record%h2o(1:record%record_no)), MAXVAL(record%h2o(1:record%record_no))


	CALL amplitude_resolution(record, AMPLITUDE_RESOLUTION_FLAG, flux_record%amplitude_resolution_par)

	CALL dropouts(record, DROPOUTS_FLAG, flux_record%dropouts_par, flux_record%extreme_dropouts_par)

	CALL absolute_limits(record, ABSOLUTE_LIMITS_FLAG)

!	DO i = 1,record_no
!		time(i) = real(record(i)%minute)*60+record(i)%second
!		IF (i > 1 .AND. time(i)<time(i-1)) time(i) = time(i) + 60*60
!	END DO

 
	!-------------------------------------------------------------------------------
	! Transformation to the natural wind coordinate system
	! (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 4)
	!-------------------------------------------------------------------------------

	record%Uz(1:record%record_no) = record%Uz(1:record%record_no) - b0
	u_bar = SUM(record%Ux(1:record%record_no))/REAL(record%record_no)
	v_bar = SUM(record%Uy(1:record%record_no))/REAL(record%record_no)
	w_bar = SUM(record%Uz(1:record%record_no))/REAL(record%record_no)

!	WRITE(*,*) SUM(record%Ux(1:record%record_no)), record%record_no


	h2o_bar = SUM(record%h2o(1:record%record_no))/REAL(record%record_no)
	co2_bar = SUM(record%co2(1:record%record_no))/REAL(record%record_no)
	T_bar = SUM(record%Ts(1:record%record_no))/REAL(record%record_no)

	CE = u_bar/SQRT(u_bar*u_bar+v_bar*v_bar)
	SE = v_bar/SQRT(u_bar*u_bar+v_bar*v_bar)

	eta = ACOS(CE)/Pi*180
	IF (SE<0) eta = 360 - eta
	eta = MODULO(360.0 - eta + CSAT3_AZIMUTH, 360.0)

	U1 = (/u_bar,v_bar,w_bar/)

!	WRITE(*,*) U1

	CALL cross(unit_k,U1,unit_j)
	unit_j = unit_j/SQRT(SUM(unit_j*unit_j))
	CALL cross(unit_j,unit_k,unit_i)

!	WRITE(*,"('i_vector = [',F,',',F,',',F,']')") unit_i(1),unit_i(2),unit_i(3)

	CALL detrend(record%time,record%Ux, U_trend, perturb%Ux, record%record_no)
	CALL detrend(record%time,record%Uy, V_trend, perturb%Uy, record%record_no)
	CALL detrend(record%time,record%Uz, W_trend, perturb%Uz, record%record_no)
	CALL detrend(record%time,record%co2, CO2_trend, perturb%co2, record%record_no)
	CALL detrend(record%time,record%h2o, H2O_trend, perturb%h2o, record%record_no)
	CALL detrend(record%time,record%Ts, T_trend, perturb%Ts, record%record_no)

	CALL higher_moment(perturb,HIGHER_MOMENT_FLAG, flux_record%higher_moment_par1, flux_record%higher_moment_par2)

	CALL discontinuities(record, DISCONTINUITIES_FLAG, flux_record%discontinuities_par1, flux_record%discontinuities_par2)

	CALL nonstationary(record%Ux, record%Uy, record%record_no, NONSTATIONARITY_FLAG, flux_record%nonstationarity_par)

!	uw = unit_i(1)*unit_k(1)*(U_trend*U_trend + 2*U_trend*u_perturb + u_perturb*u_perturb)
!	uw = uw + unit_i(2)*unit_k(2)*(V_trend*V_trend + 2*V_trend*v_perturb + v_perturb*v_perturb)
!	uw = uw + (unit_i(1)*unit_k(2)+unit_i(2)*unit_k(1))*(U_trend*V_trend + U_trend*v_perturb + u_perturb*V_trend + u_perturb*v_perturb)
!	uw = uw + (unit_i(1)*unit_k(3)+unit_i(3)*unit_k(1))*(U_trend*W_trend + U_trend*w_perturb + u_perturb*W_trend + u_perturb*w_perturb)
!	uw = uw + (unit_i(2)*unit_k(3)+unit_i(3)*unit_k(2))*(V_trend*W_trend + V_trend*w_perturb + v_perturb*W_trend + v_perturb*w_perturb)

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

	flux_record%TA = T_bar
	flux_record%WD = eta
	flux_record%WS = SQRT(SUM(unit_i*U1)*SUM(unit_i*U1)+SUM(unit_j*U1)*SUM(unit_j*U1))
	flux_record%VWS = w_bar
	flux_record%FC = 1000/44*(uco2_bar*unit_k(1)+vco2_bar*unit_k(2)+wco2_bar*unit_k(3))
	flux_record%H = rho*c_air*(uT_bar*unit_k(1)+vT_bar*unit_k(2)+wT_bar*unit_k(3))
	flux_record%LE = rho*Lv/1.2/1000*(uh2o_bar*unit_k(1)+vh2o_bar*unit_k(2)+wh2o_bar*unit_k(3))
	flux_record%CO2 = co2_bar
	flux_record%H2O = h2o_bar
	flux_record%U_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(1)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(1)*(2**2)+DROPOUTS_FLAG(1)*(2**3)+ABSOLUTE_LIMITS_FLAG(1)*(2**4)+HIGHER_MOMENT_FLAG(1)*(2**5)+DISCONTINUITIES_FLAG(1)*(2**6)+ NONSTATIONARITY_FLAG*(2**7)

	flux_record%W_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(2)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(2)*(2**2)+DROPOUTS_FLAG(2)*(2**3)+ABSOLUTE_LIMITS_FLAG(2)*(2**4)+HIGHER_MOMENT_FLAG(2)*(2**5)+DISCONTINUITIES_FLAG(2)*(2**6)+ NONSTATIONARITY_FLAG*(2**7)


	flux_record%T_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(3)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(3)*(2**2)+DROPOUTS_FLAG(3)*(2**3)+ABSOLUTE_LIMITS_FLAG(3)*(2**4)+HIGHER_MOMENT_FLAG(3)*(2**5)+DISCONTINUITIES_FLAG(3)*(2**6)+ NONSTATIONARITY_FLAG*(2**7)

	flux_record%H2O_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(4)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(4)*(2**2)+DROPOUTS_FLAG(4)*(2**3)+ABSOLUTE_LIMITS_FLAG(4)*(2**4)+HIGHER_MOMENT_FLAG(4)*(2**5)+DISCONTINUITIES_FLAG(4)*(2**6)+ NONSTATIONARITY_FLAG*(2**7)

	flux_record%CO2_FLAG = INSTRUMENT_FLAG + SPIKES_FLAG(5)*(2**1)+AMPLITUDE_RESOLUTION_FLAG(5)*(2**2)+DROPOUTS_FLAG(5)*(2**3)+ABSOLUTE_LIMITS_FLAG(5)*(2**4)+HIGHER_MOMENT_FLAG(5)*(2**5)+DISCONTINUITIES_FLAG(5)*(2**6)+ NONSTATIONARITY_FLAG*(2**7)


	CALL write_to_file(outfid, flux_record)

	! Flags:
	! Bit 1:
	! Bit 2: Spikes (2)
	! Bit 3: Amplitude resolution (4)
	! Bit 4: Dropouts (8)
	! Bit 5: Absolute limits (16)
	! Bit 6: Higher-moment statistics (32)
	! Bit 7: Discontinuities (64)
	! Bit 8: Nonstationarity of the horizontal wind (128)

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
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL :: m
	REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
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
	
	!Reduce augmented matrix to upper traingular form
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
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END SUBROUTINE FINDinv

SUBROUTINE doy(year,month,day, jday)

	IMPLICIT NONE
	INTEGER :: year,month,day, jday
	INTEGER :: eom(12)
	INTEGER :: daysum,i

	IF ((mod(year,4)==0 .AND. mod(year,100)/=0) .OR. mod(year,400)==0) THEN
		eom=(/31,29,31,30,31,30,31,31,30,31,30,31/)
	ELSE
		eom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
	END IF

	IF (month==1) THEN
	    jday=day
	ELSE
	    daysum=0
		DO i = 1, month-1
			daysum=daysum+eom(i)
		END DO
		jday=daysum+day
	END IF

END SUBROUTINE doy

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

SUBROUTINE write_to_file(fid, flux_record)

	USE fluxtype

	IMPLICIT NONE
	INTEGER :: fid
	CHARACTER(15) :: year_char, doy_char, hrmin_char, dtime_char
	CHARACTER(15) :: ust_char, t_bar_char, eta_char, WS_char, FC_char, H_char, LE_char, co2_bar_char, h2o_bar_char, vws_char
	CHARACTER(15) :: u_flag_char, w_flag_char, t_flag_char, h2o_flag_char, co2_flag_char

	TYPE(fluxdata) :: flux_record

	WRITE(year_char, "(I)") flux_record%YEAR
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

	WRITE(fid,"(A)") trim(adjustl(year_char))//","//trim(adjustl(doy_char))//","//trim(adjustl(hrmin_char))//","//trim(adjustl(dtime_char))//","//&
		trim(adjustl(ust_char))//","//trim(adjustl(t_bar_char))//","//trim(adjustl(eta_char))//","//trim(adjustl(WS_char))//","//trim(adjustl(vws_char))//","//&
		trim(adjustl(FC_char))//","//trim(adjustl(H_char))//","//trim(adjustl(LE_char))//","//trim(adjustl(co2_bar_char))//","//trim(adjustl(h2o_bar_char))//","//&
		trim(adjustl(u_flag_char))//","//trim(adjustl(w_flag_char))//","//trim(adjustl(t_flag_char))//","//trim(adjustl(h2o_flag_char))//","//&
		trim(adjustl(co2_flag_char))

	WRITE(*,"(A)") trim(adjustl(year_char))//","//trim(adjustl(doy_char))//","//trim(adjustl(hrmin_char))//","//trim(adjustl(dtime_char))//","//&
		trim(adjustl(ust_char))//","//trim(adjustl(t_bar_char))//","//trim(adjustl(eta_char))//","//trim(adjustl(WS_char))//","//trim(adjustl(vws_char))//","//&
		trim(adjustl(FC_char))//","//trim(adjustl(H_char))//","//trim(adjustl(LE_char))//","//trim(adjustl(co2_bar_char))//","//trim(adjustl(h2o_bar_char))//","//&
		trim(adjustl(u_flag_char))//","//trim(adjustl(w_flag_char))//","//trim(adjustl(t_flag_char))//","//trim(adjustl(h2o_flag_char))//","//&
		trim(adjustl(co2_flag_char))
	WRITE(*,*) flux_record%instrument_par, flux_record%spikes_par, flux_record%amplitude_resolution_par, flux_record%dropouts_par, flux_record%extreme_dropouts_par, &
		flux_record%nonstationarity_par, flux_record%higher_moment_par1, flux_record%higher_moment_par2, flux_record%discontinuities_par1, &
		flux_record%discontinuities_par2

END SUBROUTINE write_to_file

!SUBROUTINE write_to_file1(fid, flux_record)

!	USE fluxtype

!	IMPLICIT NONE
!	INTEGER :: fid
!	CHARACTER(15) :: year_char, doy_char, hrmin_char, dtime_char
!	CHARACTER(15) :: ust_char, t_bar_char, eta_char, WS_char, FC_char, H_char, LE_char, co2_bar_char, h2o_bar_char, vws_char, flag_char
!	CHARACTER(15) :: instrument_par_char, spikes_par_char, amplitude_resolution_par_char, dropouts_par_char, extreme_dropouts_par_char, nonstationarity_par_char
!	CHARACTER(15) :: higher_moment_par1_char, higher_moment_par2_char, discontinuities_par1_char, discontinuities_par2_char

!	TYPE(fluxdata) :: flux_record

!	WRITE(year_char, "(I)") flux_record%YEAR
!	WRITE(doy_char, "(I)") flux_record%DOY
!	WRITE(hrmin_char, "(I4.4)") flux_record%HRMIN
!	WRITE(dtime_char, "(F10.5)") flux_record%DTIME
!	WRITE(ust_char, "(F10.3)") flux_record%UST
!	WRITE(t_bar_char, "(F10.2)") flux_record%TA
!	WRITE(eta_char, "(F10.3)") flux_record%WD
!	WRITE(WS_char, "(F10.3)") flux_record%WS
!	WRITE(FC_char, "(F10.3)") flux_record%FC
!	WRITE(H_char, "(F10.3)") flux_record%H
!	WRITE(LE_char, "(F10.3)") flux_record%LE
!	WRITE(co2_bar_char, "(F10.3)") flux_record%CO2
!	WRITE(h2o_bar_char, "(F10.3)") flux_record%H2O
!	WRITE(vws_char, "(F10.3)") flux_record%VWS
!	WRITE(flag_char, "(I)") flux_record%FLAG

!	WRITE(instrument_par_char, "(F10.3)") flux_record%instrument_par
!	WRITE(spikes_par_char, "(F10.3)") flux_record%spikes_par
!	WRITE(amplitude_resolution_par_char, "(F10.3)") flux_record%amplitude_resolution_par
!	WRITE(dropouts_par_char, "(F10.3)") flux_record%dropouts_par
!	WRITE(extreme_dropouts_par_char, "(F10.3)") flux_record%extreme_dropouts_par
!	WRITE(nonstationarity_par_char, "(F10.3)") flux_record%nonstationarity_par
!	WRITE(higher_moment_par1_char, "(F10.3)") flux_record%higher_moment_par1
!	WRITE(higher_moment_par2_char, "(F10.3)") flux_record%higher_moment_par2
!	WRITE(discontinuities_par1_char, "(F10.3)") flux_record%discontinuities_par1
!	WRITE(discontinuities_par2_char, "(F10.3)") flux_record%discontinuities_par2

!	WRITE(fid,"(A)") trim(adjustl(year_char))//","//trim(adjustl(doy_char))//","//trim(adjustl(hrmin_char))//","//trim(adjustl(dtime_char))//","//&
!		trim(adjustl(ust_char))//","//trim(adjustl(t_bar_char))//","//trim(adjustl(eta_char))//","//trim(adjustl(WS_char))//","//trim(adjustl(vws_char))//","//&
!		trim(adjustl(FC_char))//","//trim(adjustl(H_char))//","//trim(adjustl(LE_char))//","//trim(adjustl(co2_bar_char))//","//trim(adjustl(h2o_bar_char))//","//&
!		trim(adjustl(flag_char))//","//&
!		trim(adjustl(instrument_par_char))//","//trim(adjustl(spikes_par_char))//","//trim(adjustl(amplitude_resolution_par_char))//","//&
!		trim(adjustl(dropouts_par_char))//","//trim(adjustl(extreme_dropouts_par_char))//","//trim(adjustl(nonstationarity_par_char))//","//&
!		trim(adjustl(higher_moment_par1_char))//","//trim(adjustl(higher_moment_par2_char))//","//trim(adjustl(discontinuities_par1_char))//","//&
!		trim(adjustl(discontinuities_par2_char))
!	WRITE(*,"(A)") trim(adjustl(year_char))//","//trim(adjustl(doy_char))//","//trim(adjustl(hrmin_char))//","//trim(adjustl(dtime_char))//","//&
!		trim(adjustl(ust_char))//","//trim(adjustl(t_bar_char))//","//trim(adjustl(eta_char))//","//trim(adjustl(WS_char))//","//trim(adjustl(vws_char))//","//&
!		trim(adjustl(FC_char))//","//trim(adjustl(H_char))//","//trim(adjustl(LE_char))//","//trim(adjustl(co2_bar_char))//","//trim(adjustl(h2o_bar_char))//","//&
!		trim(adjustl(flag_char))//","//&
!		trim(adjustl(instrument_par_char))//","//trim(adjustl(spikes_par_char))//","//trim(adjustl(amplitude_resolution_par_char))//","//&
!		trim(adjustl(dropouts_par_char))//","//trim(adjustl(extreme_dropouts_par_char))//","//trim(adjustl(nonstationarity_par_char))//","//&
!		trim(adjustl(higher_moment_par1_char))//","//trim(adjustl(higher_moment_par2_char))//","//trim(adjustl(discontinuities_par1_char))//","//&
!		trim(adjustl(discontinuities_par2_char))

!END SUBROUTINE write_to_file1
