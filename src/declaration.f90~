MODULE ECtype

!-------------------------------------------------------------------------------
! Define eddy-covariance record data type
!-------------------------------------------------------------------------------

	TYPE :: ECdata
		INTEGER :: record_no
		REAL :: Ux(20000), Uy(20000), Uz(20000), Ts(20000), co2(20000), h2o(20000), time(20000)
		INTEGER :: diag(20000)
		REAL :: h2o_hmp(20000), T_hmp(20000)
	END TYPE ECdata

END MODULE ECtype

MODULE fluxtype

!-------------------------------------------------------------------------------
! Define surface flux data type
!-------------------------------------------------------------------------------

	TYPE :: fluxdata
		INTEGER :: YEAR, MONTH, MDAY, DOY, HRMIN
		REAL :: DTIME, UST, TA, WD, WS, VWS, FC, H, LE, CO2, H2O
		INTEGER :: U_FLAG, W_FLAG, T_FLAG, H2O_FLAG, CO2_FLAG
		REAL :: instrument_par, nonstationarity_par, rn_par
		REAL :: U_spikes_par, U_amplitude_resolution_par, U_dropouts_par, U_extreme_dropouts_par
		REAL :: U_skewness_par, U_kurtosis_par, U_haar_mean_par, U_haar_variance_par
		REAL :: W_spikes_par, W_amplitude_resolution_par, W_dropouts_par, W_extreme_dropouts_par
		REAL :: W_skewness_par, W_kurtosis_par, W_haar_mean_par, W_haar_variance_par
		REAL :: T_spikes_par, T_amplitude_resolution_par, T_dropouts_par, T_extreme_dropouts_par
		REAL :: T_skewness_par, T_kurtosis_par, T_haar_mean_par, T_haar_variance_par
		REAL :: H2O_spikes_par, H2O_amplitude_resolution_par, H2O_dropouts_par, H2O_extreme_dropouts_par
		REAL :: H2O_skewness_par, H2O_kurtosis_par, H2O_haar_mean_par, H2O_haar_variance_par
		REAL :: CO2_spikes_par, CO2_amplitude_resolution_par, CO2_dropouts_par, CO2_extreme_dropouts_par
		REAL :: CO2_skewness_par, CO2_kurtosis_par, CO2_haar_mean_par, CO2_haar_variance_par
	END TYPE fluxdata

END MODULE fluxtype

MODULE Pressure

!-------------------------------------------------------------------------------
! Define pressure record data type
!-------------------------------------------------------------------------------

	TYPE :: Pressuredata
		INTEGER :: time(100000)
		REAL :: pressure(100000)
	END TYPE Pressuredata

END MODULE Pressure
