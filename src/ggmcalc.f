 !-------------------------------------------------------------------------------------------------
 !	GGMCalc, A program to compute some components of gravity field using Global Geopotential Models
 !
 !	Copyright (C) 2010-2018 Siamak Moazezi
 !	
 !	GGMCalc is free software: you can redistribute it and/or modify
 !	it under the terms of the GNU General Public License as published by
 !	the Free Software Foundation, either version 3 of the License, or
 !	(at your option) any later version.
 !
 !	GGMCalc is distributed in the hope that it will be useful,
 !	but WITHOUT ANY WARRANTY; without even the implied warranty of
 !	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 !	GNU General Public License for more details.
 !
 !	You should have received a copy of the GNU General Public License
 !	along with GGMCalc.  If not, see <http://www.gnu.org/licenses/>.
 !
 !	Contact info: http://www.sourceforge.net/projects/xgravity
 !	              Siamak Moazezi <s.moazezi@srbiau.ac.ir>
 !-------------------------------------------------------------------------------------------------

 use nrtype
 use Legendre_mod
 use C_2n_mod
 use gamma_mod
 use W_mod
 use U_mod
 use undulation_mod
 use height_anomaly_mod
 use gravity_disturbance_mod
 use gravity_anomaly_mod
 use ProgressBar_mod
 use duration
 implicit none
 integer(i4b) :: nmin=0, nmax, nmax_p, n, m, ii, jj, nn
 real(longdp) :: GM_ellipsoid, a, b, f_reciprocal, omega, U_0, m_ellipsoid, gamma_0
 real(longdp) :: GM, R, c, s, sigma_c, sigma_s
 real(longdp), allocatable, dimension(:) :: C_bar_nm
 real(longdp), allocatable, dimension(:,:) :: p_bar_nm, dp_bar_nm
 integer(i4b) :: Iteration
 real(longdp) :: N_Precision, zeta_Precision, Latitude, Longitude, h, depth, ice, rho, N_O, N_E
 real(longdp), allocatable, dimension(:) :: Latitude_p, Latitude_p_geodetic, Longitude_p, h_p, Undulation_p, Undulation_Simple_Assumed_Density_p, Undulation_Assumed_Density_p, Height_Anomaly_p, Gravity_Disturbance_p, Classical_Gravity_Anomaly_p, Classical_Gravity_Anomaly_2_p, Molodensky_Gravity_Anomaly_p
 real(longdp), allocatable, dimension(:) :: depth_p, ice_p
 character*(1) :: temp
 character*(50) str1, str2, product_type, modelname, errors, norm, tide_system, coefficient_kind, ellipsoidname
 character*(256) buffer, label, productionname, outputname, outputlogname
 integer(i4b) :: pos, Number_p
 logical ellipsoidal_true, T_B_IT_true, writing

 integer :: days1, hours1, minutes1, days2, hours2, minutes2
 real(longdp) :: seconds1, seconds2

! nproc = omp_get_num_procs()
! write(*,*) ' Number of available threads ',nproc
! write(*,*) " Insert number of threads "
! read(*,*) nproc
! write(*,*)' Threads used ', nproc
! call omp_set_num_threads(nproc)

 !----------------------------------------------------------------------------------------------------
 !---------------------------- Reading of Reference Ellipsoid Parameters -----------------------------
 !----------------------------------------------------------------------------------------------------
 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters . . . '

 open(12, file='ellipsoid.dat')
 !----------------------------------------------------------------------------------------------------
 read(12, *) GM_ellipsoid
 read(12, *) a
 read(12, *) b
 read(12, *) f_reciprocal
 read(12, *) omega
 read(12, *) U_0
 read(12, *) ellipsoidname
 !----------------------------------------------------------------------------------------------------
 close(12)
 !----------------------------------------------------------------------------------------------------
 m_ellipsoid = (omega**2.0_longdp * a**2.0_longdp *b) / GM_ellipsoid
 gamma_0 = GM_ellipsoid / a**2.0_longdp
 !----------------------------------------------------------------------------------------------------
 !------------------------- End of Reading of Reference Ellipsoid Parameters -------------------------
 !----------------------------------------------------------------------------------------------------


 !----------------------------------------------------------------------------------------------------
 !----------------------------- Reading of Geopotential Model Parameters -----------------------------
 !----------------------------------------------------------------------------------------------------

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters . . . '

 product_type = 'Unknown'
 modelname = 'Unknown'
 errors = 'Unknown'
 norm = 'Unknown'
 tide_system = 'Unknown'
 open(12, file='coeff.dat')
 do
	read(12, '(A)') buffer
	pos = scan(buffer, ' ')
	if (pos /= 0) then
		label = buffer(1:pos)
		buffer = adjustl(buffer(pos+1:))
		select case (label)
		case ('product_type')
			read(buffer, *) str2
			product_type = trim(str2)
		case ('modelname')
			read(buffer, *) str2
			modelname = trim(str2)
		case ('earth_gravity_constant')
			read(buffer, *) str2
			read(str2, *) GM
		case ('radius')
			read(buffer, *) str2
			read(str2, *) R
		case ('max_degree')
			read(buffer, *) str2
			read(str2, *) nmax
		case ('errors')
			read(buffer, *) str2
			errors = trim(str2)
		case ('norm')
			read(buffer, *) str2
			norm = trim(str2)
		case ('tide_system')
			read(buffer, *) str2
			tide_system = trim(str2)
		!key    L    M         C                  S           sigma C    sigma S
		case ('end_of_head')
			exit
		end select
	end if
 end do

 !----------------------------------------------------------------------------------------------------
 allocate(C_bar_nm(nmax*(nmax+2)+1))
 do ii = 0, nmax
	do jj = 0, ii !nmax
		if (jj == 0) then
			C_bar_nm((nmax - jj) * (nmax - jj + 1) + nmax - ii + 1) = 0.0_longdp
		else
			nn = (nmax - jj) * (nmax - jj + 1) + 2*(nmax - ii) + 1
			C_bar_nm(nn) = 0.0_longdp
			C_bar_nm(nn+1) = 0.0_longdp
		end if
	end do
 end do
 !----------------------------------------------------------------------------------------------------
 100 continue
 if (index(errors, 'no') /= 0) then
	read(12, *, end=200) coefficient_kind, n, m, c, s
 else
	read(12, *, end=200) coefficient_kind, n, m, c, s, sigma_c, sigma_s
 end if

 if (n < nmin .or. n > nmax) goto 100
 if (index(coefficient_kind, 'gfc') == 0) goto 100

 if (m == 0) then
	C_bar_nm((nmax - m) * (nmax - m + 1) + nmax - n + 1) = c
 else
	nn = (nmax - m) * (nmax - m + 1) + 2*(nmax - n) + 1
	C_bar_nm(nn) = c
	C_bar_nm(nn+1) = s
 end if
 goto 100
 200 continue
 !----------------------------------------------------------------------------------------------------
 close(12)
 !----------------------------------------------------------------------------------------------------
 !-------------------------- End of Reading of Geopotential Model Parameters -------------------------
 !----------------------------------------------------------------------------------------------------



 !----------------------------------------------------------------------------------------------------
 !---------------------------------- Reading of Points Coordinates -----------------------------------
 !----------------------------------------------------------------------------------------------------

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points . . . '

 open(12, file='input.dat')

 nmax_p = 0
 301 continue
 nmax_p = nmax_p + 1
 read(12, *, end=401) temp
 goto 301
 401 continue
 
 allocate(latitude_p(nmax_p))
 allocate(latitude_p_geodetic(nmax_p))
 allocate(longitude_p(nmax_p))
 allocate(h_p(nmax_p))
 allocate(Undulation_p(nmax_p))
 allocate(Undulation_Simple_Assumed_Density_p(nmax_p))
 allocate(Undulation_Assumed_Density_p(nmax_p))
 allocate(Height_Anomaly_p(nmax_p))
 allocate(Gravity_Disturbance_p(nmax_p))
 allocate(Classical_Gravity_Anomaly_p(nmax_p))
 allocate(Classical_Gravity_Anomaly_2_p(nmax_p))
 allocate(Molodensky_Gravity_Anomaly_p(nmax_p))
 !----------------------------------------------------------------------------------------------------
 do ii = 1, nmax_p
	latitude_p(ii) = 0._longdp
	latitude_p_geodetic(ii) = 0._longdp
	longitude_p(ii) = 0._longdp
	h_p(ii) = 0.0_longdp
	Undulation_p(ii) = 0.0_longdp
	Undulation_Simple_Assumed_Density_p(ii) = 0.0_longdp
	Undulation_Assumed_Density_p(ii) = 0.0_longdp
	Height_Anomaly_p(ii) = 0.0_longdp
	Gravity_Disturbance_p(ii) = 0.0_longdp
	Classical_Gravity_Anomaly_p(ii) = 0.0_longdp
	Classical_Gravity_Anomaly_2_p(ii) = 0.0_longdp
	Molodensky_Gravity_Anomaly_p(ii) = 0.0_longdp
 end do
 !----------------------------------------------------------------------------------------------------

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points . . . '

 write(6, *) 'Is Height of Points Ellipsoidal(T) or Orthometric(F):'
 read(5, *) ellipsoidal_true

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points . . . '
 if (ellipsoidal_true) then 
	write(6, *) 'Ellipsoidal height expected!'
 else
	write(6, *) 'Orthometric height expected!'
 end if

 write(6, *) 'Are Topography, Bathymetry, and Ice Thickness Available Separately (T/F):'
 read(5, *) T_B_IT_true

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points . . . '
 if (ellipsoidal_true) then 
	write(6, *) 'Ellipsoidal height expected!'
 else
	write(6, *) 'Orthometric height expected!'
 end if
 if (T_B_IT_true) then
	write(6, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'

	allocate(depth_p(nmax_p))
	allocate(ice_p(nmax_p))

	do ii = 1, nmax_p
		depth_p(ii) = 0.0_longdp
		ice_p(ii) = 0.0_longdp
	end do

 end if

 write(6, *) 'Please Enter Required Iterations for Calculation of the Undulations and the Height Anomalies (0 for Entering of Precisions):'
 read(5, *) Iteration
 if (Iteration == 0) then
	write(*, *) char(27), '[2J'
	write(*, *) char(27), '[1;1H'

	write(6, *) 'Read Ellipsoid Parameters --> Done! '
	write(6, *) 'Read Geopotential Model Parameters --> Done! '
	write(6, *) 'Read Data Points . . . '
	if (ellipsoidal_true) then 
		write(6, *) 'Ellipsoidal height expected!'
	else
		write(6, *) 'Orthometric height expected!'
	end if
	if (T_B_IT_true) then
		write(6, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'
	end if
	write(6, *) 'Number of Iteration(s) = ', Iteration
    write(6, *) 'Please Enter Required Precision for Calculation of the Undulations (m):'
    read(5, *) N_Precision

	write(*, *) char(27), '[2J'
	write(*, *) char(27), '[1;1H'

	write(6, *) 'Read Ellipsoid Parameters --> Done! '
	write(6, *) 'Read Geopotential Model Parameters --> Done! '
	write(6, *) 'Read Data Points . . . '
	if (ellipsoidal_true) then 
		write(6, *) 'Ellipsoidal height expected!'
	else
		write(6, *) 'Orthometric height expected!'
	end if
	if (T_B_IT_true) then
		write(6, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'
	end if
	write(6, *) 'Number of Iteration(s) = ', Iteration
	write(6, fmt="(' Undulation Precision (m) =     ', f7.4)") N_Precision
    write(6, *) 'Please Enter Required Precision for Calculation of the Height Anomalies (m):'
    read(5, *) Zeta_Precision
 else
    N_Precision = 0.0_longdp
    Zeta_Precision = 0.0_longdp
 end if

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points . . . '
 if (ellipsoidal_true) then 
	write(6, *) 'Ellipsoidal height expected!'
 else
	write(6, *) 'Orthometric height expected!'
 end if
 if (T_B_IT_true) then
	write(6, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'
 end if
 write(6, *) 'Number of Iteration(s) = ', Iteration
 if (Iteration == 0) then
	write(6, fmt="(' Undulation Precision (m) =     ', f7.4)") N_Precision
	write(6, fmt="(' Height Anomaly Precision (m) = ', f7.4)") zeta_Precision
 end if

 write(6, *) 'Please Enter Production Name:'
 read(5, *) productionname

 close(12)
 open(12, file='input.dat')

 ii = 0
 300 continue
 if (T_B_IT_true) then
	!------------============= Reading Topography, Bathymetry, and Ice Thickness =============------------
	read(12, *, end=400) latitude, longitude, h, depth, ice
	!------------============= Reading Topography, Bathymetry, and Ice Thickness =============------------
 else
	read(12, *, end=400) latitude, longitude, h
 end if

 ii = ii + 1

 Latitude_p_geodetic(ii) = Latitude
 Latitude_p(ii) = Latitude
 Longitude_p(ii) = Longitude
 h_p(ii) = h

 if (T_B_IT_true) then
	depth_p(ii) = depth
	ice_p(ii) = ice
 end if

 goto 300
 400 continue
 nmax_p = ii
 ii = 0
 !write(6, *) nmax_p
 !pause
 !----------------------------------------------------------------------------------------------------
 close(12)
 !----------------------------------------------------------------------------------------------------
 !------------------------------- End of Reading of Points Coordinates -------------------------------
 !----------------------------------------------------------------------------------------------------

 allocate(p_bar_nm(nmax+1, nmax+1))
 do ii = 0, nmax
	do jj = 0, nmax
		p_bar_nm(ii+1, jj+1) = 0.0_longdp
	end do
 end do


 allocate(dp_bar_nm(nmax+1, nmax+1))
 do ii = 0, nmax
	do jj = 0, nmax
		dp_bar_nm(ii+1, jj+1) = 0.0_longdp
	end do
 end do

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points --> Done! '
 if (ellipsoidal_true) then 
	write(6, *) 'Ellipsoidal height expected!'
 else
	write(6, *) 'Orthometric height expected!'
 end if
 if (T_B_IT_true) then
	write(6, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'
 end if
 write(6, *) 'Number of Iteration(s) = ', Iteration
 if (Iteration == 0) then
	write(6, fmt="(' Undulation Precision (m) =     ', f7.4)") N_Precision
	write(6, fmt="(' Height Anomaly Precision (m) = ', f7.4)") zeta_Precision
 end if

 write(6, *) ''

 write(6, *) 'Number of Data Points = ', nmax_p

 write(6, *) ''

 write(6, *) 'Ellipsoid Parameters:'
 write(6, *) '----------------------------------------------------------------------'
 write(6, *) 'Reference Ellipsoid Name:           ', trim(adjustl(ellipsoidname))
 write(6, *) 'Ellipsoid Gravity Constant =        ', GM_ellipsoid
 write(6, *) 'Semimajor Axis of the Ellipsoid =   ', a
 write(6, *) 'Semiminor Axis of the Ellipsoid =   ', b
 write(6, *) 'Reciprocal Flattening =             ', f_reciprocal
 write(6, *) 'Angular Velocity of the Earth =     ', omega
 write(6, *) 'Normal Potential at the Ellipsoid = ', U_0
 write(6, *) 'm = ω**2 a**2 b/(GM) =              ', m_ellipsoid
 write(6, *) 'Normal Gravity at the Equator =     ', gamma_0
 write(6, *) '----------------------------------------------------------------------'

 write(6, *) ''

 write(6, *) 'Geopotential Model Parameters:'
 write(6, *) '----------------------------------------------------------------------'
 write(6, *) 'Product Type:             ', trim(adjustl(product_type))
 write(6, *) 'Model Name:               ', trim(adjustl(modelname))
 write(6, *) 'Earth Gravity Constant =  ', GM
 write(6, *) 'Radius of the Earth =     ', R
 write(6, *) 'Maximum Degree of Model = ', nmax
 write(6, *) 'Errors Type of Model:     ', trim(adjustl(errors))
 write(6, *) 'Normalization Type:       ', trim(adjustl(norm))
 write(6, *) 'Tide System:              ', trim(adjustl(tide_system))
 write(6, *) '----------------------------------------------------------------------'

 write(6, *) ''
 write(6, *) 'Production Name: ', productionname
 write(6, *) ''

 write(6, *) 'Computing . . .'

 Number_p = 0
 call progress(Number_p, nmax_p) ! generate the progress bar.

 writing = .false.

 call exec_time(days1, hours1, minutes1, seconds1)

!!$OMP PARALLEL num_threads(200)
!!$OMP DO private(ii)
 do ii = 1, nmax_p
!	if (writing .eqv. .false.) then
!		writing = .true.
		call progress(Number_p, nmax_p) ! generate the progress bar.
!		writing = .false.
!	end if

	call Legendre(p_bar_nm, dp_bar_nm, 10, h_p(ii), Latitude_p(ii), longitude_p(ii), a, b)

	if (h_p(ii) >= 0.0_longdp) then
		rho = 2670.0_longdp
	else
		rho = 1030.0_longdp
	end if

	Undulation_p(ii) = undulation(latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm, Iteration, N_Precision)

	if (ellipsoidal_true) then 
		N_O = 0.0_longdp
		if (T_B_IT_true) then
			N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * (h_p(ii) - Undulation_p(ii))**2.0_longdp - 2670.0_longdp * (depth_p(ii) + Undulation_p(ii))**2.0_longdp + 1030.0_longdp * (depth_p(ii) + Undulation_p(ii))**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
			N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * (h_p(ii) - N_E)**2.0_longdp - 2670.0_longdp * (depth_p(ii) + N_E)**2.0_longdp + 1030.0_longdp * (depth_p(ii) + N_E)**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
		else
			N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * (h_p(ii) - Undulation_p(ii))**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
			N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * (h_p(ii) - N_E)**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
		end if
	else
		if (T_B_IT_true) then
			N_O = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * h_p(ii)**2.0_longdp - 2670.0_longdp * depth_p(ii)**2.0_longdp + 1030.0_longdp * depth_p(ii)**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
		else
			N_O = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * h_p(ii)**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
		end if
		N_E = 0.0_longdp
	end if

	Undulation_Simple_Assumed_Density_p(ii) = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * (h_p(ii) - N_E)**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
	Height_Anomaly_p(ii) = height_anomaly(h_p(ii) + N_O, latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm, Iteration, zeta_Precision)
	Gravity_Disturbance_p(ii) = gravity_disturbance(h_p(ii) + N_O, latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100000.0_longdp
	Classical_Gravity_Anomaly_p(ii) = gravity_anomaly_classic_first_approximation(Undulation_Simple_Assumed_Density_p(ii), latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100000.0_longdp
	Classical_Gravity_Anomaly_2_p(ii) = gravity_anomaly_classic_second_approximation(h_p(ii) + N_O, Undulation_Simple_Assumed_Density_p(ii), latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100000.0_longdp
	Molodensky_Gravity_Anomaly_p(ii) = gravity_anomaly_molodensky(h_p(ii) + N_O, Height_Anomaly_p(ii), latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100000.0_longdp

	if (T_B_IT_true) then
		Undulation_Assumed_Density_p(ii) = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * (h_p(ii) - N_E)**2.0_longdp - 2670.0_longdp * (depth_p(ii) + N_E)**2.0_longdp + 1030.0_longdp * (depth_p(ii) + N_E)**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
	end if

	Number_p = Number_p + 1
 end do
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

 call exec_time(days1, hours1, minutes1, seconds1)
 write(6, *) 'calculation time: days = ', days1, 'hours = ', hours1, 'minutes = ', minutes1, 'seconds = ', seconds1
 call exec_time(days2, hours2, minutes2, seconds2)

 call progress_close

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points --> Done! '
 if (ellipsoidal_true) then 
	write(6, *) 'Ellipsoidal height expected!'
 else
	write(6, *) 'Orthometric height expected!'
 end if
 if (T_B_IT_true) then
	write(6, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'
 end if
 write(6, *) 'Number of Iteration(s) = ', Iteration
 if (Iteration == 0) then
	write(6, fmt="(' Undulation Precision (m) =     ', f7.4)") N_Precision
	write(6, fmt="(' Height Anomaly Precision (m) = ', f7.4)") zeta_Precision
 end if

 write(6, *) ''

 write(6, *) 'Number of Data Points = ', nmax_p

 write(6, *) ''

 write(6, *) 'Ellipsoid Parameters:'
 write(6, *) '----------------------------------------------------------------------'
 write(6, *) 'Reference Ellipsoid Name:           ', trim(adjustl(ellipsoidname))
 write(6, *) 'Ellipsoid Gravity Constant =        ', GM_ellipsoid
 write(6, *) 'Semimajor Axis of the Ellipsoid =   ', a
 write(6, *) 'Semiminor Axis of the Ellipsoid =   ', b
 write(6, *) 'Reciprocal Flattening =             ', f_reciprocal
 write(6, *) 'Angular Velocity of the Earth =     ', omega
 write(6, *) 'Normal Potential at the Ellipsoid = ', U_0
 write(6, *) 'm = ω**2 a**2 b/(GM) =              ', m_ellipsoid
 write(6, *) 'Normal Gravity at the Equator =     ', gamma_0
 write(6, *) '----------------------------------------------------------------------'

 write(6, *) ''

 write(6, *) 'Geopotential Model Parameters:'
 write(6, *) '----------------------------------------------------------------------'
 write(6, *) 'Product Type:             ', trim(adjustl(product_type))
 write(6, *) 'Model Name:               ', trim(adjustl(modelname))
 write(6, *) 'Earth Gravity Constant =  ', GM
 write(6, *) 'Radius of the Earth =     ', R
 write(6, *) 'Maximum Degree of Model = ', nmax
 write(6, *) 'Errors Type of Model:     ', trim(adjustl(errors))
 write(6, *) 'Normalization Type:       ', trim(adjustl(norm))
 write(6, *) 'Tide System:              ', trim(adjustl(tide_system))
 write(6, *) '----------------------------------------------------------------------'

 write(6, *) ''
 write(6, *) 'Production Name: ', productionname
 write(6, *) ''

 write(6, *) 'Computing --> Done!'
 write(6, *) 'Saving Results . . .'
 write(6, *) 'calculation time: days = ', days1, 'hours = ', hours1, 'minutes = ', minutes1, 'seconds = ', seconds1

 write(outputname, *) 'output_', trim(modelname), '_', trim(productionname), '.dat'
 write(outputlogname, *) 'output_', trim(modelname), '_', trim(productionname), '.log'
 open(12, file=outputname)

 if (T_B_IT_true) then
	write(12,*) 'Latitude(degree) Longitude(degree) Height(m) Water_Depth(m) Ice_Thickness(m) Undulation(m) Undulation_with_Simple_Assumed_Density(m) Undulation_with_Assumed_Density(m) Heigh_Anomaly(m) Gravity_Disturbance(mGal) Classical_Gravity_Anomaly(mGal) Classical_Gravity_Anomaly_2(mGal) Molodensky_Gravity_Anomaly(mGal)'
 else
	write(12,*) 'Latitude(degree) Longitude(degree) Height(m) Undulation(m) Undulation_with_Simple_Assumed_Density(m) Heigh_Anomaly(m) Gravity_Disturbance(mGal) Classical_Gravity_Anomaly(mGal) Classical_Gravity_Anomaly_2(mGal) Molodensky_Gravity_Anomaly(mGal)'
 endif

 do ii = 1, nmax_p
	call progress(ii, nmax_p) ! generate the progress bar.

	if (T_B_IT_true) then
		write(12, *) latitude_p_geodetic(ii), longitude_p(ii), h_p(ii), depth_p(ii), ice_p(ii), Undulation_p(ii), Undulation_Simple_Assumed_Density_p(ii), Undulation_Assumed_Density_p(ii), Height_Anomaly_p(ii), Gravity_Disturbance_p(ii), Classical_Gravity_Anomaly_p(ii), Classical_Gravity_Anomaly_2_p(ii), Molodensky_Gravity_Anomaly_p(ii)
	else
		write(12, *) latitude_p_geodetic(ii), longitude_p(ii), h_p(ii), Undulation_p(ii), Undulation_Simple_Assumed_Density_p(ii), Height_Anomaly_p(ii), Gravity_Disturbance_p(ii), Classical_Gravity_Anomaly_p(ii), Classical_Gravity_Anomaly_2_p(ii), Molodensky_Gravity_Anomaly_p(ii)
	end if
 end do

 call progress_close

 close(12)

 write(*, *) char(27), '[2J'
 write(*, *) char(27), '[1;1H'

 write(6, *) 'Read Ellipsoid Parameters --> Done! '
 write(6, *) 'Read Geopotential Model Parameters --> Done! '
 write(6, *) 'Read Data Points --> Done! '
 if (ellipsoidal_true) then 
	write(6, *) 'Ellipsoidal height expected!'
 else
	write(6, *) 'Orthometric height expected!'
 end if
 if (T_B_IT_true) then
	write(6, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'
 end if
 write(6, *) 'Number of Iteration(s) = ', Iteration
 if (Iteration == 0) then
	write(6, fmt="(' Undulation Precision (m) =     ', f7.4)") N_Precision
	write(6, fmt="(' Height Anomaly Precision (m) = ', f7.4)") zeta_Precision
 end if

 write(6, *) ''

 write(6, *) 'Number of Data Points = ', nmax_p

 write(6, *) ''

 write(6, *) 'Ellipsoid Parameters:'
 write(6, *) '----------------------------------------------------------------------'
 write(6, *) 'Reference Ellipsoid Name:           ', trim(adjustl(ellipsoidname))
 write(6, *) 'Ellipsoid Gravity Constant =        ', GM_ellipsoid
 write(6, *) 'Semimajor Axis of the Ellipsoid =   ', a
 write(6, *) 'Semiminor Axis of the Ellipsoid =   ', b
 write(6, *) 'Reciprocal Flattening =             ', f_reciprocal
 write(6, *) 'Angular Velocity of the Earth =     ', omega
 write(6, *) 'Normal Potential at the Ellipsoid = ', U_0
 write(6, *) 'm = ω**2 a**2 b/(GM) =              ', m_ellipsoid
 write(6, *) 'Normal Gravity at the Equator =     ', gamma_0
 write(6, *) '----------------------------------------------------------------------'

 write(6, *) ''

 write(6, *) 'Geopotential Model Parameters:'
 write(6, *) '----------------------------------------------------------------------'
 write(6, *) 'Product Type:             ', trim(adjustl(product_type))
 write(6, *) 'Model Name:               ', trim(adjustl(modelname))
 write(6, *) 'Earth Gravity Constant =  ', GM
 write(6, *) 'Radius of the Earth =     ', R
 write(6, *) 'Maximum Degree of Model = ', nmax
 write(6, *) 'Errors Type of Model:     ', trim(adjustl(errors))
 write(6, *) 'Normalization Type:       ', trim(adjustl(norm))
 write(6, *) 'Tide System:              ', trim(adjustl(tide_system))
 write(6, *) '----------------------------------------------------------------------'

 write(6, *) ''
 write(6, *) 'Production Name: ', productionname
 write(6, *) ''

 write(6, *) 'Computing --> Done!'
 write(6, *) 'Saving Results --> Done!'

 write(6, 4141) days1, hours1, minutes1, seconds1
 2121 format(' Computing time: ', I3, ' days ', I3, ' hours ', I3, ' minutes ', F7.3, ' seconds')
 call exec_time(days2, hours2, minutes2, seconds2)
 write(6, 5151) days2, hours2, minutes2, seconds2
 3131 format(' Saving time:    ', I3, ' days ', I3, ' hours ', I3, ' minutes ', F7.3, ' seconds')

 open(12, file=outputlogname)

 write(12, *) 'Read Ellipsoid Parameters --> Done! '
 write(12, *) 'Read Geopotential Model Parameters --> Done! '
 write(12, *) 'Read Data Points --> Done! '
 if (ellipsoidal_true) then 
	write(12, *) 'Ellipsoidal height expected!'
 else
	write(12, *) 'Orthometric height expected!'
 end if
 if (T_B_IT_true) then
	write(12, *) 'Topography, Bathymetry, and Ice Thickness is Available Separately!'
 end if
 write(12, *) 'Number of Iteration(s) = ', Iteration
 if (Iteration == 0) then
	write(12, fmt="(' Undulation Precision (m) =     ', f7.4)") N_Precision
	write(12, fmt="(' Height Anomaly Precision (m) = ', f7.4)") zeta_Precision
 end if

 write(12, *) ''

 write(12, *) 'Number of Data Points = ', nmax_p

 write(12, *) ''

 write(12, *) 'Ellipsoid Parameters:'
 write(12, *) '----------------------------------------------------------------------'
 write(12, *) 'Reference Ellipsoid Name:           ', trim(adjustl(ellipsoidname))
 write(12, *) 'Ellipsoid Gravity Constant =        ', GM_ellipsoid
 write(12, *) 'Semimajor Axis of the Ellipsoid =   ', a
 write(12, *) 'Semiminor Axis of the Ellipsoid =   ', b
 write(12, *) 'Reciprocal Flattening =             ', f_reciprocal
 write(12, *) 'Angular Velocity of the Earth =     ', omega
 write(12, *) 'Normal Potential at the Ellipsoid = ', U_0
 write(12, *) 'm = ω**2 a**2 b/(GM) =              ', m_ellipsoid
 write(12, *) 'Normal Gravity at the Equator =     ', gamma_0
 write(12, *) '----------------------------------------------------------------------'

 write(12, *) ''

 write(12, *) 'Geopotential Model Parameters:'
 write(12, *) '----------------------------------------------------------------------'
 write(12, *) 'Product Type:             ', trim(adjustl(product_type))
 write(12, *) 'Model Name:               ', trim(adjustl(modelname))
 write(12, *) 'Earth Gravity Constant =  ', GM
 write(12, *) 'Radius of the Earth =     ', R
 write(12, *) 'Maximum Degree of Model = ', nmax
 write(12, *) 'Errors Type of Model:     ', trim(adjustl(errors))
 write(12, *) 'Normalization Type:       ', trim(adjustl(norm))
 write(12, *) 'Tide System:              ', trim(adjustl(tide_system))
 write(12, *) '----------------------------------------------------------------------'

 write(12, *) ''
 write(12, *) 'Production Name: ', trim(adjustl(productionname))
 write(12, *) ''

 write(12, *) 'Computing --> Done!'
 write(12, *) 'Saving Results --> Done!'

 write(12, 4141) days1, hours1, minutes1, seconds1
 4141 format(' Computing time: ', I3, ' days ', I3, ' hours ', I3, ' minutes ', F7.3, ' seconds')
 write(12, 5151) days2, hours2, minutes2, seconds2
 5151 format(' Saving time:    ', I3, ' days ', I3, ' hours ', I3, ' minutes ', F7.3, ' seconds')

 close(12)

 deallocate(C_bar_nm)
 deallocate(latitude_p)
 deallocate(latitude_p_geodetic)
 deallocate(longitude_p)
 deallocate(h_p)
 deallocate(Undulation_p)
 deallocate(Undulation_Simple_Assumed_Density_p)
 deallocate(Undulation_Assumed_Density_p)
 deallocate(Height_Anomaly_p)
 deallocate(Classical_Gravity_Anomaly_p)
 deallocate(Classical_Gravity_Anomaly_2_p)
 deallocate(Molodensky_Gravity_Anomaly_p)
 deallocate(p_bar_nm)
 deallocate(dp_bar_nm)

 if (T_B_IT_true) then
	deallocate(depth_p)
	deallocate(ice_p)
 end if

 end
