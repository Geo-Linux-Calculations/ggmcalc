 !-------------------------------------------------------------------------------------------------
 !	Copyright (C) 2010-2018 Siamak Moazezi
 !	
 !	This file is part of GGMCalc.
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

 module W_mod
   implicit none
 contains

 subroutine W(W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp, output_status, h, latitude, longitude, a, b, GM, R, C_bar_nm, nmax)
 !Subroutine to calculate eqs. (4, 16) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                                   the geoid undulation and the height anomaly using the iteration method,
 !                                                   and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !OUTPUTs:
 !W_temp, dW_dr_temp, dW_dtheta_temp = Real potential and its derivatives
 !
 !INPUTs:
 !output_status = determination of output kind
 !h             = ellipsoidal height of the studying point
 !latitude      = latitude of the studying point
 !a             = semi-major axis of the used ellipsoid
 !b             = semi-minor axis of the used ellipsoid
 !GM            = geocentric gravitational constant of the Earth
 !R             = mean radius of the Earth
 !C_bar_nm      = the coefficient of the fully normalized spherical harmonics of the model
 !nmax          = maximum degrees of the geopotential model
 !
 !Reference:
 !	Barthelms F., 2009. Definition of functional of the geopotential and their calculation from spherical harmonic models,
 !                      Scientific Technical Report STR09/02, Helmholtz Zentrum, Potsdam.
 !  Tscherning, C.C., Poder, K. (1982) Some Geodetic Applications of Clenshaw Summation, Bollettino di Geodesia e Scienze
 !          Affini 4: 351-364
 !
	use nrtype
	use coordinates_mod
	use Legendre_mod
	implicit none
	real(longdp), intent(in) :: h, latitude, longitude, a, b
	real(longdp), intent(in) :: GM, R
	real(longdp), dimension(:), intent(in) :: C_bar_nm
	real(longdp), dimension(2*nmax+2) :: s_m_c, ds_m_dr_c, s_dot_m_c, ds_m_dlambda_c
	integer(i4b), intent(in) :: nmax
	integer, intent(in) :: output_status
	real(longdp), intent(out) :: W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp
	real(longdp) :: W_n_temp, dW_dr_n_temp, dW_dtheta_n_temp, dW_dlambda_n_temp
	real(longdp) :: rr=0.0_longdp, latitude_geocentric
	real(longdp) :: t, u, t_u, q, q2, tq, uq, a_m
	real(longdp) :: cos_lambda_2, cos_m_lambda, cos_m_1_lambda, cos_m_2_lambda, sin_m_lambda, sin_m_1_lambda, sin_m_2_lambda
	real(longdp) :: s_m, ds_m_dr, s_dot_m, ds_m_dlambda
	integer :: m, output_status_temp
	logical :: W_status, dW_dr_status, dW_dtheta_status, dW_dlambda_status

	W_n_temp=0.0_longdp
	dW_dr_n_temp=0.0_longdp
	dW_dtheta_n_temp=0.0_longdp
	dW_dlambda_n_temp=0.0_longdp

	W_temp=0.0_longdp
	dW_dr_temp=0.0_longdp
	dW_dtheta_temp=0.0_longdp
	dW_dlambda_temp=0.0_longdp

	W_status=.false.
	dW_dr_status=.false.
	dW_dtheta_status=.false.
	dW_dlambda_status=.false.

	output_status_temp = output_status

	if (output_status_temp >= 8) then
		dW_dlambda_status = .true.
		output_status_temp = output_status_temp - 8
	end if
	if (output_status_temp >= 4) then
		dW_dtheta_status = .true.
		output_status_temp = output_status_temp - 4
	end if
	if (output_status_temp >= 2) then
		dW_dr_status = .true.
		output_status_temp = output_status_temp - 2
	end if
	if (output_status_temp == 1) then
		W_status = .true.
	end if

	call Geodetic_to_Geocentric(rr, latitude_geocentric, latitude, longitude, h, a, b)

	q = R / rr
	q2 = q * q
	t = sin(Radian(latitude_geocentric))
	u = cos(Radian(latitude_geocentric))
	t_u = t / u
	tq = t * q
	uq = u * q
	cos_lambda_2 = 2.0d0 * cos(Radian(longitude))
	
	if (W_status) then
		!$OMP PARALLEL num_threads(200)
		!$OMP DO private(m)
		do m = nmax, 0, -1
			call Clenshaw_s_m(s_m_c, tq, q2, cos_lambda_2, m, nmax, C_bar_nm)
		end do
		!$OMP END DO NOWAIT
		!$OMP END PARALLEL

		m = nmax + 1
		cos_m_2_lambda = cos(m*Radian(longitude))
		sin_m_2_lambda = sin(m*Radian(longitude))

		m = nmax
		cos_m_1_lambda = cos(m*Radian(longitude))
		sin_m_1_lambda = sin(m*Radian(longitude))
		s_m = s_m_c(2*m)*cos_m_1_lambda + s_m_c(2*m+1)*sin_m_1_lambda

		!!$OMP PARALLEL num_threads(4)
		!!$OMP DO private(m)
		do m = nmax-1, 1, -1
			cos_m_lambda = cos_lambda_2 * cos_m_1_lambda - cos_m_2_lambda
			sin_m_lambda = cos_lambda_2 * sin_m_1_lambda - sin_m_2_lambda

			cos_m_2_lambda = cos_m_1_lambda
			cos_m_1_lambda = cos_m_lambda
			sin_m_2_lambda = sin_m_1_lambda
			sin_m_1_lambda = sin_m_lambda
			a_m = sqrt((2.0_longdp*m + 3.0_longdp) / (2.0_longdp*m + 2))
			s_m = a_m * uq * s_m + s_m_c(2*m)*cos_m_lambda + s_m_c(2*m+1)*sin_m_lambda
		end do
		!!$OMP END DO NOWAIT
		!!$OMP END PARALLEL

		m = 0
		a_m = sqrt(3.0_longdp)
		s_m = a_m * uq * s_m + s_m_c(m+1)

		W_temp = (GM/rr) * s_m
	end if

	if (dW_dr_status) then
		!$OMP PARALLEL num_threads(200)
		!$OMP DO private(m)
		do m = nmax, 0, -1
			call Clenshaw_ds_m_dr(ds_m_dr_c, t, t_u, u, tq, q, q2, cos_lambda_2, m, nmax, C_bar_nm)
		end do
		!$OMP END DO NOWAIT
		!$OMP END PARALLEL

		m = nmax + 1
		cos_m_2_lambda = cos(m*Radian(longitude))
		sin_m_2_lambda = sin(m*Radian(longitude))

		m = nmax
		cos_m_1_lambda = cos(m*Radian(longitude))
		sin_m_1_lambda = sin(m*Radian(longitude))
		ds_m_dr = ds_m_dr_c(2*m)*cos_m_1_lambda + ds_m_dr_c(2*m+1)*sin_m_1_lambda

		!!$OMP PARALLEL num_threads(4)
		!!$OMP DO private(m)
		do m = nmax-1, 1, -1
			cos_m_lambda = cos_lambda_2 * cos_m_1_lambda - cos_m_2_lambda
			sin_m_lambda = cos_lambda_2 * sin_m_1_lambda - sin_m_2_lambda

			cos_m_2_lambda = cos_m_1_lambda
			cos_m_1_lambda = cos_m_lambda
			sin_m_2_lambda = sin_m_1_lambda
			sin_m_1_lambda = sin_m_lambda
			!a_m = sqrt((2.0_longdp*(m+1) + 1.0_longdp) / (2.0_longdp*(m+1)))
			a_m = sqrt((2.0_longdp*m + 3.0_longdp) / (2.0_longdp*m + 2))
			!s_m = a_m * uq * s_m + s_dot_m_c(2*m)*cos_m_lambda + s_dot_m_c(2*m+1)*sin_m_lambda
			ds_m_dr = a_m * uq * ds_m_dr + ds_m_dr_c(2*m)*cos_m_lambda + ds_m_dr_c(2*m+1)*sin_m_lambda
		end do
		!!$OMP END DO NOWAIT
		!!$OMP END PARALLEL

		m = 0
		a_m = sqrt(3.0_longdp)
		ds_m_dr = a_m * uq * ds_m_dr + ds_m_dr_c(1)

		dW_dr_temp = (GM / rr**2.0_longdp) * ds_m_dr
	end if

	if (dW_dtheta_status) then
		!$OMP PARALLEL num_threads(200)
		!$OMP DO private(m)
		 do m = nmax, 0, -1
			call Clenshaw_s_dot_m(s_dot_m_c, t, t_u, u, tq, q, q2, cos_lambda_2, m, nmax, C_bar_nm)
		 end do
		!$OMP END DO NOWAIT
		!$OMP END PARALLEL

		m = nmax + 1
		cos_m_2_lambda = cos(m*Radian(longitude))
		sin_m_2_lambda = sin(m*Radian(longitude))

		m = nmax
		cos_m_1_lambda = cos(m*Radian(longitude))
		sin_m_1_lambda = sin(m*Radian(longitude))
		s_dot_m = s_dot_m_c(2*m)*cos_m_1_lambda + s_dot_m_c(2*m+1)*sin_m_1_lambda

		!!$OMP PARALLEL num_threads(4)
		!!$OMP DO private(m)
		do m = nmax-1, 1, -1
			cos_m_lambda = cos_lambda_2 * cos_m_1_lambda - cos_m_2_lambda
			sin_m_lambda = cos_lambda_2 * sin_m_1_lambda - sin_m_2_lambda

			cos_m_2_lambda = cos_m_1_lambda
			cos_m_1_lambda = cos_m_lambda
			sin_m_2_lambda = sin_m_1_lambda
			sin_m_1_lambda = sin_m_lambda
			a_m = sqrt((2.0_longdp*m + 3.0_longdp) / (2.0_longdp*m + 2))
			s_dot_m = a_m * uq * s_dot_m + s_dot_m_c(2*m)*cos_m_lambda + s_dot_m_c(2*m+1)*sin_m_lambda
		end do
		!!$OMP END DO NOWAIT
		!!$OMP END PARALLEL

		m = 0
		a_m = sqrt(3.0_longdp)
		s_dot_m = a_m * uq * s_dot_m + s_dot_m_c(1)!*cos_m_lambda !+ s_m_s(m+1)*sin_m_lambda

		dW_dtheta_temp = (GM/rr) * (1.0_longdp + s_dot_m)
	end if

	if (dW_dlambda_status) then
		!$OMP PARALLEL num_threads(200)
		!$OMP DO private(m)
		do m = nmax, 0, -1
			call Clenshaw_ds_m_dlambda(ds_m_dlambda_c, tq, q2, cos_lambda_2, m, nmax, C_bar_nm)
		end do
		!$OMP END DO NOWAIT
		!$OMP END PARALLEL

		m = nmax + 1
		cos_m_2_lambda = cos(m*Radian(longitude))
		sin_m_2_lambda = sin(m*Radian(longitude))

		m = nmax
		cos_m_1_lambda = cos(m*Radian(longitude))
		sin_m_1_lambda = sin(m*Radian(longitude))
		ds_m_dlambda = ds_m_dlambda_c(2*m)*sin_m_1_lambda - ds_m_dlambda_c(2*m+1)*cos_m_1_lambda

		!!$OMP PARALLEL num_threads(4)
		!!$OMP DO private(m)
		do m = nmax-1, 1, -1
			cos_m_lambda = cos_lambda_2 * cos_m_1_lambda - cos_m_2_lambda
			sin_m_lambda = cos_lambda_2 * sin_m_1_lambda - sin_m_2_lambda

			cos_m_2_lambda = cos_m_1_lambda
			cos_m_1_lambda = cos_m_lambda
			sin_m_2_lambda = sin_m_1_lambda
			sin_m_1_lambda = sin_m_lambda
			a_m = sqrt((2.0_longdp*m + 3.0_longdp) / (2.0_longdp*m + 2))
			ds_m_dlambda = a_m * uq * ds_m_dlambda + ds_m_dlambda_c(2*m)*sin_m_lambda - ds_m_dlambda_c(2*m+1)*cos_m_lambda
		end do
		!!$OMP END DO NOWAIT
		!!$OMP END PARALLEL

		m = 0
		a_m = sqrt(3.0_longdp)
		ds_m_dlambda = a_m * uq * ds_m_dlambda !+ ds_m_dlambda_c(1)!*cos_m_lambda !+ s_m_s(m+1)*sin_m_lambda

		dW_dlambda_temp = (GM/rr) * ds_m_dlambda
	end if
 end subroutine W
 end module W_mod
