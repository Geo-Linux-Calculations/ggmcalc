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

 module Legendre_mod
   implicit none
 contains
 subroutine Legendre(p_bar_nm, dp_bar_nm, nmax, h, latitude_geodetic, longitude, a, b)
 !Subroutine to calculate eqs. (10, 11, 14, 15) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                                            the geoid undulation and the height anomaly using the iteration method,
 !                                                            and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !OUTPUTs:
 !p_bar_nm  = Legendre function values for mentioned latitude
 !dp_bar_nm = derivatives of Legendre function values for mentioned latitude
 !
 !INPUTs:
 !nmax              = maximum degrees of the geopotential model
 !h                 = ellipsoidal height of the studying point
 !latitude_geodetic = latitude of the studying point
 !longitude         = longitude of the studying point
 !a                 = semi-major axis of the used ellipsoid
 !b                 = semi-minor axis of the used ellipsoid
 !
 !Reference:
 !	Holmes, S.A., Featherstone,W.E., 2002. A Unified Approach to the Clenshaw Summation and the Recursive Computation
 !          of Very High Degree and Order Normalised Associated Legendre Functions, Journal of Geodesy, 76(5), pp.279-299
 !
	use nrtype
	use coordinates_mod
	implicit none
	real(longdp), dimension(:,:), intent(out) :: p_bar_nm
	real(longdp), dimension(:,:), intent(out) :: dp_bar_nm
	integer(i4b), intent(in) :: nmax
	integer(i4b) :: ii, jj
	real(longdp), intent(in) :: h, latitude_geodetic, longitude, a, b
	real(longdp) :: rr, latitude
	real(longdp) :: colatitude, colatitude_rad, sin_colatitude_rad, cos_colatitude_rad, n, m, a_nm, b_nm, f_nm
	do ii = 0, nmax
		do jj = 0, nmax
			p_bar_nm(ii+1, jj+1) = 0.0_longdp
		end do
	end do
	do ii = 0, nmax
		do jj = 0, nmax
			dp_bar_nm(ii+1, jj+1) = 0.0_longdp
		end do
	end do

	call Geodetic_to_Geocentric(rr, latitude, latitude_geodetic, longitude, h, a, b)

!----------------------------------------------------------------------------------------------------
!----------------------------------------- Legendre Function ----------------------------------------
!----------------------------------------------------------------------------------------------------
	colatitude = 90.0_longdp - latitude
	colatitude_rad = Radian(colatitude)
	sin_colatitude_rad = sin(colatitude_rad)  !u = sin(colatitude_rad)
	cos_colatitude_rad = cos(colatitude_rad)  !t = cos(colatitude_rad)
	p_bar_nm(0+1, 0+1) = 1.0_longdp
	p_bar_nm(1+1, 1+1) = sqrt(3.0_longdp) * sin_colatitude_rad

	do ii = 2, nmax
		m = real(ii)
		p_bar_nm(ii+1, ii+1) = sin_colatitude_rad * sqrt((2.0_longdp*m + 1.0_longdp) / (2.0_longdp*m)) * p_bar_nm(ii-1+1, ii-1+1)
	end do

	do ii = 0, nmax-1
		n = real(ii) + 1.0_longdp
		m = real(ii)
		a_nm = sqrt(((2.0_longdp*n - 1.0_longdp) * (2.0_longdp*n + 1.0_longdp)) / ((n - m) * (n + m)))
		p_bar_nm(ii+1+1, ii+1) = a_nm * cos_colatitude_rad * p_bar_nm(ii+1, ii+1)
	end do

!$OMP PARALLEL num_threads(200)
!$OMP DO private(ii,jj,m,n,a_nm,b_nm)
	do jj = 0, nmax-2
		do ii = jj+2, nmax
			n = real(ii)
			m = real(jj)
			a_nm = sqrt(((2.0_longdp*n - 1.0_longdp) * (2.0_longdp*n + 1.0_longdp)) / ((n - m) * (n + m)))
			b_nm = sqrt(((2.0_longdp*n + 1.0_longdp) * (n + m - 1.0_longdp) * (n - m - 1.0_longdp)) / ((n - m) * (n + m) * (2.0_longdp*n - 3.0_longdp)))
			p_bar_nm(ii+1, jj+1) = a_nm * cos_colatitude_rad * p_bar_nm(ii-1+1, jj+1) - b_nm * p_bar_nm(ii-2+1, jj+1)
		end do
	end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!----------------------------------------------------------------------------------------------------
!----------------------------------------- Legendre Function ----------------------------------------
!----------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------
!------------------------------- First Derivation of Legendre Function ------------------------------
!----------------------------------------------------------------------------------------------------
!$OMP PARALLEL num_threads(200)
!$OMP DO private(ii,m)
	do ii = 0, nmax
		m = real(ii)
		dp_bar_nm(ii+1, ii+1) = m * (cos_colatitude_rad / sin_colatitude_rad) * p_bar_nm(ii+1, ii+1)
	end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$OMP PARALLEL num_threads(200)
!$OMP DO private(ii,jj,m,n,f_nm)
	do jj = 0, nmax-1
		do ii = jj+1, nmax
			n = real(ii)
			m = real(jj)
			f_nm = sqrt(((n**2.0_longdp - m**2.0_longdp) * (2.0_longdp*n + 1.0_longdp)) / (2.0_longdp*n - 1.0_longdp))
			dp_bar_nm(ii+1, jj+1) = (1.0_longdp/sin_colatitude_rad) * (n * cos_colatitude_rad * p_bar_nm(ii+1, jj+1) - f_nm * p_bar_nm(ii-1+1, jj+1))
		end do
	end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!----------------------------------------------------------------------------------------------------
!------------------------------- First Derivation of Legendre Function ------------------------------
!----------------------------------------------------------------------------------------------------
 end subroutine Legendre

 subroutine Clenshaw_s_m(s_m_c, tq, q2, cos_lambda_2, m, nmax, C_bar_nm)
 !Subroutine to calculate eq. (12) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                               the geoid undulation and the height anomaly using the iteration method,
 !                                               and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !OUTPUTs:
 !s_m_c  = Clenshaw summation of the fully normalized associated Legendre’s function for constant order m
 !
 !INPUTs:
 !tq           = sin(latitude) * (R / r)
 !q2           = (R / r)**2
 !cos_lambda_2 = 2*cos(longitude)
 !m            = order
 !nmax         = maximum degrees of the geopotential model
 !C_bar_nm     = the coefficient of the fully normalized spherical harmonics of the model
 !
 !Reference:
 !	Holmes, S.A., Featherstone,W.E., 2002. A Unified Approach to the Clenshaw Summation and the Recursive Computation
 !          of Very High Degree and Order Normalised Associated Legendre Functions, Journal of Geodesy, 76(5), pp.279-299
 !	Tscherning, C.C., Poder, K. (1982) Some Geodetic Applications of Clenshaw Summation, Bollettino di Geodesia e Scienze
 !          Affini 4: 351-364
 !
	use nrtype
	implicit none
	real(longdp), intent(out), dimension(:) :: s_m_c
	real(longdp), intent(in) :: tq, q2, cos_lambda_2
	integer, intent(in) :: m, nmax
	real(longdp), intent(in), dimension(:) :: C_bar_nm
	real(longdp) :: a_lm, b_lm
	real(longdp) :: s_mm_c, s_mm_s, s_mm_c_pre_1, s_mm_s_pre_1, s_mm_c_pre_2, s_mm_s_pre_2, s_m, s_m_pre
	real(longdp) :: cos_m_lambda, cos_m_1_lambda, cos_m_2_lambda, sin_m_lambda, sin_m_1_lambda, sin_m_2_lambda
	integer :: l, n

	if (m == nmax) then
		s_m_c(2*m) = C_bar_nm(1)
		s_m_c(2*m+1) = C_bar_nm(2)
	else if (m == nmax-1) then
		s_mm_c_pre_1 = C_bar_nm(3)
		s_mm_s_pre_1 = C_bar_nm(4)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * tq
		s_mm_c = a_lm * s_mm_c_pre_1 + C_bar_nm(5)
		s_mm_s = a_lm * s_mm_s_pre_1 + C_bar_nm(6)

		s_m_c(2*m) = s_mm_c
		s_m_c(2*m+1) = s_mm_s
	else if (m <= nmax-2 .and. m >= 1) then
		s_mm_c_pre_2 = C_bar_nm((nmax - m) * (nmax - m + 1) + 1)
		s_mm_s_pre_2 = C_bar_nm((nmax - m) * (nmax - m + 1) + 2)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * tq
		s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 3)
		s_mm_s_pre_1 = a_lm * s_mm_s_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 4)

		do l = nmax-2, m, -1
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1 - m) * (l+1 + m))) * tq
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + m + 1.0_longdp) * (l - m + 1.0_longdp)) / ((l+2 - m) * (l+2 + m) * (2.0_longdp*l + 1.0_longdp))) * q2
			s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 1)
			s_mm_s = a_lm * s_mm_s_pre_1 - b_lm * s_mm_s_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 2)
			s_mm_c_pre_2 = s_mm_c_pre_1
			s_mm_s_pre_2 = s_mm_s_pre_1
			s_mm_c_pre_1 = s_mm_c
			s_mm_s_pre_1 = s_mm_s
		end do

		s_m_c(2*m) = s_mm_c
		s_m_c(2*m+1) = s_mm_s
	else if (m == 0) then
		s_mm_c_pre_2 = C_bar_nm(nmax * (nmax + 1) + 1)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / (nmax * nmax)) * tq
		s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + C_bar_nm(nmax * (nmax + 1) + 2)

		do l = nmax-2, m, -1
			!n = l*l + 1
			!n = nmax * (nmax + 1) + nmax - l + 1
			!a_lm = sqrt(((2.0_longdp*(l+1) - 1.0_longdp) * (2.0_longdp*(l+1) + 1.0_longdp)) / ((l+1) * (l+1)))
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1) * (l+1))) * tq
			!b_lm = sqrt(((2.0_longdp*(l+2) + 1.0_longdp) * ((l+2) - 1.0_longdp) * ((l+2) - 1.0_longdp)) / ((l+2) * (l+2) * (2.0_longdp*(l+2) - 3.0_longdp)))
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + 1.0_longdp) * (l + 1.0_longdp)) / ((l+2) * (l+2) * (2.0_longdp*l + 1.0_longdp))) * q2
			s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + C_bar_nm(nmax * (nmax + 1) + nmax - l + 1)
			s_mm_c_pre_2 = s_mm_c_pre_1
			s_mm_c_pre_1 = s_mm_c
		end do

		s_m_c(1) = s_mm_c
	end if
 end subroutine Clenshaw_s_m


 subroutine Clenshaw_s_dot_m(s_m_c, t, t_u, u, tq, q, q2, cos_lambda_2, m, nmax, C_bar_nm)
 !Subroutine to calculate derivative (with respect to latitude) of
 !eq. (12) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                       the geoid undulation and the height anomaly using the iteration method,
 !                       and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !OUTPUT:
 !s_m_c  = Clenshaw summation of derivative of the fully normalized associated Legendre’s function respect to the latitude for constant order m
 !
 !INPUTs:
 !t            = sin(latitude)
 !t_u          = sin(latitude) / cos(latitude)
 !u            = cos(latitude)
 !tq           = sin(latitude) * (R / r)
 !q            = (R / r)
 !q2           = (R / r)**2
 !cos_lambda_2 = 2*cos(longitude)
 !m            = order
 !nmax         = maximum degrees of the geopotential model
 !C_bar_nm     = the coefficient of the fully normalized spherical harmonics of the model
 !
 !Reference:
 !	Holmes, S.A., Featherstone,W.E., 2002. A Unified Approach to the Clenshaw Summation and the Recursive Computation
 !          of Very High Degree and Order Normalised Associated Legendre Functions, Journal of Geodesy, 76(5), pp.279-299
 !	Tscherning, C.C., Poder, K. (1982) Some Geodetic Applications of Clenshaw Summation, Bollettino di Geodesia e Scienze
 !          Affini 4: 351-364
 !
	use nrtype
	implicit none
	real(longdp), intent(out), dimension(:) :: s_m_c
	real(longdp), intent(in) :: t, t_u, u, tq, q, q2, cos_lambda_2
	integer, intent(in) :: m, nmax
	real(longdp), intent(in), dimension(:) :: C_bar_nm
	real(longdp) :: a_lm, b_lm
	real(longdp) :: s_mm_c, s_mm_s, s_mm_c_pre_1, s_mm_s_pre_1, s_mm_c_pre_2, s_mm_s_pre_2, s_m, s_m_pre
	real(longdp) :: s_dot_mm_c, s_dot_mm_s, s_dot_mm_c_pre_1, s_dot_mm_s_pre_1, s_dot_mm_c_pre_2, s_dot_mm_s_pre_2, s_dot_m, s_dot_m_pre
	real(longdp) :: cos_m_lambda, cos_m_1_lambda, cos_m_2_lambda, sin_m_lambda, sin_m_1_lambda, sin_m_2_lambda
	integer :: l, n

	if (m == nmax) then
		s_m_c(2*m) = m * t_u * C_bar_nm(1)
		s_m_c(2*m+1) = m * t_u * C_bar_nm(2)
	else if (m == nmax-1) then
		s_mm_c_pre_1 = C_bar_nm(3)
		s_mm_s_pre_1 = C_bar_nm(4)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * q
		s_dot_mm_c = a_lm * s_mm_c_pre_1
		s_dot_mm_s = a_lm * s_mm_s_pre_1
		a_lm = a_lm * t
		s_mm_c = a_lm * s_mm_c_pre_1 + C_bar_nm(5)
		s_mm_s = a_lm * s_mm_s_pre_1 + C_bar_nm(6)

		s_m_c(2*m) = m * t_u * s_mm_c - u * s_dot_mm_c
		s_m_c(2*m+1) = m * t_u * s_mm_s - u * s_dot_mm_s
	else if (m <= nmax-2 .and. m >= 1) then
		s_mm_c_pre_2 = C_bar_nm((nmax - m) * (nmax - m + 1) + 1)
		s_mm_s_pre_2 = C_bar_nm((nmax - m) * (nmax - m + 1) + 2)
		s_dot_mm_c_pre_2 = 0.0_longdp !s_mm_c_pre_2
		s_dot_mm_s_pre_2 = 0.0_longdp !s_mm_s_pre_2
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * q
		s_dot_mm_c_pre_1 = a_lm * s_mm_c_pre_2
		s_dot_mm_s_pre_1 = a_lm * s_mm_s_pre_2
		a_lm = a_lm * t
		s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 3)
		s_mm_s_pre_1 = a_lm * s_mm_s_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 4)

		do l = nmax-2, m, -1
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1 - m) * (l+1 + m))) * q
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + m + 1.0_longdp) * (l - m + 1.0_longdp)) / ((l+2 - m) * (l+2 + m) * (2.0_longdp*l + 1.0_longdp))) * q2
			s_dot_mm_c = a_lm * (s_dot_mm_c_pre_1 * t + s_mm_c_pre_1) - b_lm * s_dot_mm_c_pre_2
			s_dot_mm_s = a_lm * (s_dot_mm_s_pre_1 * t + s_mm_s_pre_1) - b_lm * s_dot_mm_s_pre_2
			a_lm = a_lm * t
			s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 1)
			s_mm_s = a_lm * s_mm_s_pre_1 - b_lm * s_mm_s_pre_2 + C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 2)

			s_mm_c_pre_2 = s_mm_c_pre_1
			s_mm_s_pre_2 = s_mm_s_pre_1
			s_mm_c_pre_1 = s_mm_c
			s_mm_s_pre_1 = s_mm_s
			s_dot_mm_c_pre_2 = s_dot_mm_c_pre_1
			s_dot_mm_s_pre_2 = s_dot_mm_s_pre_1
			s_dot_mm_c_pre_1 = s_dot_mm_c
			s_dot_mm_s_pre_1 = s_dot_mm_s
		end do

		s_m_c(2*m) = m * t_u * s_mm_c - u * s_dot_mm_c
		s_m_c(2*m+1) = m * t_u * s_mm_s - u * s_dot_mm_s
	else if (m == 0) then
		s_mm_c_pre_2 = C_bar_nm(nmax * (nmax + 1) + 1)
		s_dot_mm_c_pre_2 = 0.0_longdp !s_mm_c_pre_2
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / (nmax * nmax)) * q
		s_dot_mm_c_pre_1 = a_lm * s_mm_c_pre_2
		a_lm = a_lm * t
		s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + C_bar_nm(nmax * (nmax + 1) + 2)

		do l = nmax-2, m, -1
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1) * (l+1))) * q
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + 1.0_longdp) * (l + 1.0_longdp)) / ((l+2) * (l+2) * (2.0_longdp*l + 1.0_longdp))) * q2
			s_dot_mm_c = a_lm * (s_dot_mm_c_pre_1 * t + s_mm_c_pre_1) - b_lm * s_dot_mm_c_pre_2
			a_lm = a_lm * t
			s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + C_bar_nm(nmax * (nmax + 1) + nmax - l + 1)

			s_mm_c_pre_2 = s_mm_c_pre_1
			s_mm_c_pre_1 = s_mm_c
			s_dot_mm_c_pre_2 = s_dot_mm_c_pre_1
			s_dot_mm_c_pre_1 = s_dot_mm_c
		end do

		s_m_c(1) = - u * s_dot_mm_c
	end if
 end subroutine Clenshaw_s_dot_m

 subroutine Clenshaw_ds_m_dlambda(s_m_c, tq, q2, cos_lambda_2, m, nmax, C_bar_nm)
 !Subroutine to calculate derivative (with respect to longitude) of
 !eq. (12) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                       the geoid undulation and the height anomaly using the iteration method,
 !                       and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !OUTPUT:
 !s_m_c  = Clenshaw summation of derivative of the fully normalized associated Legendre’s function respect to the longitude for constant order m
 !
 !INPUTs:
 !tq           = sin(latitude) * (R / r)
 !q2           = (R / r)**2
 !cos_lambda_2 = 2*cos(longitude)
 !m            = order
 !nmax         = maximum degrees of the geopotential model
 !C_bar_nm     = the coefficient of the fully normalized spherical harmonics of the model
 !
 !Reference:
 !	Holmes, S.A., Featherstone,W.E., 2002. A Unified Approach to the Clenshaw Summation and the Recursive Computation
 !          of Very High Degree and Order Normalised Associated Legendre Functions, Journal of Geodesy, 76(5), pp.279-299
 !	Tscherning, C.C., Poder, K. (1982) Some Geodetic Applications of Clenshaw Summation, Bollettino di Geodesia e Scienze
 !          Affini 4: 351-364
 !
	use nrtype
	implicit none
	real(longdp), intent(out), dimension(:) :: s_m_c
	real(longdp), intent(in) :: tq, q2, cos_lambda_2
	integer, intent(in) :: m, nmax
	real(longdp), intent(in), dimension(:) :: C_bar_nm
	real(longdp) :: a_lm, b_lm
	real(longdp) :: s_mm_c, s_mm_s, s_mm_c_pre_1, s_mm_s_pre_1, s_mm_c_pre_2, s_mm_s_pre_2, s_m, s_m_pre
	real(longdp) :: cos_m_lambda, cos_m_1_lambda, cos_m_2_lambda, sin_m_lambda, sin_m_1_lambda, sin_m_2_lambda
	integer :: l, n

	if (m == nmax) then
		s_m_c(2*m) = m*C_bar_nm(1)
		s_m_c(2*m+1) = m*C_bar_nm(2)
	else if (m == nmax-1) then
		s_mm_c_pre_1 = m*C_bar_nm(3)
		s_mm_s_pre_1 = m*C_bar_nm(4)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * tq
		s_mm_c = a_lm * s_mm_c_pre_1 + m*C_bar_nm(5)
		s_mm_s = a_lm * s_mm_s_pre_1 + m*C_bar_nm(6)

		s_m_c(2*m) = s_mm_c
		s_m_c(2*m+1) = s_mm_s
	else if (m <= nmax-2 .and. m >= 1) then
		s_mm_c_pre_2 = m*C_bar_nm((nmax - m) * (nmax - m + 1) + 1)
		s_mm_s_pre_2 = m*C_bar_nm((nmax - m) * (nmax - m + 1) + 2)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * tq
		s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + m*C_bar_nm((nmax - m) * (nmax - m + 1) + 3)
		s_mm_s_pre_1 = a_lm * s_mm_s_pre_2 + m*C_bar_nm((nmax - m) * (nmax - m + 1) + 4)

		do l = nmax-2, m, -1
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1 - m) * (l+1 + m))) * tq
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + m + 1.0_longdp) * (l - m + 1.0_longdp)) / ((l+2 - m) * (l+2 + m) * (2.0_longdp*l + 1.0_longdp))) * q2
			s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + m*C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 1)
			s_mm_s = a_lm * s_mm_s_pre_1 - b_lm * s_mm_s_pre_2 + m*C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 2)
			s_mm_c_pre_2 = s_mm_c_pre_1
			s_mm_s_pre_2 = s_mm_s_pre_1
			s_mm_c_pre_1 = s_mm_c
			s_mm_s_pre_1 = s_mm_s
		end do

		s_m_c(2*m) = s_mm_c
		s_m_c(2*m+1) = s_mm_s
	else if (m == 0) then
		s_mm_c_pre_2 = C_bar_nm(nmax * (nmax + 1) + 1)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / (nmax * nmax)) * tq
		s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 !+ m*C_bar_nm(nmax * (nmax + 1) + 2)

		do l = nmax-2, m, -1
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1) * (l+1))) * tq
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + 1.0_longdp) * (l + 1.0_longdp)) / ((l+2) * (l+2) * (2.0_longdp*l + 1.0_longdp))) * q2
			s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 !+ m*C_bar_nm(nmax * (nmax + 1) + nmax - l + 1)
			s_mm_c_pre_2 = s_mm_c_pre_1
			s_mm_c_pre_1 = s_mm_c
		end do

		s_m_c(1) = s_mm_c
	end if
 end subroutine Clenshaw_ds_m_dlambda

 subroutine Clenshaw_ds_m_dr(s_m_c, t, t_u, u, tq, q, q2, cos_lambda_2, m, nmax, C_bar_nm)
 !Subroutine to calculate derivative (with respect to radius) of
 !eq. (12) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                       the geoid undulation and the height anomaly using the iteration method,
 !                       and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !OUTPUT:
 !s_m_c  = Clenshaw summation of derivative of the fully normalized associated Legendre’s function respect to the radius for constant order m
 !
 !INPUTs:
 !t            = sin(latitude)
 !t_u          = sin(latitude) / cos(latitude)
 !u            = cos(latitude)
 !tq           = sin(latitude) * (R / r)
 !q            = (R / r)
 !q2           = (R / r)**2
 !cos_lambda_2 = 2*cos(longitude)
 !m            = order
 !nmax         = maximum degrees of the geopotential model
 !C_bar_nm     = the coefficient of the fully normalized spherical harmonics of the model
 !
 !Reference:
 !	Holmes, S.A., Featherstone,W.E., 2002. A Unified Approach to the Clenshaw Summation and the Recursive Computation
 !          of Very High Degree and Order Normalised Associated Legendre Functions, Journal of Geodesy, 76(5), pp.279-299
 !	Tscherning, C.C., Poder, K. (1982) Some Geodetic Applications of Clenshaw Summation, Bollettino di Geodesia e Scienze
 !          Affini 4: 351-364
 !
	use nrtype
	implicit none
	real(longdp), intent(out), dimension(:) :: s_m_c
	real(longdp), intent(in) :: t, t_u, u, tq, q, q2, cos_lambda_2
	integer, intent(in) :: m, nmax
	real(longdp), intent(in), dimension(:) :: C_bar_nm
	real(longdp) :: a_lm, b_lm
	real(longdp) :: s_mm_c, s_mm_s, s_mm_c_pre_1, s_mm_s_pre_1, s_mm_c_pre_2, s_mm_s_pre_2, s_m, s_m_pre
	real(longdp) :: ds_mm_dr_c, ds_mm_dr_s, ds_mm_dr_c_pre_1, ds_mm_dr_s_pre_1, ds_mm_dr_c_pre_2, ds_mm_dr_s_pre_2, ds_m_dr, ds_m_dr_pre
	real(longdp) :: cos_m_lambda, cos_m_1_lambda, cos_m_2_lambda, sin_m_lambda, sin_m_1_lambda, sin_m_2_lambda
	integer :: l, n

	if (m == nmax) then
		s_m_c(2*m) = -(nmax + 1) * C_bar_nm(1)
		s_m_c(2*m+1) = -(nmax + 1) * C_bar_nm(2)
	else if (m == nmax-1) then
		ds_mm_dr_c_pre_1 = -(nmax + 1) * C_bar_nm(3)
		ds_mm_dr_s_pre_1 = -(nmax + 1) * C_bar_nm(4)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * tq
		ds_mm_dr_c = a_lm * ds_mm_dr_c_pre_1 - nmax * C_bar_nm(5)
		ds_mm_dr_s = a_lm * ds_mm_dr_s_pre_1 - nmax * C_bar_nm(6)

		s_m_c(2*m) = ds_mm_dr_c
		s_m_c(2*m+1) = ds_mm_dr_s
	else if (m <= nmax-2 .and. m >= 1) then
		ds_mm_dr_c_pre_2 = -(nmax + 1) * C_bar_nm((nmax - m) * (nmax - m + 1) + 1)
		ds_mm_dr_s_pre_2 = -(nmax + 1) * C_bar_nm((nmax - m) * (nmax - m + 1) + 2)

		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / ((nmax - m) * (nmax + m))) * tq
		ds_mm_dr_c_pre_1 = a_lm * ds_mm_dr_c_pre_2 - nmax * C_bar_nm((nmax - m) * (nmax - m + 1) + 3)
		ds_mm_dr_s_pre_1 = a_lm * ds_mm_dr_s_pre_2 - nmax * C_bar_nm((nmax - m) * (nmax - m + 1) + 4)

		do l = nmax-2, m, -1
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1 - m) * (l+1 + m))) * tq
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + m + 1.0_longdp) * (l - m + 1.0_longdp)) / ((l+2 - m) * (l+2 + m) * (2.0_longdp*l + 1.0_longdp))) * q2
			ds_mm_dr_c = a_lm * ds_mm_dr_c_pre_1 - b_lm * ds_mm_dr_c_pre_2 - (l + 1) * C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 1)
			ds_mm_dr_s = a_lm * ds_mm_dr_s_pre_1 - b_lm * ds_mm_dr_s_pre_2 - (l + 1) * C_bar_nm((nmax - m) * (nmax - m + 1) + 2*(nmax - l) + 2)

			ds_mm_dr_c_pre_2 = ds_mm_dr_c_pre_1
			ds_mm_dr_s_pre_2 = ds_mm_dr_s_pre_1
			ds_mm_dr_c_pre_1 = ds_mm_dr_c
			ds_mm_dr_s_pre_1 = ds_mm_dr_s
		end do

		s_m_c(2*m) = ds_mm_dr_c
		s_m_c(2*m+1) = ds_mm_dr_s
	else if (m == 0) then
		ds_mm_dr_c_pre_2 = -(nmax + 1) * C_bar_nm(nmax * (nmax + 1) + 1)
		a_lm = sqrt(((2.0_longdp*nmax - 1.0_longdp) * (2.0_longdp*nmax + 1.0_longdp)) / (nmax * nmax)) * tq
		ds_mm_dr_c_pre_1 = a_lm * ds_mm_dr_c_pre_2 - nmax * C_bar_nm(nmax * (nmax + 1) + 2)

		do l = nmax-2, m, -1
			a_lm = sqrt(((2.0_longdp*l + 1.0_longdp) * (2.0_longdp*l + 3.0_longdp)) / ((l+1) * (l+1))) * tq
			b_lm = sqrt(((2.0_longdp*l + 5.0_longdp) * (l + 1.0_longdp) * (l + 1.0_longdp)) / ((l+2) * (l+2) * (2.0_longdp*l + 1.0_longdp))) * q2
			ds_mm_dr_c = a_lm * ds_mm_dr_c_pre_1 - b_lm * ds_mm_dr_c_pre_2 - (l + 1) * C_bar_nm(nmax * (nmax + 1) + nmax - l + 1)

			ds_mm_dr_c_pre_2 = ds_mm_dr_c_pre_1
			ds_mm_dr_c_pre_1 = ds_mm_dr_c
		end do

		s_m_c(1) = ds_mm_dr_c
	end if
 end subroutine Clenshaw_ds_m_dr
 end module Legendre_mod
