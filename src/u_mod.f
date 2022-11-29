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

 module U_mod
   implicit none
 contains
 subroutine U(U_temp, dU_dr_temp, dU_dtheta_temp, output_status, h, latitude, longitude, GM_ellipsoid, a, b, f_reciprocal, omega, nmax, p_bar_nm, dp_bar_nm)
 !Subroutine to calculate eq. (5) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                              the geoid undulation and the height anomaly using the iteration method,
 !                                              and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !OUTPUTs:
 !U_temp, dU_dr_temp, dU_dtheta_temp = Normal potential and its derivatives
 !
 !INPUTs:
 !output_status  = determination of output kind
 !h              = ellipsoidal height of the studying point
 !latitude       = latitude of the studying point
 !GM_ellipsoid   = geocentric gravitational constant of the ellipsoid
 !a              = semi-major axis of the used ellipsoid
 !b              = semi-minor axis of the used ellipsoid
 !f_reciprocal   = reciprocal flattening of the ellipsoid
 !omega          = angular velocity of the earth
 !nmax           = maximum degrees of the geopotential model
 !p_bar_nm       = Legendre function values for mentioned latitude
 !dp_bar_nm      = derivatives of Legendre function values for mentioned latitude
 !
 !Reference:
 !	Barthelms F., 2009. Definition of functional of the geopotential and their calculation from spherical harmonic models,
 !                      Scientific Technical Report STR09/02, Helmholtz Zentrum, Potsdam.
 !
	use nrtype
	use coordinates_mod
	use C_2n_mod
	implicit none
	real(longdp), intent(in) :: h, latitude, longitude
	real(longdp), intent(in) :: GM_ellipsoid, a, b, f_reciprocal, omega
	integer(i4b), intent(in) :: nmax
	real(longdp) :: R, rr=0.0_longdp, latitude_geocentric
	integer, intent(in) :: output_status
	real(longdp), intent(out) :: U_temp, dU_dr_temp, dU_dtheta_temp
	real(longdp), dimension(:,:), intent(in) :: p_bar_nm, dp_bar_nm
	integer(i4b) :: ii=0
	integer :: output_status_temp
	logical :: U_status, dU_dr_status, dU_dtheta_status
	real(dp) :: J2, J4, J6, J8, J10, C_2, C_4, C_6, C_8, C_10

	U_temp = 0.0_longdp
	dU_dr_temp = 0.0_longdp
	dU_dtheta_temp = 0.0_longdp

	U_status=.false.
	dU_dr_status=.false.
	dU_dtheta_status=.false.

	output_status_temp = output_status

	if (output_status_temp >= 4) then
		dU_dtheta_status = .true.
		output_status_temp = output_status_temp - 4
	end if
	if (output_status_temp >= 2) then
		dU_dr_status = .true.
		output_status_temp = output_status_temp - 2
	end if
	if (output_status_temp == 1) then
		U_status = .true.
	end if

	R = a
	call Geodetic_to_Geocentric(rr, latitude_geocentric, latitude, longitude, h, a, b)

	J2=0.108262982131D-2
	J4=-.237091120053D-05
	J6=0.608346498882D-8
	J8=-0.142681087920D-10
	J10=0.121439275882D-13
	C_2 = -J2/DSQRT(5.D0)
	C_4 = -J4/3.0D0
	C_6 = -J6/DSQRT(13.D0)
	C_8 = -J8/DSQRT(17.D0)
	C_10 = -J10/DSQRT(21.D0)
	if (U_status) then
		U_temp = (GM_ellipsoid/rr) * (1.0_longdp + (R / rr)**2.0_longdp * C_2 * p_bar_nm(2+1, 0+1) + (R / rr)**4.0_longdp * C_4 * p_bar_nm(4+1, 0+1) + (R / rr)**6.0_longdp * C_6 * p_bar_nm(6+1, 0+1) + (R / rr)**8.0_longdp * C_8 * p_bar_nm(8+1, 0+1) + (R / rr)**10.0_longdp * C_10 * p_bar_nm(10+1, 0+1))
	end if
	if (dU_dr_status) then
		dU_dr_temp = GM_ellipsoid * (-1.0_longdp / rr**2.0_longdp + (-real(2 + 1)) * (R**2.0_longdp / rr**4.0_longdp) * C_2 * p_bar_nm(2+1, 0+1) + (-real(4 + 1)) * (R**4.0_longdp / rr**6.0_longdp) * C_4 * p_bar_nm(4+1, 0+1) + (-real(6 + 1)) * (R**6.0_longdp / rr**8.0_longdp) * C_6 * p_bar_nm(6+1, 0+1) + (-real(8 + 1)) * (R**8.0_longdp / rr**10.0_longdp) * C_8 * p_bar_nm(8+1, 0+1) + (-real(10 + 1)) * (R**10.0_longdp / rr**12.0_longdp) * C_10 * p_bar_nm(10+1, 0+1))
	end if
	if (dU_dtheta_status) then
		dU_dtheta_temp = (GM_ellipsoid/rr) * (1.0_longdp + (R / rr)**2.0_longdp * C_2 * dp_bar_nm(2+1, 0+1) + (R / rr)**4.0_longdp * C_4 * dp_bar_nm(4+1, 0+1) + (R / rr)**6.0_longdp * C_6 * dp_bar_nm(6+1, 0+1) + (R / rr)**8.0_longdp * C_8 * dp_bar_nm(8+1, 0+1) + (R / rr)**10.0_longdp * C_10 * dp_bar_nm(10+1, 0+1))
	end if
 end subroutine U
 end module U_mod
