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

 module gravity_anomaly_mod
   implicit none
 contains
 function gravity_anomaly_classic_first_approximation(N, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) result(ddelta_g)
 !Function to calculate eq. (34) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                             the geoid undulation and the height anomaly using the iteration method,
 !                                             and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !N            = geoidal undulation
 !latitude     = latitude of the studying point
 !longitude    = longitude of the studying point
 !GM           = geocentric gravitational constant of the Earth
 !R            = mean radius of the Earth
 !GM_ellipsoid = geocentric gravitational constant of the ellipsoid
 !a            = semi-major axis of the used ellipsoid
 !b            = semi-minor axis of the used ellipsoid
 !f_reciprocal = reciprocal flattening of the ellipsoid
 !omega        = angular velocity of the Earth
 !C_bar_nm     = the coefficient of the fully normalized spherical harmonics of the model
 !nmax         = maximum degrees of the geopotential model
 !p_bar_nm     = Legendre function values for mentioned latitude
 !dp_bar_nm    = derivatives of Legendre function values for mentioned latitude
 !
 !Reference:
 !	Barthelms F., 2009. Definition of functional of the geopotential and their calculation from spherical harmonic models,
 !                      Scientific Technical Report STR09/02, Helmholtz Zentrum, Potsdam.
 !
	use nrtype
	use gravity_disturbance_mod
	use gamma_mod
	implicit none
	real(longdp), intent(in) :: N, latitude, longitude
	real(longdp), intent(in) :: GM, R, GM_ellipsoid, a, b, f_reciprocal, omega
	real(longdp), dimension(:), intent(in) :: C_bar_nm
	integer(i4b), intent(in) :: nmax
	real(longdp), dimension(:,:), intent(in) :: p_bar_nm, dp_bar_nm
	real(longdp) :: ddelta_g

	ddelta_g = gravity_disturbance(N, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) + N * Dgamma_Dh(latitude, GM_ellipsoid, a, b, f_reciprocal, omega)

 end function gravity_anomaly_classic_first_approximation

 function gravity_anomaly_classic_second_approximation(h, N, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) result(ddelta_g)
 !Function to calculate eq. (35) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                             the geoid undulation and the height anomaly using the iteration method,
 !                                             and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !h            = ellipsoidal height of the studying point
 !N            = geoidal undulation
 !latitude     = latitude of the studying point
 !longitude    = longitude of the studying point
 !GM           = geocentric gravitational constant of the Earth
 !R            = mean radius of the Earth
 !GM_ellipsoid = geocentric gravitational constant of the ellipsoid
 !a            = semi-major axis of the used ellipsoid
 !b            = semi-minor axis of the used ellipsoid
 !f_reciprocal = reciprocal flattening of the ellipsoid
 !omega        = angular velocity of the Earth
 !C_bar_nm     = the coefficient of the fully normalized spherical harmonics of the model
 !nmax         = maximum degrees of the geopotential model
 !p_bar_nm     = Legendre function values for mentioned latitude
 !dp_bar_nm    = derivatives of Legendre function values for mentioned latitude
 !
 !Reference:
 !	Barthelms F., 2009. Definition of functional of the geopotential and their calculation from spherical harmonic models,
 !                      Scientific Technical Report STR09/02, Helmholtz Zentrum, Potsdam.
 !
	use nrtype
	use gravity_disturbance_mod
	use gamma_mod
	implicit none
	real(longdp), intent(in) :: h, N, latitude, longitude
	real(longdp), intent(in) :: GM, R, GM_ellipsoid, a, b, f_reciprocal, omega
	real(longdp), dimension(:), intent(in) :: C_bar_nm
	integer(i4b), intent(in) :: nmax
	real(longdp), dimension(:,:), intent(in) :: p_bar_nm, dp_bar_nm
	real(longdp) :: ddelta_g

	ddelta_g = gravity_disturbance(h, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) - (h - N) * Dgamma_Dh(latitude, GM_ellipsoid, a, b, f_reciprocal, omega) + gamma_h(latitude, h, GM_ellipsoid, a, b, f_reciprocal, omega) - gamma_h(latitude, 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega)

 end function gravity_anomaly_classic_second_approximation

 function gravity_anomaly_molodensky(h, zeta, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) result(ddelta_g)
 !Function to calculate eq. (37) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                             the geoid undulation and the height anomaly using the iteration method,
 !                                             and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !h            = ellipsoidal height of the studying point
 !zeta         = height anomaly
 !latitude     = latitude of the studying point
 !longitude    = longitude of the studying point
 !GM           = geocentric gravitational constant of the Earth
 !R            = mean radius of the Earth
 !GM_ellipsoid = geocentric gravitational constant of the ellipsoid
 !a            = semi-major axis of the used ellipsoid
 !b            = semi-minor axis of the used ellipsoid
 !f_reciprocal = reciprocal flattening of the ellipsoid
 !omega        = angular velocity of the Earth
 !C_bar_nm     = the coefficient of the fully normalized spherical harmonics of the model
 !nmax         = maximum degrees of the geopotential model
 !p_bar_nm     = Legendre function values for mentioned latitude
 !dp_bar_nm    = derivatives of Legendre function values for mentioned latitude
 !
 !Reference:
 !	Barthelms F., 2009. Definition of functional of the geopotential and their calculation from spherical harmonic models,
 !                      Scientific Technical Report STR09/02, Helmholtz Zentrum, Potsdam.
 !
	use nrtype
	use gravity_disturbance_mod
	use gamma_mod
	implicit none
	real(longdp), intent(in) :: h, zeta, latitude, longitude
	real(longdp), intent(in) :: GM, R, GM_ellipsoid, a, b, f_reciprocal, omega
	!real(longdp), dimension(:,:), intent(in) :: C_nm, S_nm
	real(longdp), dimension(:), intent(in) :: C_bar_nm
	integer(i4b), intent(in) :: nmax
	real(longdp), dimension(:,:), intent(in) :: p_bar_nm, dp_bar_nm
	real(longdp) :: ddelta_g

	ddelta_g = gravity_disturbance(h, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) + zeta * Dgamma_Dh(latitude, GM_ellipsoid, a, b, f_reciprocal, omega)

 end function gravity_anomaly_molodensky
 end module gravity_anomaly_mod
