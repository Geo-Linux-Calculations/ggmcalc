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

 module gravity_disturbance_mod
   implicit none
 contains
 function gravity_disturbance(h, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) result(delta_g_h)
 !Function to calculate gravity disturbance used in eqs. (34, 36) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                                                              the geoid undulation and the height anomaly using the iteration method,
 !                                                                              and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !h            = ellipsoidal height of the studying point
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
	use W_mod
	use U_mod
	use gamma_mod
	use coordinates_mod
	implicit none
	real(longdp), intent(in) :: h, latitude, longitude
	real(longdp), intent(in) :: GM, R, GM_ellipsoid, a, b, f_reciprocal, omega
	real(longdp), dimension(:), intent(in) :: C_bar_nm
	integer(i4b), intent(in) :: nmax
	real(longdp), dimension(:,:), intent(in) :: p_bar_nm, dp_bar_nm
	real(longdp) :: W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp
	real(longdp) :: U_temp, dU_dr_temp, dU_dtheta_temp
	real(longdp) :: rr, latitude_geocentric
	real(longdp) :: delta_g_h

	delta_g_h = 0.0_longdp

	call Geodetic_to_Geocentric(rr, latitude_geocentric, latitude, longitude, h, a, b)
	
	!call W(W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp, 14, h, latitude, longitude, a, b, GM, R, C_bar_nm, nmax)
	!call U(U_temp, dU_dr_temp, dU_dtheta_temp, 6, h, latitude, longitude, GM_ellipsoid, a, b, f_reciprocal, omega, nmax, p_bar_nm, dp_bar_nm)
	call W(W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp, 2, h, latitude, longitude, a, b, GM, R, C_bar_nm, nmax)
	call U(U_temp, dU_dr_temp, dU_dtheta_temp, 2, h, latitude, longitude, GM_ellipsoid, a, b, f_reciprocal, omega, nmax, p_bar_nm, dp_bar_nm)

	!delta_g_h = sqrt((dW_dr_temp + omega**2*rr*cos(Radian(latitude))**2)**2 + (1.0_longdp/rr * (dW_dtheta_temp - omega**2*rr**2*cos(Radian(latitude))*sin(Radian(latitude))))**2 + (1.0_longdp/(rr*cos(Radian(latitude))) * dW_dlambda_temp)**2) - sqrt((dU_dr_temp + omega**2*rr*cos(Radian(latitude))**2)**2 + (1.0_longdp/rr * (dU_dtheta_temp - omega**2*rr**2*cos(Radian(latitude))*sin(Radian(latitude))))**2)
	delta_g_h = -(dW_dr_temp - dU_dr_temp)

 end function gravity_disturbance
 end module gravity_disturbance_mod
