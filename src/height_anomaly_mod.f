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

 module height_anomaly_mod
   implicit none
 contains
 function height_anomaly(h, latitude, longitude, GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm, Iteration, zeta_Precision) result(zeta)
 !Function to calculate eq. (28) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                             the geoid undulation and the height anomaly using the iteration method,
 !                                             and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !h              = ellipsoidal height of the studying point
 !latitude       = latitude of the studying point
 !GM             = geocentric gravitational constant of the Earth
 !R              = mean radius of the Earth
 !GM_ellipsoid   = geocentric gravitational constant of the ellipsoid
 !a              = semi-major axis of the used ellipsoid
 !b              = semi-minor axis of the used ellipsoid
 !f_reciprocal   = reciprocal flattening of the ellipsoid
 !omega          = angular velocity of the earth
 !C_bar_nm       = the coefficient of the fully normalized spherical harmonics of the model
 !nmax           = maximum degrees of the geopotential model
 !p_bar_nm       = Legendre function values for mentioned latitude
 !dp_bar_nm      = derivatives of Legendre function values for mentioned latitude
 !zeta_precision = arbitrary precision for calculating the height anomaly
 !
 !Reference:
 !	Barthelms F., 2009. Definition of functional of the geopotential and their calculation from spherical harmonic models,
 !                      Scientific Technical Report STR09/02, Helmholtz Zentrum, Potsdam.
 !
	use nrtype
	use W_mod
	use U_mod
	use gamma_mod
	implicit none
	real(longdp), intent(in) :: h, latitude, longitude
	real(longdp), intent(in) :: GM, R, GM_ellipsoid, a, b, f_reciprocal, omega
	real(longdp), dimension(:), intent(in) :: C_bar_nm
	integer(i4b), intent(in) :: nmax
	real(longdp), dimension(:,:), intent(in) :: p_bar_nm, dp_bar_nm
	real(longdp) :: W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp
	real(longdp) :: U_temp, dU_dr_temp, dU_dtheta_temp
	integer(i4b), intent(in) :: Iteration
	real(longdp), intent(in) :: zeta_Precision
	real(longdp) :: zeta, zeta_1, zeta_2
	integer(i4b) :: ii
	zeta = 0.0_longdp
	zeta_1 = 0.0_longdp
	zeta_2 = 0.0_longdp
	call W(W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp, 1, h, latitude, longitude, a, b, GM, R, C_bar_nm, nmax)
	call U(U_temp, dU_dr_temp, dU_dtheta_temp, 1, h, latitude, longitude, GM_ellipsoid, a, b, f_reciprocal, omega, nmax, p_bar_nm, dp_bar_nm)
	zeta_1 = (W_temp - U_temp) / gamma_h(latitude, h, GM_ellipsoid, a, b, f_reciprocal, omega)
	if (Iteration > 0) then
		if (Iteration == 1) then
			zeta = zeta_1
			return
		end if
		do ii = 2, Iteration
			call W(W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp, 1, h, latitude, longitude, a, b, GM, R, C_bar_nm, nmax)
			call U(U_temp, dU_dr_temp, dU_dtheta_temp, 1, h-zeta_1, latitude, longitude, GM_ellipsoid, a, b, f_reciprocal, omega, nmax, p_bar_nm, dp_bar_nm)
			zeta_1 = zeta_1 + (W_temp - U_temp) / gamma_h(latitude, h-zeta_1, GM_ellipsoid, a, b, f_reciprocal, omega)
		end do
		zeta = zeta_1
		return
	else
		do
			call W(W_temp, dW_dr_temp, dW_dtheta_temp, dW_dlambda_temp, 1, h, latitude, longitude, a, b, GM, R, C_bar_nm, nmax)
			call U(U_temp, dU_dr_temp, dU_dtheta_temp, 1, h-zeta_1, latitude, longitude, GM_ellipsoid, a, b, f_reciprocal, omega, nmax, p_bar_nm, dp_bar_nm)
			zeta_2 = zeta_1 + (W_temp - U_temp) / gamma_h(latitude, h-zeta_1, GM_ellipsoid, a, b, f_reciprocal, omega)
			if (abs(zeta_2 - zeta_1) <= zeta_Precision) exit
			zeta_1 = zeta_2
		end do
		zeta = zeta_2
	end if
 end function height_anomaly
 end module height_anomaly_mod
