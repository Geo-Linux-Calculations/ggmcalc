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

 module gamma_mod
   implicit none
 contains
 function gamma_h(latitude, h, GM_ellipsoid, a, b, f_reciprocal, omega) result(gamma_temp)
 !Function to calculate eq. (6) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                            the geoid undulation and the height anomaly using the iteration method,
 !                                            and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !	Revised by Siamak Moazezi 2017
 !
 !latitude     = latitude of the studying point
 !GM_ellipsoid = geocentric gravitational constant of the ellipsoid
 !a            = semi-major axis of the used ellipsoid
 !b            = semi-minor axis of the used ellipsoid
 !f_reciprocal = reciprocal flattening of the ellipsoid
 !omega        = angular velocity of the earth
 !
 !Reference:
 !	Hofmann-Wellenhof, B., Moritz, H., 2006. Physical geodesy, 2nd, Corrected Ed., Springer, Wien New York.
 !
	use nrtype
	use coordinates_mod
	implicit none
	real(longdp), intent(in) :: latitude, h, GM_ellipsoid, a, b, f_reciprocal, omega
	real(longdp) :: latitude_rad, gamma_temp, m_ellipsoid, f, E, e_first, e_second, gamma_a, gamma_b, f_2, f_4, gamma_0
	latitude_rad = radian(latitude)

	f = 1.0_longdp / f_reciprocal
	E = sqrt(a**2.0_longdp - b**2.0_longdp)
	e_second = E / b
	m_ellipsoid = (omega**2.0_longdp * a**2.0_longdp *b) / GM_ellipsoid

	f_2 = -f + (5.0_longdp / 2.0_longdp)*m_ellipsoid + (1.0_longdp / 2.0_longdp)*f**2.0_longdp - (26.0_longdp / 7.0_longdp)*f*m_ellipsoid + (15.0_longdp / 4.0_longdp)*m_ellipsoid**2.0_longdp
	f_4 = -(1.0_longdp / 2.0_longdp)*f**2.0_longdp + (5.0_longdp / 2.0_longdp)*f*m_ellipsoid

	!equation 2-186
	gamma_a = (GM_ellipsoid / (a * b)) * (1.0_longdp - (3.0_longdp / 2.0_longdp)*m_ellipsoid - (3.0_longdp / 14.0_longdp)*e_second**2.0_longdp*m_ellipsoid)
	!equation 2-199
	gamma_0 = gamma_a * (1.0_longdp + f_2*sin(latitude_rad)**2.0_longdp + f_4*sin(latitude_rad)**4.0_longdp)

	!equation 2-215
	gamma_temp = gamma_0 * (1.0_longdp - (2.0_longdp/a) * (1.0_longdp + f + m_ellipsoid - 2.0_longdp*f*sin(latitude_rad)**2.0_longdp) * h + 3.0_longdp/a**2.0_longdp * h**2.0_longdp)
 end function gamma_h

 function Dgamma_Dh(latitude, GM_ellipsoid, a, b, f_reciprocal, omega) result(Dgamma_Dh_temp)
 !Function to calculate derivative (with respect to h) of
 !eq. (6) in paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                  the geoid undulation and the height anomaly using the iteration method,
 !                  and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !	Revised by Siamak Moazezi 2017
 !
 !latitude     = latitude of the studying point
 !GM_ellipsoid = geocentric gravitational constant of the ellipsoid
 !a            = Semi-major axis of the used ellipsoid
 !b            = Semi-minor axis of the used ellipsoid
 !f_reciprocal = reciprocal flattening of the ellipsoid
 !omega        = Angular velocity of the Earth
 !
 !Reference:
 !	Hofmann-Wellenhof, B., Moritz, H., 2006. Physical geodesy, 2nd, Corrected Ed., Springer, Wien New York.
 !
	use nrtype
	use coordinates_mod
	implicit none
	real(longdp), intent(in) :: latitude, GM_ellipsoid, a, b, f_reciprocal, omega
	real(longdp) :: latitude_rad, Dgamma_Dh_temp, m_ellipsoid, f, E, e_first, e_second, gamma_a, gamma_b, f_2, f_4, gamma_0
	latitude_rad = radian(latitude)

	f = 1.0_longdp / f_reciprocal
	E = sqrt(a**2.0_longdp - b**2.0_longdp)
	e_second = E / b
	m_ellipsoid = (omega**2.0_longdp * a**2.0_longdp *b) / GM_ellipsoid

	f_2 = -f + (5.0_longdp / 2.0_longdp)*m_ellipsoid + (1.0_longdp / 2.0_longdp)*f**2.0_longdp - (26.0_longdp / 7.0_longdp)*f*m_ellipsoid + (15.0_longdp / 4.0_longdp)*m_ellipsoid**2.0_longdp
	f_4 = -(1.0_longdp / 2.0_longdp)*f**2.0_longdp + (5.0_longdp / 2.0_longdp)*f*m_ellipsoid

	gamma_a = (GM_ellipsoid / (a * b)) * (1.0_longdp - (3.0_longdp / 2.0_longdp)*m_ellipsoid - (3.0_longdp / 14.0_longdp)*e_second**2.0_longdp*m_ellipsoid)
	gamma_0 = gamma_a * (1.0_longdp + f_2*sin(latitude_rad)**2.0_longdp + f_4*sin(latitude_rad)**4.0_longdp)

	Dgamma_Dh_temp = ((-2.0_longdp * gamma_0) / a) * (1.0_longdp + f + m_ellipsoid - 2.0_longdp*f*sin(latitude_rad)**2.0_longdp)
 end function Dgamma_Dh

 end module gamma_mod
