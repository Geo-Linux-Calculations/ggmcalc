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

 module C_2n_mod
   implicit none
 contains
 function C_2n(n, a, b, m_ellipsoid) result(C_2n_temp)
 !Function to calculate eq. (3) in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !                                            the geoid undulation and the height anomaly using the iteration method,
 !                                            and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !
 !n           = index of coefficient C_2n
 !a           = semi-major axis of the used ellipsoid
 !b           = semi-minor axis of the used ellipsoid
 !m_ellipsoid = centrifugal force at equator / gravity at equator
 !
 !Reference:
 !	Moritz, H., 1980. Geodetic Reference System 1980, Journal of Geodesy, 54(3), pp.395-405.
 !
	use nrtype
	implicit none
	integer(i4b), intent(in) :: n
	real(longdp), intent(in) :: a, b, m_ellipsoid
	real(longdp) :: C_2n_temp, E=0.0_longdp, e_1=0.0_longdp, e_2=0.0_longdp, q_0=0.0_longdp, CA_ME2=0.0_longdp
	E = sqrt(a**2.0_longdp - b**2.0_longdp)
	e_1 = E / a
	e_2 = E / b
	q_0 = 0.5_longdp * ((1.0_longdp + 3.0_longdp*(b**2.0_longdp / E**2.0_longdp)) * atan(E/b) - 3.0_longdp*(b/E))
	CA_ME2 = (1.0_longdp/3.0_longdp) * (1.0_longdp - (2.0_longdp/15.0_longdp) * ((m_ellipsoid*e_2)/q_0))
	C_2n_temp = (-1.0_longdp)**real(n) * ((3.0_longdp*e_1**(2.0_longdp*real(n)))/((2.0_longdp*real(n) + 1.0_longdp) * (2.0_longdp*real(n) + 3.0_longdp))) * (1.0_longdp - real(n) + 5.0_longdp*real(n) * CA_ME2)
	C_2n_temp = (1.0_longdp / sqrt(2.0_longdp*n + 1.0_longdp)) * C_2n_temp
 end function C_2n
 end module C_2n_mod
