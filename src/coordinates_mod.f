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

 module coordinates_mod
   implicit none
 contains
 !Module to convert some coordinates parameters
 !in the paper "Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of
 !              the geoid undulation and the height anomaly using the iteration method,
 !              and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136."
 !	Written by Siamak Moazezi 2010
 !

 function Radian(Deg) result(Rad)
	use nrtype
	implicit none
	real(longdp), intent(in) :: Deg
	real(longdp) :: Rad

	Rad = (PI_D / 180.0_longdp) * Deg
 end function Radian

 function Degree(Rad) result(Deg)
	use nrtype
	implicit none
	real(longdp), intent(in) :: Rad
	real(longdp) :: Deg

	Deg = (180.0_longdp / PI_D) * Rad
 end function Degree

 subroutine Geodetic_to_Geocentric(rr, latitude_geocentric, latitude_geodetic, longitude, h, a, b)
	use nrtype
	implicit none
	real(longdp), intent(in) :: latitude_geodetic, longitude, h, a, b
	real(longdp), intent(out) :: rr, latitude_geocentric
	real(longdp) :: latitude_geocentric_rad, latitude_geodetic_rad, longitude_rad
	real(longdp) :: E, e_first, N, X, Y, Z

	latitude_geodetic_rad = radian(latitude_geodetic)
	longitude_rad = radian(longitude)

	E = sqrt(a**2.0_longdp - b**2.0_longdp)
	e_first = E / a
	N = a / sqrt(1.0_longdp - e_first**2.0_longdp*sin(latitude_geodetic_rad)**2.0_longdp)

	X = (N + h) * cos(latitude_geodetic_rad) * cos(longitude_rad)
	Y = (N + h) * cos(latitude_geodetic_rad) * sin(longitude_rad)
	Z = (N * (1.0_longdp - e_first**2.0_longdp) + h) * sin(latitude_geodetic_rad)

	rr =sqrt(X**2.0_longdp + Y**2.0_longdp + Z**2.0_longdp) 
	latitude_geocentric_rad = atan(Z / sqrt(X**2.0_longdp + Y**2.0_longdp))
	latitude_geocentric = Degree(latitude_geocentric_rad)

 end subroutine Geodetic_to_Geocentric

 end module coordinates_mod
