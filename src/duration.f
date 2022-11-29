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

 module duration

 use date_sub

 implicit none

 contains

 subroutine exec_time(days, hours, minutes, seconds)
	use nrtype
	implicit none
	integer :: test_1, test_2
	integer,dimension(8) :: values
	integer :: days, days1, days2
	real(longdp) :: seconds, seconds1, seconds2
	integer :: hours, minutes
	save test_1, test_2, days1, seconds1

	days = 0
	hours = 0
	minutes = 0
	seconds = 0.0d0

	if (test_1.ne.-317.and.test_2.ne.751) then
		call date_and_time(VALUES=values)
		days1 = jd(values(1), values(2), values(3))
		seconds1 = 3600 * values(5) + 60 * values(6) + values(7) + 1d-3*values(8)

		test_1 = -317
		test_2 = 751
	else
		call date_and_time(VALUES=values)
		days2 = jd(values(1), values(2), values(3))
		seconds2 = 3600 * values(5) + 60 * values(6) + values(7) + 1d-3*values(8)

		seconds = seconds2 - seconds1
		if (seconds >= 0.0d0) then
			days = days2 - days1
		else
			days = days2 - days1 - 1
			seconds = seconds + 86400.0d0
		end if

		hours = int(seconds) / 3600
		seconds = mod(seconds, 3600.0d0)
		minutes = int(seconds) / 60
		seconds = mod(seconds, 60.0d0)

		test_1 = 0
		test_2 = 0
	end if
 end subroutine exec_time

 end module duration
