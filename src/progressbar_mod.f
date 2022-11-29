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

 module ProgressBar_mod
	 implicit none
 contains

 subroutine progress(first, total)
 implicit none
	integer(kind=4)::first, total,k
	character*15 first_str, total_str
	character(len=60)::bar="[                                                  ] ??????%"
	character(len=65)::text="                                                                 "
	character,parameter :: esc = char(27)
	write (first_str, '(I15)') first
	write (total_str, '(I15)') total
	write(unit=bar(54:59),fmt="(f6.2)") 100.0*real(first)/real(total)
	if (first == 1) then
		write(text, *) '    (', trim(adjustl(first_str)), ' point of ', trim(adjustl(total_str)), ' points completed)'
	else
		write(text, *) '    (', trim(adjustl(first_str)), ' points of ', trim(adjustl(total_str)), ' points completed)'
	end if
	do k=1, int(50.0*real(first)/real(total))
		bar(1+k:1+k)='='!char(254)
	end do
	write(unit=6,fmt="(a1,a60)",advance='no') char(13), bar
	write(6,"(a65)",advance='no') text
	return
 end subroutine progress

 subroutine progress_close
	write(6, *) ''
 end subroutine progress_close
 end module ProgressBar_mod
