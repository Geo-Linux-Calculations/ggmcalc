`ORG.GLC:`
![GitHub release (latest by date)](https://img.shields.io/github/v/release/Geo-Linux-Calculations/ggmcalc)
![GitHub Release Date](https://img.shields.io/github/release-date/Geo-Linux-Calculations/ggmcalc)
![GitHub repo size](https://img.shields.io/github/repo-size/Geo-Linux-Calculations/ggmcalc)
![GitHub all releases](https://img.shields.io/github/downloads/Geo-Linux-Calculations/ggmcalc/total)
![GitHub](https://img.shields.io/github/license/Geo-Linux-Calculations/ggmcalc)  

# GGMCalc

### About

GGMCalc is a program to compute some components of the gravity field, using Global Geopotential Models, including
1. Undulation
2. Heigh anomaly
3. Gravity disturbance
4. Classical gravity anomaly
5. Molodensky gravity anomaly

### Reference

Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of the geoid undulation and the height anomaly using the iteration method, and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136. https://doi.org/10.1007/s12145-012-0102-2

### Building

`gfortran` and "GNU Make" should be installed on the OS. Using `make` command one can build `ggmcalc`.

### Usage

One can run the program by `./ggmcalc` command under `src` directory.

The global geopotential model (GGM) file `dat/coeff.dat` is used by the software. One can download GGM files from ICGEM website at <http://icgem.gfz-potsdam.de/ICGEM/> and rename the coefficients file to `coeff.dat`.

The ellipsoid parameters file "ellipsoid.dat" is used by the software.
- `dat/ellpsoid.dat.GRS80` contains parameters of GRS80 ellipsoid.
- `dat/ellpsoid.dat.WGS84` contains parameters of WGS84 ellipsoid.
To use a particular ellipsoid, rename its file to `dat/ellipsoid.dat`.

The `dat/input.dat` file contains the coordinates of the calculating points.

Descriptions of the functions and the subroutines including correspoding references can be found in the source codes.

### License

Copyright (C) 2010-2018 Siamak Moazezi

This file is part of GGMCalc.

GGMCalc is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GGMCalc is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with GGMCalc. If not, see <http://www.gnu.org/licenses/>.

Contact info:
* http://www.sourceforge.net/projects/xgravity
* Siamak Moazezi <s.moazezi@srbiau.ac.ir>
