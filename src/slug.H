/*********************************************************************
Copyright (C) 2014 Robert da Silva, Michele Fumagalli, Mark Krumholz
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

////////////////////////////////////////////////////////////////////////
// slug Header file
// This header introduces common definitions and includes used through
// SLUG.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_H_
#define _slug_H_

#ifdef ENABLE_MPI
#   include "mpi.h"
#endif

// Declare all slug classes here
class slug_PDF;
class slug_PDF_delta;
class slug_PDF_exponential;
class slug_PDF_lognormal;
class slug_PDF_normal;
class slug_PDF_powerlaw;
class slug_PDF_schechter;
class slug_PDF_segment;
class slug_cluster;
class slug_extinction;
class slug_filter;
class slug_filter_set;
class slug_galaxy;
class slug_line;
class slug_line_list;
class slug_parmParser;
class slug_sim;
class slug_specsyn;
class slug_specsyn_hillier;
class slug_specsyn_kurucz;
class slug_specsyn_pauldrach;
class slug_specsyn_planck;
class slug_specsyn_powr;
class slug_specsyn_sb99;
class slug_specsyn_sb99hruv;
class slug_tracks;
class slug_yields;
class slug_yields_agb;
class slug_yields_multiple;
class slug_yields_snii;

// Enum for WR star types
enum WRtype { NONE, WC, WN };

// Enum for output modes
enum outputMode { ASCII, BINARY
#ifdef ENABLE_FITS
		  , FITS
#endif
};

// Struct to represent basic data on stars
typedef struct { 
  double mass;       // Mass at birth
  double birth_time; // Formation time
  double death_time; // Time when star dies
} slug_star;

// Struct to hold a lot of data about a star
typedef struct { 
  double logM;       // Log_10 current mass in Msun
  double logL;       // Log_10 luminoisty in Lsun
  double logR;       // Log_10 radius in Rsun
  double logTeff;    // Log_10 Teff in K
  double logg;       // Log_10 g in cm/s^2
#if 0
  double Mdot;       // Mass loss rate in Msun/yr
  double vWind;      // Wind velocity in km/s
#endif
  WRtype WR;         // Type of WR star
} slug_stardata;

#endif
// _slug_H_
