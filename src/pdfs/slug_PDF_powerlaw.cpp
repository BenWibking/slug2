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

#include "slug_PDF_powerlaw.H"
#include <cmath>
#include <vector>

using namespace boost;
using namespace boost::random;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_PDF_powerlaw::
slug_PDF_powerlaw(double sMin_, double sMax_, double sSlope_, 
		  rng_type *rng_, slug_ostreams& ostreams_)
  : slug_PDF_segment(sMin_, sMax_, rng_, ostreams_) 
{
  std::vector<double> tokenVals(1, sSlope_);
  initialize(tokenVals);
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_PDF_powerlaw::~slug_PDF_powerlaw() { 
  delete unidist;
}


////////////////////////////////////////////////////////////////////////
// Initializer
////////////////////////////////////////////////////////////////////////
void
slug_PDF_powerlaw::initialize(const std::vector<double> &tokenVals) {

  //Clean up the distribution object if we already have one
  if (initialised==true && variable_seg==true)
  {
    delete unidist; 
  }

  // Save the slope
  segSlope = tokenVals[0];

  // Build a uniform distribution object with the specified parameters
  uniform_01<> uni01;
  unidist =
    new variate_generator<rng_type&, uniform_01<> >(*rng, uni01);

  // Set the norm, min, max, and expectation values
  if (segSlope != -1.0) {
    norm = (segSlope + 1.0) /
      (pow(segMax, segSlope+1.0) - pow(segMin, segSlope+1.0));
  } else {
    norm = 1.0 / log(segMax/segMin);
  }
  segMinVal = norm * pow(segMin, segSlope);
  segMaxVal = norm * pow(segMax, segSlope);
  expectVal = expectationVal(segMin, segMax);
}


////////////////////////////////////////////////////////////////////////
// Value evaluated at a point
////////////////////////////////////////////////////////////////////////
double
slug_PDF_powerlaw::operator() (double x) {
  if ((x < segMin) || (x >= segMax)) return 0.0;
  return norm * pow(x, segSlope);
}


////////////////////////////////////////////////////////////////////////
// Function to return the expectation value over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_powerlaw::expectationVal(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  if ((segSlope != -1.0) && (segSlope != -2.0)) {
    return (segSlope+1.0) / (segSlope+2.0) * 
      (pow(b1, 2.0+segSlope) - pow(a1, 2.0+segSlope)) /
      (pow(b1, 1.0+segSlope) - pow(a1, 1.0+segSlope));
  } else if (segSlope == -1.0) {
    return (b1 - a1) / log(b1/a1);
  } else {
    return a1 * b1 * log(b1/a1) / (b1 - a1);
  }
}


////////////////////////////////////////////////////////////////////////
// Function to return the integral over a finite interval
////////////////////////////////////////////////////////////////////////
double 
slug_PDF_powerlaw::integral(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  if (segSlope != -1.0) {
    return norm / (segSlope + 1.0) *
      (pow(b1, segSlope+1.0) - pow(a1, segSlope+1.0));
  } else {
    return norm * log(b1/a1);
  }
}


////////////////////////////////////////////////////////////////////////
// Draw a mass from a powerlaw with minimum and maximum
////////////////////////////////////////////////////////////////////////
double
slug_PDF_powerlaw::draw(double a, double b) {
  double a1 = a < segMin ? segMin : a;
  double b1 = b > segMax ? segMax : b;
  double val, udev;

  // Draw a uniform deviate
  udev = (*unidist)();

  // Transform to powerlaw distribution
  if (segSlope != -1.0) {
    val = pow(udev*pow(b1, segSlope+1.0) + 
	      (1.0-udev)*pow(a1, segSlope+1.0), 
	      1.0/(segSlope+1.0));
  } else {
    val = pow(b1, udev) / pow(a1, udev-1.0);
  }
  return val;
}



