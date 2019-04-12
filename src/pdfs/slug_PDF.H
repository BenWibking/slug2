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
// This file defines the slug_PDF class, which is a generic holder for
// probability distribution functions made up of one or more
// continuous segments
////////////////////////////////////////////////////////////////////////

#ifndef _slug_PDF_H_
#define _slug_PDF_H_

#include <boost/filesystem.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include "../slug_IO.H"
#include "slug_PDF_segment.H"

// Enum for sampling methods
enum samplingMethod { NO_METHOD, STOP_NEAREST, STOP_BEFORE, STOP_AFTER, 
		      STOP_50, NUMBER, POISSON, SORTED_SAMPLING };

class slug_PDF {

public:
  // Create from PDF descriptor file and using specified random number
  // generator. The optional parameter normalized specifies whether
  // this PDF is normalized to have an integral of unity or not.
  slug_PDF(const char *fname, rng_type *my_rng, slug_ostreams& ostreams_,
	   bool is_normalized=true);

  // Create a PDF that wraps a single segment
  slug_PDF(slug_PDF_segment *new_seg, rng_type *my_rng,
	   slug_ostreams& ostreams_,
	   double normalization = 1.0);

  // Destructor
  ~slug_PDF();

  // Return the min, max, and stochstic limit min and max of the PDF
  double get_xMin() const { return xMin; }
  double get_xMax() const { return xMax; }
  double get_xStochMin() const { return xStochMin; }
  double get_xStochMax() const { return xStochMax; }

  // Get and change sampling method
  void setMethod(const samplingMethod new_method) { method = new_method; }
  const samplingMethod& getMethod() const { return method; }

  // Change the normalization
  void setNorm(double new_norm);

  // Add a segment
  void add_segment(double x_min, double x_max, double wgt, 
		   std::string type_name);

  // Methods to intereact with the limits for stochasticity
  bool has_stoch_lim() const { return range_restrict; }
  void remove_stoch_lim();
  void set_stoch_lim(double x_stoch_min);
  void set_stoch_lim(double x_stoch_min, double x_stoch_max);

  // Expectation value and integral; both defined for the full PDF,
  // for a specified range within the PDF, and only over the part of
  // the PDF that is being treated stochastically
  double expectationVal() const { return expectVal; }
  double expectationVal_restricted() const 
  { return expectVal_restrict; }
  double expectationVal(double a, double b) const;
  double integral() const { return PDFintegral; }
  double integral_restricted() const { return PDFintegral_restrict; }
  double integral(double a, double b) const;

  // Fraction of the PDF by number and by mass that is in the
  // restricted range that we're treating stochstically; also a
  // generic mass fraction over any range
  double num_frac_restrict() const { 
    return PDFintegral_restrict / PDFintegral;
  }
  double mass_frac_restrict() const {
    return expectVal_restrict * PDFintegral_restrict / 
      (expectVal * PDFintegral);
  }
  double mass_frac(double a, double b) const {
    return expectationVal(a, b) * integral(a, b) / 
      (expectVal * PDFintegral);
  }

  // Functions to draw a single value, either over full PDF or a
  // limited range
  double draw() const;             // Draw from full range
  double draw_restricted() const;  // Draw from stochastically-limited range
  double draw(double a, double b) const; // Draw from specified range
  std::vector<double> draw(double a, double b, unsigned long n) const; 
                               // Draw n samples from specified range

  // Function to draw a population with the goal of reaching a certain
  // sum; return value is the sum of the actual population drawn, and
  // the population drawn is returned in pop. Note that the population
  // returned will only be drawn from the stochastically-limited
  // range, and any input target will be adjusted down by the mass
  // fraction that is in the stochastic range.
  double drawPopulation(double target, std::vector<double>& pop) const;

  // Operator to return the PDF evaluated at a particular x or set of x's
  double operator()(const double x) const;
  std::vector<double> operator()(const std::vector<double> x) const;
  
  //Functions to check for (and initialise) variable segments 
  bool init_vsegs();                          //Check & initialise
  bool check_for_vsegs(){return vsegcheck;}   //Check only
  
  //Function to draw a value from each variable segment
  std::vector<double> vseg_draw();
  
  //Clean up variable segments
  void cleanup();   

protected:
  slug_ostreams& ostreams;
	
private:
  // Functions to parse basic and advanced mode PDF files
  void parseBasic(std::ifstream& PDFFile,
		  std::vector<std::string> firstline,
		  int& lineCount);
  void parseAdvanced(std::ifstream& PDFFile, int& lineCount);

  // Data
  // Vector of segments in the PDF
  std::vector<slug_PDF_segment *> segments;  // Segments in the PDF
  std::vector<slug_PDF_segment *> seg_restricted; // Segments for
	  				          // restricted draws

  // Are any segments variable?
  bool vsegcheck = false;	

  // Pointer to the rng
  rng_type *rng;

  // The distributions that will decide which segment to draw from
  boost::variate_generator<rng_type&, 
			   boost::random::discrete_distribution<> > *disc;
  boost::variate_generator<rng_type&, 
			   boost::random::discrete_distribution<> > *disc_restricted;

  // A 50/50 coin toss generator; used for the STOP_50 sampling method
  boost::variate_generator<rng_type&, 
			   boost::random::uniform_smallint<> > *coin;

  // Function invoked if parser fails to parse a file
  [[ noreturn ]] void parseError(int lineCount, std::string line = "", 
				 std::string message = "");

  // Function invoked if the parser reaches EOF without getting all
  // the data it needed
  [[ noreturn ]] void eofError(std::string message = "");

  // Name of the PDF file
  std::string PDFFileName;

  // Weights of the various segments
  std::vector<double> weights;
  std::vector<double> weights_restricted;

  // Lower and upper limits of the PDF
  double xMin, xMax;

  // Lower and upper limits for the stochstic treatment
  double xStochMin, xStochMax;

  // Expectation value and integral under PDF, and under the
  // stochastically-treated part of it only
  double expectVal, expectVal_restrict;
  double PDFintegral, PDFintegral_restrict;

  // Sampling method
  samplingMethod method;

  // Flags for whether this is a normalized PDF, and whether there are
  // limits to the stochastic range
  bool normalized;
  bool range_restrict;
};

#endif
// _slug_PDF_H_
