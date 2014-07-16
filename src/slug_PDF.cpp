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

#include <cstdlib>
#include <cmath>
#include <limits>
#include <string>
#include "slug_PDF.H"
#include "slug_PDF_lognormal.H"
#include "slug_PDF_normal.H"
#include "slug_PDF_powerlaw.H"
#include "slug_PDF_schechter.H"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/poisson_distribution.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace boost::random;

#define BIG (numeric_limits<double>::max())

////////////////////////////////////////////////////////////////////////
// Constructor from a single segment
////////////////////////////////////////////////////////////////////////
slug_PDF::slug_PDF(slug_PDF_segment *new_seg, rng_type *my_rng,
		   double normalization) :
  disc(NULL), disc_restricted(NULL) {

  // Set pointer to the rng
  rng = my_rng;

  // Set the sampling method to its default
  method = STOP_NEAREST;

  // Set up the normalization
  if (normalization==1.0) {
    normalized = true;
    PDFintegral = PDFintegral_restrict = 1.0;
  } else {
    normalized = false;
    PDFintegral = PDFintegral_restrict = normalization;
  }

  // Set the weights
  weights.push_back(PDFintegral);

  // Set the expectation value
  expectVal = expectVal_restrict = new_seg->expectationVal();

  // Store the limits
  xMin = xStochMin = new_seg->sMin();
  xMax = xStochMax = new_seg->sMax();
  range_restrict = false;

  // Store the segment
  segments.push_back(new_seg);
}

////////////////////////////////////////////////////////////////////////
// Constructor from file
////////////////////////////////////////////////////////////////////////
slug_PDF::slug_PDF(const char *PDF, rng_type *my_rng,
		   bool is_normalized) :
  disc(NULL), disc_restricted(NULL) {

  // Set pointer to the rng
  rng = my_rng;

  // Set the sampling method to its default value
  method = STOP_NEAREST;

  // Record if this PDF is assumed to be normalized to unity or not
  normalized = is_normalized;

  // Try to open the PDF file; search for it relative to SLUG_DIR
  // environment variable first if that is set, and then search
  // relative to current directory
  ifstream PDFFile;
  char *slug_dir = getenv("SLUG_DIR");
  path PDFpath(PDF), PDFfullPath;
  if (slug_dir != NULL) {
    // Try opening relative to SLUG_DIR
    PDFfullPath = path(slug_dir) / PDFpath;
    PDFFile.open(PDFfullPath.c_str());
  }
  if (PDFFile.is_open()) {
    PDFpath = PDFfullPath;
  } else {
    // Try opening relative to current path
    PDFFile.open(PDFpath.c_str());
  }
  if (!PDFFile.is_open()) {
    // Couldn't open file, so bail out
    cerr << "slug error: unable to open PDF file " 
	 << PDFpath.string();
    if (slug_dir != NULL)
      cerr << " or " << PDFfullPath.string();
    cerr << endl;
    exit(1);
  }

  // Save name of PDF file
  PDFFileName = PDFpath.string();

  // We've successfully opened the PDF file. Read the first
  // non-whitespace line.
  string line, linecopy;
  int lineCount = 0;
  while (!PDFFile.eof()) {
    getline(PDFFile, line);
    linecopy = line;
    lineCount++;
    trim(line);
    if ((line.compare(0, 1, "#") != 0) && 
	(line.length() != 0)) break;
  }

  // Tokenize the line by whitespace and commas
  vector<string> tokens;
  split(tokens, line, is_any_of("\t ,"), token_compress_on);

  // Look at the first token
  to_lower(tokens[0]);
  if (tokens[0].compare("breakpoints")==0) {
    // First token is "breakpoints", so pass to basic mode parser
    parseBasic(PDFFile, tokens, lineCount);
  } else if (tokens[0].compare("advanced")==0) {
    // First token is "advanced". Make sure there's no extraneous junk
    // on this line, and then call advanced mode parser
    if (tokens.size() > 1) {
      if (tokens[1].compare(0, 1, "#") != 0) {
	parseError(lineCount, linecopy, 
		   "Expected: 'breakpoints' or 'advanced'");
      }
    }
    parseAdvanced(PDFFile, lineCount);
  } else {
    // First word is not breakpoints or advanced, so bail out
    parseError(lineCount, linecopy, 
	       "Expected: 'breakpoints' or 'advanced'");
  }

  // Close file
  PDFFile.close();

  // Set up the discrete distribution picker
  if (segments.size() > 1) {
    discrete_distribution<> dist(weights.begin(), weights.end());
    disc = new variate_generator<rng_type&,
				 discrete_distribution <> >(*rng, dist);
  } else {
    disc = NULL;
  }

  // Set up the 50-50 generator
  uniform_smallint<> udist(0, 1);
  coin = new variate_generator<rng_type&,
			       uniform_smallint <> >(*rng, udist);

  // Compute the expectation value
  expectVal = 0.0;
  for (unsigned int i=0; i<segments.size(); i++)
    expectVal += weights[i]*segments[i]->expectationVal();
  expectVal = expectVal / PDFintegral;

  // Get the lower and upper limits of the PDF
  xMin = BIG;
  xMax = -BIG;
  for (unsigned int i=0; i<segments.size(); i++) {
    xMin = xMin < segments[i]->sMin() ? xMin : segments[i]->sMin();
    xMax = xMax > segments[i]->sMax() ? xMax : segments[i]->sMax();
  }

  // Initialize range restriction parameters
  xStochMin = xMin;
  xStochMax = xMax;
  range_restrict = false;
  expectVal_restrict = expectVal;
  PDFintegral_restrict = PDFintegral;
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_PDF::~slug_PDF() { 
  for (unsigned int i=0; i<segments.size(); i++)
    delete segments[i];
  if (disc_restricted != NULL)
    delete disc_restricted;
  if (disc != NULL)
    delete disc;
  delete coin;
}

////////////////////////////////////////////////////////////////////////
// Function to return integral over a finite range
////////////////////////////////////////////////////////////////////////
double
slug_PDF::integral(double a, double b) {
  double val = 0.0;
  for (unsigned int i=0; i<segments.size(); i++) {
    if (a >= segments[i]->sMax()) continue;
    if (b <= segments[i]->sMin()) continue;
    val += weights[i] * segments[i]->integral(a, b);
  }
  return val;
}

////////////////////////////////////////////////////////////////////////
// Function to return expectation value over a finite range
////////////////////////////////////////////////////////////////////////
double
slug_PDF::expectationVal(double a, double b) {
  double num = 0.0;
  double denom = 0.0;
  for (unsigned int i=0; i<segments.size(); i++) {
    if (a >= segments[i]->sMax()) continue;
    if (b <= segments[i]->sMin()) continue;
    double segInt = segments[i]->integral(a, b);
    num += weights[i] * segments[i]->expectationVal(a,b) * segInt;
    denom += weights[i] * segInt;
  }
  return num/denom;
}


////////////////////////////////////////////////////////////////////////
// Operator to return the PDF evaluated at a particular value x
////////////////////////////////////////////////////////////////////////
double
slug_PDF::operator() (const double x) {
  double val = 0.0;
  for (unsigned int i=0; i<segments.size(); i++)
    val += weights[i] * (*segments[i])(x);
  return val;
}


////////////////////////////////////////////////////////////////////////
// Method to restrict the range of the PDF from which we sample when
// producing populations
////////////////////////////////////////////////////////////////////////
void
slug_PDF::set_stoch_lim(double x_stoch_min) {
  set_stoch_lim(x_stoch_min, xMax);
}

void
slug_PDF::set_stoch_lim(double x_stoch_min, double x_stoch_max) {

  // Step 1: if we currently have a different range restriction set,
  // clear it
  if (range_restrict) {
    seg_restricted.resize(0);
    weights_restricted.resize(0);
    delete disc_restricted;
  }
  
  // Step 2: set flags and store data for new range restriction
  xStochMin = x_stoch_min;
  xStochMax = x_stoch_max;
  range_restrict = true;

  // Step 3: generate a restricted segment list
  vector<double> wgt_temp;
  for (unsigned int i=0; i<segments.size(); i++) {
    if (xStochMin >= segments[i]->sMax()) continue; // Segment is out of range
    if (xStochMax <= segments[i]->sMin()) continue; // Segment is out of range
    if (xStochMin <= segments[i]->sMin() && xStochMax >= segments[i]->sMax()) {
      // Segment entirely in range, so just copy it into the
      // restricted segment list
      seg_restricted.push_back(segments[i]);
      weights_restricted.push_back(weights[i]);
    } else {
      // Segment partly in range; push onto segment list, and compute
      // reduced weight
      seg_restricted.push_back(segments[i]);
      weights_restricted.
	push_back(weights[i] * 
		  segments[i]->integral(xStochMin, xStochMax));
    }
  }

  // Step 4: set up a new discrete picker with the new weights
  if (seg_restricted.size() > 1) {
    discrete_distribution<> dist(weights_restricted.begin(), 
				 weights_restricted.end());
    disc_restricted = new 
      variate_generator<rng_type&,
			discrete_distribution <> >(*rng, dist);
  } else {
    disc_restricted = NULL;
  }

  // Step 5: get the integral and expection value for the resticted
  // region
  PDFintegral_restrict = 0.0;
  for (unsigned int i=0; i<weights_restricted.size(); i++)
    PDFintegral_restrict += weights_restricted[i];
  expectVal_restrict = 0.0;
  for (unsigned int i=0; i<seg_restricted.size(); i++)
    expectVal_restrict += weights_restricted[i] * 
      seg_restricted[i]->expectationVal(xStochMin, xStochMax);
  expectVal_restrict = expectVal_restrict / PDFintegral_restrict;
}


////////////////////////////////////////////////////////////////////////
// Method to remove stochastic range restrictions
////////////////////////////////////////////////////////////////////////
void
slug_PDF::remove_stoch_lim() {
  if (range_restrict) {
    seg_restricted.resize(0);
    weights_restricted.resize(0);
    delete disc_restricted;
    disc_restricted = NULL;
    xStochMin = xMin;
    xStochMax = xMax;
    PDFintegral_restrict = PDFintegral;
    expectVal_restrict = expectVal;
  }
  range_restrict = false;
}


////////////////////////////////////////////////////////////////////////
// Draw function
////////////////////////////////////////////////////////////////////////
double
slug_PDF::draw() {

  // First decide which segment to draw from
  int segNum;
  if (segments.size() > 1) {
    segNum = (*disc)();
  } else {
    segNum = 0;
  }

  // Draw from that segment and return
  return(segments[segNum]->draw());
}


////////////////////////////////////////////////////////////////////////
// Draw function with stochastic range restriction
////////////////////////////////////////////////////////////////////////
double
slug_PDF::draw_restricted() {

  // If we have no restrictions, just call the regular draw function
  if (!range_restrict) return draw();

  // First decide which segment to draw from
  int segNum;
  if (seg_restricted.size() > 1) {
    segNum = (*disc_restricted)();
  } else {
    segNum = 0;
  }

  // Draw from that segment and return
  return(seg_restricted[segNum]->draw(xStochMin, xStochMax));
}


////////////////////////////////////////////////////////////////////////
// Draw function over limited range
////////////////////////////////////////////////////////////////////////
double
slug_PDF::draw(double a, double b) {

  // If there's just one segment, this is trival: just draw from it
  // with a restricted range and return the result
  if (segments.size() == 1) {
    return segments[0]->draw(a, b);
  }

  // If we're here, we need to construct temporary list of segments
  // and weights with the restricted range we've been given
  vector<slug_PDF_segment *> seg_temp;
  vector<double> wgt_temp;
  for (unsigned int i=0; i<segments.size(); i++) {
    if (a >= segments[i]->sMax()) continue; // Segment is out of range
    if (b <= segments[i]->sMin()) continue; // Segment is out of range
    if (a <= segments[i]->sMin() && b >= segments[i]->sMax()) {
      // Segment entirely in range, just copy weight
      seg_temp.push_back(segments[i]);
      wgt_temp.push_back(weights[i]);
    } else {
      // Segment partly in range; store segment and compute new weight
      seg_temp.push_back(segments[i]);
      wgt_temp.push_back(weights[i] * segments[i]->integral(a, b));
    }
  }

  // Create a new discrete distribution generator from the new weights
  discrete_distribution<> dist(wgt_temp.begin(), wgt_temp.end());
  variate_generator<rng_type&,
		    discrete_distribution <> > disc_temp(*rng, dist);

  // Draw from discrete generator
  int segNum;
  if (seg_temp.size() > 1) {
    segNum = disc_temp();
  } else {
    segNum = 0;
  }

  // Draw from that segment and return
  return seg_temp[segNum]->draw(a, b);
}

////////////////////////////////////////////////////////////////////////
// Draw population function
////////////////////////////////////////////////////////////////////////
double
slug_PDF::drawPopulation(double target, vector<double>& pop) {
  double sum = 0.0;

  // If we're only using stochasticity over a limited range, reduce
  // the target value by the fraction of the PDF that is being treated
  // stochastically
  if (range_restrict) target *= mass_frac_restrict();

  // Procedure depends on sampling method
  if ((method == STOP_NEAREST) || (method == STOP_BEFORE) ||
      (method == STOP_AFTER) || (method == STOP_50)) {

    // Sampling methods based on mass instead of number

    // Draw until we exceed target
    while (sum < target) {
      pop.push_back(draw_restricted());
      sum += pop[pop.size()-1];
    }

    // Decide what to do based on sampling method
    if (method == STOP_BEFORE) {

      // Stop before: always drop last element
      sum -= pop.back();
      pop.pop_back();

    } else if (method == STOP_NEAREST) {

      // Stop nearest: drop last element if that reduces the error
      double sum_minus = sum - pop.back();
      if (sum - target > target - sum_minus) {
	pop.pop_back();
	sum = sum_minus;
      }

    } else if (method == STOP_50) {

      // Stop 50: flip a coin to decide whether to keep or not
      if ((*coin)() == 0) {
	sum -= pop[pop.size()-1];
	pop.pop_back();
      }

    }

  } else {

    // Sampling methods based on number

    if (method == NUMBER) {

      // Draw exactly the expected number of stars
      int nExpect = (int) round(target/expectVal);
      for (int i=0; i<nExpect; i++) {
	pop.push_back(draw_restricted());
	sum += pop[pop.size()-1];
      }

    } else if (method == POISSON) {

      // Draw from a Poisson distribution to get number of stars
      double nExpect = target/expectVal;
      poisson_distribution<> pdist(nExpect);
      variate_generator<rng_type&,
			poisson_distribution <> > poisson(*rng, pdist);
      int nStar = poisson();
      for (int i=0; i<nStar; i++) {
	pop.push_back(draw_restricted());
	sum += pop[pop.size()-1];
      }

    } else if (method == SORTED_SAMPLING) {

      // Step 1: draw expected number of stars repeatedly until target
      // mass is exceeded
      while (sum < target) {
	int nExpect = (int) round((target-sum)/expectVal);
	if (nExpect == 0) nExpect = 1;
	for (int i=0; i<nExpect; i++) {
	  pop.push_back(draw_restricted());
	  sum += pop[pop.size()-1];
	}
      }

      // Step 2: remove the most massive star using a stop_nearest
      // policy
      sort(pop.begin(), pop.end());
      double sum_minus = sum - pop.back();
      if (sum - target > target - sum_minus) {
	pop.pop_back();
	sum = sum_minus;
      }
    }

  }

  return(sum);
}

////////////////////////////////////////////////////////////////////////
// Basic parser
////////////////////////////////////////////////////////////////////////
void
slug_PDF::parseBasic(ifstream& PDFFile, vector<string> firstline,
		     int& lineCount) {

  // First token of first line is breakpoints; make sure that we have
  // at least tokens total on that line
  if (firstline.size() < 3)
    parseError(lineCount, "", "Need at least two breakpoints");

  // Read the breakpoints
  unsigned int nbreak = firstline.size() - 1;
  unsigned int nsegment = nbreak - 1;
  double *breakpoints = new double[nbreak];
  for (unsigned int i=0; i<nbreak; i++) {
    try {
      breakpoints[i] = lexical_cast<double>(firstline[i+1]);
    } catch (const bad_lexical_cast& ia) {
      // If we're here, a type conversion failed
      parseError(lineCount, "", 
		 "Expected: 'breakpoints M1 M2 M3 ... MN'");
    }
  }

  // Now read the segments associated with those breakpoints. Each
  // segment should be formatted as
  // segment
  // type TYPE
  // var1 VALUE
  // var2 VALUE
  // ...
  // Where the names of var1, var2, etc. depend on the type of
  // segment. Lines of the form
  // method METHOD
  // are also allowed.
  unsigned int bptr = 0;
  bool inSegment = false;
  string line, linecopy;
  vector<string> tokens;
  while (!PDFFile.eof()) {

    // Get a line and trim leading whitespace
    getline(PDFFile, line);
    lineCount++;
    linecopy = line;
    trim(line);

    // Skip comment and blank lines
    if (line.length() == 0) continue;
    if (line.compare(0, 1, "#") == 0) continue;

    // Split line into tokens, and lowercase the first one
    split(tokens, line, is_any_of("\t ,"), token_compress_on);
    to_lower(tokens[0]);

    // Action depends on whether we're in a segment
    if (inSegment) {

      // We are reading a segment

      // Make sure this is a type specification
      if (tokens[0].compare("type") != 0)
	parseError(lineCount, linecopy, "Expected: 'type TYPE'");

      // Make sure there's no extraneous junk on the line
      if (tokens.size() > 2) {
	if (tokens[2].compare(0, 1, "#") != 0) {
	  parseError(lineCount, linecopy, "Expected: 'type TYPE'");
	}
      }

      // Read the segment type, and call the appropriate constructor
      to_lower(tokens[1]);
      slug_PDF_segment *seg = NULL;
      if (tokens[1].compare("lognormal")==0) {
	slug_PDF_lognormal *new_seg = 
	  new slug_PDF_lognormal(breakpoints[bptr], breakpoints[bptr+1]);
	seg = (slug_PDF_segment *) new_seg;
      } else if (tokens[1].compare("normal")==0) {
	slug_PDF_normal *new_seg = 
	  new slug_PDF_normal(breakpoints[bptr], breakpoints[bptr+1]);
	seg = (slug_PDF_segment *) new_seg;
      } else if (tokens[1].compare("powerlaw")==0) {
	slug_PDF_powerlaw *new_seg = 
	  new slug_PDF_powerlaw(breakpoints[bptr], breakpoints[bptr+1]);
	seg = (slug_PDF_segment *) new_seg;
      } else if (tokens[1].compare("schechter")==0) {
	slug_PDF_schechter *new_seg = 
	  new slug_PDF_schechter(breakpoints[bptr], breakpoints[bptr+1]);
	seg = (slug_PDF_segment *) new_seg;
      } else {
	string errStr("Unknown segment type ");
	errStr += tokens[1];
	parseError(lineCount, linecopy, errStr);
      }

      // Call the parser for the segment we just created to get
      // whatever data it needs
      string errMsg;
      parseStatus stat = seg->parse(PDFFile, lineCount, errMsg, *rng);
      if (stat == PARSE_ERROR)
	parseError(lineCount, "", errMsg);
      else if (stat == EOF_ERROR)
	eofError(errMsg);

      // Push this segment onto the segment vector
      segments.push_back(seg);
      weights.push_back(1.0);   // We'll fix the weights below
      bptr++;

      // We're done with this segment
      inSegment = false;

    } else {

      // If we're not in the middle of a segment, two things are
      // acceptable:
      //    segment
      // and 
      //    method METHOD
      // Check that we got one of these.
      if (tokens[0].compare("method") == 0) {

	// This is a method line, so make sure we have exactly 1
	// more non-comment token, and then read the method
	if (tokens.size() > 2) {
	  if (tokens[2].compare(0, 1, "#") != 0) {
	    parseError(lineCount, linecopy, "Expected: 'method METHOD'");
	  }
	}
	to_lower(tokens[1]);
	if (tokens[1] == "stop_nearest") {
	  method = STOP_NEAREST;
	} else if (tokens[1] == "stop_before") {
	  method = STOP_BEFORE;
	} else if (tokens[1] == "stop_after") {
	  method = STOP_AFTER;
	} else if (tokens[1] == "stop_50") {
	  method = STOP_50;
	} else if (tokens[1] == "number") {
	  method = NUMBER;
	} else if (tokens[1] == "poisson") {
	  method = POISSON;
	} else if (tokens[1] == "sorted_sampling") {
	  method = SORTED_SAMPLING;
	} else {
	  string errStr("Unknown sampling method ");
	  errStr += tokens[1];
	  parseError(lineCount, linecopy, errStr);
	}

      } else if (tokens[0].compare("segment") == 0) {

	// This is the start of a segment. Make sure there's no extra
	// junk on the line
	if (tokens.size() > 1) {
	  if (tokens[1].compare(0, 1, "#") != 0) {
	    parseError(lineCount, linecopy, "Expected: 'segment'");
	  }
	}

	// We've found a valid segment line, so set a flag
	inSegment = true;

	// Make sure we aren't trying to read too many segments
	if (bptr >= nsegment)
	  parseError(lineCount, linecopy,
		     "number of segments must equal number of breakpoints - 1");
      } else {

	// This is not a valid line
	parseError(lineCount, linecopy, 
		   "Expected 'segment' or 'method METHOD'");

      }
    }
  }

  // Make sure we have the right number of segments; if not, throw
  // error
  if (segments.size() != nsegment) {
    string errStr = "Expected " + lexical_cast<string>(nsegment) +
      " segments, found " + lexical_cast<string>(segments.size());
    parseError(lineCount, "", errStr);
  }

  // Now figure out the correct weights on all segments in order to
  // make them continuous across the breakpoints
  double cumWeight = weights[0];
  for (unsigned int i=1; i<nsegment; i++) {
    weights[i] = weights[i-1] * 
      segments[i-1]->sMaxVal() / segments[i]->sMinVal();
    cumWeight += weights[i];
  }

  // If normalizing to unity, normalize; if not, store integral under
  // function
  if (normalized) {
    for (unsigned int i=0; i<nsegment; i++)
      weights[i] /= cumWeight;
    PDFintegral = 1.0;
  } else {
    PDFintegral = cumWeight;
  }
}


////////////////////////////////////////////////////////////////////////
// Advanced parser
////////////////////////////////////////////////////////////////////////
void
slug_PDF::parseAdvanced(ifstream& PDFFile, int& lineCount) {

  // Advanced files are formatted as a series of segments, each of
  // which follows the format
  // segment
  // type TYPE
  // weight WEIGHT
  // var1 VALUE
  // var2 VALUE
  // ...
  // Where the names of var1, var2, etc. depend on the type of segment
  string line, linecopy;
  vector<string> tokens;
  bool inSegment = false;
  while (!PDFFile.eof()) {

    // Get a line and trim leading whitespace
    getline(PDFFile, line);
    lineCount++;
    linecopy = line;
    trim(line);

    // Skip comment and blank lines
    if (line.length() == 0) continue;
    if (line.compare(0, 1, "#") == 0) continue;

    // Split line into tokens, and lowercase the first one
    split(tokens, line, is_any_of("\t ,"), token_compress_on);
    to_lower(tokens[0]);

    // Action depends on whether we're in a segment
    if (inSegment) {

      // We are reading a segment, so figure out what kind

      // Make sure this is a type specification
      if (tokens[0].compare("type") != 0)
	parseError(lineCount, linecopy, "Expected: 'type TYPE'");

      // Make sure there's no extraneous junk on the line
      if (tokens.size() > 2) {
	if (tokens[2].compare(0, 1, "#") != 0) {
	  parseError(lineCount, linecopy, "Expected: 'type TYPE'");
	}
      }

      // Read the segment type, and call the appropriate constructor
      to_lower(tokens[1]);
      slug_PDF_segment *seg = NULL;
      if (tokens[1].compare("lognormal")==0) {
	slug_PDF_lognormal *new_seg = new slug_PDF_lognormal;
	seg = (slug_PDF_segment *) &new_seg;
      } else if (tokens[1].compare("normal")==0) {
	slug_PDF_normal *new_seg = new slug_PDF_normal;
	seg = (slug_PDF_segment *) &new_seg;
      } else if (tokens[1].compare("powerlaw")==0) {
	slug_PDF_powerlaw *new_seg = new slug_PDF_powerlaw;
	seg = (slug_PDF_segment *) &new_seg;
      } else if (tokens[1].compare("schechter")==0) {
	slug_PDF_schechter *new_seg = new slug_PDF_schechter;
	seg = (slug_PDF_segment *) &new_seg;
      } else {
	string errStr("Unknown segment type ");
	errStr += tokens[1];
	parseError(lineCount, linecopy, errStr);
      }

      // Call the parser for the segment we just created to get
      // whatever data it needs
      string errMsg;
      double wgt;
      parseStatus stat = seg->parse(PDFFile, lineCount, errMsg, *rng,
				    &wgt);
      if (stat == PARSE_ERROR)
	parseError(lineCount, "", errMsg);
      else if (stat == EOF_ERROR)
	eofError(errMsg);

      // Push this segment onto the segment vector
      segments.push_back(seg);
      weights.push_back(wgt);

    } else {

      // If we're not in the middle of a segment, two things are
      // acceptable:
      //    segment
      // and 
      //    method METHOD
      // Check that we got one of these.
      if (tokens[0].compare("method") == 0) {

	// This is a method line, so make sure we have exactly 1
	// more non-comment token, and then read the method
	if (tokens.size() > 2) {
	  if (tokens[2].compare(0, 1, "#") != 0) {
	    parseError(lineCount, linecopy, "Expected: 'method METHOD'");
	  }
	}
	to_lower(tokens[1]);
	if (tokens[1] == "stop_nearest") {
	  method = STOP_NEAREST;
	} else if (tokens[1] == "stop_before") {
	  method = STOP_BEFORE;
	} else if (tokens[1] == "stop_after") {
	  method = STOP_AFTER;
	} else if (tokens[1] == "stop_50") {
	  method = STOP_50;
	} else if (tokens[1] == "number") {
	  method = NUMBER;
	} else if (tokens[1] == "sorted_sampling") {
	  method = SORTED_SAMPLING;
	} else {
	  string errStr("Unknown sampling method ");
	  errStr += tokens[1];
	  parseError(lineCount, linecopy, errStr);
	}

      } else if (tokens[0].compare("segment") == 0) {

	// This is the start of a segment. Make sure there's no extra
	// junk on the line
	if (tokens.size() > 1) {
	  if (tokens[1].compare(0, 1, "#") != 0) {
	    parseError(lineCount, linecopy, "Expected: 'segment'");
	  }
	}

	// We've found a valid segment line, so set a flag
	inSegment = true;

      } else {

	// This is not a valid line
	parseError(lineCount, linecopy, 
		   "Expected 'segment' or 'method METHOD'");

      }
    }
  }

  // Make sure we got at least one segment. If not, bail out.
  if (segments.size()==0)
    eofError("Expected to find at least 1 segment.");

  // Normalize segment weights if this is a normalized PDF; store
  // integrated value if not
  double cumWeight = weights[0];
  for (unsigned int i=1; i<segments.size(); i++) {
    cumWeight += weights[i];
  }
  if (normalized) {
    for (unsigned int i=0; i<segments.size(); i++)
      weights[i] /= cumWeight;
    PDFintegral = 1.0;
  } else {
    PDFintegral = cumWeight;
  }
}


////////////////////////////////////////////////////////////////////////
// Parsing error handler
////////////////////////////////////////////////////////////////////////
void
  slug_PDF::parseError(int lineCount, string line, string message) {
  cerr << "slug error: parsing error in file " 
       << PDFFileName 
       << " on line " << lineCount;
  if (line.length() > 0) 
    cerr << ": " << endl << line << endl;
  else
    cerr << endl;
  if (message.length() > 0)
    cerr << message << endl;
  exit(1);
}


////////////////////////////////////////////////////////////////////////
// Unexpected EOF error handler
////////////////////////////////////////////////////////////////////////
void
slug_PDF::eofError(string message) {
  cerr << "slug error: unxepctedly reached end of PDF file "
       << PDFFileName << endl;
  if (message.length() > 0)
    cerr << message << endl;
  exit(1);
}
