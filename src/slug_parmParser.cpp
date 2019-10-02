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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include "constants.H"
#include "slug_parmParser.H"
#include "slug_IO.H"
#include "tracks/slug_tracks.H"
#ifdef ENABLE_FITS
extern "C" {
#include "fitsio.h"
}
#endif

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////

slug_parmParser::slug_parmParser(int argc, char **argv,
				 slug_ostreams &ostreams_
#ifdef ENABLE_MPI
				 , MPI_Comm comm_
#endif
				 ) :
  ostreams(ostreams_)
#ifdef ENABLE_MPI
  , comm(comm_)
#endif
{
#ifdef ENABLE_MPI
  // Get rank if working in MPI mode
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);
#endif

  // First make sure we have the right number of arguments; if not,
  // print error and exit with error
  if (argc == 1 || argc > 3) {
    ostreams.slug_err_one << "expected 1 or 2 arguments" << std::endl;
    printUsage();
    exit(1);
  }

  // Grab arguments
  string arg1(argv[1]);
  string arg2;
  if (argc == 3) arg2 = argv[2];
  
  // If we got "-h" or "--help", print usage message and then exit
  // normally
  if (!arg1.compare("-h") || !arg1.compare("--help") ||
      !arg2.compare("-h") || !arg2.compare("--help")) {
    printUsage();
    exit(0);
  }

  // See if either of our arguments was "-r" or "--restart",
  // indicating that this is a restart
  bool restart;
  if (!arg1.compare("-r") || !arg1.compare("--restart") ||
      !arg2.compare("-r") || !arg2.compare("--restart")) {
    restart = true;
  } else {
    restart = false;
  }

  // Get parameter file name
  string paramFileName;
  if (argc == 2) paramFileName = arg1;
  else {
    if (!arg1.compare("-r") || !arg1.compare("--restart")) {
      paramFileName = arg2;
    } else if (!arg2.compare("-r") || !arg2.compare("--restart")) {
      paramFileName = arg1;
    } else {
      ostreams.slug_err_one
	<< "unable to parse command line" << std::endl;
      printUsage();
      exit(1);
    }
  }

  // Start by setting all parameters to their default values
  setDefaults();

  // Try to open parameter file, and exit with error message if we
  // can't
  std::ifstream paramFile;
  paramFile.open(paramFileName.c_str(), ios::in);
  if (!paramFile.is_open()) {
    ostreams.slug_err_one << "unable to open file " 
			  << paramFileName << endl;
    exit(1);
  }

  // Parse parameter file
  parseFile(paramFile);

  // Close file
  paramFile.close();

  // Check that all parameters are set to valid values
  checkParams();

  // If this is a restart, parse the restart files to figure out how
  // many completed trials they contain, and to set the checkpoint
  // counter
  if (restart) restartSetup();
}


////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////

slug_parmParser::~slug_parmParser() { }


////////////////////////////////////////////////////////////////////////
// Method to print a usage message
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::printUsage() {
  ostreams.slug_out_one << "Usage: slug slug.param" << std::endl;
  ostreams.slug_out_one << "       slug [-r or --restart] slug.param"
			<< std::endl;
  ostreams.slug_out_one << "       slug [-h or --help]" << std::endl;
}


////////////////////////////////////////////////////////////////////////
// Method to initialize all keywords to default values; this also has
// the important side effect of specifying which type each keyword is
// supposed to have, information we will use to parse the parameter
// file.
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::setDefaults() {

  // Checkpointing data
  checkpoint_ctr = checkpoint_trials = 0;

  // Basic data
  keywords["model_name"] = "SLUG_DEF";
  keywords["out_dir"] = "";
  keywords["verbosity"] = 1;
  
  // Control flow parameters
  keywords["sim_type"] = galaxy_type;
  keywords["n_trials"] = 1;
  keywords["checkpoint_interval"] = 0;
  keywords["log_time"] = false;
  keywords["time_step"] = -constants::big;
  keywords["start_time"] = -constants::big;
  keywords["end_time"] = -constants::big;
  keywords["output_times"] = vector<double>();
  keywords["sfr"] = -constants::big;
  keywords["sfh"] = "";
  keywords["cluster_mass"] = -constants::big;
  keywords["redshift"] = 0.0;
  keywords["rng_offset"] = 0;
  keywords["save_seed"] = false;
  keywords["read_seed"] = false;
  keywords["rng_seed_file"] = "slug_rng_seed.txt";

  // Output control parameters
  keywords["out_cluster"] = 1;
  keywords["out_cluster_phot"] = 1;
  keywords["out_cluster_spec"] = 1;
  keywords["out_cluster_sn"] = 1;
  keywords["out_cluster_yield"] = 1;
  keywords["out_cluster_ew"] = 0;
#ifdef WINDS_ON
  keywords["out_cluster_winds"] = 1;
#endif
  keywords["out_integrated"] = 1;
  keywords["out_integrated_phot"] = 1;
  keywords["out_integrated_spec"] = 1;
  keywords["out_integrated_sn"] = 1;
  keywords["out_integrated_yield"] = 1;
  keywords["output_mode"] = ASCII;

  // Stellar model parameters
  path p = path("lib");
  keywords["imf"] = (p / path("imf") / path("chabrier.imf")).string();
  keywords["cmf"] = (p / path("cmf") / path("slug_default.cmf")).string();
  keywords["clf"] = (p / path("clf") / path("slug_default.clf")).string();
  keywords["tracks"] = GENEVA_2013_VVCRIT_00;
  keywords["track_dir"] = (p / path("tracks")).string();
  keywords["atmospheres"] = (p / path("atmospheres")).string();
  keywords["specsyn_mode"] = SB99;
  keywords["clust_frac"] = 1.0;
  keywords["min_stoch_mass"] = 0.0;
  keywords["metallicity"] = -constants::big;

  // Extinction keywords
  keywords["a_v"] = 0.0;
  keywords["extinction_curve"] = (p / path("extinction") /
				  path("SB_ATT_SLUG.dat")).string();
  keywords["nebular_extinction_factor"] = 1.0;

  // Nebular keywords
  keywords["compute_nebular"] = 1;
  keywords["atomic_data"] = (p / path("atomic")).string();
  keywords["nebular_no_metals"] = 0;
  keywords["nebular_den"] = 1.0e2;
  keywords["nebular_temp"] = -1.0;
  keywords["nebular_logu"] = -3.0;
  keywords["nebular_phi"] = 0.73;

  // Photometric filter keywords
  keywords["phot_bands"] = vector<string>();
  keywords["filters"] = (p / path("filters")).string();
  keywords["phot_mode"] = L_NU;

  // Yield keywords
  keywords["yield_dir"] = (p / path("yields")).string();
  keywords["yield_mode"] = SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14;
  keywords["no_decay_isotopes"] = false;
  keywords["isotopes_included"] = 0;

  // Line keywords
  keywords["line_dir"] = (p / path("lines")).string();
  keywords["spectral_lines"] = vector<string>();

  // Miscellaneous keywords
  keywords["lamers_loss"] = 0;
  keywords["lamers_t4"] = 1.9e8;
  keywords["lamers_gamma"] = 0.65;
}

////////////////////////////////////////////////////////////////////////
// Method to parse an open parameter file
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::parseFile(std::ifstream &paramFile) {
  string line;
  while (!(paramFile.eof())) {

    // Read a line
    getline(paramFile, line);

    // Trim whitespace
    trim(line);

    // If there is nothing left after trimming, or the first remaining
    // character is # indicating a comment, then continue to next line
    if (line.length() == 0) continue;
    if (line.compare(0, 1 ,"#") == 0) continue;

    // Break string up into whitespace-separated tokens; save original
    // in case we need it to print error message
    string linecopy(line);
    vector<string> tokens;
    split(tokens, linecopy, is_any_of("\t "), token_compress_on);

    // Make sure we have at least two tokens; if not, print error
    // message and exit
    if (tokens.size() < 2) parseError(line);

    // If we're here, line is in valid format; check the keyword to
    // make sure that it is known
    to_lower(tokens[0]);
    map<const std::string, slug_parmParser_keyword>::iterator it =
      keywords.find(tokens[0]);
    if (it == keywords.end()) {
      ostreams.slug_err_one << "unknown keyword "
			    << tokens[0]
			    << " on line: "
			    << line << std::endl;
      bailout(1);
    }

    // Parse keyword value depending on type; we also need to handle
    // various special cases where keywords can have multiple types,
    // or where the type of the input doesn't match the variable type
    // we eventually want to store
    unsigned int nTokExpected = 2;
    try {

      // Special cases: these are keywords where the entry the user
      // provides is not of the same "type" as the data, or where
      // multiple "types" are acceptable in the parameter file
      if (!(tokens[0].compare("sim_type"))) {

	// Acceptable values are "cluster" or "galaxy"
	to_lower(tokens[1]);
	if (tokens[1].compare("cluster") == 0) {
	  keywords["sim_type"] = cluster_type;
	  keywords["out_integrated"]
	    = keywords["out_integrated_phot"]
	    = keywords["out_integrated_spec"]
	    = keywords["out_integrated_phot"]
	    = keywords["out_integrated_sn"] = false;
	}
	else if (tokens[1].compare("galaxy") == 0)
	  keywords["sim_type"] = galaxy_type;
	else {
	  ostreams.slug_err_one << "unknown sim_type: " << line << std::endl;
	  bailout(1);
	}
	
      } else if (!(tokens[0].compare("time_step"))) {

	// Acceptable values are a number of the name of a PDF file
	try {
	  keywords["time_step"] = lexical_cast<double>(tokens[1]);
	} catch (const bad_lexical_cast& ia) {
	  keywords["time_step"] = tokens[1];
	  filepaths["output_time_dist"] = path(tokens[1]);
	}
	
      } else if (!(tokens[0].compare("output_mode"))) {

	// Acceptable values are 'ascii', 'binary' or 'bin', or 'fits'
	to_lower(tokens[1]);
	if (tokens[1].compare("ascii") == 0)
	  keywords["output_mode"] = ASCII;
	else if (tokens[1].compare("binary") == 0 ||
		 tokens[1].compare("bin") == 0)
	  keywords["output_mode"] = BINARY;
#ifdef ENABLE_FITS
	else if (tokens[1].compare("fits") == 0)
	  keywords["output_mode"] = FITS;
#endif
	else {
	  ostreams.slug_err_one << "unknown output_mode: "
				<< line << std::endl;
	  bailout(1);
	}
	
      } else if (!(tokens[0].compare("sfr"))) {

	// Acceptable values are a number, the string 'sfh', or the
	// name of a PDF file
	to_lower(tokens[1]);
	if (tokens[1].compare("sfh") == 0) {
	  keywords["constant_sfr"] = false;
	  keywords["random_sfr"] = false;
	} else {
	  try {
	    // See if the SFR is a number, indicating a constant SFR
	    keywords["sfr"] = lexical_cast<double>(tokens[1]);
	    keywords["constant_sfr"] = true;
	    keywords["random_sfr"] = false;
	  } catch (const bad_lexical_cast& ia) {
	    // SFR is not a number or "sfh", so interpret it as the
	    // name of a PDF file
	    keywords["constant_sfr"] = false;
	    keywords["random_sfr"] = true;
	    filepaths["sfr"] = path(tokens[1]);
	  }
	}

      } else if (!(tokens[0].compare("cluster_mass"))) {

	// Acceptable values are the string 'cmf' or a numerical
	// value
	to_lower(tokens[1]);
	if (tokens[1].compare("cmf") == 0) {
	  keywords["cluster_mass"] = tokens[1];
	  keywords["random_cluster_mass"] = true;
	} else {
	  keywords["cluster_mass"] = lexical_cast<double>(tokens[1]);
	  keywords["random_cluster_mass"] = false;
	}

      } else if (!(tokens[0].compare("tracks"))) {
	
	// User can specify either one of the known track sets, or can
	// give the file name of a valid track file
	string track_val = tokens[1];
	to_lower(track_val);
	if (!track_val.compare("geneva_2013_vvcrit_00")) {
	  keywords["tracks"] = GENEVA_2013_VVCRIT_00;
	} else if (!track_val.compare("geneva_2013_vvcrit_40")) {
	  keywords["tracks"] = GENEVA_2013_VVCRIT_40;
	} else if (!track_val.compare("geneva_mdot_std")) {
	  keywords["tracks"] = GENEVA_MDOT_STD;
	} else if (!track_val.compare("geneva_mdot_enhanced")) {
	  keywords["tracks"] = GENEVA_MDOT_ENHANCED;
	} else if (!track_val.compare("padova_tpagb_yes")) {
	  keywords["tracks"] = PADOVA_TPAGB_YES;
	} else if (!track_val.compare("padova_tpagb_no")) {
	  keywords["tracks"] = PADOVA_TPAGB_NO;
	} else if (!track_val.compare("mist_2016_vvcrit_00")) {
	  keywords["tracks"] = MIST_2016_VVCRIT_00;
	} else if (!track_val.compare("mist_2016_vvcrit_40")) {
	  keywords["tracks"] = MIST_2016_VVCRIT_40;
	} else if (!track_val.compare("parsec_v1.2")) {
	  keywords["tracks"] = PARSEC_V1_2;
	} else {
	  // Track is not a special value, so interpret as a file name
	  keywords["tracks"] = NO_TRACK_SET;
	  filepaths["tracks"] = path(tokens[1]);
	}
	
      } else if (!(tokens[0].compare("specsyn_mode"))) {

	// Spectral synthesis mode, specified as one of several
	// strings
	to_lower(tokens[1]);
	if (tokens[1].compare("planck") == 0)
	  keywords["specsyn_mode"] = PLANCK;
	else if (tokens[1].compare("kurucz") == 0)
	  keywords["specsyn_mode"] = KURUCZ;
	else if (tokens[1].compare("kurucz+hillier") == 0)
	  keywords["specsyn_mode"] = KURUCZ_HILLIER;
	else if (tokens[1].compare("kurucz+pauldrach") == 0)
	  keywords["specsyn_mode"] = KURUCZ_PAULDRACH;
	else if (tokens[1].compare("sb99") == 0)
	  keywords["specsyn_mode"] = SB99;
	else if (tokens[1].compare("sb99hruv") == 0)
	  keywords["specsyn_mode"] = SB99_HRUV;
	else {
	  ostreams.slug_err_one << "unknown specsyn_mode: "
				<< line << std::endl;
	  bailout(1);
	}

      } else if (!(tokens[0].compare("phot_mode"))) {

	// Photometry mode: acceptable values are l_nu, l_lambda,
	// ab, vega, and stmag
	to_lower(tokens[1]);
	if (tokens[1].compare("l_nu") == 0)
	  keywords["phot_mode"] = L_NU;
	else if (tokens[1].compare("l_lambda") == 0)
	  keywords["phot_mode"] = L_LAMBDA;
	else if (tokens[1].compare("ab") == 0)
	  keywords["phot_mode"] = AB;
	else if (tokens[1].compare("stmag") == 0)
	  keywords["phot_mode"] = STMAG;
	else if (tokens[1].compare("vega") == 0)
	  keywords["phot_mode"] = VEGA;
	else {
	  ostreams.slug_err_one << "unknown output_mode: "
				<< line << std::endl;
	  bailout(1);
	}
	
      } else if (!(tokens[0].compare("a_v"))) {

	// Extinction can be specified as a number or a PDF file name
	try {
	  // See if the A_V is a number, indicating a constant A_V
	  keywords["A_V"] = lexical_cast<double>(tokens[1]);
	  keywords["constant_A_V"] = true;
	  if (query<double>("A_V") > 0.0)
	    keywords["use_extinct"] = true;
	} catch (const bad_lexical_cast& ia) {
	  // A_V is not a number, so assume it is a distribution file
	  // name
	  keywords["A_V"] = tokens[1];
	  filepaths["A_V"] = path(tokens[1]);
	  keywords["constant_A_V"] = false;
	  keywords["use_extinct"] = true;
	}

      } else if (!(tokens[0].compare("nebular_extinction_factor"))) {

	// Nebular extinction factor: can be a number or a PDF
	try {
	  // See if A_V,neb/A_V,star is a number
	  keywords["nebular_extincton_factor"]
	    = lexical_cast<double>(tokens[1]);
	  keywords["constant_nebular_extinction_factor"] = true;
	  if (query<double>("nebular_extinction_factor") != 1.0)
	    keywords["use_nebular_extinction"] = true;
	  else
	    keywords["use_nebular_extinction"] = false;
	} catch (const bad_lexical_cast& ia) {
	  // Not a number, so assume it is a distribution file name
	  keywords["nebular_extinction_factor"] = tokens[1];
	  keywords["constant_nebular_extinction_factor"] = false;
	  keywords["use_nebular_extinction"] = true;
	  filepaths["nebular_extinction_factor"] = path(tokens[1]);
	}
	
      } else if (!(tokens[0].compare("yield_mode"))) {

	// Yield mode: list of known combinations of sources for
	// calcualting yields
	to_lower(tokens[1]);
	if (tokens[1].compare("sukhbold16") == 0)
	  keywords["yield_mode"] = SNII_SUKHBOLD16;
	else if (tokens[1].compare("karakas16+doherty14") == 0)
	  keywords["yield_mode"] = AGB_KARAKAS16_DOHERTY14;
	else if (tokens[1].compare("sukhbold16+karakas16+doherty14") == 0)
	  keywords["yield_mode"] = SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14;
	else {
	  ostreams.slug_err_one << "unknown yield_mode: "
				<< line << std::endl;
	  bailout(1);
	}
	
      } else if (!(tokens[0].compare("isotopes_included"))) {

	// Acceptable values are "union" and "intersection"
	to_lower(tokens[1]);
	if (tokens[1].compare("union") == 0)
	  keywords["isotopes_included"] = 1;
	else if (tokens[1].compare("intersection") == 0)
	  keywords["isotopes_included"] = 0;
	else {
	  ostreams.slug_err_one << "unknown isotopes_included value: "
				<< line << std::endl;
	  bailout(1);
	}
      }
      
      // The rest of these are generic cases where we process by type
      else if (int *val = boost::get<int>(&(it->second))) {

	it->second = lexical_cast<int>(tokens[1]);
	
      } else if (double *val = boost::get<double>(&(it->second))) {

	it->second = lexical_cast<double>(tokens[1]);

      } else if (string *val = boost::get<string>(&(it->second))) {

	it->second = tokens[1];

      } else if (vector<double> *val =
		 boost::get<vector<double> >(&(it->second))) {

	// Count tokens
	nTokExpected = 1;

	// For this key type, we don't know in advance how many
	// entries there will be, so parse them one at a time; store
	// them in a temporary vector, then copy them when we're done
	vector<double> t;
	for (vector<string>::size_type tokPtr = 1;
	     tokPtr < tokens.size(); tokPtr++) {

	  // Check if this is a comment; if so, stop iterating; if
	  // not, increment the number of tokens expected
	  if ((tokens[tokPtr]).compare(0, 1, "#") == 0) break;
	  nTokExpected++;

	  // This is not a comment; break up by commas
	  vector<string> tokTmp;
	  split(tokTmp, tokens[tokPtr], is_any_of(", "),
		token_compress_on);

	  // Store in temporary vector
	  for (vector<string>::size_type i = 0; i < tokTmp.size(); i++) {
	    if (tokTmp[i].length() == 0) continue;
	    t.push_back(lexical_cast<double>(tokTmp[i]));
	  }
	}

	// Store final result
	it->second = t;

      } else if (vector<string> *val =
		 boost::get<vector<string> >(&(it->second))) {

	// Count tokens
	nTokExpected = 1;

	// For this key type, we don't know in advance how many
	// entries there will be, so parse them one at a time; store
	// them in a temporary vector, then copy them when we're done
	vector<string> t;
	for (vector<string>::size_type tokPtr = 1;
	     tokPtr < tokens.size(); tokPtr++) {

	  // Check if this is a comment; if so, stop iterating; if
	  // not, increment the number of tokens expected
	  if ((tokens[tokPtr]).compare(0, 1, "#") == 0) break;
	  nTokExpected++;

	  // This is not a comment; break up by commas
	  vector<string> tokTmp;
	  split(tokTmp, tokens[tokPtr], is_any_of(", "),
		token_compress_on);

	  // Store in temporary vector
	  for (vector<string>::size_type i = 0; i < tokTmp.size(); i++) {
	    if (tokTmp[i].length() == 0) continue;
	    t.push_back(tokTmp[i]);
	  }
	}

	// Store final result
	it->second = t;
      }
      
    } catch (const bad_lexical_cast& ia) {
      // If we're here, a type conversion failed
      (void) ia; // No-op to suppress compiler warning
      parseError(line);
    }

    // If we have more than the expected number of tokens, make sure
    // that the extra ones start with #, indicating a comment
    if (tokens.size() > nTokExpected) {
      if (tokens[nTokExpected].compare(0, 1, "#")) parseError(line);
    }
  }

  // Set various "derived" keywords: these are keywords that are not
  // allowed in the input file, but are derived from ones that are
  // allowed, and which we set here to default values if they have not
  // already been set
  if (keywords.find("random_cluster_mass") == keywords.end())
    keywords["random_cluster_mass"] = true;
  if (keywords.find("constant_sfr") == keywords.end())
    keywords["constant_sfr"] = true;
  if (keywords.find("random_sfr") == keywords.end())
    keywords["random_sfr"] = false;
  if (keywords.find("use_extinct") == keywords.end())
    keywords["use_extinct"] = false;

  // Initialise file paths from keywords
  filepaths["rng_seed_file"] = path(query<string>("rng_seed_file"));
  filepaths["track_dir"] = path(query<string>("track_dir"));
  filepaths["filter_dir"] = path(query<string>("filters"));
  filepaths["yield_dir"] = path(query<string>("yield_dir"));
  filepaths["line_dir"] = path(query<string>("line_dir"));
  filepaths["imf"] = path(query<string>("imf"));
  filepaths["cmf"] = path(query<string>("cmf"));
  filepaths["clf"] = path(query<string>("clf"));
  filepaths["sfh"] = path(query<string>("sfh"));
  filepaths["atmospheres"] = path(query<string>("atmospheres"));
  filepaths["extinction_curve"] = path(query<string>("extinction_curve"));
  filepaths["atomic_data"] = path(query<string>("atomic_data"));
}


////////////////////////////////////////////////////////////////////////
// Methods to throw error and exit
////////////////////////////////////////////////////////////////////////

[[noreturn]]
void
slug_parmParser::parseError(string line) {
  ostreams.slug_err_one << "unable to parse line: "
			<< line << std::endl;
  bailout(1);
}

[[noreturn]]
void
slug_parmParser::valueError(string line) {
  ostreams.slug_err_one << "bad parameter value: "
			<< line << std::endl;
  exit(1);
}


////////////////////////////////////////////////////////////////////////
// Method to check the validity of parameters, and exit if any are
// invalid
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::checkParams() {

  // Make sure parameters have acceptable values
  if (query<int>("verbosity") > 2) valueError("verbosity must be 0, 1, or 2");
  if (query<int>("n_trials") < 1) valueError("n_trials must be >= 1");

  // Output times must be set in one of 4 valid ways: user must
  // specify either (1) output_times as a list of positive, strictly
  // increasing real numbers, (2) time_step as a string that gives the
  // name of a PDF file, (3) if log_time is false, all of start_time,
  // end_time, and time_step must be positive real numbers, (4) if
  // log_time is false, end_time and time_step must be positive real
  // numbers. Make sure we have one of these valid
  // configurations. Also, store the keyword "random_output_times" for
  // future convenience, as this is easier to query than seeing if
  // time_step returns a string or a double.
  if (query<vector<double> >("output_times").size() > 0) {
    keywords["random_output_times"] = false;
    vector<double> output_times = query<vector<double> >("output_times");
    for (vector<double>::size_type i=0; i<output_times.size()-1; i++) {
      if (output_times[i] >= output_times[i+1] ||
	  output_times[i] <= 0.0) {
	valueError("output_times must be a strictly increasing series of positive real numbers");
      }
    }
  } else if (string *val = boost::get<string>(&keywords["time_step"])) {
    keywords["random_output_times"] = true;
  } else {
    keywords["random_output_times"] = false;
    if (!query<int>("log_time") &&
	query<double>("start_time") == -constants::big)
      keywords["start_time"] = query<double>("time_step");
    if (query<double>("start_time") < 0)
      valueError("start_time must be set to a positive real number");
    if (query<double>("end_time") < 0)
      valueError("end_time must be set to a positive real number");
    if (query<double>("time_step") < 0)
      valueError("time_step must be set to a positive real number or a PDF file name");
  }

  // SFR can be set by specifying a positive constant SFR, a PDF file,
  // or a star formation history file; make sure one of these is done
  // correctly
  if (query<int>("sim_type") == galaxy_type) {
    if (!query<int>("constant_sfr") &&
	!query<int>("random_sfr") &&
	query<string>("sfh").length() == 0) {
      valueError("SFH requested, but no SFH file specified");
    } else if (query<int>("constant_sfr")) {
      if (query<double>("sfr") <= 0.0)
	valueError("SFR must be a positive number");
    }
  }

  // If we have been told to use a constant cluster mass, make sure it
  // is positive
  if (query<int>("sim_type") == cluster_type &&
      !query<int>("random_cluster_mass")) {
    if (query<double>("cluster_mass") < 0)
      valueError("cluster_mass must be either cmf or a number > 0 for cluster sim");
  }

  // Make sure clustering fraction is 0 - 1
  if ((query<double>("clust_frac") < 0 ||
       query<double>("clust_frac") > 1) &&
      query<int>("sim_type") == galaxy_type) {
    valueError("clust_frac must be in the range [0,1]");
  }

  // Make sure we have a valid, non-contradictory metallicity
  if (query<int>("tracks") == NO_TRACK_SET &&
      query<double>("metallicity") != -constants::big) {
    ostreams.slug_warn_one << "metallicity set but tracks specified "
			   << "by file name; metallicity will be set "
			   << "to value corresponding to track file, and "
			   << "input metallicity will be ignored"
			   << std::endl;
  }
  if (query<double>("metallicity") < 0) {
    if (query<double>("metallicity") == -constants::big) {
      keywords["metallicity"] = 1.0;
    } else {
      valueError("metallicity must be >= 0");
    }
  }

  // Make sure that nebular phi parameter is physically allowed
  if (query<double>("nebular_phi") < 0 ||
      query<double>("nebular_phi") > 1) {
    valueError("nebular_phi must be in the range [0,1]");
  }

  // Make sure that, if we've been told ot use the Lamers mass loss
  // model, we are in the case where it works
  if (query<int>("lamers_loss") &&
      (!query<int>("random_cluster_mass") ||
       !query<int>("random_output_times"))) {
    valueError("Lamers (2005) cluster evaporation model only avaialable with random cluster masses and output times");
  }
  if (query<int>("lamers_loss") && query<double>("lamers_t4") <= 0.0) {
    valueError("Lamers (2005) t4 parameter must be > 0");
  }
  if (query<int>("lamers_loss") && query<double>("lamers_gamma") <= 0.0) {
    valueError("Lamers (2005) gamma parameter must be > 0");
  }

  // Make sure some output has been requested
  if (!query<int>("out_cluster") &&
      !query<int>("out_cluster_phot") &&
      !query<int>("out_cluster_spec") &&
      !query<int>("out_cluster_yield") &&
      !query<int>("out_cluster_sn") &&
      !query<int>("out_cluster_ew") &&
#ifdef WINDS_ON
      !query<int>("out_cluster_winds") &&
#endif
      !query<int>("out_integrated") &&
      !query<int>("out_integrated_phot") &&
      !query<int>("out_integrated_spec") &&
      !query<int>("out_integrated_yield") &&
      !query<int>("out_integrated_sn")) {  
    valueError("no output requested");
  }

  // Make sure that, if photometry has been requested, photometric
  // bands have been specified
  if ((query<int>("out_cluster_phot") ||
       query<int>("out_integrated_phot")) && 
      (query<vector<string> >("phot_bands").size() == 0)) {
    valueError("photometry requested, but no photometric bands specified");
  }

  // Make sure that, if equivalent widths have been requested, lines
  // have been specified
  vector<string> spectral_lines = query<vector<string> >("spectral_lines");
  if (query<int>("out_cluster_ew") && spectral_lines.size()==0) {
    valueError("equivalent widths requested, but no lines specified");
  }

  // Make sure that equivalent widths are requested only in the modes
  // that support them
  if (query<int>("out_cluster_ew") &&
      query<int>("specsyn_mode") != SB99_HRUV) {
    valueError("equivalent widths only available with sb99hruv spectral synthesis mode");
  }
  if (query<int>("out_cluster_ew") &&
      query<int>("sim_type") == galaxy_type) {
    valueError("equivalent widths not yet supported in galaxy simulations");
  }
  if (query<int>("out_cluster_ew") &&
      (query<int>("output_mode") == ASCII ||
       query<int>("output_mode") == BINARY)) {
    valueError("equivalent widths not yet supported in ASCII or BINARY output modes");
  }

  // Make sure filter names are unique; if not, eliminate duplicates
  // and spit out a warning
  vector<string> phot_bands = query<vector<string> >("phot_bands");
  vector<vector<double>::size_type> duplicates;
  for (vector<double>::size_type i=0; i<phot_bands.size(); i++)
    for (vector<double>::size_type j=i+1; j<phot_bands.size(); j++)
      if (phot_bands[i] == phot_bands[j]) duplicates.push_back(j);
  vector<vector<double>::size_type>::reverse_iterator 
    rit = duplicates.rbegin();
  for ( ; rit != duplicates.rend(); ++rit) {
    vector<double>::size_type i = *rit;
    ostringstream ss;
    ss << "ignoring duplicate photometric band " << phot_bands[i];
    ostreams.slug_warn_one << ss.str() << std::endl;
    phot_bands.erase(phot_bands.begin() + i);
  }
  keywords["phot_bands"] = phot_bands;

  // Make sure lines are unique; if not, eliminate duplicates
  // and spit out a warning
  duplicates.clear();
  for (vector<double>::size_type i=0; i<spectral_lines.size(); i++) {
    for (vector<double>::size_type j=i+1; j<spectral_lines.size(); j++) {
      if (spectral_lines[i] == spectral_lines[j]) {
        duplicates.push_back(j);
      }
    }
    rit = duplicates.rbegin();
  }
  for ( ; rit != duplicates.rend(); ++rit) {
    vector<double>::size_type i = *rit;
    ostringstream ss;
    ss << "ignoring duplicate spectral line " << spectral_lines[i];
    ostreams.slug_warn_one << ss.str() << std::endl;
    spectral_lines.erase(spectral_lines.begin() + i);
  }
  keywords["spectral_lines"] = spectral_lines;

  // See if the SLUG_DIR environment variable is set, for use in
  // setting up default paths. If not, set it to current working
  // directory.
  char *slug_dir_ptr = getenv("SLUG_DIR");
  string slug_dir;
  if (slug_dir_ptr != NULL)
    slug_dir = slug_dir_ptr;
  if (slug_dir.length() == 0) slug_dir = current_path().string();
  path slug_path(slug_dir);

  // If any of the file paths we have been given are relative paths, look
  // for a file of that name in the current working directory, and, if
  // one doesn't exist, prepend the environment variable SLUG_DIR
  for (map<const string, path>::iterator it=filepaths.begin();
       it != filepaths.end(); ++it) {
    if (!it->second.is_absolute()) {
      if (!exists(it->second)) {
	it->second = (slug_path / it->second);
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Method to parse restarts
////////////////////////////////////////////////////////////////////////
void slug_parmParser::restartSetup() {

  // Get list of output file types
  vector<string> outtypes;
  if (query<int>("out_integrated") &&
      query<int>("sim_type") == galaxy_type)
    outtypes.push_back("_integrated_prop");
  if (query<int>("out_integrated_spec") &&
      query<int>("sim_type") == galaxy_type)
    outtypes.push_back("_integrated_spec");
  if (query<int>("out_integrated_phot") &&
      query<int>("sim_type") == galaxy_type)
    outtypes.push_back("_integrated_phot");
  if (query<int>("out_integrated_yield") &&
      query<int>("sim_type") == galaxy_type)
    outtypes.push_back("_integrated_yield");
  if (query<int>("out_integrated_sn") &&
      query<int>("sim_type") == galaxy_type)
    outtypes.push_back("_integrated_sn");
  if (query<int>("out_cluster"))
    outtypes.push_back("_cluster_prop");
  if (query<int>("out_cluster_spec"))
    outtypes.push_back("_cluster_spec");
  if (query<int>("out_cluster_phot"))
    outtypes.push_back("_cluster_phot");
  if (query<int>("out_cluster_sn"))
    outtypes.push_back("_cluster_sn");
  if (query<int>("out_cluster_yield"))
    outtypes.push_back("_cluster_yield");
  if (query<int>("out_cluster_ew"))
    outtypes.push_back("_cluster_ew");
#ifdef WINDS_ON
  if (query<int>("out_cluster_winds"))
    outtypes.push_back("_cluster_winds");
#endif
  
  // Get parallel rank indicator
  string par_str;
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    par_str = ss.str();
  }
#endif

  // Get extension
  string ext;
  if (query<int>("output_mode") == ASCII) ext = ".txt";
  else if (query<int>("output_mode") == BINARY) ext = ".bin";
#ifdef ENABLE_FITS
  else if (query<int>("output_mode") == FITS) ext = ".fits";
#endif

  // Now loop through checkpoints, seeing which ones exist, and
  // counting the trials they contain; if this is an MPI calculation,
  // we need to loop over processor numbers
  vector<unsigned int> trials_ctr(outtypes.size());
  while (true) {

    // Loop over file types
    bool checkpoint_valid = true;
    for (vector<int>::size_type i=0; i<outtypes.size(); i++) {

      // Construct checkpoint file name
      ostringstream ss;
      ss << "_chk" << setfill('0') << setw(4) << checkpoint_ctr;
      string fname = query<string>("model_name")
	+ par_str + ss.str() + outtypes[i] + ext;
      path full_path = path(query<string>("out_dir")) / fname;

      // Try to open file and read number of trials from it; bail if
      // any of this fails
      if (query<int>("output_mode") == ASCII) {
	
	// Try to open
	std::ifstream checkpoint_file;
	checkpoint_file.open(full_path.c_str(), ios::in);
	if (!checkpoint_file.is_open()) {
	  checkpoint_valid = false;
	  break;
	}
	
	// Try to read a line
	if (checkpoint_file.eof()) {
	  checkpoint_valid = false;
	  checkpoint_file.close();
	  break;
	}
	string line;
	getline(checkpoint_file, line);

	// Parse the line to get number of trials in file
	vector<string> tokens;
	split(tokens, line, is_any_of("="), token_compress_on);
	if (tokens.size() != 2) {
	  checkpoint_valid = false;
	  break;
	}
	trim(tokens[0]);
	if (tokens[0].compare("N_Trials") != 0) {
	  checkpoint_valid = false;
	  break;
	}
	trim(tokens[1]);
	trials_ctr[i] = lexical_cast<unsigned int>(tokens[1]);

	// Close
	checkpoint_file.close();
	
      } else if (query<int>("output_mode") == BINARY) {
	
	// Try to open
	std::ifstream checkpoint_file;
	checkpoint_file.open(full_path.c_str(), ios::in | ios::binary);
	if (!checkpoint_file.is_open()) {
	  checkpoint_valid = false;
	  break;
	}
	
	// Try to read an unsigned int from file
	if (checkpoint_file.eof()) {
	  checkpoint_valid = false;
	  break;
	}
	unsigned int trials_in_file;
	checkpoint_file.read((char *) &trials_in_file,
			     sizeof trials_in_file);
	if (!(checkpoint_file.good())) {
	  checkpoint_valid = false;
	  checkpoint_file.close();
	  break;
	}
	trials_ctr[i] = trials_in_file;

	// Close
	checkpoint_file.close();

      }
#ifdef ENABLE_FITS
      else if (query<int>("output_mode") == FITS) {

	// Try to open file
	fitsfile *checkpoint_file;
	int fits_status = 0;
	fits_open_table(&checkpoint_file, full_path.c_str(),
			READONLY, &fits_status);
	if (fits_status != 0) {
	  checkpoint_valid = false;
	  break;
	}

	// Read the number of trials
	unsigned int trials_in_file;
	char comment[100];
	fits_read_key(checkpoint_file, TUINT, "N_Trials",
		      &trials_in_file, comment, &fits_status);
	if (fits_status != 0) {
	  checkpoint_valid = false;
	  fits_close_file(checkpoint_file, &fits_status);
	  break;
	}
        trials_ctr[i] = trials_in_file;

	// Close
	fits_close_file(checkpoint_file, &fits_status);
      }
#endif
    }

    // If we failed at any point, bail out
    if (!checkpoint_valid) break;

    // Check to make sure that all checkpoint files say they have the
    // same number of trials
    for (vector<int>::size_type i=1; i<trials_ctr.size(); i++) {
      if (trials_ctr[i] != trials_ctr[0]) {
	checkpoint_valid = false;
	break;
      }
    }

    // If we made it to here, this is a valid checkpoint. Add the
    // number of trials it contains to our running count, and
    // increment the checkpoint counter
    checkpoint_ctr++;
    checkpoint_trials += trials_ctr[0];
  }

  // If we are in MPI mode, we need to sum the number of completed
  // trials and files over all processors
#ifdef ENABLE_MPI
  unsigned int global_files, global_trials;
  MPI_Allreduce(&checkpoint_ctr, &global_files, 1, MPI_UNSIGNED_LONG,
		MPI_SUM, comm);
  MPI_Allreduce(&checkpoint_trials, &global_trials, 1, MPI_UNSIGNED_LONG,
		MPI_SUM, comm);
  checkpoint_trials = global_trials;
#endif

  // Print if verbose
  if (query<int>("verbosity") > 0) {
    ostreams.slug_out_one
      << "found "
#ifdef ENABLE_MPI
      << global_files
#else
      << checkpoint_ctr
#endif
      << " checkpoint "
      << "files containing " << checkpoint_trials << " trials"
      << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////
// Method to write parameters to a file
////////////////////////////////////////////////////////////////////////

void
slug_parmParser::writeParams() const {

#ifdef ENABLE_MPI
  if (rank != 0) return;
#endif

  // Form output file name
  string fname(query<string>("model_name") + "_summary.txt");
  path full_path(query<string>("out_dir"));
  full_path /= fname;

  // Open file for output
  std::ofstream paramFile;
  paramFile.open(full_path.c_str(), ios::out);
  if (!paramFile.is_open()) {
    ostreams.slug_err_one
      << "unable to open parameter summmary file " 
      << full_path.string() << std::endl;
    bailout(1);
  }

  // Write parameters to file
  paramFile << "SLUG WAS RUN WITH THE FOLLOWING PARAMETERS" << endl;
  paramFile << "model_name           " << query<string>("model_name") << endl;
  paramFile << "out_dir              " << query<string>("out_dir") << endl;
  paramFile << "sim_type             ";
  if (query<int>("sim_type") == galaxy_type)
    paramFile << "galaxy" << endl;
  else
    paramFile << "cluster" << endl;
  paramFile << "n_trials             " << query<int>("n_trials") << endl;
  if (query<vector<double> >("output_times").size() > 0) {
    paramFile << "output_times        ";
    vector<double> output_times = query<vector<double> >("output_times");
    for (vector<double>::size_type i=0; i<output_times.size(); i++)
      paramFile << " " << output_times[i];
    paramFile << endl;
  } else if (query<int>("random_output_times")) {
    paramFile << "time_step            "
	      << query<string>("time_step") << endl;
  } else {
    paramFile << "time_step            " << query<double>("time_step") << endl;
    paramFile << "start_time           " << query<double>("start_time") << endl;
    paramFile << "end_time             " << query<double>("end_time") << endl;
  }
  if (query<int>("sim_type") == galaxy_type) {
    if (query<int>("constant_sfr"))
      paramFile << "SFR                  " << query<double>("sfr") << endl;
    else if (query<int>("random_sfr"))
      paramFile << "SFR                  " << fpath("sfr") << endl;
    else
      paramFile << "SFH                  " << fpath("sfh") << endl;
  }
  paramFile << "IMF                  " << fpath("imf") << endl;
  if (query<int>("sim_type") == cluster_type ||
      query<int>("random_cluster_mass"))
    paramFile << "CMF                  " << fpath("cmf") << endl;
  if (query<int>("sim_type") == cluster_type &&
      !query<int>("random_cluster_mass"))
    paramFile << "cluster_mass         " << query<double>("cluster_mass")
	      << endl;
  paramFile << "CLF                  " << fpath("clf") << endl;
  switch (static_cast<trackSet>(query<int>("tracks"))) {
  case NO_TRACK_SET: {
    paramFile << "tracks               " << fpath("tracks") << endl;
    break;
  }
  case GENEVA_2013_VVCRIT_00: {
    paramFile << "tracks               geneva_2013_vvcrit_00" << endl;
    break;
  }
  case GENEVA_2013_VVCRIT_40: {
    paramFile << "tracks               geneva_2013_vvcrit_40" << endl;
    break;
  }
  case GENEVA_MDOT_STD: {
    paramFile << "tracks               geneva_mdot_std" << endl;
    break;
  }
  case GENEVA_MDOT_ENHANCED: {
    paramFile << "tracks               geneva_mdot_enhanced" << endl;
    break;
  }
  case PADOVA_TPAGB_YES: {
    paramFile << "tracks               padova_tp_agb_yes" << endl;
    break;
  }
  case PADOVA_TPAGB_NO: {
    paramFile << "tracks               padova_tp_agb_no" << endl;
    break;
  }
  case MIST_2016_VVCRIT_00: {
    paramFile << "tracks               mist_2016_vvcrit_00" << endl;
    break;
  }
  case MIST_2016_VVCRIT_40: {
    paramFile << "tracks               mist_2016_vvcrit_40" << endl;
    break;
  }
  case PARSEC_V1_2: {
    paramFile << "tracks               parsec_v1.2" << endl;
    break;
  }
  }
  paramFile << "atmospheres          " << fpath("atmospheres") << endl;
  paramFile << "yield_dir            " << fpath("yield_dir") << endl;
  paramFile << "min_stoch_mass       " << query<double>("min_stoch_mass")
	    << endl;
  paramFile << "redshift             " << query<double>("redshift") << endl;
  if (query<double>("metallicity") != -constants::big)
    paramFile << "metallicity          "
	      << query<double>("metallicity") << endl;
  paramFile << "specsyn_mode         ";
  if (query<int>("specsyn_mode") == PLANCK) {
    paramFile << "planck" << endl;
  } else if (query<int>("specsyn_mode") == KURUCZ) {
    paramFile << "kurucz" << endl;
  } else if (query<int>("specsyn_mode") == KURUCZ_HILLIER) {
    paramFile << "kurucz+hillier" << endl;
  } else if (query<int>("specsyn_mode") == SB99) {
    paramFile << "sb99" << endl;
  } else if (query<int>("specsyn_mode") == SB99_HRUV) {    
    paramFile << "sb99hruv" << endl; 
  }  
  paramFile << "yield_mode           ";
  switch (query<int>("yield_mode")) {
  case SNII_SUKHBOLD16: {
    paramFile << "SNII: Sukhbold+16; AGB: none" << endl;
    break;
  }
  case AGB_KARAKAS16_DOHERTY14: {
    paramFile << "SNII: none; AGB: Karakas & Lugaro 2016 + Doherty+ 2014" << endl;
    break;
  }
  case SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14: {
    paramFile << "SNII: Sukhbold+16; AGB: Karakas & Lugaro 2016 + Doherty+ 2014" << endl;
    break;
  }
  }
  if (query<int>("use_extinct")) {
    if (!query<int>("constant_A_V")) {
      paramFile << "A_V                  " << fpath("A_V") << endl;
      paramFile << "extinction_curve     " << fpath("extinction_curve") << endl;
    } else {
      paramFile << "A_V                  " << query<double>("A_V") << endl;
      paramFile << "extinction_curve     " << fpath("extinction_curve") << endl;
    }
  }
  if (query<int>("compute_nebular")) {
    paramFile << "nebular_emission     " << "yes" << endl;
    paramFile << "nebular_density      " << query<double>("nebular_den")
	      << endl;
    paramFile << "nebular_temperature  " << query<double>("nebular_temp")
	      << endl;
    paramFile << "nebular_phi          " << query<double>("nebular_phi")
	      << endl;
    paramFile << "nebular_logU         " << query<double>("nebular_logu")
	      << endl;
    paramFile << "nebular_no_metals    " << query<int>("nebular_no_metals")
	      << endl;
  } else {
    paramFile << "nebular_emission     " << "no" << endl;
  }
  if (query<int>("lamers_loss")) {
    paramFile << "lamers_evaporation  " << "yes" << endl;
    paramFile << "lamers_t4           " << query<double>("lamers_t4")
	      << endl;
    paramFile << "lamers_gamma        " << query<double>("lamers_gamma")
	      << endl;
  }
  if (query<int>("sim_type") == galaxy_type)
    paramFile << "clust_frac           " << query<double>("clust_frac") << endl;
  if (query<int>("out_cluster_phot") || query<int>("out_integrated_phot")) {
    paramFile << "phot_mode            ";
    switch (query<int>("phot_mode")) {
    case L_NU:     { paramFile << "L_nu" << endl; break; }
    case L_LAMBDA: { paramFile << "L_lambda" << endl; break; }
    case AB:       { paramFile << "AB" << endl; break; }
    case STMAG:    { paramFile << "STMAG" << endl; break; }
    case VEGA:     { paramFile << "Vega" << endl; break; }
    }
  }
  paramFile << "out_cluster          "
	    << query<int>("out_cluster") << endl;
  paramFile << "out_cluster_phot     "
	    << query<int>("out_cluster_phot") << endl;
  paramFile << "out_cluster_spec     "
	    << query<int>("out_cluster_spec") << endl;
  paramFile << "out_cluster_sn       "
	    << query<int>("out_cluster_sn") << endl;
  paramFile << "out_cluster_yield    "
	    << query<int>("out_cluster_yield") << endl;
  paramFile << "out_cluster_ew       "
	    << query<int>("out_cluster_ew")   << endl;
#ifdef WINDS_ON
  paramFile << "out_cluster_winds    "
	    << query<int>("out_cluster_winds")   << endl;
#endif
  if (query<int>("sim_type") == galaxy_type) {
    paramFile << "out_integrated       "
	      << query<int>("out_integrated") << endl;
    paramFile << "out_integrated_phot  "
	      << query<int>("out_integrated_phot") << endl;
    paramFile << "out_integrated_spec  "
	      << query<int>("out_integrated_spec") << endl;
    paramFile << "out_integrated_sn    "
	      << query<int>("out_integrated_sn") << endl;
    paramFile << "out_integrated_yield "
	      << query<int>("out_integrated_yield") << endl;
  }
  const vector<string>& phot_bands = query<vector<string> >("phot_bands");
  if (phot_bands.size() > 0) {
    paramFile << "phot_bands           ";
    for (vector<string>::size_type i=0; i<phot_bands.size(); i++) {
      paramFile << phot_bands[i];
      if (i < phot_bands.size()-1) paramFile << ", ";
    }
    paramFile << endl;
  }
  const vector<string>& spectral_lines =
    query<vector<string> >("spectral_lines");
  if (spectral_lines.size() > 0) {
    paramFile << "spectral_lines        ";
    for (vector<string>::size_type i=0; i<spectral_lines.size(); i++) {
      paramFile << spectral_lines[i];
      if (i < spectral_lines.size()-1) paramFile << ", ";
    }
  }
  if (query<int>("output_mode") == BINARY)
    paramFile << "output_mode          binary" << endl;
  else if (query<int>("output_mode") == ASCII)
    paramFile << "output_mode          ascii" << endl;
#ifdef ENABLE_FITS
  else
    paramFile << "output_mode          fits" << endl;
#endif

  // Close
  paramFile.close();
}




