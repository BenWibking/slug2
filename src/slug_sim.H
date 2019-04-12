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
// slug_sim class
//
// This defines the slug_sim class, which is the main driver class for
// running simulations. It is responsible for initializing
// simulations, running them, and managing the output.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_sim_H_
#define _slug_sim_H_

#include "slug.H"
#include "slug_IO.H"
#include "slug_cluster.H"
#include "slug_galaxy.H"
#include "slug_nebular.H"
#include "slug_parmParser.H"
#include "pdfs/slug_PDF.H"
#include "tracks/slug_tracks.H"
#include "yields/slug_yields.H"
#include <vector>
#ifdef ENABLE_FITS
extern "C" {
#   include "fitsio.h"
}
#endif
#ifdef ENABLE_MPI
#   include "mpi.h"
#endif

enum simType { galaxy_type, cluster_type };

class slug_sim {

public:

  // Constructor
  slug_sim(const slug_parmParser& pp_, slug_ostreams &ostreams_
#ifdef ENABLE_MPI
	   , MPI_Comm comm_
#endif
	   );

  // Destructor
  ~slug_sim();

  // Class to run a galaxy simulation
  void galaxy_sim();

  // Class to run a cluster simulation
  void cluster_sim();

private:

  // Functions to open output files and write headers of the various
  // output files
  void open_output(slug_output_files &outfiles, int chknum = -1);
  void close_output(slug_output_files &outfiles, int chknum = -1,
		    unsigned int ntrials = 0);
  void open_integrated_prop(slug_output_files &outfiles, int chknum = -1);
  void open_cluster_prop(slug_output_files &outfiles, int chknum = -1);
  void open_integrated_spec(slug_output_files &outfiles, int chknum = -1);
  void open_cluster_spec(slug_output_files &outfiles, int chknum = -1);
  void open_integrated_phot(slug_output_files &outfiles, int chknum = -1);
  void open_cluster_phot(slug_output_files &outfiles, int chknum = -1);
  void open_integrated_yield(slug_output_files &outfiles, int chknum = -1);
  void open_cluster_yield(slug_output_files &outfiles, int chknum = -1);
  void open_integrated_sn(slug_output_files &outfiles, int chknum = -1);
  void open_cluster_sn(slug_output_files &outfiles, int chknum = -1);
  void open_cluster_ew(slug_output_files &outfiles, int chknum = -1);
  
  // Function to write separators between trials to files
  void write_separator(std::ofstream& file, 
		       const unsigned int width = 80);

  // Private data to be used in the simulations
  const slug_parmParser &pp;  // Parameter parser
  rng_type *rng;              // Random number generator
  slug_tracks *tracks;        // Stellar evolution tracks
  slug_PDF *imf;              // Stellar IMF
  slug_PDF *cmf;              // Cluster mass function
  slug_PDF *clf;              // Cluster lifetime function
  slug_PDF *sfh;              // Star formation history
  slug_PDF *sfr_pdf;          // PDF of constant SFRs
  slug_PDF *out_time_pdf;     // PDF of output times
  slug_specsyn *specsyn;      // Spectral synthesizer
  slug_filter_set *filters;   // Photometric filters
  slug_line_list *lines;      // Line list
  slug_extinction *extinct;   // Extinction class
  slug_nebular *nebular;      // Nebular emission class
  slug_yields *yields;        // Class for calculating stellar yields
  slug_cluster *cluster;      // A single star cluster
  slug_galaxy *galaxy;        // A single galaxy
  simType sim_type;           // Simulation type
  outputMode out_mode;        // Output mode
  int checkpoint_ctr;         // Checkpoint counter
  bool is_imf_var = false;          //Does the IMF contain variable segments?
  std::vector<double> outTimes;     // Output times
  std::vector<double> imf_vpdraws;  //Variable parameter draws from the IMF

  // Classes to handle I/O
  slug_ostreams &ostreams;

#ifdef ENABLE_MPI
  // Information needed for MPI runs
  MPI_Comm comm;
  int rank;
#endif
};
  
#endif
// _slug_sim_H_
