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
#include "slug_MPI.H"
#include "slug_IO.H"
#include "slug_sim.H"
#include "pdfs/slug_PDF_powerlaw.H"
#include "specsyn/slug_specsyn_hillier.H"
#include "specsyn/slug_specsyn_kurucz.H"
#include "specsyn/slug_specsyn_pauldrach.H"
#include "specsyn/slug_specsyn_planck.H"
#include "specsyn/slug_specsyn_sb99.H"
#include "specsyn/slug_specsyn_sb99hruv.H"
#include "tracks/slug_tracks_parsec.H"
#include "tracks/slug_tracks_mist.H"
#include "tracks/slug_tracks_sb99.H"
#include "yields/slug_yields_multiple.H"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "fcntl.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_sim::slug_sim(const slug_parmParser& pp_, slug_ostreams &ostreams_
#ifdef ENABLE_MPI
		   , MPI_Comm comm_
#endif
		   ) : pp(pp_),
		       ostreams(ostreams_)
#ifdef ENABLE_MPI
		       , comm(comm_)
#endif
{

  // If running in MPI mode, record our rank
#ifdef ENABLE_MPI
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);
  else rank = 0;
#endif
  
  // Either read a random seed from a file, or generate one
  unsigned int seed;
  if (pp.query<int>("read_seed")) {

    // Read seed from file; for MPI runs, only root processor does this
#ifdef ENABLE_MPI
    if (rank == 0) {
#endif
      std::ifstream seed_file;
      string seed_file_name(pp.fpath("rng_seed_file"));
      if (pp.query<int>("rng_offset") != 0) {
	stringstream ss;
	ss << seed_file_name << "_off_" << pp.query<int>("rng_offset");
	seed_file_name = ss.str();
      }
      seed_file.open(seed_file_name.c_str(), ios::in);
    
      // Check if the file exists
      if (!seed_file) {
	ostreams.slug_err_one << "can't open RNG seed file "
			      << seed_file_name << endl;
	exit(1);
      }
      seed_file >> seed;
      seed_file.close();
#ifdef ENABLE_MPI
    }

    // On MPI runs, the root processor broadcasts the seed value to
    // all other processors, which then add an offset equal to their
    // rank
    if (comm != MPI_COMM_NULL)
      MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, comm);
    seed += rank;
#endif
      
  } else {

    // Get a random see from /dev/random if possible; this will give
    // us good seeds even in MPI mode
    int fn;
    bool rand_set = false;
    fn = open("/dev/random", O_RDONLY);
    if (fn != -1) {
      rand_set = (read(fn, &seed, 4) == 4); // True if read succeeds
      close(fn);
    }
    if (!rand_set) {
      // Failed to set from /dev/random; seed using system time instead.
      // If we're running in MPI mode, we need to avoid having the
      // system time be the same for all MPI threads. This code hashes
      // the sytem time with the process ID; it comes from
      // https://arxiv.org/pdf/1005.4117
#ifdef ENABLE_MPI
      long s = time(0);
      long pid = getpid();
      seed = static_cast<unsigned int>
	(abs(((s*181)*((pid-83)*359))%104729));
#else
      seed = static_cast<unsigned int>(time(0));
#endif
      
    }
    // Add offset if requested; this probably isn't necessary if
    // /dev/random worked, but do it anyway in case it failed.
    seed += pp.query<int>("rng_offset");

    // Save the rng seed if requested; only root rank does this on MPI
    // jobs
#ifdef ENABLE_MPI
    if (rank == 0) {
#endif
      if (pp.query<int>("save_seed")) {
	std::ofstream seed_file;
	string seed_file_name(pp.fpath("rng_seed_file"));
	if (pp.query<int>("rng_offset") != 0) {
	  stringstream ss;
	  ss << seed_file_name << "_off_" << pp.query<int>("rng_offset");
	  seed_file_name = ss.str();
	}
	seed_file.open(seed_file_name.c_str(), ios::out);
	seed_file << seed;
	seed_file.close();
      }
    }
#ifdef ENABLE_MPI
  }
#endif

  // Set up the random number generator
  rng = new rng_type(seed);

  // Warm up the rng by drawing 1000 random numbers. This is another
  // safety measure to avoid getting correlated sequences of random
  // numbers if we're running in parallel.
  boost::random::uniform_int_distribution<> six_sided_die(1,6);
  for (int i=0; i<1000; i++) six_sided_die(*rng);
  
  // Set up the time stepping; these can be specified in multiple
  // ways, which we have to check for
  if (pp.query<vector<double> >("output_times").size() > 0) {
    // User has specified exact output times
    outTimes = pp.query<vector<double> >("output_times");
    out_time_pdf = nullptr;
  } else if (pp.query<int>("random_output_times")) {
    // User has specified a PDF of output times
    out_time_pdf = new slug_PDF(pp.fpath("output_time_dist"),
				rng, ostreams);
  } else {
    // User has specified start time, end time, and time step
    double t = pp.query<double>("start_time");
    double dt = pp.query<double>("time_step");
    double end_time = pp.query<double>("end_time");
    while (t <= end_time) {
      outTimes.push_back(t);
      if (pp.query<int>("log_time")) t *= pow(10., dt);
      else t += dt;
    }
    out_time_pdf = nullptr;
  }

  // Set up the photometric filters
  const vector<string>& phot_bands = pp.query<vector<string> >("phot_bands");
  if (phot_bands.size() > 0) {
    if (pp.query<int>("verbosity") > 1)
      ostreams.slug_out_one << "reading filters" << std::endl;
    filters = new
      slug_filter_set(phot_bands, 
		      pp.fpath("filter_dir"), 
		      static_cast<photMode>(pp.query<int>("phot_mode")),
		      ostreams,
		      pp.fpath("atmospheres"));
  } else {
    filters = nullptr;
  }

  // Read the tracks
  if (pp.query<int>("verbosity") > 1)
    ostreams.slug_out_one << "reading tracks" << std::endl;
  switch (static_cast<trackSet>(pp.query<int>("tracks"))) {
  case NO_TRACK_SET: {
    // User has manually specified the file name; decide if it is a
    // starburst99 or MIST file based on its extension; note that
    // MIST files are only usable if we have FITS capability
    string fname(pp.fpath("tracks"));
    size_t pos = fname.find_last_of(".");
    bool mist_ext;
    if (pos == string::npos) {
      mist_ext = false;
    } else if (fname.substr(pos).compare(".gz") == 0) {
      mist_ext = true;
    } else {
      mist_ext = false;
    }
    if (mist_ext) {
#ifdef ENABLE_FITS
      tracks = (slug_tracks *)
	new slug_tracks_mist(pp.fpath("tracks"), ostreams);
#else
      ostreams.slug_err_one
        << "MIST gzip'ed track files only available if "
        << "SLUG was compiled with ENABLE_FITS option"
        << endl;
      exit(1);
#endif
    } else {
      tracks = (slug_tracks *)
	new slug_tracks_sb99(pp.fpath("tracks"), ostreams);
    }
    break;
  }
  case GENEVA_2013_VVCRIT_00:
  case GENEVA_2013_VVCRIT_40:
  case GENEVA_MDOT_STD:
  case GENEVA_MDOT_ENHANCED:
  case PADOVA_TPAGB_YES:
  case PADOVA_TPAGB_NO: {
    // These are starburst99 track sets
    tracks = (slug_tracks *)
      new slug_tracks_sb99(static_cast<trackSet>(pp.query<int>("tracks")),
			   pp.query<double>("metallicity"),
			   pp.fpath("track_dir"), ostreams);
    break;
  }
#ifdef ENABLE_FITS
  case MIST_2016_VVCRIT_00:
  case MIST_2016_VVCRIT_40: {
    // These are MIST track sets
    tracks = (slug_tracks *)
      new slug_tracks_mist(static_cast<trackSet>(pp.query<int>("tracks")),
			   pp.query<double>("metallicity"),
			   pp.fpath("track_dir"), ostreams);
    break;
  }
  case PARSEC_V1_2: {
    // PARSEC tracks
    tracks = (slug_tracks *)
      new slug_tracks_parsec(pp.query<double>("metallicity"),
			     pp.fpath("track_dir"), ostreams);
    break;
  }
#endif
  }

  // If we're computing yields, set up yield tables
  if (pp.query<int>("out_cluster_yield") ||
      pp.query<int>("out_integrated_yield") ||
      pp.query<int>("out_cluster_sn") ||
      pp.query<int>("out_integrated_sn")) {
    if (pp.query<int>("verbosity") > 1)
      ostreams.slug_out_one << "reading yield tables" << std::endl;
    yields = (slug_yields *) new
      slug_yields_multiple(pp.fpath("yield_dir"),
			   static_cast<yieldMode>(pp.query<int>("yield_mode")),
			   ((slug_tracks_2d *) tracks)->get_metallicity(),
			   ostreams,
			   pp.query<int>("no_decay_isotopes"),
			   pp.query<int>("isotopes_included"));
  } else {
    yields = nullptr;
  }
  
  // Set up the IMF, including the limts on its stochasticity
  imf = new slug_PDF(pp.fpath("imf"), rng, ostreams);
  imf->set_stoch_lim(pp.query<double>("min_stoch_mass"));

  // Compare IMF and tracks, and issue warning if IMF extends outside
  // range of tracks
  if (imf->get_xMin() < tracks->min_mass()*(1.0-1.0e-10)) {
    ostreams.slug_warn_one
      << "minimum IMF mass " << imf->get_xMin() 
      << " Msun < minimum evolution track mass " << tracks->min_mass()
      << " Msun. Calculation will proceed, but stars with mass "
      << imf->get_xMin() << " Msun to " << tracks->min_mass()
      << " Msun will be treated as having zero luminosity." << std::endl;
  }
  if (imf->get_xMax() > tracks->max_mass()*(1.0+1.0e-10))
    ostreams.slug_warn_one
      << "maximum IMF mass " << imf->get_xMax() 
      << " Msun > maximum evolution track mass " << tracks->max_mass()
      << " Msun. Calculation will proceed, but stars with mass "
      << tracks->max_mass() << " Msun to " << imf->get_xMax()
      << " Msun will be treated as having zero luminosity." << std::endl;
  
  // Ditto for IMF and yield table
  if (yields) {
    if (imf->get_xMin() < yields->min_mass()*(1.0-1.0e-10))
      ostreams.slug_warn_one
	<< "minimum IMF mass " << imf->get_xMin() 
	<< " Msun < minimum yield table mass " << yields->min_mass()
	<< " Msun. Calculation will proceed, but stars with mass "
	<< imf->get_xMin() << " Msun to " << yields->min_mass()
	<< " Msun will be treated as having zero yield." << endl;
    if (imf->get_xMax() > yields->max_mass()*(1.0+1.0e-10))
      ostreams.slug_warn_one
	<< "maximum IMF mass " << imf->get_xMax() 
	<< " Msun > maximum yield table mass " << yields->max_mass()
	<< " Msun. Calculation will proceed, but stars with mass "
	<< yields->max_mass() << " Msun to " << imf->get_xMax()
	<< " Msun will be treated as having zero yield." << endl;
  }
  
  // Check for variable segments & initialise them
  is_imf_var = imf->init_vsegs();	
	
  // Check for variable segments & have a first draw
  if (is_imf_var == true) {
    // Draw new values for variable parameters
    // Update IMF segments and recompute weights 
    imf_vpdraws = imf->vseg_draw();
    // Reset range restrictions
    imf->set_stoch_lim(pp.query<double>("min_stoch_mass"));
  }

  // Record whether we're running a galaxy simulation or a cluster
  // simulation
  sim_type = static_cast<simType>(pp.query<int>("sim_type"));

  // Set the cluster lifetime function
  clf = new slug_PDF(pp.fpath("clf"), rng, ostreams);

  // Set the cluster mass function
  if (sim_type == galaxy_type || pp.query<int>("random_cluster_mass"))
    cmf = new slug_PDF(pp.fpath("cmf"), rng, ostreams);
  else
    cmf = nullptr;

  // Set the star formation history
  if (sim_type == cluster_type) {
    sfh = nullptr;
    sfr_pdf = nullptr;
  } else {
    if (pp.query<int>("constant_sfr")) {
      // SFR is constant, so create a powerlaw segment of slope 0 with
      // the correct normalization
      slug_PDF_powerlaw *sfh_segment = 
	new slug_PDF_powerlaw(0.0, outTimes.back(), 0.0, rng, ostreams);
      sfh = new slug_PDF(sfh_segment, rng, ostreams,
			 outTimes.back()*pp.query<double>("sfr"));
      sfr_pdf = nullptr;
    } else if (pp.query<int>("random_sfr")) {
      // SFR is to be drawn from a PDF, so read the PDF, and
      // initialize the SFH from it
      sfr_pdf = new slug_PDF(pp.fpath("sfr"), rng, ostreams, false);
      slug_PDF_powerlaw *sfh_segment = 
	new slug_PDF_powerlaw(0.0, outTimes.back(), 0.0, rng, ostreams);
      sfh = new slug_PDF(sfh_segment, rng, ostreams,
			 outTimes.back()*sfr_pdf->draw());
    } else {
      // SFR is not constant, so read SFH from file
      sfh = new slug_PDF(pp.fpath("sfh"), rng, ostreams, false);
      sfr_pdf = nullptr;
    }
  }

  // Initialize the spectral synthesizer
  if (pp.query<int>("verbosity") > 1)
    ostreams.slug_out_one << "reading atmospheres" << std::endl;
  double z = pp.query<double>("redshift");
  if (pp.query<int>("specsyn_mode") == PLANCK) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_planck(tracks, imf, sfh, ostreams, z));
  } else if (pp.query<int>("specsyn_mode") == KURUCZ) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_kurucz(pp.fpath("atmospheres"), tracks, 
			       imf, sfh, ostreams, z));
  } else if (pp.query<int>("specsyn_mode") == KURUCZ_HILLIER) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_hillier(pp.fpath("atmospheres"), tracks, 
				imf, sfh, ostreams, z));
  } else if (pp.query<int>("specsyn_mode") == KURUCZ_PAULDRACH) {
    specsyn = static_cast<slug_specsyn *>
      (new slug_specsyn_pauldrach(pp.fpath("atmospheres"), tracks, 
				  imf, sfh, ostreams, z));
  } else if (pp.query<int>("specsyn_mode") == SB99) {
    specsyn = static_cast<slug_specsyn *> 
      (new slug_specsyn_sb99(pp.fpath("atmospheres"), tracks,
			     imf, sfh, ostreams, z));
  }
  else if  (pp.query<int>("specsyn_mode") == SB99_HRUV) {
    specsyn = static_cast<slug_specsyn *> 
      (new slug_specsyn_sb99hruv(pp.fpath("atmospheres"), tracks,
				 imf, sfh, ostreams, z));   
  }

  // Load line list for equivalent width calculations
  if (pp.query<int>("out_cluster_ew")) {
    if (pp.query<int>("verbosity") > 1) {
      ostreams.slug_out_one << "reading line list" << std::endl;
    }    
    lines = new slug_line_list(pp.query<vector<string> >("spectral_lines"),
			       pp.fpath("line_dir"),
			       ostreams);
  } else {
    lines = nullptr;
  }
  
  // If using nebular emission, initialize the computation of that
  if (pp.query<int>("compute_nebular")) {
    if (pp.query<int>("tracks") == NO_TRACK_SET) {
      // User has manually specified a file name
      nebular = new slug_nebular(pp.fpath("atomic_data"),
				 specsyn->lambda(true),
				 pp.fpath("tracks"),
				 ostreams,
				 pp.query<double>("nebular_den"),
				 pp.query<double>("nebular_temp"),
				 pp.query<double>("nebular_logu"),
				 pp.query<double>("nebular_phi"),
				 pp.query<double>("redshift"),
				 pp.query<int>("nebular_no_metals"));
    } else {
      // User has specified a track set
      nebular
	= new slug_nebular(pp.fpath("atomic_data"),
			   specsyn->lambda(true),
			   ((slug_tracks_2d *) tracks)->get_metallicity(),
			   ((slug_tracks_2d *) tracks)->trackset_filenames(),
			   ((slug_tracks_2d *) tracks)->trackset_metallicities(),
			   ((slug_tracks_2d *) tracks)->trackset_Z_int_meth(),
			   ostreams,
			   pp.query<double>("nebular_den"),
			   pp.query<double>("nebular_temp"),
			   pp.query<double>("nebular_logu"),
			   pp.query<double>("nebular_phi"),
			   pp.query<double>("redshift"),
			   pp.query<int>("nebular_no_metals"));
    }
  } else {
    nebular = nullptr;
  }

  // If using extinction, initialize the extinction curve
  if (pp.query<int>("use_extinct")) {
    if (nebular != nullptr)
      extinct = new slug_extinction(pp, specsyn->lambda(true),
				    nebular->lambda(), rng, ostreams);
    else
      extinct = new slug_extinction(pp, specsyn->lambda(true), rng,
				    ostreams);
  } else {
    extinct = nullptr;
  }

  // Initialize either a galaxy or a single cluster, depending on
  // which type of simulation we're running
  if (sim_type == galaxy_type) {
    galaxy = new slug_galaxy(pp, imf, cmf, clf, sfh, tracks, 
			     specsyn, filters, extinct, nebular, 
			     yields, ostreams);
    cluster = nullptr;
  } else {
    double cluster_mass;
    if (pp.query<int>("random_cluster_mass"))
      cluster_mass = cmf->draw();
    else
      cluster_mass = pp.query<double>("cluster_mass");
    cluster = new slug_cluster(0, cluster_mass, 0.0, imf,
    			       tracks, specsyn, filters,
    			       extinct, nebular, yields,
			       lines, ostreams, clf);
    galaxy = nullptr;
  }

  // Record the output mode and set the checkpoint counter
  out_mode = (outputMode) pp.query<int>("output_mode");
  if (pp.query<int>("checkpoint_interval") == 0)
    checkpoint_ctr = -1; // Indicate no checkpointing
  else
    checkpoint_ctr = pp.get_checkpoint_ctr();  
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_sim::~slug_sim() {

  // Clean up variable parameter storage
  imf_vpdraws.clear();

  // Delete the various objects we created
  if (galaxy != nullptr) delete galaxy;
  if (cluster != nullptr) delete cluster;
  if (specsyn != nullptr) delete specsyn;
  if (sfh != nullptr) delete sfh;
  if (clf != nullptr) delete clf;
  if (cmf != nullptr) delete cmf;
  if (imf != nullptr) delete imf;
  if (tracks != nullptr) delete tracks;
  if (yields != nullptr) delete yields;
  if (rng != nullptr) delete rng;
  if (filters != nullptr) delete filters;
  if (out_time_pdf != nullptr) delete out_time_pdf;
  if (sfr_pdf != nullptr) delete sfr_pdf;
  if (extinct != nullptr) delete extinct;
  if (nebular != nullptr) delete nebular;
  if (lines != nullptr) delete lines;
}


////////////////////////////////////////////////////////////////////////
// Methods to open and close output files
////////////////////////////////////////////////////////////////////////
void slug_sim::open_output(slug_output_files &outfiles, int chknum) {
  if (sim_type == galaxy_type && pp.query<int>("out_integrated"))
    open_integrated_prop(outfiles, chknum);
  if (sim_type == galaxy_type && pp.query<int>("out_integrated_spec")) 
    open_integrated_spec(outfiles, chknum);
  if (sim_type == galaxy_type && pp.query<int>("out_integrated_phot")) 
    open_integrated_phot(outfiles, chknum);
  if (sim_type == galaxy_type && pp.query<int>("out_integrated_yield"))
    open_integrated_yield(outfiles, chknum);
  if (sim_type == galaxy_type && pp.query<int>("out_integrated_sn"))
    open_integrated_sn(outfiles, chknum);
  if (pp.query<int>("out_cluster")) open_cluster_prop(outfiles, chknum);
  if (pp.query<int>("out_cluster_phot")) open_cluster_phot(outfiles, chknum);
  if (pp.query<int>("out_cluster_spec")) open_cluster_spec(outfiles, chknum);
  if (pp.query<int>("out_cluster_yield")) open_cluster_yield(outfiles, chknum);
  if (pp.query<int>("out_cluster_sn")) open_cluster_sn(outfiles, chknum);
  if (pp.query<int>("out_cluster_ew")) open_cluster_ew(outfiles, chknum);
  outfiles.is_open = true;
}

void slug_sim::close_output(slug_output_files &outfiles,
			    int checkpoint_ctr, unsigned int ntrials) {
  
  // If we are closing a checkpoint, edit the number of trials it
  // contains
  if (checkpoint_ctr >= 0) {

    // Figure out which files are open
    vector<std::ofstream *> open_files;
    if (outfiles.int_prop_file.is_open())
      open_files.push_back(&(outfiles.int_prop_file));
    if (outfiles.int_spec_file.is_open())
      open_files.push_back(&(outfiles.int_spec_file));
    if (outfiles.int_phot_file.is_open())
      open_files.push_back(&(outfiles.int_phot_file));
    if (outfiles.int_yield_file.is_open())
      open_files.push_back(&(outfiles.int_yield_file));
    if (outfiles.int_sn_file.is_open())
      open_files.push_back(&(outfiles.int_sn_file));
    if (outfiles.cluster_prop_file.is_open())
      open_files.push_back(&(outfiles.cluster_prop_file));
    if (outfiles.cluster_spec_file.is_open())
      open_files.push_back(&(outfiles.cluster_spec_file));
    if (outfiles.cluster_phot_file.is_open())
      open_files.push_back(&(outfiles.cluster_phot_file));
    if (outfiles.cluster_yield_file.is_open())
      open_files.push_back(&(outfiles.cluster_yield_file));
    if (outfiles.cluster_sn_file.is_open())
      open_files.push_back(&(outfiles.cluster_sn_file));

    // Fix number of trials in file
    for (vector<std::ifstream>::size_type i=0; i<open_files.size(); i++) {
      open_files[i]->seekp(ios_base::beg);
      if (out_mode == ASCII) {
	// ASCII mode; replace the first line
	(*open_files[i]) << "N_Trials = " << setfill('0')
			 << setw(14) << ntrials << endl;
      } else if (out_mode == BINARY) {
	// Binary mode; replace the first uint in file
	open_files[i]->write((char *) &ntrials, sizeof ntrials);
      }
    }
  }

  // Close open files
  if (outfiles.int_prop_file.is_open())
    outfiles.int_prop_file.close();
  if (outfiles.cluster_prop_file.is_open())
    outfiles.cluster_prop_file.close();
  if (outfiles.int_spec_file.is_open())
    outfiles.int_spec_file.close();
  if (outfiles.cluster_spec_file.is_open())
    outfiles.cluster_spec_file.close();
  if (outfiles.int_phot_file.is_open())
    outfiles.int_phot_file.close();
  if (outfiles.cluster_phot_file.is_open())
    outfiles.cluster_phot_file.close();
  if (outfiles.int_yield_file.is_open())
    outfiles.int_yield_file.close();
  if (outfiles.cluster_yield_file.is_open())
    outfiles.cluster_yield_file.close();
  if (outfiles.int_sn_file.is_open())
    outfiles.int_sn_file.close();
  if (outfiles.cluster_sn_file.is_open())
    outfiles.cluster_sn_file.close();

#ifdef ENABLE_FITS
  // Reset number of trials for FITS files
  if (checkpoint_ctr >= 0) {

    // Figure out which files are open
    vector<fitsfile *> open_files;
    if (outfiles.int_prop_fits != nullptr)
      open_files.push_back(outfiles.int_prop_fits);
    if (outfiles.int_spec_fits != nullptr)
      open_files.push_back(outfiles.int_spec_fits);
    if (outfiles.int_phot_fits != nullptr)
      open_files.push_back(outfiles.int_phot_fits);
    if (outfiles.int_yield_fits != nullptr)
      open_files.push_back(outfiles.int_yield_fits);
    if (outfiles.cluster_prop_fits != nullptr)
      open_files.push_back(outfiles.cluster_prop_fits);
    if (outfiles.cluster_spec_fits != nullptr)
      open_files.push_back(outfiles.cluster_spec_fits);
    if (outfiles.cluster_phot_fits != nullptr)
      open_files.push_back(outfiles.cluster_phot_fits);
    if (outfiles.cluster_yield_fits != nullptr)
      open_files.push_back(outfiles.cluster_yield_fits);
    if (outfiles.cluster_ew_fits != nullptr)
      open_files.push_back(outfiles.cluster_ew_fits);

    // Fix number of trials in file
    for (vector<std::ifstream>::size_type i=0; i<open_files.size(); i++) {
      int fits_status = 0;
      fits_movabs_hdu(open_files[i], 2, NULL, &fits_status);
      fits_update_key(open_files[i], TUINT, "N_Trials", 
		      &ntrials, "Number of trials in file",
		      &fits_status);
    }
  }
 
  // Close FITS files
  int fits_status = 0;
  if (out_mode == FITS) {
    if (outfiles.int_prop_fits != nullptr)
      fits_close_file(outfiles.int_prop_fits, &fits_status);
    if (outfiles.cluster_prop_fits != nullptr)
      fits_close_file(outfiles.cluster_prop_fits, &fits_status);
    if (outfiles.int_spec_fits != nullptr)
      fits_close_file(outfiles.int_spec_fits, &fits_status);
    if (outfiles.cluster_spec_fits != nullptr)
      fits_close_file(outfiles.cluster_spec_fits, &fits_status);
    if (outfiles.int_phot_fits != nullptr)
      fits_close_file(outfiles.int_phot_fits, &fits_status);
    if (outfiles.cluster_phot_fits != nullptr)
      fits_close_file(outfiles.cluster_phot_fits, &fits_status);
    if (outfiles.int_yield_fits != nullptr)
      fits_close_file(outfiles.int_yield_fits, &fits_status);
    if (outfiles.cluster_yield_fits != nullptr)
      fits_close_file(outfiles.cluster_yield_fits, &fits_status);
    if (outfiles.int_sn_fits != nullptr)
      fits_close_file(outfiles.int_sn_fits, &fits_status);
    if (outfiles.cluster_sn_fits != nullptr)
      fits_close_file(outfiles.cluster_sn_fits, &fits_status);
    if (outfiles.cluster_ew_fits != nullptr)
      fits_close_file(outfiles.cluster_ew_fits, &fits_status);
    outfiles.int_prop_fits = outfiles.cluster_prop_fits
      = outfiles.int_spec_fits = outfiles.cluster_spec_fits
      = outfiles.int_phot_fits = outfiles.cluster_phot_fits
      = outfiles.int_sn_fits = outfiles.cluster_sn_fits 
      = outfiles.int_yield_fits = outfiles.cluster_yield_fits 
      = outfiles.cluster_ew_fits = nullptr;
  }
#endif

  // Flag that files are closed
  outfiles.is_open = false;
}

////////////////////////////////////////////////////////////////////////
// Method to run a galaxy simulation
////////////////////////////////////////////////////////////////////////
void slug_sim::galaxy_sim() {

  // Prepare to count trials; in serial mode this is trivial, while in
  // MPI mode we need to set up a remote access window so that we can
  // synchronize trial counts across processes
  unsigned long trials_to_do = pp.query<int>("n_trials");
  unsigned long trial_ctr = pp.get_checkpoint_trials();
  unsigned long trial_ctr_loc = 0; // Counts trials on this processor
  unsigned long trial_ctr_last = 1; // Trial counter at last write
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  unsigned long trial_ctr_buf = trial_ctr; // Buffer to hold global trial counter
  MPI_Win win;
  if (comm != MPI_COMM_NULL) {
    if (rank == 0) {
      // Root process holds the accumulated number of trials
      MPI_Win_create(&trial_ctr_buf, sizeof(trial_ctr_buf),
		     sizeof(trial_ctr_buf), MPI_INFO_NULL, comm, &win);
    } else {
      // All other processes don't need a remote access window
      MPI_Win_create(NULL, 0, sizeof(int), MPI_INFO_NULL, comm, &win);
    }
  }
#endif

  // Initialize output file information
  slug_output_files outfiles;
  vector<string>::size_type nphot =
    pp.query<vector<string> >("phot_bands").size();
  
  // Main loop
  while (true) {

    // Increment the trial counters; under MPI, the global counter
    // update is done as a Fetch_and_op operation. Note that
    // Fetch_and_op returns the value of the buffer before the
    // addition is done, so we still need to increment the trial
    // counter
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
    if (comm != MPI_COMM_NULL) {
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
      unsigned long i=1;
      MPI_Fetch_and_op(&i, &trial_ctr, MPI_UNSIGNED_LONG, 0, 0,
		       MPI_SUM, win);
      MPI_Win_unlock(0, win);
    }
#endif
    trial_ctr++;
    trial_ctr_loc++;

    // If global trial counter exceeds number of trials we're supposed
    // to do, exit the loop
    if (trial_ctr > trials_to_do) break;

    // Figure out if we need to open a new output file
    bool open_new_output = false;
    if (trial_ctr_loc == 1) {
      // Always open new file on first trial
      open_new_output = true;
    } else if (pp.query<int>("checkpoint_interval") != 0) {
      if ((trial_ctr_loc-1) % pp.query<int>("checkpoint_interval") == 0) {
	// Create a new checkpoint at periodic intervals if we are
	// checkpointing
	open_new_output = true;
      }
    }

    // If we are opening new outputs, do so now
    if (open_new_output) {
      // Write status message
      if (pp.query<int>("verbosity") > 0 &&
	  checkpoint_ctr > 0 &&
	  trial_ctr_loc != 1)
	ostreams.slug_out << "finalizing checkpoint "
			  << checkpoint_ctr << std::endl;
      // Close old output if it is open
      if (outfiles.is_open)
	close_output(outfiles, checkpoint_ctr-1,
		     trial_ctr_loc - trial_ctr_last);
      // Open new checkpoint files
      open_output(outfiles, checkpoint_ctr);
      // Increment checkpoint counter, and update number of trials
      // written counter, if we are checkpointing
      if (checkpoint_ctr >= 0) {
	checkpoint_ctr++;
	trial_ctr_last = trial_ctr_loc;
      }
    }

    // If sufficiently verbose, print status
    if (pp.query<int>("verbosity") > 0)
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
      ostreams.slug_out << " starting trial " << trial_ctr << " of "
			<< trials_to_do << endl;
#else
      ostreams.slug_out << "starting trial " << trial_ctr << " of "
			<< trials_to_do << endl;
#endif

    // Check for variable segments
    if (is_imf_var == true) {
      // Draw new values for variable parameters
      // Update IMF segments and recompute weights 
      imf_vpdraws = imf->vseg_draw();
      // Reset range restrictions
      imf->set_stoch_lim(pp.query<double>("get_min_stoch_mass"));
    }

    // Reset the galaxy
    galaxy->reset();

    // Write trial separator to ASCII files if operating in ASCII
    // mode
    if ((out_mode == ASCII) && !open_new_output) {
	
      // Skip if we just opened new files
      if (pp.query<int>("checkpoint_interval") > 0) {
	if (((trial_ctr_loc-1) % pp.query<int>("checkpoint_interval") == 0) &&
	    (trial_ctr_loc != 1))
	  break;
      }
      
      if (pp.query<int>("out_integrated")) { 
        // Write separators
        int ncol = 9*14-3;      
	if (is_imf_var==true) ncol += (imf_vpdraws.size())*14;
	write_separator(outfiles.int_prop_file, ncol);
      }

      unsigned int extrafields = 0;       //Store IMFs in prop file
      unsigned int nfield = 1;
      if (nebular != nullptr) nfield++;
      if (extinct != nullptr) {
	nfield++;
	extrafields++;
	if (nebular != nullptr) nfield++;
      }
      // Add on IMF fields
      if (is_imf_var==true) extrafields += (imf_vpdraws.size());
      if (pp.query<int>("out_integrated_spec")) 
	write_separator(outfiles.int_spec_file, (2+nfield)*14-3);
      if (pp.query<int>("out_integrated_phot"))
	write_separator(outfiles.int_phot_file,
			(1+nfield*nphot)*21-3);
      if (pp.query<int>("out_integrated_sn"))
	write_separator(outfiles.int_sn_file, 14*3-3);
      if (pp.query<int>("out_integrated_yield"))
	write_separator(outfiles.int_yield_file, 14*5-3);
      if (pp.query<int>("out_cluster"))
	write_separator(outfiles.cluster_prop_file, (10+extrafields)*14-3);
      if (pp.query<int>("out_cluster_spec"))
	write_separator(outfiles.cluster_spec_file, (4+nfield)*14-3);
      if (pp.query<int>("out_cluster_phot"))
	write_separator(outfiles.cluster_phot_file,
			(2+nfield*nphot)*21-3);
      if (pp.query<int>("out_cluster_sn"))
	write_separator(outfiles.cluster_sn_file, 14*4-3);
      if (pp.query<int>("out_cluster_yield"))
	write_separator(outfiles.cluster_yield_file, 14*6-3);
    }

    // If the output time is randomly changing, draw a new output time
    // for this trial
    if (pp.query<int>("random_output_times")) {
      outTimes.resize(0);
      outTimes.push_back(out_time_pdf->draw());
    }

    // If the SFR is randomly changing, draw a new SFR for this trial
    if (pp.query<int>("random_sfr")) {
      double sfr = sfr_pdf->draw();
      sfh->setNorm(outTimes.back()*sfr);
    }

    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // Flag if we should delete clusters on this pass
      bool del_cluster = (j==outTimes.size()-1) &&
	(!pp.query<int>("out_cluster_spec")) && (!pp.query<int>("out_cluster_phot")) &&
	(!pp.query<int>("out_cluster_yield"));

      // If sufficiently verbose, print status
      if (pp.query<int>("verbosity") > 1)
	ostreams.slug_out << "  trial " << trial_ctr << ", advance to time " 
			  << outTimes[j] << std::endl;
      
      // Advance to next time
      galaxy->advance(outTimes[j]);

      // Write physical properties if requested
      if (pp.query<int>("out_integrated")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_prop(outfiles.int_prop_file, out_mode,
					trial_ctr_loc, imf_vpdraws);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_prop(outfiles.int_prop_fits, trial_ctr_loc,
					imf_vpdraws);
	}
#endif
      }
      if (pp.query<int>("out_cluster")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_prop(outfiles.cluster_prop_file, out_mode,
				     trial_ctr_loc, imf_vpdraws);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_prop(outfiles.cluster_prop_fits,
				     trial_ctr_loc, imf_vpdraws);
	}
#endif
      }

      // Write yield if requested
      if (pp.query<int>("out_integrated_yield")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_yield(outfiles.int_yield_file, out_mode,
					 trial_ctr_loc, del_cluster);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_yield(outfiles.int_yield_fits,
					 trial_ctr_loc, del_cluster);
	}
#endif
      }
      if (pp.query<int>("out_cluster_yield")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_yield(outfiles.cluster_yield_file, out_mode,
				      trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_yield(outfiles.cluster_yield_fits,
				      trial_ctr_loc);
	}
#endif
      }

      // Write supernova count if requested
      if (pp.query<int>("out_integrated_sn")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_sn(outfiles.int_sn_file, out_mode,
					 trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_sn(outfiles.int_sn_fits,
					 trial_ctr_loc);
	}
#endif
      }
      if (pp.query<int>("out_cluster_sn")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_sn(outfiles.cluster_sn_file, out_mode,
				   trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_sn(outfiles.cluster_sn_fits,
				   trial_ctr_loc);
	}
#endif
      }

      // Write spectra if requested
      if (pp.query<int>("out_integrated_spec")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_spec(outfiles.int_spec_file, out_mode,
					trial_ctr_loc, del_cluster);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_spec(outfiles.int_spec_fits, trial_ctr_loc,
					del_cluster);
	}
#endif
      }
      if (pp.query<int>("out_cluster_spec")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_spec(outfiles.cluster_spec_file, out_mode,
				     trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_spec(outfiles.cluster_spec_fits,
				     trial_ctr_loc);
	}
#endif
      }
      
      // Write photometry if requested
      if (pp.query<int>("out_integrated_phot")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_integrated_phot(outfiles.int_phot_file, out_mode,
					trial_ctr_loc, del_cluster);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_integrated_phot(outfiles.int_phot_fits,
					trial_ctr_loc, del_cluster);
	}
#endif
      }
      if (pp.query<int>("out_cluster_phot")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  galaxy->write_cluster_phot(outfiles.cluster_phot_file, out_mode,
				     trial_ctr_loc);
#ifdef ENABLE_FITS
	} else {
	  galaxy->write_cluster_phot(outfiles.cluster_phot_fits,
				     trial_ctr_loc);
	}
#endif
      }
    }
  }
  
  // Close last output file
  if ((pp.query<int>("verbosity") > 0) && (trial_ctr_loc - trial_ctr_last > 0) &&
      (checkpoint_ctr > 0))
    ostreams.slug_out << "finalizing checkpoint "
		      << checkpoint_ctr << std::endl;
  close_output(outfiles, checkpoint_ctr, trial_ctr_loc - trial_ctr_last);
  
  // Free MPI window
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) MPI_Win_free(&win);
#endif
  
  // Clean up vector of variable pdfs
  if (is_imf_var == true) {
    imf->cleanup();
  }
}


////////////////////////////////////////////////////////////////////////
// Method to run a cluster simulation
////////////////////////////////////////////////////////////////////////
void slug_sim::cluster_sim() {

  // Prepare to count trials; in serial mode this is trivial, while in
  // MPI mode we need to set up a remote access window so that we can
  // synchronize trial counts across processes
  unsigned long trials_to_do = pp.query<int>("n_trials");
  unsigned long trial_ctr = pp.get_checkpoint_trials();
  unsigned long trial_ctr_loc = 0; // Counts trials on this processor
  unsigned long trial_ctr_last = 1; // Trial counter at last write
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  unsigned long trial_ctr_buf = trial_ctr; // Buffer to hold global trial counter
  MPI_Win win;
  if (comm != MPI_COMM_NULL) {
    if (rank == 0) {
      // Root process holds the accumulated number of trials
      MPI_Win_create(&trial_ctr_buf, sizeof(trial_ctr_buf),
		     sizeof(trial_ctr_buf), MPI_INFO_NULL, comm, &win);
    } else {
      // All other processes don't need a remote access window
      MPI_Win_create(NULL, 0, sizeof(int), MPI_INFO_NULL, comm, &win);
    }
  }
#endif

  // Initialize output file information
  slug_output_files outfiles;
  outfiles.is_open = false;
#ifdef ENABLE_FITS
  outfiles.int_prop_fits = outfiles.int_spec_fits = outfiles.int_phot_fits
    = outfiles.int_yield_fits = outfiles.cluster_prop_fits
    = outfiles.cluster_spec_fits = outfiles.cluster_phot_fits
    = outfiles.cluster_yield_fits = outfiles.cluster_ew_fits= nullptr;
#endif
  vector<string>::size_type nphot =
    pp.query<vector<string> >("phot_bands").size();
  
  // Loop over trials
  unsigned long id = 0;
  while (true) {

    // Increment the global trial counter
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
    // Method for MPI: add 1 to the global RMA variable trial_ctr
    if (comm != MPI_COMM_NULL) {
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
      unsigned long i=1;
      MPI_Fetch_and_op(&i, &trial_ctr, MPI_UNSIGNED_LONG, 0, 0,
		       MPI_SUM, win);
      MPI_Win_unlock(0, win);
    }
#endif
    trial_ctr++;

    // Update local trial counter
    trial_ctr_loc++;

    // If global trial counter exceeds number of trials we're supposed
    // to do, exit the loop
    if (trial_ctr > trials_to_do) break;

    // Figure out if we need to open a new output file
    bool open_new_output = false;
    if (trial_ctr_loc == 1) {
      // Always open new file on first trial
      open_new_output = true;
    } else if (pp.query<int>("checkpoint_interval") != 0) {
      if ((trial_ctr_loc-1) % pp.query<int>("checkpoint_interval") == 0) {
	// Create a new checkpoint at periodic intervals if we are
	// checkpointing
	open_new_output = true;
      }
    }

    // Open a new output file if we need to
    if (open_new_output) {
      // Write status message
      if (pp.query<int>("verbosity") > 0 && checkpoint_ctr > 0 &&
	  trial_ctr_loc != 1)
	ostreams.slug_out << "finalizing checkpoint "
			  << checkpoint_ctr << std::endl;
      // Close old output if it is open
      if (outfiles.is_open) {
	close_output(outfiles, checkpoint_ctr-1,
		     trial_ctr_loc - trial_ctr_last);
      }
      // Open new output files
      open_output(outfiles, checkpoint_ctr);
      // Increment checkpoint counter, and update number of trials
      // written counter, if we are checkpointing
      if (checkpoint_ctr >= 0) {
	checkpoint_ctr++;
	trial_ctr_last = trial_ctr_loc;
      }
    }
	
    // If sufficiently verbose, print status
    if (pp.query<int>("verbosity") > 0)
      ostreams.slug_out << "starting trial " << trial_ctr << " of "
			<< trials_to_do << std::endl;
    
    // If the output time is randomly changing, draw a new output time
    // for this trial
    if (pp.query<int>("random_output_times")) {
      outTimes.resize(0);
      outTimes.push_back(out_time_pdf->draw());
    }

    // Check for variable segments
    if (is_imf_var == true) {
      // Draw new values for variable parameters
      // Update IMF segments and recompute weights 
      imf_vpdraws = imf->vseg_draw();
      // Reset range restrictions
      imf->set_stoch_lim(pp.query<double>("min_stoch_mass"));
    }
    
    // Reset the cluster if the mass is constant, destroy it and build
    // a new one if not
    if (pp.query<int>("random_cluster_mass")) {
      if (cluster) {
	id = cluster->get_id();
	delete cluster;
	cluster = nullptr;
      }
      double m_cl = cmf->draw();
      if (pp.query<int>("lamers_loss")) {
	double lgamma = pp.query<double>("lamers_gamma");
	double t4 = pp.query<double>("lamers_t4");
	double fac = lgamma*pow(1.0e4/m_cl, lgamma)*outTimes.back()/t4;
	if (fac > 1.0) {
	  // Cluster gone completely, end trial
	  if (pp.query<int>("verbosity") > 1) {
	    ostreams.slug_out << "  Lamers evaporation removed cluster"
			      << ", ending trial" << endl;
	  }
	  continue;
	}  
	m_cl *= pow(1.0 - fac, 1.0/lgamma);
      }
      cluster = new slug_cluster(id+1, m_cl, 0.0, imf,
				 tracks, specsyn, filters,
				 extinct, nebular, yields, lines, ostreams, clf);
    } else {
      cluster->reset();
    }

    // Write trial separator to ASCII files if operating in ASCII
    // mode
    if ((out_mode == ASCII) && !open_new_output) {
      if (pp.query<int>("out_cluster")) {
	int ncol = 10*14-3;
	if (pp.query<int>("use_extinct")) ncol += 14;
	if (pp.query<int>("use_extinct") &&
	    pp.query<int>("use_nebular_extinction")) ncol += 14;
	if (is_imf_var==true) ncol += (imf_vpdraws.size())*14;
	write_separator(outfiles.cluster_prop_file, ncol);
      }
      if (pp.query<int>("out_cluster_spec"))
	write_separator(outfiles.cluster_spec_file, 4*14-3);
      if (pp.query<int>("out_cluster_phot"))
	write_separator(outfiles.cluster_phot_file, (2+nphot)*21-3);
      if (pp.query<int>("out_cluster_yield"))
	write_separator(outfiles.cluster_yield_file, 5*14-3);
      if (pp.query<int>("out_cluster_sn"))
	write_separator(outfiles.cluster_sn_file, 4*14-3);
    }
    
    // Loop over time steps
    for (unsigned int j=0; j<outTimes.size(); j++) {

      // If sufficiently verbose, print status
      if (pp.query<int>("verbosity") > 1)
	ostreams.slug_out << "  trial " << trial_ctr
			  << ", advance to time " 
			  << outTimes[j] << std::endl;
      
      // Advance to next time
      cluster->advance(outTimes[j]);

      // See if cluster has disrupted; if so, terminate this iteration
      if (cluster->disrupted()) {
	if (pp.query<int>("verbosity") > 1)
	  ostreams.slug_out << "  cluster disrupted, terminating trial"
			    << std::endl;
	break;
      }

      // Write physical properties if requested
      if (pp.query<int>("out_cluster")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_prop(outfiles.cluster_prop_file, out_mode,
			      trial_ctr_loc, true,
			      imf_vpdraws);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_prop(outfiles.cluster_prop_fits, trial_ctr_loc,
			      imf_vpdraws);
	}
#endif
      }

      // Write spectrum if requested
      if (pp.query<int>("out_cluster_spec")) { 
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_spectrum(outfiles.cluster_spec_file, out_mode,
				  trial_ctr_loc, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_spectrum(outfiles.cluster_spec_fits, trial_ctr_loc);
	}
#endif
      }
#ifdef ENABLE_FITS    
    // Write equivalent width if requested
    if (pp.query<int>("out_cluster_ew")) {
      cluster->write_ew(outfiles.cluster_ew_fits, trial_ctr_loc);
    }
#endif

      // Write photometry if requested
      if (pp.query<int>("out_cluster_phot")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_photometry(outfiles.cluster_phot_file, out_mode,
				    trial_ctr_loc, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_photometry(outfiles.cluster_phot_fits, trial_ctr_loc);
	}
#endif
      }

      // Write yields if requested
      if (pp.query<int>("out_cluster_yield")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_yield(outfiles.cluster_yield_file, out_mode,
			       trial_ctr_loc, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_yield(outfiles.cluster_yield_fits, trial_ctr_loc);
	}
#endif
      }

      // Write supernova counts if requested
      if (pp.query<int>("out_cluster_sn")) {
#ifdef ENABLE_FITS
	if (out_mode != FITS) {
#endif
	  cluster->write_sn(outfiles.cluster_sn_file, out_mode,
			    trial_ctr_loc, true);
#ifdef ENABLE_FITS
	} else {
	  cluster->write_sn(outfiles.cluster_sn_fits, trial_ctr_loc);
	}
#endif
      }
    }
  }

  // Close last output file
  if ((pp.query<int>("verbosity") > 0) &&
      (trial_ctr_loc - trial_ctr_last > 0) &&
      (checkpoint_ctr > 0))
    ostreams.slug_out << "finalizing checkpoint "
		      << checkpoint_ctr << std::endl;
  close_output(outfiles, checkpoint_ctr, trial_ctr_loc - trial_ctr_last);

  // Free MPI window
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) MPI_Win_free(&win);
#endif

  // Clean up vector of variable pdfs
  if (is_imf_var == true) {
    imf->cleanup();
  }
}


////////////////////////////////////////////////////////////////////////
// Open integrated properties file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_prop(slug_output_files &outfiles,
				    int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_integrated_prop";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.int_prop_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    outfiles.int_prop_file.open(full_path.c_str(), ios::out | ios::binary);
  } 
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.int_prop_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err
	<< "unable to open integrated properties file "
	<< full_path.string()
	<< "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.int_prop_file.is_open()) {
      ostreams.slug_err << "unable to open intergrated properties file " 
			<< full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif 

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.int_prop_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    outfiles.int_prop_file << setw(14) << left << "Time"
		  << setw(14) << left << "TargetMass"
		  << setw(14) << left << "ActualMass"
		  << setw(14) << left << "LiveMass"
		  << setw(14) << left << "StellarMass"
		  << setw(14) << left << "ClusterMass"
		  << setw(14) << left << "NumClusters"
		  << setw(14) << left << "NumDisClust"
		  << setw(14) << left << "NumFldStar";

    //If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var) {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++) {
        outfiles.int_prop_file << setw(14) << left 
		      << "VP"+ std::to_string(static_cast<long long>(p));
      }
    }
    outfiles.int_prop_file << endl;

    // Write unit line
    outfiles.int_prop_file << setw(14) << left << "(yr)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << "(Msun)"
		  << setw(14) << left << ""
		  << setw(14) << left << ""
		  << setw(14) << left << "";
    if (is_imf_var) {
      //Loop over the variable parameters
      for (vector<double>::size_type p=0; p<imf_vpdraws.size(); p++) {
        outfiles.int_prop_file << setw(14) << left << "";
      }
    }      
    outfiles.int_prop_file << endl;
		  		  		  
    // Line of underscores
    outfiles.int_prop_file << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------";
    if (is_imf_var) {
      for (vector<double>::size_type p=0; p<imf_vpdraws.size(); p++) {
        outfiles.int_prop_file << setw(14) << left << "-----------";
      }
    }      
    outfiles.int_prop_file << endl;
  
  } else if (out_mode == BINARY) {

    // State number of trials expected to be in this file
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.int_prop_file.write((char *) &ntrial, sizeof ntrial);
    }

    // File starts with an integer for the number of variable parameters
    if (is_imf_var == true) {
        int nvps = imf_vpdraws.size();
        outfiles.int_prop_file.write((char *) &nvps, sizeof nvps);
    } else {    
      int nvps = 0;
      outfiles.int_prop_file.write((char *) &nvps, sizeof nvps);    
    }      
  }
  
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    int ncol = 10;    
    vector<string> ttype_str = 
      { "Trial", "Time", "TargetMass", "ActualMass", "LiveMass", 
	"StellarMass", "ClusterMass", "NumClusters", "NumDisClust", 
	"NumFldStar" };
    vector<string> tform_str = 
      { "1K", "1D", "1D", "1D", "1D", "1D", "1D", "1K", "1K", "1K" };
    vector<string> tunit_str = 
      { "Msun", "Msun", "Msun", "Msun", "Msun", "Msun", "", "", "", "" };
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var) {
      // Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size(); p++) {
        ttype_str.push_back("VP"+std::to_string(static_cast<long long>(p)));
        tform_str.push_back("1D");
        tunit_str.push_back("");
        ncol++; 
      }
    }

    // Make things we can pass to c
    /*
    char *ttype[ncol], *tform[ncol], *tunit[ncol];
    for (int i=0; i<ncol; i++)  {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    */
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }

    // Create FITS table
    int fits_status = 0;
    fits_create_tbl(outfiles.int_prop_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, nullptr, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.int_prop_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster properties file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_prop(slug_output_files &outfiles,
				 int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_cluster_prop";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.cluster_prop_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    outfiles.cluster_prop_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.cluster_prop_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err
	<< "unable to open cluster properties file "
	<< full_path.string()
	<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.cluster_prop_file.is_open()) {
      ostreams.slug_err << "unable to open cluster properties file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.cluster_prop_file << "N_Trials = " << setfill('0')
			<< setw(14) << 0 << setfill(' ') << endl;
    outfiles.cluster_prop_file << setw(14) << left << "UniqueID"
		      << setw(14) << left << "Time"
		      << setw(14) << left << "FormTime"
		      << setw(14) << left << "Lifetime"
		      << setw(14) << left << "TargetMass"
		      << setw(14) << left << "BirthMass"
		      << setw(14) << left << "LiveMass"
		      << setw(14) << left << "StellarMass"
		      << setw(14) << left << "NumStar"
		      << setw(14) << left << "MaxStarMass";
    if (extinct != nullptr) {
      outfiles.cluster_prop_file << setw(14) << left << "A_V";
      if (extinct->excess_neb_extinct()) 
	outfiles.cluster_prop_file << setw(14) << left << "A_Vneb";
    }
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
#ifndef __INTEL_COMPILER
        outfiles.cluster_prop_file << setw(14) << left << "VP"+ std::to_string(p);
#else
	// This is needed to fix a deficiency in the intel implementation of
	// the C++11 STL, which fails to define std::to_string(int)
	outfiles.cluster_prop_file << setw(14) << left << "VP" +
	  std::to_string(static_cast<long long>(p));
#endif
      }
    }
            
    outfiles.cluster_prop_file << endl;
    outfiles.cluster_prop_file << setw(14) << left << ""
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << "(Msun)"
		      << setw(14) << left << ""
		      << setw(14) << left << "(Msun)";
		      
    if (extinct != nullptr) {
      outfiles.cluster_prop_file << setw(14) << left << "(mag)";
      if (extinct->excess_neb_extinct()) 
	outfiles.cluster_prop_file << setw(14) << left << "(mag)";
    }
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
        outfiles.cluster_prop_file << setw(14) << left << "";
      }
    
    }      
      
    outfiles.cluster_prop_file << endl;
    outfiles.cluster_prop_file << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------";
    if (extinct != nullptr) {
      outfiles.cluster_prop_file << setw(14) << left << "-----------";
      if (extinct->excess_neb_extinct()) 
	outfiles.cluster_prop_file << setw(14) << left << "-----------";
    }
    
    // If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
        outfiles.cluster_prop_file << setw(14) << left << "-----------";
      }
    }
      
    outfiles.cluster_prop_file << endl;
  } else if (out_mode == BINARY) {

    // State number of trials in this file
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.cluster_prop_file.write((char *) &ntrial, sizeof ntrial);
    }

    // File starts with a bit indicating whether we're using extinction
    bool use_extinct = (extinct != nullptr);
    outfiles.cluster_prop_file.write((char *) &use_extinct, sizeof use_extinct);
    bool use_neb_excess = false;
    if (use_extinct) {
      if (extinct->excess_neb_extinct()) use_neb_excess = true;
    }
    outfiles.cluster_prop_file.write((char *) &use_neb_excess,
			    sizeof use_neb_excess);

    // File then contains an integer number of variable parameters
    if (is_imf_var == true)
    {
      int nvps = imf_vpdraws.size();
      outfiles.cluster_prop_file.write((char *) &nvps, sizeof nvps);
    }
    else
    {    
      int nvps = 0;
      outfiles.cluster_prop_file.write((char *) &nvps, sizeof nvps);    
    }
    
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    int ncol = 11;
    vector<string> ttype_str = 
      { "Trial", "UniqueID", "Time", "FormTime", "Lifetime",
	"TargetMass", "BirthMass", "LiveMass", "StellarMass",
	"NumStar", "MaxStarMass" };
    vector<string> tform_str = 
      { "1K", "1K", "1D", "1D", "1D", "1D", "1D", "1D", "1D", 
	"1K", "1D" };
    vector<string> tunit_str = 
      { "", "", "yr", "yr", "yr", "Msun", "Msun", "Msun", "Msun", 
	"", "Msun" };
    if (extinct != nullptr) {
      ttype_str.push_back("A_V");
      tform_str.push_back("1D");
      tunit_str.push_back("mag");
      ncol++;
      if (extinct->excess_neb_extinct()) {
	ttype_str.push_back("A_Vneb");
	tform_str.push_back("1D");
	tunit_str.push_back("mag");
	ncol++;
      }
    }
    
    //If we have a variable IMF, write headers for the variable parameters
    if (is_imf_var == true)
    {
      //Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imf_vpdraws.size();p++)
      {
#ifndef __INTEL_COMPILER
        ttype_str.push_back("VP"+std::to_string(p));
#else
	// This is needed to fix a deficiency in the intel implementation of
	// the C++11 STL, which fails to define std::to_string(int)
        ttype_str.push_back("VP"+std::to_string(static_cast<long long>(p)));
#endif
        tform_str.push_back("1D");
        tunit_str.push_back("");
        ncol++; 
      }    
    } 
    
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }

    // Create the table
    int fits_status = 0;
    fits_create_tbl(outfiles.cluster_prop_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, nullptr, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.cluster_prop_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated spectra file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_spec(slug_output_files &outfiles,
				    int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_integrated_spec";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.int_spec_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      outfiles.int_spec_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.int_spec_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err << "unable to open integrated spectrum file "
			<< full_path.string()
			<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.int_spec_file.is_open()) {
      ostreams.slug_err << "unable to open intergrated spectrum file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.int_spec_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    outfiles.int_spec_file << setw(14) << left << "Time"
		  << setw(14) << left << "Wavelength"
		  << setw(14) << left << "L_lambda";
    if (nebular != nullptr)
      outfiles.int_spec_file << setw(14) << left << "L_l_neb";
    if (extinct != nullptr) {
      outfiles.int_spec_file << setw(14) << left << "L_lambda_ex";
      if (nebular != nullptr)
	outfiles.int_spec_file << setw(14) << left << "L_l_neb_ex";
    }
    outfiles.int_spec_file << endl;
    outfiles.int_spec_file << setw(14) << left << "(yr)"
		  << setw(14) << left << "(Angstrom)"
		  << setw(14) << left << "(erg/s/A)";
    if (nebular != nullptr)
      outfiles.int_spec_file << setw(14) << left << "(erg/s/A)";
    if (extinct != nullptr) {
      outfiles.int_spec_file << setw(14) << left << "(erg/s/A)";
      if (nebular != nullptr)
	outfiles.int_spec_file << setw(14) << left << "(erg/s/A)";
    }
    outfiles.int_spec_file << endl;
    outfiles.int_spec_file << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------"
		  << setw(14) << left << "-----------";
    if (nebular != nullptr)
      outfiles.int_spec_file << setw(14) << left << "-----------";
    if (extinct != nullptr) {
      outfiles.int_spec_file << setw(14) << left << "-----------";
      if (nebular != nullptr)
	outfiles.int_spec_file << setw(14) << left << "-----------";
    }
    outfiles.int_spec_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.int_spec_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with two bits indicating whether we're using
    // nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    outfiles.int_spec_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    outfiles.int_spec_file.write((char *) &use_extinct, sizeof use_extinct);
    // Write list of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    outfiles.int_spec_file.write((char *) &nl, sizeof nl);
    outfiles.int_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
    // Write list of nebular wavelengths if using nebular emission
    if (use_nebular) {
      vector<double> lambda_neb = nebular->lambda();
      vector<double>::size_type nl_neb = lambda_neb.size();
      outfiles.int_spec_file.write((char *) &nl_neb, sizeof nl_neb);
      outfiles.int_spec_file.write((char *) &(lambda_neb[0]),
			  nl_neb*sizeof(double));
    }
    // Write list of extincted wavelengths if using extinction
    if (use_extinct) {
      vector<double> lambda_ext = extinct->lambda();
      vector<double>::size_type nl_ext = lambda_ext.size();
      outfiles.int_spec_file.write((char *) &nl_ext, sizeof nl_ext);
      outfiles.int_spec_file.write((char *) &(lambda_ext[0]),
			  nl_ext*sizeof(double));
    }
    // Write list of nebular and extincted wavelengths if using both
    // nebular emission and extinction
    if (use_nebular && use_extinct) {
      vector<double> lambda_neb_ext = extinct->lambda_neb();
      vector<double>::size_type nl_neb_ext = lambda_neb_ext.size();
      outfiles.int_spec_file.write((char *) &nl_neb_ext, sizeof nl_neb_ext);
      outfiles.int_spec_file.write((char *) &(lambda_neb_ext[0]),
			  nl_neb_ext*sizeof(double));
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {

    // In FITS mode, write the wavelength information in the first HDU
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    vector<double> lambda_neb, lambda_ext, lambda_neb_ext;
    vector<double>::size_type nl_ext = 0, nl_neb = 0, nl_neb_ext = 0;
    vector<string> ttype_str = { "Wavelength" };
    vector<string> tform_str;
    tform_str.push_back(lexical_cast<string>(nl) + "D");
    vector<string> tunit_str = { "Angstrom" };
    int ncol=1;
    if (nebular != nullptr) {
      lambda_neb = nebular->lambda();
      nl_neb = lambda_neb.size();
      ttype_str.push_back("Wavelength_neb");
      tform_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr) {
      lambda_ext = extinct->lambda();
      nl_ext = lambda_ext.size();
      ttype_str.push_back("Wavelength_ex");
      tform_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr && nebular != nullptr) {
      lambda_neb_ext = extinct->lambda_neb();
      nl_neb_ext = lambda_neb_ext.size();
      ttype_str.push_back("Wavelength_neb_ex");
      tform_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char wl_name[] = "Wavelength";

    // Create the table
    fits_create_tbl(outfiles.int_spec_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, wl_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.int_spec_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
    
    // Write wavelength data to table
    fits_write_col(outfiles.int_spec_fits, TDOUBLE, 1, 1, 1, nl,
		   lambda.data(), &fits_status);
    int col=2;
    if (nebular != nullptr) {
      fits_write_col(outfiles.int_spec_fits, TDOUBLE, col, 1, 1, nl_neb,
		     lambda_neb.data(), &fits_status);
      col++;
    }
    if (extinct != nullptr) {
      fits_write_col(outfiles.int_spec_fits, TDOUBLE, col, 1, 1, nl_ext,
		     lambda_ext.data(), &fits_status);
      col++;
    }
    if (nebular != nullptr && extinct != nullptr) {
      fits_write_col(outfiles.int_spec_fits, TDOUBLE, col, 1, 1, nl_neb_ext,
		     lambda_neb_ext.data(), &fits_status);
      col++;
    }

    // Create a new table to hold the computed spectra.
    char spec_name[] = "Spectra";
    vector<string> ttype2_str = { "Trial", "Time", "L_lambda" };
    vector<string> tform2_str = { "1K", "1D", "" };
    tform2_str[2] = lexical_cast<string>(nl) + "D";
    vector<string> tunit2_str = { "", "yr", "erg/s/A" };
    ncol = 3;
    if (nebular != nullptr) {
      ttype2_str.push_back("L_lambda_neb");
      tform2_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    if (extinct != nullptr) {
      ttype2_str.push_back("L_lambda_ex");
      tform2_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    if (nebular != nullptr && extinct != nullptr) {
      ttype2_str.push_back("L_lambda_neb_ex");
      tform2_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    fits_create_tbl(outfiles.int_spec_fits, BINARY_TBL, 0, ncol,
		    ttype2, tform2, tunit2, spec_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster spectra file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_spec(slug_output_files &outfiles,
				 int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_cluster_spec";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.cluster_spec_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      outfiles.cluster_spec_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.cluster_spec_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err << "unable to open cluster spectrum file "
			<< full_path.string()
			<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.cluster_spec_file.is_open()) {
      ostreams.slug_err << "unable to open cluster spectrum file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.cluster_spec_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    outfiles.cluster_spec_file << setw(14) << left << "UniqueID"
		      << setw(14) << left << "Time"
		      << setw(14) << left << "Wavelength"
		      << setw(14) << left << "L_lambda";
    if (nebular != nullptr)
      outfiles.cluster_spec_file << setw(14) << left << "L_l_neb";
    if (extinct != nullptr) {
      outfiles.cluster_spec_file << setw(14) << left << "L_lambda_ex";
      if (nebular != nullptr)
	outfiles.cluster_spec_file << setw(14) << left << "L_l_neb_ex";
    }
    outfiles.cluster_spec_file << endl;
    outfiles.cluster_spec_file << setw(14) << left << ""
		      << setw(14) << left << "(yr)"
		      << setw(14) << left << "(Angstrom)"
		      << setw(14) << left << "(erg/s/A)";
    if (nebular != nullptr)
      outfiles.cluster_spec_file << setw(14) << left << "(erg/s/A)";
    if (extinct != nullptr) {
      outfiles.cluster_spec_file << setw(14) << left << "(erg/s/A)";
      if (nebular != nullptr)
	outfiles.cluster_spec_file << setw(14) << left << "(erg/s/A)";
    }
    outfiles.cluster_spec_file << endl;
    outfiles.cluster_spec_file << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------"
		      << setw(14) << left << "-----------";
    if (nebular != nullptr)
      outfiles.cluster_spec_file << setw(14) << left << "-----------";
    if (extinct != nullptr) {
      outfiles.cluster_spec_file << setw(14) << left << "-----------";
      if (nebular != nullptr)
	outfiles.cluster_spec_file << setw(14) << left << "-----------";
    }
    outfiles.cluster_spec_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.cluster_spec_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with two bits indicating whether we're using
    // nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    outfiles.cluster_spec_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    outfiles.cluster_spec_file.write((char *) &use_extinct, sizeof use_extinct);
    // List of wavelengths
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    outfiles.cluster_spec_file.write((char *) &nl, sizeof nl);
    outfiles.cluster_spec_file.write((char *) &(lambda[0]), nl*sizeof(double));
    // List of nebular wavelengths
    if (use_nebular) {
      vector<double> lambda_neb = nebular->lambda();
      vector<double>::size_type nl_neb = lambda_neb.size();
      outfiles.cluster_spec_file.write((char *) &nl_neb, sizeof nl_neb);
      outfiles.cluster_spec_file.write((char *) &(lambda_neb[0]),
			      nl_neb*sizeof(double));
    }
    // List of extincted wavelengths
    if (use_extinct) {
      vector<double> lambda_ext = extinct->lambda();
      vector<double>::size_type nl_ext = lambda_ext.size();
      outfiles.cluster_spec_file.write((char *) &nl_ext, sizeof nl_ext);
      outfiles.cluster_spec_file.write((char *) &(lambda_ext[0]),
			      nl_ext*sizeof(double));
    }
    // List of nebular extincted wavelengths
    if (use_extinct && use_nebular) {
      vector<double> lambda_neb_ext = extinct->lambda_neb();
      vector<double>::size_type nl_neb_ext = lambda_neb_ext.size();
      outfiles.cluster_spec_file.write((char *) &nl_neb_ext, sizeof nl_neb_ext);
      outfiles.cluster_spec_file.write((char *) &(lambda_neb_ext[0]),
			      nl_neb_ext*sizeof(double));
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {

    // In FITS mode, write the wavelength information in the first HDU
    vector<double> lambda = specsyn->lambda();
    vector<double>::size_type nl = lambda.size();
    vector<double> lambda_neb, lambda_ext, lambda_neb_ext;
    vector<double>::size_type nl_neb = 0, nl_ext = 0, nl_neb_ext = 0;
    vector<string> ttype_str = { "Wavelength" };
    vector<string> tform_str;
    tform_str.push_back(lexical_cast<string>(nl) + "D");
    vector<string> tunit_str = { "Angstrom" };
    
    //For rectified spectrum
    vector<double> lambda_rs;
    vector<double>::size_type nl_rs = 0;
    
    int ncol=1;
    if (nebular != nullptr) {
      lambda_neb = nebular->lambda();
      nl_neb = lambda_neb.size();
      ttype_str.push_back("Wavelength_neb");
      tform_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr) {
      lambda_ext = extinct->lambda();
      nl_ext = lambda_ext.size();
      ttype_str.push_back("Wavelength_ex");
      tform_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    if (extinct != nullptr && nebular != nullptr) {
      lambda_neb_ext = extinct->lambda_neb();
      nl_neb_ext = lambda_neb_ext.size();
      ttype_str.push_back("Wavelength_neb_ex");
      tform_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;
    }
    // If there is a rectified spectrum present, make space for this
    if (specsyn->get_rectify())
    {
      lambda_rs = specsyn->get_recspec_wl();
      nl_rs = lambda_rs.size();
      ttype_str.push_back("Wavelength_Rectified");
      tform_str.push_back(lexical_cast<string>(nl_rs) + "D");
      tunit_str.push_back("Angstrom");
      ncol++;      
    }
    
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char wl_name[] = "Wavelength";

    // Create the table
    fits_create_tbl(outfiles.cluster_spec_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, wl_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.cluster_spec_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
    
    // Write wavelength data to table
    fits_write_col(outfiles.cluster_spec_fits, TDOUBLE, 1, 1, 1, nl,
		   lambda.data(), &fits_status);
    int col=2;
    if (nebular != nullptr) {
      fits_write_col(outfiles.cluster_spec_fits, TDOUBLE, col, 1, 1, nl_neb,
		     lambda_neb.data(), &fits_status);
      col++;
    }
    if (extinct != nullptr) {
      fits_write_col(outfiles.cluster_spec_fits, TDOUBLE, col, 1, 1, nl_ext,
		     lambda_ext.data(), &fits_status);
      col++;
    }
    if (nebular != nullptr && extinct != nullptr) {
      fits_write_col(outfiles.cluster_spec_fits, TDOUBLE, col, 1, 1, nl_neb_ext,
		     lambda_neb_ext.data(), &fits_status);
      col++;
    }

    if (specsyn->get_rectify())
    {
      fits_write_col(outfiles.cluster_spec_fits, TDOUBLE, col, 1, 1, nl_rs,
		     lambda_rs.data(), &fits_status);
      col++;
    }


    // Create a new table to hold the computed spectra
    char spec_name[] = "Spectra";
    vector<string> ttype2_str = 
      { "Trial", "UniqueID", "Time", "L_lambda" };
    vector<string> tform2_str = { "1K", "1K", "1D", "" };
    tform2_str[3] = lexical_cast<string>(nl) + "D";
    vector<string> tunit2_str = { "", "", "yr", "erg/s/A" };
    ncol = 4;
    if (nebular != nullptr) {
      ttype2_str.push_back("L_lambda_neb");
      tform2_str.push_back(lexical_cast<string>(nl_neb) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
    }
    if (extinct != nullptr) {
      ttype2_str.push_back("L_lambda_ex");
      tform2_str.push_back(lexical_cast<string>(nl_ext) + "D");
      tunit2_str.push_back("erg/s/A");
      ncol++;
      if (nebular != nullptr) {
	ttype2_str.push_back("L_lambda_neb_ex");
	tform2_str.push_back(lexical_cast<string>(nl_neb_ext) + "D");
	tunit2_str.push_back("erg/s/A");
	ncol++;
      }
    }
    
    //If there is a rectified spectrum present, make space for this
    if (specsyn->get_rectify())
    {
      ttype2_str.push_back("Rectified_Spec");
      tform2_str.push_back(lexical_cast<string>(nl_rs) + "D");
      tunit2_str.push_back(" ");
      ncol++;     
    }

    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    fits_create_tbl(outfiles.cluster_spec_fits, BINARY_TBL, 0, ncol,
		    ttype2, tform2, tunit2, spec_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated photometry file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_phot(slug_output_files &outfiles,
				    int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_integrated_phot";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.int_phot_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      outfiles.int_phot_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.int_phot_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err << "unable to open integrated photometry file "
			<< full_path.string()
			<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.int_phot_file.is_open()) {
      ostreams.slug_err << "unable to open intergrated photometry file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Grab the names and units of the photometric filters
  const vector<string> filter_names = filters->get_filter_names();
  const vector<string> filter_units = filters->get_filter_units();

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.int_phot_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    outfiles.int_phot_file << setw(21) << left << "Time";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.int_phot_file << setw(21) << left << filter_names[i];
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.int_phot_file << setw(21) << left << filter_names[i]+"_n";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.int_phot_file << setw(21) << left << filter_names[i]+"_ex";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  outfiles.int_phot_file << setw(21) << left << filter_names[i]+"_nex";
      }
    }
    outfiles.int_phot_file << endl;
    outfiles.int_phot_file << setw(21) << left << "(yr)";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.int_phot_file << setw(21) << left 
		    << "(" + filter_units[i] + ")";
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.int_phot_file << setw(21) << left 
		      << "(" + filter_units[i] + ")";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.int_phot_file << setw(21) << left 
		      << "(" + filter_units[i] + ")";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  outfiles.int_phot_file << setw(21) << left 
			<< "(" + filter_units[i] + ")";
      }
    }
    outfiles.int_phot_file << endl;
    outfiles.int_phot_file << setw(21) << left << "------------------";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.int_phot_file << setw(21) << left << "------------------";
    if (nebular != nullptr)
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.int_phot_file << setw(21) << left << "------------------";
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.int_phot_file << setw(21) << left << "------------------";
      if (nebular != nullptr)
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  outfiles.int_phot_file << setw(21) << left << "------------------";
    }
    outfiles.int_phot_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.int_phot_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with the number of filters and then the list
    // of filters names and units in ASCII; rest is binary
    outfiles.int_phot_file << filter_names.size() << endl;
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.int_phot_file << filter_names[i] << " " 
		    << filter_units[i] << endl;
    // First piece of binary information: two bits saying whether we're
    // using nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    outfiles.int_phot_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    outfiles.int_phot_file.write((char *) &use_extinct, sizeof use_extinct);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    vector<char *> ttype, tform, tunit;
    // First column is trial number, second column is time
    vector<string> ttype_str = { "Trial", "Time" };
    vector<string> form_str = { "1K", "1D" };
    vector<string> unit_str = { "", "yr" };
    for (int i=0; i<2; i++) {
      ttype.push_back(const_cast<char*>(ttype_str[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[i].c_str()));
      tunit.push_back(const_cast<char*>(unit_str[i].c_str()));
    }
    // Next n columns are filters
    for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
      ttype.push_back(const_cast<char*>(filter_names[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[1].c_str()));
      tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
    }
    int ncol = 2+filter_names.size();
    // If using nebular emission and/or extinction, add columns for
    // filters with those effects included
    vector<string> filter_names_neb(filter_names.size());
    vector<string> filter_names_ex(filter_names.size());
    vector<string> filter_names_neb_ex(filter_names.size());
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names_neb.size(); i++) {
	filter_names_neb[i] = filter_names[i] + "_neb";
	ttype.push_back(const_cast<char*>(filter_names_neb[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[1].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names_neb.size();
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	filter_names_ex[i] = filter_names[i] + "_ex";
	ttype.push_back(const_cast<char*>(filter_names_ex[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[1].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names_ex.size();
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names_neb_ex.size(); 
	     i++) {
	  filter_names_neb_ex[i] = filter_names[i] + "_neb_ex";
	  ttype.push_back(const_cast<char*>(filter_names_neb_ex[i].c_str()));
	  tform.push_back(const_cast<char*>(form_str[1].c_str()));
	  tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
	}
	ncol += filter_names_neb.size();
      }
    }
    // Create table
    int fits_status = 0;
    fits_create_tbl(outfiles.int_phot_fits, BINARY_TBL, 0, ncol,
		    ttype.data(), tform.data(), tunit.data(), nullptr, 
		    &fits_status);
    
    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.int_phot_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster photometry file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_phot(slug_output_files &outfiles,
				 int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_cluster_phot";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.cluster_phot_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
      full_path /= fname;
      outfiles.cluster_phot_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.cluster_phot_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err << "unable to open cluster photometry file "
			<< full_path.string()
			<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.cluster_phot_file.is_open()) {
      ostreams.slug_err << "unable to open cluster photometry file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Grab the names and units of the photometric filters
  const vector<string> filter_names = filters->get_filter_names();
  const vector<string> filter_units = filters->get_filter_units();

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.cluster_phot_file << "N_Trials = " << setfill('0')
		    << setw(14) << 0 << setfill(' ') << endl;
    outfiles.cluster_phot_file << setw(21) << left << "UniqueID"
		      << setw(21) << left << "Time";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.cluster_phot_file << setw(21) << left << filter_names[i];
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.cluster_phot_file << setw(21) << left << filter_names[i]+"_n";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.cluster_phot_file << setw(21) << left << filter_names[i]+"_ex";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  outfiles.cluster_phot_file << setw(21) << left << filter_names[i]+"_nex";
      }
    }
    outfiles.cluster_phot_file << endl;
    outfiles.cluster_phot_file << setw(21) << left << ""
		      << setw(21) << left << "(yr)";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.cluster_phot_file << setw(21) << left 
		    << "(" + filter_units[i] + ")";
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.cluster_phot_file << setw(21) << left 
			  << "(" + filter_units[i] + ")";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.cluster_phot_file << setw(21) << left 
			  << "(" + filter_units[i] + ")";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  outfiles.cluster_phot_file << setw(21) << left 
			    << "(" + filter_units[i] + ")";
      }
    }
    outfiles.cluster_phot_file << endl;
    outfiles.cluster_phot_file << setw(21) << left << "------------------"
		      << setw(21) << left << "------------------";
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.cluster_phot_file << setw(21) << left << "------------------";
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.cluster_phot_file << setw(21) << left << "------------------";
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	outfiles.cluster_phot_file << setw(21) << left << "------------------";
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++)
	  outfiles.cluster_phot_file << setw(21) << left << "------------------";
      }
    }
    outfiles.cluster_phot_file << endl;
  } else if (out_mode == BINARY) {
    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.cluster_phot_file.write((char *) &ntrial, sizeof ntrial);
    }
    // File starts with the number of filters and then the list
    // of filters names and units in ASCII; rest is binary
    outfiles.cluster_phot_file << filter_names.size() << endl;
    for (vector<string>::size_type i=0; i<filter_names.size(); i++)
      outfiles.cluster_phot_file << filter_names[i] << " " 
		    << filter_units[i] << endl;
    // First pieces of binary information: two bits saying whether we're
    // using nebular emission and extinction
    bool use_nebular = (nebular != nullptr);
    outfiles.cluster_phot_file.write((char *) &use_nebular, sizeof use_nebular);
    bool use_extinct = (extinct != nullptr);
    outfiles.cluster_phot_file.write((char *) &use_extinct, sizeof use_extinct);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    vector<char *> ttype, tform, tunit;
    // Columns are trial, uniqueID, time
    vector<string> ttype_str = { "Trial", "UniqueID", "Time" };
    vector<string> form_str = { "1K", "1K", "1D" };
    vector<string> unit_str = { "", "", "yr" };
    for (int i=0; i<3; i++) {
      ttype.push_back(const_cast<char*>(ttype_str[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[i].c_str()));
      tunit.push_back(const_cast<char*>(unit_str[i].c_str()));
    }
    // Next n columns are filters
    for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
      ttype.push_back(const_cast<char*>(filter_names[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[2].c_str()));
      tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
    }
    int ncol = 3+filter_names.size();
    // If using extinction and/or nebular emission, add columns for
    // filters with those effect included
    vector<string> filter_names_neb(filter_names.size());
    vector<string> filter_names_ex(filter_names.size());
    vector<string> filter_names_neb_ex(filter_names.size());
    if (nebular != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	filter_names_neb[i] = filter_names[i] + "_neb";
	ttype.push_back(const_cast<char*>(filter_names_neb[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[2].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names.size();
    }
    if (extinct != nullptr) {
      for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	filter_names_ex[i] = filter_names[i] + "_ex";
	ttype.push_back(const_cast<char*>(filter_names_ex[i].c_str()));
	tform.push_back(const_cast<char*>(form_str[2].c_str()));
	tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
      }
      ncol += filter_names.size();
      if (nebular != nullptr) {
	for (vector<string>::size_type i=0; i<filter_names.size(); i++) {
	  filter_names_neb_ex[i] = filter_names[i] + "_neb_ex";
	  ttype.push_back(const_cast<char*>(filter_names_neb_ex[i].c_str()));
	  tform.push_back(const_cast<char*>(form_str[2].c_str()));
	  tunit.push_back(const_cast<char*>(filter_units[i].c_str()));
	}
	ncol += filter_names.size();
      }
    }
    // Create table
    int fits_status = 0;
    fits_create_tbl(outfiles.cluster_phot_fits, BINARY_TBL, 0, ncol,
		    ttype.data(), tform.data(), tunit.data(), nullptr, 
		    &fits_status);
    
    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.cluster_phot_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster equivalent width file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_ew(slug_output_files &outfiles,
				 int chknum) 
{

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) 
  {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) 
  {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_cluster_ew";
  path full_path(pp.query<string>("out_dir"));
  
#ifdef ENABLE_FITS
  if (out_mode == FITS) 
  {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.cluster_ew_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) 
    {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err << "unable to open cluster equivalent width file "
			<< full_path.string()
			<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif


  // Grab the names of the lines
  const vector<string> line_names = lines->get_line_names();


  // Write header
#ifdef ENABLE_FITS
  if (out_mode == FITS) 
  {
    vector<char *> ttype, tform, tunit;
    // Columns are trial, uniqueID, time
    vector<string> ttype_str = { "Trial", "UniqueID", "Time" };
    vector<string> form_str = { "1K", "1K", "1D" };
    vector<string> unit_str = { "", "", "yr" };
    string ew_unit = "Angstrom";
    for (int i=0; i<3; i++) {
      ttype.push_back(const_cast<char*>(ttype_str[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[i].c_str()));
      tunit.push_back(const_cast<char*>(unit_str[i].c_str()));
    }
    // Next n columns are lines
    for (vector<string>::size_type i=0; i<line_names.size(); i++) 
    {
      ttype.push_back(const_cast<char*>(line_names[i].c_str()));
      tform.push_back(const_cast<char*>(form_str[2].c_str()));
      tunit.push_back(const_cast<char*>(ew_unit.c_str()));
    }
    int ncol = 3+line_names.size();
    
    // Create table
    int fits_status = 0;
    fits_create_tbl(outfiles.cluster_ew_fits, BINARY_TBL, 0, ncol,
		    ttype.data(), tform.data(), tunit.data(), nullptr, 
		    &fits_status);
    
    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.cluster_phot_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated supernovae file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_sn(slug_output_files &outfiles,
				  int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_integrated_sn";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.int_sn_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    outfiles.int_sn_file.open(full_path.c_str(), ios::out | ios::binary);
  } 
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.int_sn_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err
	<< "unable to open integrated supernovae file "
	<< full_path.string()
	<< "; cfitsio says: " << err_txt << endl;
      exit(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.int_sn_file.is_open()) {
      ostreams.slug_err << "unable to open intergrated supernovae file " 
			<< full_path.string() << endl;
      exit(1);
    }
#ifdef ENABLE_FITS
  }
#endif 

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.int_sn_file << "N_Trials = " << setfill('0')
			   << setw(14) << 0 << setfill(' ') << endl;
    outfiles.int_sn_file << setw(14) << left << "Time"
			 << setw(14) << left << "TotSN"
			 << setw(14) << left << "StochSN"
			 << endl;

    // Write unit line
    outfiles.int_sn_file << setw(14) << left << "(yr)"
			 << setw(14) << left << ""
			 << setw(14) << left << ""
			 << endl;
		  		  		  
    // Line of underscores
    outfiles.int_sn_file << setw(14) << left << "-----------"
			 << setw(14) << left << "-----------"
			 << setw(14) << left << "-----------"
			 << endl;
  
  } else if (out_mode == BINARY) {

    // State number of trials expected to be in this file
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.int_sn_file.write((char *) &ntrial, sizeof ntrial);
    }

  }
  
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    int ncol = 4;
    vector<string> ttype_str = 
      { "Trial", "Time", "TotSN", "StochSN" };
    vector<string> tform_str = 
      { "1K", "1D", "1D", "1K" };
    vector<string> tunit_str = 
      { "", "yr", "", "" };
    
    // Make things we can pass to c
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }

    // Create FITS table
    int fits_status = 0;
    fits_create_tbl(outfiles.int_sn_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, nullptr, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.int_sn_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster supernova file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_sn(slug_output_files &outfiles,
			       int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_cluster_sn";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.cluster_sn_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    outfiles.cluster_sn_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.cluster_sn_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err
	<< "unable to open cluster supernovae file "
	<< full_path.string()
	<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

  // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.cluster_sn_file.is_open()) {
      ostreams.slug_err << "unable to open cluster supernovae file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.cluster_sn_file << "N_Trials = " << setfill('0')
			       << setw(14) << 0 << setfill(' ') << endl;
    outfiles.cluster_sn_file << setw(14) << left << "UniqueID"
			     << setw(14) << left << "Time"
			     << setw(14) << left << "TotSN"
			     << setw(14) << left << "StochSN"
			     << endl;
    outfiles.cluster_sn_file << setw(14) << left << ""
			     << setw(14) << left << "(yr)"
			     << setw(14) << left << ""
			     << setw(14) << left << ""
			     << endl;
    outfiles.cluster_sn_file << setw(14) << left << "-----------"
			     << setw(14) << left << "-----------"
			     << setw(14) << left << "-----------"
			     << setw(14) << left << "-----------"
			     << endl;

  } else if (out_mode == BINARY) {

    // State number of trials in this file
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.cluster_sn_file.write((char *) &ntrial, sizeof ntrial);
    }

  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    // Note: this is pretty awkward -- we have to declare a vector of
    // string, then cast them to arrays of char *, because the cfitsio
    // library wants that. Unfortunately this awkwardness is the only
    // way to avoid lots of compiler warnings.
    int ncol = 5;
    vector<string> ttype_str = 
      { "Trial", "UniqueID", "Time", "TotSN", "StochSN" };
    vector<string> tform_str = 
      { "1K", "1K", "1D", "1D", "1K" };
    vector<string> tunit_str = 
      { "", "", "yr", "", "" };
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }

    // Create the table
    int fits_status = 0;
    fits_create_tbl(outfiles.cluster_sn_fits, BINARY_TBL, 0, ncol,
		    ttype, tform, tunit, nullptr, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.cluster_sn_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open integrated yield file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_integrated_yield(slug_output_files &outfiles,
				     int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_integrated_yield";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.int_yield_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    outfiles.int_yield_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.int_yield_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err << "unable to open integrated yield file "
			<< full_path.string()
			<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

    // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.int_yield_file.is_open()) {
      ostreams.slug_err << "unable to open integrated yield file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.int_yield_file << "N_Trials = " << setfill('0')
		     << setw(14) << 0 << setfill(' ') << endl;
    outfiles.int_yield_file << setw(14) << left << "Time"
          << setw(14) << left << "Symbol"
          << setw(14) << left << "Z"
          << setw(14) << left << "A"
          << setw(14) << left << "Yield";
    outfiles.int_yield_file << endl;
    outfiles.int_yield_file << setw(14) << left << "(yr)"
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << "(Msun)";
    outfiles.int_yield_file << endl;
    outfiles.int_yield_file << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------";
    outfiles.int_yield_file << endl;
  } else if (out_mode == BINARY) {

    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.int_yield_file.write((char *) &ntrial, sizeof ntrial);
    }

    // In binary mode, we begin the file by writing out the number of
    // isotopes, then for each isotope writing out its symbol as 4
    // characters (to maintain alignment), then writing out its atomic
    // number and weight as unsigned integers
    const vector<const isotope_data *>& isotopes = yields->get_isotopes();
    vector<isotope_data>::size_type niso = isotopes.size();
    outfiles.int_yield_file.write((char *) &niso, sizeof niso);
    for (vector<isotope_data>::size_type i=0; i<niso; i++) {
      char symbol[4];
      size_t len = isotopes[i]->symbol().copy(symbol, 4);
      for (size_t j = len; j<4; j++) symbol[j] = ' ';
      outfiles.int_yield_file.write(symbol, 4);
      unsigned int Z = isotopes[i]->num();
      outfiles.int_yield_file.write((char *) &Z, sizeof Z);
      unsigned int A = isotopes[i]->wgt();
      outfiles.int_yield_file.write((char *) &A, sizeof A);
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    
    // In FITS mode, write the isotope data to the first HDU

    // Grab data
    const vector<const isotope_data *>& isotopes = yields->get_isotopes();
    char **isotope_names = new char*[isotopes.size()];
    vector<long> isotope_Z(isotopes.size());
    vector<long> isotope_A(isotopes.size());
    for (vector<long>::size_type i=0; i<isotopes.size(); i++) {
      isotope_names[i] = new char[3];
      sprintf(isotope_names[i], "%3s", isotopes[i]->symbol().c_str());
      isotope_Z[i] = isotopes[i]->num();
      isotope_A[i] = isotopes[i]->wgt();
    }

    // Set up binary table
    vector<string> ttype_str = { "Name", "Z", "A" };
    vector<string> tform_str = { "3A", "1K", "1K" };
    vector<string> tunit_str = { "", "", "" };
    int ncol = 3;
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char table_name[] = "Isotope List";

    // Create the table
    fits_create_tbl(outfiles.int_yield_fits, BINARY_TBL, 0, ncol,
        ttype, tform, tunit, table_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.int_yield_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }

    // Write isotope data to table
    fits_write_col(outfiles.int_yield_fits, TSTRING, 1, 1, 1, isotopes.size(),
       isotope_names, &fits_status);
    fits_write_col(outfiles.int_yield_fits, TLONG, 2, 1, 1, isotopes.size(),
       isotope_Z.data(), &fits_status);
    fits_write_col(outfiles.int_yield_fits, TLONG, 3, 1, 1, isotopes.size(),
       isotope_A.data(), &fits_status);
    for (vector<long>::size_type i=0; i<isotopes.size(); i++)
      delete isotope_names[i];
    delete[] isotope_names;
    
    // Creat a new table to hold the computed yields
    vector<string> ttype2_str = 
      { "Trial", "Time", "Yield" };
    vector<string> tform2_str = 
      { "1K", "1D" };
    tform2_str.push_back(lexical_cast<string>(isotopes.size()) + "D");
    vector<string> tunit2_str = 
      { "", "yr", "Msun" };

    ncol=3;
    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (vector<double>::size_type i=0; i<(unsigned int)ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    char yield_table_name[] = "Yield";

    // Create the table
    fits_create_tbl(outfiles.int_yield_fits, BINARY_TBL, 0, ncol,
        ttype2, tform2, tunit2, yield_table_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Open cluster yield file and write its header
////////////////////////////////////////////////////////////////////////
void slug_sim::open_cluster_yield(slug_output_files &outfiles,
				  int chknum) {

  // Construct file name and path
  string fname(pp.query<string>("model_name"));
#if defined(ENABLE_MPI) && !(MPI_VERSION == 1 || MPI_VERSION == 2)
  if (comm != MPI_COMM_NULL) {
    ostringstream ss;
    ss << "_" << setfill('0') << setw(4) << rank;
    fname += ss.str();
  }
#endif
  if (chknum >= 0) {
    ostringstream ss;
    ss << "_chk" << setfill('0') << setw(4) << chknum;
    fname += ss.str();
  }
  fname += "_cluster_yield";
  path full_path(pp.query<string>("out_dir"));
  if (out_mode == ASCII) {
    fname += ".txt";
    full_path /= fname;
    outfiles.cluster_yield_file.open(full_path.c_str(), ios::out);
  } else if (out_mode == BINARY) {
    fname += ".bin";
    full_path /= fname;
    outfiles.cluster_yield_file.open(full_path.c_str(), ios::out | ios::binary);
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    fname += ".fits";
    full_path /= fname;
    string fname_tmp = "!" + full_path.string();
    int fits_status = 0;
    fits_create_file(&outfiles.cluster_yield_fits, fname_tmp.c_str(), &fits_status);
    if (fits_status) {
      char err_txt[80] = "";
      fits_read_errmsg(err_txt);
      ostreams.slug_err << "unable to open integrated yield file "
			<< full_path.string()
			<< "; cfitsio says: " << err_txt << endl;
      bailout(1);
    }
  }
#endif

    // Make sure file is open
#ifdef ENABLE_FITS
  if (out_mode != FITS) {
#endif
    if (!outfiles.cluster_yield_file.is_open()) {
      ostreams.slug_err << "unable to open cluster yield file " 
			<< full_path.string() << endl;
      bailout(1);
    }
#ifdef ENABLE_FITS
  }
#endif

  // Write header
  if (out_mode == ASCII) {
    if (chknum >= 0)
      outfiles.cluster_yield_file << "N_Trials = " << setfill('0')
			 << setw(14) << 0 << setfill(' ') << endl;
    outfiles.cluster_yield_file << setw(14) << left << "UniqueID"
          << setw(14) << left << "Time"
          << setw(14) << left << "Symbol"
          << setw(14) << left << "Z"
          << setw(14) << left << "A"
          << setw(14) << left << "Yield";
    outfiles.cluster_yield_file << endl;
    outfiles.cluster_yield_file << setw(14) << left << ""
          << setw(14) << left << "(yr)"
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << ""
          << setw(14) << left << "(Msun)";
    outfiles.cluster_yield_file << endl;
    outfiles.cluster_yield_file << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------"
          << setw(14) << left << "-----------";
    outfiles.cluster_yield_file << endl;
  } else if (out_mode == BINARY) {

    // State number of trials contained in checkpoints
    if (chknum >= 0) {
      unsigned int ntrial = 0;
      outfiles.cluster_yield_file.write((char *) &ntrial, sizeof ntrial);
    }
    
    // In binary mode, we begin the file by writing out the number of
    // isotopes, then for each isotope writing out its symbol as 4
    // characters (to maintain alignment), then writing out its atomic
    // number and weight as unsigned integers
    const vector<const isotope_data *>& isotopes = yields->get_isotopes();
    vector<isotope_data>::size_type niso = isotopes.size();
    outfiles.cluster_yield_file.write((char *) &niso, sizeof niso);
    for (vector<isotope_data>::size_type i=0; i<niso; i++) {
      char symbol[4];
      size_t len = isotopes[i]->symbol().copy(symbol, 4);
      for (size_t j = len; j<4; j++) symbol[j] = ' ';
      outfiles.cluster_yield_file.write(symbol, 4);
      unsigned int Z = isotopes[i]->num();
      outfiles.cluster_yield_file.write((char *) &Z, sizeof Z);
      unsigned int A = isotopes[i]->wgt();
      outfiles.cluster_yield_file.write((char *) &A, sizeof A);
    }
  }
#ifdef ENABLE_FITS
  else if (out_mode == FITS) {
    
    // In FITS mode, write the isotope data to the first HDU

    // Grab data
    const vector<const isotope_data *>& isotopes = yields->get_isotopes();
    char **isotope_names = new char*[isotopes.size()];
    vector<long> isotope_Z(isotopes.size());
    vector<long> isotope_A(isotopes.size());
    for (vector<long>::size_type i=0; i<isotopes.size(); i++) {
      isotope_names[i] = new char[3];
      sprintf(isotope_names[i], "%3s", isotopes[i]->symbol().c_str());
      isotope_Z[i] = isotopes[i]->num();
      isotope_A[i] = isotopes[i]->wgt();
    }

    // Set up binary table
    vector<string> ttype_str = { "Name", "Z", "A" };
    vector<string> tform_str = { "3A", "1K", "1K" };
    vector<string> tunit_str = { "", "", "" };
    int ncol = 3;
    char **ttype = new char *[ncol];
    char **tform = new char *[ncol];
    char **tunit = new char *[ncol];
    for (int i=0; i<ncol; i++) {
      ttype[i] = const_cast<char*>(ttype_str[i].c_str());
      tform[i] = const_cast<char*>(tform_str[i].c_str());
      tunit[i] = const_cast<char*>(tunit_str[i].c_str());
    }
    int fits_status = 0;
    char table_name[] = "Isotope List";

    // Create the table
    fits_create_tbl(outfiles.cluster_yield_fits, BINARY_TBL, 0, ncol,
        ttype, tform, tunit, table_name, &fits_status);
    delete [] ttype;
    delete [] tform;
    delete [] tunit;

    // Add a keyword stating the expected number of trials if this is
    // a checkpoint file
    if (chknum >= 0) {
      unsigned int ntrials = 0;
      fits_write_key(outfiles.cluster_yield_fits, TUINT, "N_Trials", 
		     &ntrials, "Number of trials in file",
		     &fits_status);
    }
    
    // Write isotope data to table
    fits_write_col(outfiles.cluster_yield_fits, TSTRING, 1, 1, 1, isotopes.size(),
       isotope_names, &fits_status);
    fits_write_col(outfiles.cluster_yield_fits, TLONG, 2, 1, 1, isotopes.size(),
       isotope_Z.data(), &fits_status);
    fits_write_col(outfiles.cluster_yield_fits, TLONG, 3, 1, 1, isotopes.size(),
       isotope_A.data(), &fits_status);
    for (vector<long>::size_type i=0; i<isotopes.size(); i++)
      delete isotope_names[i];
    delete[] isotope_names;
    
    // Creat a new table to hold the computed yields
    vector<string> ttype2_str = 
      { "Trial", "UniqueID", "Time", "Yield" };
    vector<string> tform2_str = 
      { "1K", "1K", "1D" };
    tform2_str.push_back(lexical_cast<string>(isotopes.size()) + "D");
    vector<string> tunit2_str = 
      { "", "", "yr", "Msun" };

    ncol=4;
    char **ttype2 = new char *[ncol];
    char **tform2 = new char *[ncol];
    char **tunit2 = new char *[ncol];
    for (vector<double>::size_type i=0; i<(unsigned int) ncol; i++) {
      ttype2[i] = const_cast<char*>(ttype2_str[i].c_str());
      tform2[i] = const_cast<char*>(tform2_str[i].c_str());
      tunit2[i] = const_cast<char*>(tunit2_str[i].c_str());
    }
    char yield_table_name[] = "Yield";

    // Create the table
    fits_create_tbl(outfiles.cluster_yield_fits, BINARY_TBL, 0, ncol,
        ttype2, tform2, tunit2, yield_table_name, &fits_status);
    delete [] ttype2;
    delete [] tform2;
    delete [] tunit2;
  }
#endif
}


////////////////////////////////////////////////////////////////////////
// Write out a separator
////////////////////////////////////////////////////////////////////////
void slug_sim::write_separator(std::ofstream& file, 
			       const unsigned int width) {
  string sep;
  for (unsigned int i=0; i<width; i++) sep += "-";
  file << sep << endl;
}
