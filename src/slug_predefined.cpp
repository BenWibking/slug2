#include "slug_predefined.H"
#include "specsyn/slug_specsyn_hillier.H"
#include "specsyn/slug_specsyn_kurucz.H"
#include "specsyn/slug_specsyn_pauldrach.H"
#include "specsyn/slug_specsyn_planck.H"
#include "specsyn/slug_specsyn_sb99.H"
#include "tracks/slug_tracks_sb99.H"
#include "yields/slug_yields_multiple.H"
#include "fcntl.h"
#include <cstdlib>
#include <ctime>
#include <string>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////
slug_predefined::slug_predefined() : rng(nullptr)
{ 
  // Get location of SLUG_DIR
  char *slug_dir_ptr = getenv("SLUG_DIR");
  string slug_dir;
  if (slug_dir_ptr != NULL)
    slug_dir = slug_dir_ptr;
  if (slug_dir.length() == 0) slug_dir = current_path().string();
  path slug_path(slug_dir);

  // Set default locations of data
  path lib_path = slug_path / "lib";
  imf_dir = lib_path / "imf";
  track_dir = lib_path / "tracks";
  atmos_dir = lib_path / "atmospheres";
  yield_dir = lib_path / "yields";
  filter_dir = lib_path / "filters";

  // Load IMFs
  vector<string> predefined_imfs = {
    "chabrier", "chabrier03", "kroupa", "kroupa_T100", "kroupa_T50",
    "kroupa_T80", "kroupa_sb99", "salpeter", "wk06" };
  for (vector<string>::size_type i=0; i<predefined_imfs.size(); i++)
    known_imfs[predefined_imfs[i]] = nullptr;

  // Load tracks
  vector<string> predefined_tracks = {
    "Z0020v00.txt", "Z0020v40.txt", "Z0140v00.txt", "Z0140v40.txt",
    "modc001.dat", "modc004.dat", "modc008.dat", "modc020.dat",
    "modc040.dat", "mode001.dat", "mode004.dat", "mode008.dat",
    "mode020.dat", "mode040.dat", "modp0004.dat", "modp004.dat",
    "modp020.dat", "modp050.dat", "mods0004.dat", "mods004.dat",
    "mods008.dat", "mods020.dat", "mods050.dat" };
  for (vector<string>::size_type i=0; i<predefined_tracks.size(); i++)
    known_tracks[predefined_tracks[i]] = nullptr;    

  // Load specsyns
  vector<string> predefined_specsyn = {
    "planck", "kurucz", "kurucz_hillier", "kurucz_pauldrach", "sb99" };
  for (vector<string>::size_type i=0; i<predefined_specsyn.size(); i++)
    known_specsyn[predefined_specsyn[i]] = nullptr;

  // Load filter sets
  vector<string> predefined_filter_sets = {
    "Lbol", "QH0", "QHe0", "QHe1", "Lbol+QH0", "Lbol+QH0+QHe0+QHe1" };
  for (vector<string>::size_type i=0; i<predefined_filter_sets.size(); i++)
    known_filter_sets[predefined_filter_sets[i]] = nullptr;

  // Load yields
  vector<string> predefined_yields = {
    "SNII_Sukhbold16", "SNII_Sukhbold16_nodecay" };
  for (vector<string>::size_type i=0; i<predefined_yields.size(); i++)
    known_yields[predefined_yields[i]] = nullptr;
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_predefined::~slug_predefined() {
  // Free the rng
  if (!rng) delete rng;

  // Free any IMFs we've allocated
  for (map<const string, const slug_PDF*>::iterator
	 it = known_imfs.begin();
       it != known_imfs.end(); it++) {
    if (it->second) delete it->second;
  }

  // Free any tracks we've allocated
  for (map<const string, const slug_tracks*>::iterator
	 it = known_tracks.begin();
       it != known_tracks.end(); it++) {
    if (it->second) delete it->second;
  }

  // Free any spectral synthesizers
  for (map<const string, const slug_specsyn*>::iterator
	 it = known_specsyn.begin();
       it != known_specsyn.end(); it++) {
    if (it->second) delete it->second;
  }

  // Free filter sets
  for (map<const string, const slug_filter_set*>::iterator
	 it = known_filter_sets.begin();
       it != known_filter_sets.end(); it++) {
    if (it->second) delete it->second;
  }

  // Free yields
  for (map<const string, const slug_yields*>::iterator
	 it = known_yields.begin();
       it != known_yields.end(); it++) {
    if (it->second) delete it->second;
  } 
}


////////////////////////////////////////////////////////////////////////
// Method to construct a rng
////////////////////////////////////////////////////////////////////////
void slug_predefined::build_rng() {

  // Do nothing if rng already exists
  if (rng) return;
  
  // Get a random see from /dev/urandom if possible
  unsigned int seed;
  int fn;
  bool rand_set = false;
  fn = open("/dev/urandom", O_RDONLY);
  if (fn != -1) {
    rand_set = (read(fn, &seed, 4) == 4); // True if read succeeds
    close(fn);
  }
  if (!rand_set) {
    // Failed to set from /dev/urandom; seed using system time instead.
    seed = static_cast<unsigned int>(time(0));
  }

  // Set up the random number generator
  rng = new rng_type(seed);
}

////////////////////////////////////////////////////////////////////////
// Method to construct slug_PDF objects
////////////////////////////////////////////////////////////////////////
inline const slug_PDF*
slug_predefined::build_PDF(const string& fname,
			   const double min_stoch_mass,
			   const samplingMethod method) {
  build_rng();
  slug_PDF *pdf = new slug_PDF(fname.c_str(), rng, ostreams);
  if (min_stoch_mass > 0.0) pdf->set_stoch_lim(min_stoch_mass);
  if (method != NO_METHOD) pdf->setMethod(method);
  return pdf;
}

////////////////////////////////////////////////////////////////////////
// Method to construct IMFs, tracks, filter sets
////////////////////////////////////////////////////////////////////////
inline const slug_PDF*
slug_predefined::build_IMF(const string& imfname,
			   const double min_stoch_mass,
			   const samplingMethod method) {
  path p(imfname);
  return build_PDF((imf_dir / p).string() + ".imf",
		   min_stoch_mass, method);
}

inline const slug_tracks*
slug_predefined::build_tracks(const string&trackname) {
  path p(trackname);
  // This is a quick-fix to get sb99 working again; will do a more
  // permanent patch later
  p = path("sb99") / p;
  build_rng();
  return (slug_tracks *)
    new slug_tracks_sb99((track_dir / p).string().c_str(), ostreams);
}

inline const slug_filter_set*
slug_predefined::build_filter_set(const string& filter_set_name) {
  char_separator<char> sep("+");
  tokenizer< char_separator<char> > fnames(filter_set_name, sep);
  vector<string> filter_names;
  for (tokenizer< char_separator<char> >::iterator it = fnames.begin();
       it != fnames.end(); ++it) {
    filter_names.push_back(*it);
  }
  return new slug_filter_set(filter_names, filter_dir.string().c_str(),
			     L_NU, ostreams, nullptr);
}

////////////////////////////////////////////////////////////////////////
// Methods to fetch IMFs, tracks, specsyn, yields
////////////////////////////////////////////////////////////////////////
const slug_PDF* slug_predefined::imf(const string& imfname,
				     const double min_stoch_mass,
				     const samplingMethod method) {
  map<const string, const slug_PDF*>::iterator it
    = known_imfs.find(imfname);
  if (it != known_imfs.end()) {
    if (!known_imfs[imfname])
      known_imfs[imfname] = build_IMF(imfname, min_stoch_mass, method);
    else {
      assert(known_imfs[imfname]->get_xStochMin() == min_stoch_mass);
      if (method != NO_METHOD)
	assert(known_imfs[imfname]->getMethod() == method);
    }
    return known_imfs[imfname];
  }
  else return nullptr;
}

const slug_tracks* slug_predefined::tracks(const string& trackname) {
  map<const string, const slug_tracks*>::iterator it
    = known_tracks.find(trackname);
  if (it != known_tracks.end()) {
    if (!known_tracks[trackname])
      known_tracks[trackname] = build_tracks(trackname);
    return known_tracks[trackname];
  }
  else return nullptr;
}

const slug_specsyn*
slug_predefined::specsyn(const string& specsyn_name,
			 const slug_tracks* tracks_,
			 const slug_PDF* imf_) {
  map<const string, const slug_specsyn*>::iterator it
    = known_specsyn.find(specsyn_name);
  if (it != known_specsyn.end()) {
    if (!known_specsyn[specsyn_name]) {
      if (specsyn_name == "planck") {
	known_specsyn[specsyn_name] = static_cast<slug_specsyn *>
	  (new slug_specsyn_planck(tracks_, imf_, nullptr, ostreams));
      } else if (specsyn_name == "kurucz") {
	known_specsyn[specsyn_name] = static_cast<slug_specsyn *>
	  (new slug_specsyn_kurucz(atmos_dir.string().c_str(),
				   tracks_, imf_,
				   nullptr, ostreams));
      } else if (specsyn_name == "kurucz_hillier") {
	known_specsyn[specsyn_name] = static_cast<slug_specsyn *>
	  (new slug_specsyn_hillier(atmos_dir.string().c_str(),
				    tracks_, imf_,
				    nullptr, ostreams));
      } else if (specsyn_name == "kurucz_pauldrach") {
	known_specsyn[specsyn_name] = static_cast<slug_specsyn *>
	  (new slug_specsyn_pauldrach(atmos_dir.string().c_str(),
				      tracks_, imf_,
				      nullptr, ostreams));
      } else if (specsyn_name == "sb99") {
	known_specsyn[specsyn_name] = static_cast<slug_specsyn *>
	  (new slug_specsyn_sb99(atmos_dir.string().c_str(),
				 tracks_, imf_,
				 nullptr, ostreams));
      }
    }
    return known_specsyn[specsyn_name];
  } else return nullptr;
}

const slug_filter_set*
slug_predefined::filter_set(const string& filter_set_name) {
  map<const string, const slug_filter_set*>::iterator it
    = known_filter_sets.find(filter_set_name);
  if (it != known_filter_sets.end()) {
    if (!known_filter_sets[filter_set_name])
      known_filter_sets[filter_set_name]
	= build_filter_set(filter_set_name);
    return known_filter_sets[filter_set_name];
  }
  else return nullptr;
}

const slug_yields *slug_predefined::yields(const string& yields_name) {
  map<const string, const slug_yields*>::iterator it
    = known_yields.find(yields_name);
  if (it != known_yields.end()) {
    if (!known_yields[yields_name]) {
      if (yields_name == "SNII_Sukhbold16") {
	known_yields[yields_name] = (slug_yields *)
	  new slug_yields_multiple(yield_dir.string().c_str(),
				   SNII_SUKHBOLD16, 1.0, ostreams);
      } else if (yields_name == "SNII_Sukhbold16_nodecay") {
	known_yields[yields_name] = (slug_yields *)
	  new slug_yields_multiple(yield_dir.string().c_str(),
				   SNII_SUKHBOLD16, 1.0, ostreams,
				   true);
      }
    }
    return known_yields[yields_name];
  } else return nullptr;
}

////////////////////////////////////////////////////////////////////////
// Define slug_predef
////////////////////////////////////////////////////////////////////////
slug_predefined slug_predef;
