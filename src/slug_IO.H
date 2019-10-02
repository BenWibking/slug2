/*********************************************************************
Copyright (C) 2014-7 Robert da Silva, Michele Fumagalli, Mark Krumholz
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
// This file defines classes and structures that manage IO
////////////////////////////////////////////////////////////////////////

#ifndef _slug_IO_H_
#define _slug_IO_H_

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
    typedef decltype(nullptr) nullptr_t;
}
#endif

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#ifdef ENABLE_MPI
#   include "mpi.h"
#endif
#ifdef ENABLE_FITS
extern "C" {
#   include "fitsio.h"
}
#endif
#include "slug.H"

////////////////////////////////////////////////////////////////////////
// class slug_output_files
//
// A class to hold pointers to the output files and information about
// them; created just for the convenience of passing this around
// instead of a bunch of individual file pointers. It has no methods,
// and we only define it as a class instead of a struct so that we can
// automatically initialize the fits pointers to null
////////////////////////////////////////////////////////////////////////
class slug_output_files {

public:

  slug_output_files() :
    is_open(false)
#ifdef ENABLE_FITS
    , int_prop_fits(nullptr), cluster_prop_fits(nullptr)
    , int_spec_fits(nullptr), cluster_spec_fits(nullptr)
    , int_phot_fits(nullptr), cluster_phot_fits(nullptr)
    , int_sn_fits(nullptr), cluster_sn_fits(nullptr)
    , int_yield_fits(nullptr), cluster_yield_fits(nullptr)
                             , cluster_ew_fits(nullptr)
#endif
  { }

  bool is_open;                     // Are files open
  std::ofstream int_prop_file;      // Integrated properties file
  std::ofstream cluster_prop_file;  // Cluster properties file
  std::ofstream int_spec_file;      // Integrated spectra file
  std::ofstream cluster_spec_file;  // Cluster spectra file
  std::ofstream int_phot_file;      // Integrated photometry file
  std::ofstream cluster_phot_file;  // Cluster photometry file
  std::ofstream int_sn_file;        // Integrated supernovae file
  std::ofstream cluster_sn_file;    // Cluster supernovae file
  std::ofstream int_yield_file;     // Integrated yield file
  std::ofstream cluster_yield_file; // Cluster yield file
#ifdef WINDS_ON
  std::ofstream cluster_winds_file; // Cluster winds file
#endif
#ifdef ENABLE_FITS
  // FITS file versions of the above
  fitsfile *int_prop_fits;
  fitsfile *cluster_prop_fits;
  fitsfile *int_spec_fits;
  fitsfile *cluster_spec_fits;
  fitsfile *int_phot_fits;
  fitsfile *cluster_phot_fits;
  fitsfile *int_sn_fits;
  fitsfile *cluster_sn_fits;
  fitsfile *int_yield_fits;
  fitsfile *cluster_yield_fits;
  fitsfile *cluster_ew_fits;
#ifdef WINDS_ON
  fitsfile *cluster_winds_fits;
#endif
#endif
};


////////////////////////////////////////////////////////////////////////
// class slug_prefixbuf
//
// This defines a class that appends a prefix to output to a stream
// whenever a newline is detected; the design of this class borrows
// heavily from
// http://stackoverflow.com/questions/27336335/c-cout-with-prefix
////////////////////////////////////////////////////////////////////////
class slug_prefixbuf : public std::streambuf {
  
public:
  slug_prefixbuf(std::string const& prefix, std::streambuf* sbuf)
    : prefix(prefix), sbuf(sbuf), need_prefix(true) { }

protected:
  virtual int sync();
  virtual int overflow(int c);
  
  std::string     prefix;
  std::streambuf* sbuf;
  bool            need_prefix;
};

////////////////////////////////////////////////////////////////////////
// class slug_oprefixstream
//
// This is an output stream class derived from the slug_prefixbuf
// class; it has an MPI version and a non-MPI version. The non-MPI
// version is trivial, and just defaults back to the parent class. The
// MPI version takes care of synchronizing outputs between
// processors so that they don't wind up writing over each other when
// trying to access stdout or stderr.
////////////////////////////////////////////////////////////////////////
  
#ifdef ENABLE_MPI

#define SLUG_MAX_LINE_SIZE 16384

// MPI version
class slug_oprefixstream : private virtual slug_prefixbuf,
			   public std::ostream {
public:
  slug_oprefixstream(std::string const& prefix_, std::ostream& out,
                     MPI_Comm comm_ = MPI_COMM_NULL,
                     bool only_one_proc_ = false,
                     int tag_ = 0);
  ~slug_oprefixstream() {
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) flush_communication();
  }
  // This routine forces all communications to finish, and blocks
  // until they do
  void flush_communication();
private:
  virtual int overflow(int c);
  MPI_Comm comm;         // Communicator
  int rank;              // MPI rank
  const int tag;         // Tag to use for this IO objects messages
  bool only_one_proc;    // Type of IO; from all procs or only proc 0
  std::ostream& target;  // Stream to which output eventually goes
  std::string outgoing_buf;  // Buffer for communication
  std::streamsize outgoing_ptr; // Pointer in outgoing buffer
  std::vector<MPI_Request> req; // Storage for MPI requests
  std::vector<std::streamsize> msg_size; // Sizes of messages
};

#else

// Non-MPI version; this is just a trival pass-through
class slug_oprefixstream : private virtual slug_prefixbuf,
			   public std::ostream {
public:
  slug_oprefixstream(std::string const& prefix_, std::ostream& out)
    : slug_prefixbuf(prefix_, out.rdbuf()),
      std::ios(static_cast<std::streambuf*>(this)),
      std::ostream(static_cast<std::streambuf*>(this)) { }
};
#endif


////////////////////////////////////////////////////////////////////////
// class slug_ostreams
//
// This is a class that just holds all the IO handlers we define in a
// standard run
////////////////////////////////////////////////////////////////////////
class slug_ostreams {

private:
  slug_oprefixstream *slug_out_ptr, *slug_err_ptr, *slug_warn_ptr;
#ifdef ENABLE_MPI
  slug_oprefixstream *slug_out_one_ptr, *slug_err_one_ptr,
    *slug_warn_one_ptr;
#endif

public:
  slug_oprefixstream &slug_out, &slug_err, &slug_warn,
    &slug_out_one, &slug_err_one, &slug_warn_one;

#ifdef ENABLE_MPI
  // Constructor version where we provide the MPI communicator
  slug_ostreams(MPI_Comm comm) :
    slug_out_ptr(new
      slug_oprefixstream("slug: ", std::cout, comm, false, 101)),
    slug_err_ptr(new
      slug_oprefixstream("slug: error: ", std::cerr, comm, false, 102)),
    slug_warn_ptr(new
      slug_oprefixstream("slug: warning: ", std::cerr, comm, false, 103)),
    slug_out_one_ptr(new
      slug_oprefixstream("slug: ", std::cout, comm, true)),
    slug_err_one_ptr(new
      slug_oprefixstream("slug: error: ", std::cerr, comm, true)),
    slug_warn_one_ptr(new
      slug_oprefixstream("slug: warning: ", std::cerr, comm, true)),
    slug_out(*slug_out_ptr),
    slug_err(*slug_err_ptr),
    slug_warn(*slug_warn_ptr),
    slug_out_one(*slug_out_one_ptr),
    slug_err_one(*slug_err_one_ptr),
    slug_warn_one(*slug_warn_one_ptr)
  { }    
#endif
  
  // Non-MPI version
  slug_ostreams() :
    slug_out_ptr(new
      slug_oprefixstream("slug: ", std::cout)),
    slug_err_ptr(new
      slug_oprefixstream("slug: error: ", std::cerr)),
    slug_warn_ptr(new
      slug_oprefixstream("slug: warning: ", std::cerr)),
    slug_out(*slug_out_ptr),
    slug_err(*slug_err_ptr),
    slug_warn(*slug_warn_ptr),
    slug_out_one(*slug_out_ptr),
    slug_err_one(*slug_err_ptr),
    slug_warn_one(*slug_warn_ptr)
  { }

  // Explicitly disallow the copy constructor; this object should
  // only ever be passed by reference
  slug_ostreams(const slug_ostreams& ostreams_) = delete;
    
  ~slug_ostreams() {
    delete slug_out_ptr;
    delete slug_err_ptr;
    delete slug_warn_ptr;
#ifdef ENABLE_MPI
    delete slug_out_one_ptr;
    delete slug_err_one_ptr;
    delete slug_warn_one_ptr;
#endif
  }

#ifdef ENABLE_MPI
  void flush_communication() {
    slug_out_ptr->flush_communication();
    slug_err_ptr->flush_communication();
    slug_warn_ptr->flush_communication();
    slug_out_one_ptr->flush_communication();
    slug_err_one_ptr->flush_communication();
    slug_warn_one_ptr->flush_communication();
  }
#endif
};

#endif
// _slug_IO_H_

