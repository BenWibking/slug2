# Main makefile for slug2

# Did the user tell us to use a particular machine? If so, use that.
ifdef MACHINE
     include Make.mach.$(MACHINE)
else
     # Machine not specified, so try to guess
     UNAME		= $(shell uname)
     UNAMEN		= $(shell uname -n)

     # Do we have a makefile that matches the machine name? If so, use
     # that. If not, use a generic makefile depending on the OS.
     ifeq ($(UNAMEN), hyades.ucsc.edu)
          include Make.mach.ucsc-hyades
     else ifeq ($(UNAMEN), cosma-e)
          include Make.mach.icc-cosma4
     else ifeq ($(UNAMEN), avatar)
          include Make.mach.avatar
     else ifeq ($(UNAMEN), mosura)
          include Make.mach.coala
     else ifneq (,$(findstring raijin, $(UNAMEN)))
          include Make.mach.icc-raijin
     else ifneq (,$(findstring gadi, $(UNAMEN)))
          include Make.mach.icc-gadi
     else ifeq ($(UNAME), Linux)
          include Make.mach.linux-gnu
     else ifeq ($(UNAME), Darwin)
          include Make.mach.darwin
     else
          $(info Cannot detect system type. Suggest you specify MACHINE= manually.)
          include Make.mach.generic
     endif
endif

# Set compiler
CXX		= $(MACH_CXX)
MPICXX		= $(MACH_MPICXX)

# Set optimization mode flags
CXXOPTFLAGS	= $(MACH_CXXOPTFLAGS) $(MACH_C11FLAG) -DNDEBUG \
	-DBOOST_DISABLE_ASSERTS -DHAVE_INLINE \
	-MMD -MP
LDOPTFLAGS	= $(MACH_LDOPTFLAGS) $(MACH_CXXFLAG)

# Set debug mode flags
CXXDEBFLAGS     = $(MACH_CXXDEBFLAGS) $(MACH_C11FLAG) -MMD -MP
LDDEBFLAGS	= $(MACH_LDDEBFLAGS) $(MACH_CXXFLAG)

# Read any user overrides
-include Make.config.override

# Include flags
ifdef BOOST_HDR_PATH
     INCFLAGS += -I$(BOOST_HDR_PATH)
endif
ifdef GSL_HDR_PATH
     INCFLAGS += -I$(GSL_HDR_PATH)
endif
ifdef C_INCLUDE_PATH
     INCFLAGS += -I$(subst :, -I,$(C_INCLUDE_PATH))
endif
ifdef CXX_INCLUDE_PATH
     INCFLAGS += -I$(subst :, -I,$(CXX_INCLUDE_PATH))
endif

# Link flags
LDLIBFLAGS      = -lgsl -lgslcblas -lboost_system$(MACH_BOOST_TAG) \
	-lboost_filesystem$(MACH_BOOST_TAG) \
	-lboost_regex$(MACH_BOOST_TAG)
ifdef BOOST_LIB_PATH
     LDLIBFLAGS += -L$(BOOST_LIB_PATH)
endif
ifdef GSL_LIB_PATH
     LDLIBFLAGS += -L$(GSL_LIB_PATH)
endif
ifdef LD_LIBRARY_PATH
     LDLIBFLAGS += -L$(subst :, -L,$(LD_LIBRARY_PATH))
endif
ifdef LIBRARY_PATH
     LDLIBFLAGS += -L$(subst :, -L,$(LIBRARY_PATH))
endif

# Turn FITS on or off
ifeq ($(FITS), ENABLE_FITS)
     LDLIBFLAGS += -lcfitsio
     DEFINES    += -DENABLE_FITS
endif
ifdef FITS_LIB_PATH
     LDLIBFLAGS += -L$(FITS_LIB_PATH)
endif

# Specify which version of the GSL we're using
DEFINES += -DGSLVERSION=$(GSLVERSION)

# Switch to MPI compiler if requested to compile with MPI enables
ifeq ($(MPI), ENABLE_MPI)
     CXX        =  $(MPICXX)
     DEFINES    += -DENABLE_MPI
endif

CXXFLAGS +=  $(INCFLAGS) $(DEFINES)
LDFLAGS  +=  $(LDLIBFLAGS)

# Name for executable
EXENAME		= slug

# Name for compiled library, shared or static
LIBNAME		= lib$(EXENAME)$(LIB_EXTENSION)
LIBSTATNAME	= lib$(EXENAME)$(LIBSTAT_EXTENSION)

# List of directories containing source files
VPATH		= filters interpolators lines pdfs specsyn tracks \
	utils yields

# Pointers to sources and objects
SOURCES		= $(wildcard *.cpp) $(wildcard */*.cpp)
OBJECTS		= $(SOURCES:%.cpp=%.o)
DEPS		= $(SOURCES:%.cpp=%.d)

# Default target
.PHONY: all exe debug lib lib-debug libstatic libstatic-debug compile clean

exe: CXXFLAGS += $(CXXOPTFLAGS)
exe: LDFLAGS  += $(LDOPTFLAGS)
exe: $(EXENAME) 

debug: CXXFLAGS += $(CXXDEBFLAGS)
debug: LDFLAGS  += $(LDDEBFLAGS)
debug: $(EXENAME)

lib: CXXFLAGS += $(CXXOPTFLAGS) $(CLIBFLAGS)
lib: LDFLAGS  += $(LDOPTFLAGS)
lib: $(LIBNAME)

lib-debug: CXXFLAGS += $(CXXDEBFLAGS) $(CLIBFLAGS)
lib-debug: LDFLAGS  += $(LDDEBFLAGS)
lib-debug: $(LIBNAME)

libstatic: CXXFLAGS += $(CXXOPTFLAGS) $(CLIBFLAGS)
libstatic: $(LIBSTATNAME)

libstatic-debug: CXXFLAGS += $(CXXDEBFLAGS) $(CLIBFLAGS)
libstatic-debug: $(LIBSTATNAME)

all: exe
all: lib

# Include dependencies
-include $(DEPS)

$(EXENAME): $(OBJECTS)
	$(CXX) -o $(EXENAME) $^ $(LDFLAGS) $(LDLIBFLAGS)

$(LIBNAME): $(OBJECTS)
	$(CXX) -o $(LIBNAME) $^ $(LDFLAGS) $(LDLIBFLAGS) $(DYNLIBFLAG)

$(LIBSTATNAME): $(OBJECTS)
	ar rcs $(LIBSTATNAME) $^

compile:	$(OBJECTS)

clean:
	rm -f $(EXENAME) $(OBJECTS) $(DEPS) $(LIBNAME)
