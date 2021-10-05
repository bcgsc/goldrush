bin_PROGRAMS = workpackage2

workpackage2_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

workpackage2_CPPFLAGS = -I$(top_srcdir)/workpackage2 \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)
	
workpackage2_LDADD = $(top_builddir)/Common/libcommon.a -lz \
	-lsdsl
	
workpackage2_LDFLAGS = $(OPENMP_CXXFLAGS)

workpackage2_SOURCES = workpackage2.cpp \
	MIBFConstructSupport.hpp \
	MIBFQuerySupport.hpp \
	MIBloomFilter.hpp \
	MIBFClassifier.hpp \
	Options.h Options.cpp
