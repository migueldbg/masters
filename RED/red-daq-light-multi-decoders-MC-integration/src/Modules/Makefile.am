CLEANFILES = *.d *.so *_rdict.* *_rdict_* *.rootmap
MAINTAINERCLEANFILES = Makefile.in

# EXTRA_DIST = some extra files, e.g. test-data

lib_LTLIBRARIES = libmodules.la

libmodules_la_LIBADD = \
	$(ROOTLIBS)

program_ldadd = libmodules.la $(libmodules_la_LIBADD)

AM_CPPFLAGS = 

libmodules_la_SOURCES = \
    RDCluster.cc \
	RDBaseline.cc \
	RDFind_min.cc \
	RDPulseFit.cc \
	RDPulseFitter.cc \
    RDIntegral.cc \
	RDFiltering.cc \
    RDClusterization.cc \
    MCTraces.cc

libmodules_la_headers = \
    RDCluster.hh \
	RDBaseline.hh \
	RDIntegral.hh \
	RDPulseFit.hh \
	RDPulseFitter.hh \
	peakfinder.h \
	RDFind_min.hh \
	RDFiltering.hh \
	RDClusterization.hh \
    MCTraces.hh

pkginclude_HEADERS = $(libmodules_la_headers)


bin_PROGRAMS = 

noinst_PROGRAMS =

bin_SCRIPTS = 

rootmapdir = $(libdir)
rootmap_DATA =


libmodules_la_SOURCES += libmodules_rdict.cxx
rootmap_DATA += .libs/libmodules.rootmap

if WITH_CLING

rootmap_DATA += .libs/libmodules_rdict.pcm

libmodules_rdict.cxx: $(libmodules_la_headers) modules_LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libmodules@SHLIBEXT@ -rml libmodules@SHLIBEXT@ -rmf libmodules.rootmap.tmp -c $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^/"]*/\)*\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^/"]*/\)*\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/"]*/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]*/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header  *\)\([^ ].*/\)*\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' libmodules.rootmap.tmp > libmodules.rootmap && $(RM) libmodules.rootmap.tmp

.libs/libmodules.rootmap: libmodules.la
	cp libmodules.rootmap libmodules_rdict.pcm .libs/

.libs/libmodules_rdict.pcm: .libs/libmodules.rootmap

else # not WITH_CLING

libmodules_rdict.cxx: $(libmodules_la_headers) modules_LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) -I$(includedir) $+

.libs/libmodules.rootmap: libmodules.la | $(libmodules_la_headers) modules_LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/libmodules@SHLIBEXT@ -c $|

endif # WITH_CLING
