CLEANFILES = *.d *.so *_rdict.* *_rdict_* *.rootmap
MAINTAINERCLEANFILES = Makefile.in

# EXTRA_DIST = some extra files, e.g. test-data

lib_LTLIBRARIES = liblevelone.la

liblevelone_la_LIBADD = \
	$(top_builddir)/src/PadmeRoot/libpadmeroot.la \
	$(top_builddir)/src/Modules/libmodules.la $(ROOTLIBS)

program_ldadd = liblevelone.la $(liblevelone_la_LIBADD)

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/PadmeRoot \
	-I$(top_srcdir)/src/Modules

liblevelone_la_SOURCES = \
        ADCBoard.cc  \
        ADCEvent.cc \
        ADCFile.cc \
        DBService.cc \
	RDconfig.cc \
        RootIO.cc \
	Decoder.cc \
	MCDecoder.cc \
        Analyzer.cc \
	EvHeader.cc \
	EvRaw0.cc \
	EvRec0.cc \
        VDecoder.cc \
	RunHeader.cc

liblevelone_la_headers = \
	ADCBoard.hh \
	ADCEvent.hh \
	ADCEventTags.hh \
	ADCFile.hh \
        Analyzer.hh \
	DBService.hh \
	RDconfig.h \
	RootIO.hh  \
        VDecoder.hh \
	MCDecoder.hh \
	Decoder.hh \
	EvHeader.hh \
	EvRaw0.hh \
	EvRec0.hh\
	RunHeader.hh
     
pkginclude_HEADERS = $(liblevelone_la_headers)


bin_PROGRAMS = RedLevel1 viewer viewer_online
#liveViewer viewer

noinst_PROGRAMS =

RedLevel1_SOURCES = RedLevel1.cpp
RedLevel1_LDADD = $(program_ldadd)

viewer_SOURCES = viewer.cpp
viewer_LDADD = $(program_ldadd)

viewer_online_SOURCES = viewer_online.cpp
viewer_online_LDADD = $(program_ldadd)

#liveViewer_SOURCES = liveViewer.cpp
#liveViewer_LDADD = $(program_ldadd)

bin_SCRIPTS = script/makelist.sh script/reco.sh script/makelist_multi_lns.sh script/get_last_daq_files.sh

rootmapdir = $(libdir)
rootmap_DATA =


liblevelone_la_SOURCES += liblevelone_rdict.cxx
rootmap_DATA += .libs/liblevelone.rootmap

if WITH_CLING

rootmap_DATA += .libs/liblevelone_rdict.pcm

liblevelone_rdict.cxx: $(liblevelone_la_headers) levelone_LinkDef.h
	$(ROOTCLING) -f $@.tmp -s liblevelone@SHLIBEXT@ -rml liblevelone@SHLIBEXT@ -rmf liblevelone.rootmap.tmp -c $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^/"]*/\)*\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^/"]*/\)*\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/"]*/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]*/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header  *\)\([^ ].*/\)*\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' liblevelone.rootmap.tmp > liblevelone.rootmap && $(RM) liblevelone.rootmap.tmp

.libs/liblevelone.rootmap: liblevelone.la
	cp liblevelone.rootmap liblevelone_rdict.pcm .libs/

.libs/liblevelone_rdict.pcm: .libs/liblevelone.rootmap

else # not WITH_CLING

liblevelone_rdict.cxx: $(liblevelone_la_headers) levelone_LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) -I$(includedir) $+

.libs/liblevelone.rootmap: liblevelone.la | $(liblevelone_la_headers) levelone_LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/liblevelone@SHLIBEXT@ -c $|

endif # WITH_CLING
