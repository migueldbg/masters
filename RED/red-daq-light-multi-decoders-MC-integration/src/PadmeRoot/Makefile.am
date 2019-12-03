CLEANFILES = *.d *.so *_rdict.* *_rdict_* *.rootmap
MAINTAINERCLEANFILES = Makefile.in

# EXTRA_DIST = some extra files, e.g. test-data

lib_LTLIBRARIES = libpadmeroot.la

libpadmeroot_la_LIBADD =

program_ldadd = libpadmeroot.la

libpadmeroot_la_SOURCES = \
        TADCBoard.cc \
        TADCChannel.cc \
        TADCTrigger.cc \
        TDetectorInfo.cc \
        TECalMCEvent.cc \
        TECalMCHit.cc \
        TEVetoMCEvent.cc \
        TEVetoMCHit.cc \
        THEPVetoMCEvent.cc \
        THEPVetoMCHit.cc \
        TMCEvent.cc \
        TMCVEvent.cc \
        TMCVHit.cc \
        TPadmeRun.cc \
        TPVetoMCEvent.cc \
        TPVetoMCHit.cc \
        TRawEvent.cc \
        TSACMCEvent.cc \
        TSACMCHit.cc \
        TSubDetectorInfo.cc \
        TTargetMCEvent.cc \
        TTargetMCHit.cc 

libpadmeroot_la_headers = \
        TADCBoard.hh \
        TADCChannel.hh \
        TADCTrigger.hh \
        TDetectorInfo.hh \
        TECalMCEvent.hh \
        TECalMCHit.hh \
        TEVetoMCEvent.hh \
        TEVetoMCHit.hh \
        THEPVetoMCEvent.hh \
        THEPVetoMCHit.hh \
        TMCEvent.hh \
        TMCVEvent.hh \
        TMCVHit.hh \
        TPadmeRun.hh \
        TPVetoMCEvent.hh \
        TPVetoMCHit.hh \
        TRawEvent.hh \
        TSACMCEvent.hh \
        TSACMCHit.hh \
        TSubDetectorInfo.hh \
        TTargetMCEvent.hh \
        TTargetMCHit.hh 

pkginclude_HEADERS = $(libpadmeroot_la_headers)


bin_PROGRAMS =

noinst_PROGRAMS =

rootmapdir = $(libdir)
rootmap_DATA =


libpadmeroot_la_SOURCES += libpadmeroot_rdict.cxx
rootmap_DATA += .libs/libpadmeroot.rootmap

if WITH_CLING

rootmap_DATA += .libs/libpadmeroot_rdict.pcm

libpadmeroot_rdict.cxx: $(libpadmeroot_la_headers) padmeroot_LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libpadmeroot@SHLIBEXT@ -rml libpadmeroot@SHLIBEXT@ -rmf libpadmeroot.rootmap.tmp -c $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^/"]*/\)*\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^/"]*/\)*\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/"]*/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]*/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header  *\)\([^ ].*/\)*\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' libpadmeroot.rootmap.tmp > libpadmeroot.rootmap && $(RM) libpadmeroot.rootmap.tmp

.libs/libpadmeroot.rootmap: libpadmeroot.la
	cp libpadmeroot.rootmap libpadmeroot_rdict.pcm .libs/

.libs/libpadmeroot_rdict.pcm: .libs/libpadmeroot.rootmap

else # not WITH_CLING

libpadmeroot_rdict.cxx: $(libpadmeroot_la_headers) padmeroot_LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $+

.libs/libpadmeroot.rootmap: libpadmeroot.la | $(libpadmeroot_la_headers) padmeroot_LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/libpadmeroot@SHLIBEXT@ -c $|

endif # WITH_CLING
