# Create lib from source files

#Set working path and install directory
set( WorkPath "$ENV{PWD}/.." )
set( CMAKE_INSTALL_PREFIX ${WorkPath} )

set( CMAKE_CXX_FLAGS -std=c++17 )

include_directories(${WorkPath}/src)

set(HEADER_FILES EvHeader.hh RDCluster.hh RDPulseFit.hh RedAna.h)

#---Generate dictionary
ROOT_GENERATE_DICTIONARY(G_RedDict ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(ReDANA SHARED  EvHeader.cc RDCluster.cc RDPulseFit.cc RedAna.C G_RedDict.cxx)
target_link_libraries(ReDANA ${ROOT_LIBRARIES})

# Install the library in lib/  and the include files in include/ 
install (TARGETS ReDANA DESTINATION lib )
install (FILES ${HEADER_FILES} DESTINATION include )
