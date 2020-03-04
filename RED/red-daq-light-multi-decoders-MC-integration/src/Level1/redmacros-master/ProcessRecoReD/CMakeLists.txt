cmake_minimum_required (VERSION 2.6)

project (ReD)

set(CMAKE_MACOSX_RPATH 1)

#Set working path and install directory
set( WorkPath "$ENV{PWD}/.." )
set( CMAKE_INSTALL_PREFIX ${WorkPath} )


# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(src)        #src needed by packages below

#---Create  a main program using the library

add_executable(ProcessRecoReD ProcessRecoReD.C)
target_link_libraries(ProcessRecoReD ReDANA ${ROOT_LIBRARIES})
#install (TARGETS ProcessRecoReD DESTINATION bin)
message("-- install place  =  ${WorkPath}")

get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- inc_dirs (Project.Ref)  =  ${inc_dirs}")

get_property(link_dirs DIRECTORY PROPERTY LINK_DIRECTORIES)
message("-- link_dirs (Project.Ref)  =  ${link_dirs}")

#  add_executable( PlotHists PlotHists.C )
#  target_link_libraries ( PlotHists ${ROOT_LIBRARIES} )
#  install (TARGETS PlotHists DESTINATION ${WorkPath} )


#LIST(APPEND CMAKE_PROGRAM_PATH  ${WorkPath}/bin/ )
