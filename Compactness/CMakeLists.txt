cmake_minimum_required(VERSION 3.1...3.23)

project(Compactness)

find_package(CGAL REQUIRED)

find_package(Boost OPTIONAL_COMPONENTS serialization iostreams)
include(CGAL_Boost_serialization_support)
include(CGAL_Boost_iostreams_support)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(NOT TARGET CGAL::Boost_serialization_support)
  message("NOTICE: This project requires Boost Serialization, and will not be compiled.")
endif()

if(NOT TARGET CGAL::Boost_iostreams_support)
  message("NOTICE: This project requires Boost IO Streams, and will not be compiled.")
endif()
set(Boost_BZIP2_LIBRARY_DEBUG "C:/dev/vcpkg/installed/x64-windows/debug/lib/bz2d.lib")
set(Boost_ZLIB_LIBRARY_DEBUG "C:/dev/vcpkg/installed/x64-windows/debug/lib/zlibd.lib")
set(Boost_BZIP2_LIBRARY_RELEASE "C:/dev/vcpkg/installed/x64-windows/lib/bz2.lib")
set(Boost_ZLIB_LIBRARY_RELEASE "C:/dev/vcpkg/installed/x64-windows/lib/zlib.lib")
#C:/dev/vcpkg/installed/x64-windows/debug/lib/bz2d.lib
#C:/dev/vcpkg/installed/x64-windows/debug/lib/zlibd.lib
#C:/dev/vcpkg/installed/x64-windows/lib/bz2.lib
#C:/dev/vcpkg/installed/x64-windows/lib/zlib.lib

create_single_source_cgal_program("compactness_measure.cpp")
create_single_source_cgal_program("compactness_measure_planimetric.cpp")

create_single_source_cgal_program("compactness_classifier.cpp")

create_single_source_cgal_program("compactness_classifier_eval.cpp")
target_link_libraries(compactness_classifier_eval PUBLIC CGAL::Boost_iostreams_support
                                         CGAL::Boost_serialization_support)



create_single_source_cgal_program("compactness_measure_planimetric_copie.cpp")

create_single_source_cgal_program("compactness_measure_planimetric_copie_copie.cpp")


create_single_source_cgal_program("example_ethz_random_forest.cpp")