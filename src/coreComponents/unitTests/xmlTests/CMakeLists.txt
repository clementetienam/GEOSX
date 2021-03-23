#
# Specify list of tests
#

set(gtest_geosx_tests
   testXML.cpp
   )

#
# Add gtest C++ based tests
#
set( dependencyList gtest optionparser)
set (dependencyList ${dependencyList} geosx_core)

if ( ENABLE_MPI )
  set ( dependencyList ${dependencyList} mpi )
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} "-lz"
                        )

    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/basic_input.xml
                  )
endforeach()
