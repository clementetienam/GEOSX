
#
# Specify all headers
#
set( codingUtilities_headers
     traits.hpp
     EnumStrings.hpp
     StringUtilities.hpp
     UnitTestUtilities.hpp
     Utilities.hpp
     static_if.hpp
     SFINAE_Macros.hpp
   )

#
# Specify all sources
#
set( codingUtilities_sources
     StringUtilities.cpp
   )

convert_filenames_to_full_paths(codingUtilities_sources)
convert_filenames_to_full_paths(codingUtilities_headers)

set(CORE_SOURCES ${CORE_SOURCES} ${codingUtilities_sources} PARENT_SCOPE)
set(CORE_HEADERS ${CORE_HEADERS} ${codingUtilities_headers} PARENT_SCOPE)

geosx_add_code_checks(PREFIX codingUtilities )