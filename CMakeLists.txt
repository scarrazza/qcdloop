cmake_minimum_required(VERSION 3.12...3.31)

# Disable the use of RPATHS - we probably are not
# that interested in relocatable binaries and it
# causes a lot of trouble on macOS
SET(CMAKE_SKIP_RPATH TRUE)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

project(qcdloop)
include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(VERSION 2.0.11)

# Set default CXXFLAGS but allow for environment override
# c.f. https://cmake.org/cmake/help/v3.31/envvar/CXXFLAGS.html
if (NOT CMAKE_CXX_FLAGS)
  # c.f. https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
  if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -fvisibility-inlines-hidden -fmessage-length=0 -ftree-vectorize -fstack-protector-strong -O2 -pipe -fext-numeric-literals")
  else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra")
  endif()
else()
  if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    # Ensure -fext-numeric-literals is in CMAKE_CXX_FLAGS
    string(FIND "${CMAKE_CXX_FLAGS}" "-fext-numeric-literals" _found_ext_numeric_literals)
    if (_found_ext_numeric_literals EQUAL -1)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals" CACHE STRING "Update environment CXXFLAGS" FORCE)
    endif()
  endif()
endif()

# TODO: QCDLoop is currently only well tested on GNU, so other compiler
# and linker defaults should be added as they become known
if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address ${CMAKE_CXX_FLAGS}" CACHE STRING "debug CXXFLAGS" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address" CACHE STRING "debug linker flags" FORCE)
endif()

set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(includedir ${CMAKE_INSTALL_INCLUDEDIR})
set(libdir ${CMAKE_INSTALL_LIBDIR})

configure_file(
  "${PROJECT_SOURCE_DIR}/src/qcdloop/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/qcdloop/config.h"
  )

configure_file(
  "${PROJECT_SOURCE_DIR}/src/qcdloop-config.in"
  "${PROJECT_BINARY_DIR}/src/qcdloop-config"
  )

configure_file(
  "${PROJECT_SOURCE_DIR}/src/qcdloop.pc.in"
  "${PROJECT_SOURCE_DIR}/src/qcdloop.pc"
  )
IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
set(QUADMATH_NAMES ${QUADMATH_NAMES} libquadmath.so quadmath)
find_library(QUADMATH_LIBRARY
  NAMES ${QUADMATH_NAMES}
  PATHS /usr/lib64/atlas /usr/lib/atlas
        /usr/lib64 /usr/lib /usr/local/lib64
        /usr/local/lib /usr/x86_64-linux-gnu/*
        /usr/lib/gcc/x86_64-linux-gnu/*
        /usr/lib/gcc/x86_64-redhat-linux/*
        /usr/lib/gcc/x86_64-pc-linux-gnu/*
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
)

if(QUADMATH_LIBRARY)
  set(QUADMATH_LIBRARIES ${QUADMATH_LIBRARY})
  set(QUADMATH_FOUND "YES")
else(QUADMATH_LIBRARY)
  set(QUADMATH_FOUND "NO")
endif(QUADMATH_LIBRARY)

if(QUADMATH_FOUND)
  if(NOT QUADMATH_FIND_QUIETLY)
    message(STATUS "Found QUADMATH: ${QUADMATH_LIBRARIES}")
  endif(NOT QUADMATH_FIND_QUIETLY)
else(QUADMATH_FOUND)
  if(QUADMATH_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find QuadMath")
  endif(QUADMATH_FIND_REQUIRED)
endif(QUADMATH_FOUND)
endif()
# libqcdloop configuration
include_directories(src/qcdloop src)
FILE(GLOB_RECURSE Headers "src/qcdloop/*.h")
add_library(qcdloop SHARED src/box.cc
                         src/bubble.cc
                         src/cache.cc
                         src/qcdloop.cc
                         src/tadpole.cc
                         src/tools.cc
                         src/topology.cc
                         src/triangle.cc
                         src/types.cc
                         src/wrapper.cc
                         ${Headers}
)

target_link_libraries(qcdloop ${QUADMATH_LIBRARY})
install(FILES ${CMAKE_BINARY_DIR}/src/qcdloop-config
  DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(FILES ${PROJECT_SOURCE_DIR}/src/qcdloop.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
install(DIRECTORY src/qcdloop DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS qcdloop DESTINATION ${CMAKE_INSTALL_LIBDIR})

# enable disable fortran/wrapper (to avoid name conflicts with ql1.x)
option(ENABLE_FORTRAN_WRAPPER "Enable fortran wrapper" ON)
if(ENABLE_FORTRAN_WRAPPER)
  add_definitions(-DQL_NAMES)
endif()

# qcdloop examples
option(ENABLE_EXAMPLES "Compile test examples" OFF)
if(ENABLE_EXAMPLES)
  add_executable(cache_test examples/cache_test.cc)
  target_link_libraries(cache_test qcdloop ${QUADMATH_LIBRARIES})
  set_target_properties(cache_test PROPERTIES LINK_FLAGS "-Wl,-rpath,./")

  add_executable(cmass_test examples/cmass_test.cc)
  target_link_libraries(cmass_test qcdloop ${QUADMATH_LIBRARIES})
  set_target_properties(cmass_test PROPERTIES LINK_FLAGS "-Wl,-rpath,./")

  add_executable(trigger_test examples/trigger_test.cc)
  target_link_libraries(trigger_test qcdloop ${QUADMATH_LIBRARIES})
  set_target_properties(trigger_test PROPERTIES LINK_FLAGS "-Wl,-rpath,./")
endif()
