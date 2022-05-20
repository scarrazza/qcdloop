cmake_minimum_required (VERSION 3.0.2)

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

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(VERSION "\"2.0.9\"")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=nocona -mtune=haswell -fvisibility-inlines-hidden -fmessage-length=0 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -fext-numeric-literals")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fPIC ${CMAKE_CXX_FLAGS}" CACHE STRING "debug compile flags" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=address" CACHE STRING "debug linker flags" FORCE)
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")

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

set(QUADMATH_NAMES ${QUADMATH_NAMES} libquadmath.so quadmath)
find_library(QUADMATH_LIBRARY
  NAMES ${QUADMATH_NAMES}
  PATHS /usr/lib64/atlas /usr/lib/atlas
        /usr/lib64 /usr/lib /usr/local/lib64
        /usr/local/lib /usr/x86_64-linux-gnu/*
        /usr/lib/gcc/x86_64-linux-gnu/*
        /usr/lib/gcc/x86_64-redhat-linux/*
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

target_link_libraries(qcdloop)
install(FILES ${CMAKE_BINARY_DIR}/src/qcdloop-config
  DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(FILES ${PROJECT_SOURCE_DIR}/src/qcdloop.pc DESTINATION lib/pkgconfig)
install(DIRECTORY src/qcdloop DESTINATION include)
install(TARGETS qcdloop DESTINATION lib)

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
