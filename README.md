
# OpenSubdiv Lite

## Overview

A light-weight fork of the OpenSubdiv. Contains all the topology processing logic
(Sdc/Vtr/Far layers), but none of the compute/rasterization back-ends (Osd layer).

Also includes 'Tmr' : the Topology Map Representation layer with support for topology
hashing.
(see http://www.graphics.stanford.edu/~niessner/papers/2016/4subdiv/brainerd2016efficient.pdf)

## Requirements

* Windows or Linux
* CMake 3.22
* C++ 20

## Build

1. Git clone this source tree

2. Ensure you have a recent C++ std20 compiler such as the latest Visual Studio 2022
   (run the Visual Studio updater application if unsure), Clang or GCC

3. Create a 'build' folder 
   (at the top of the repository, names such as 'build' or '_build' are recommended)

3. Use CMake to configure the build.
   * CMake GUI : make sure you set the build folder ! ("Where to build the binaries" line)

   * You can also use the command line or scripts:
     ```bash
     cmake_cmd="C:/Program\ Files/CMake/bin/cmake.exe"
     cmd="$cmake_cmd
          -G 'Visual Studio 17 2022 -A x64'
          .."
     echo $cmd
     eval $cmd
     ```

5. Build : either open the solution in MSVC or use the CMake command line builder
   ```bash
   cmake --build . --config Release --target package
   ```
**_Note:_** it is possible to generate a packaged build with both debug and release
binaries, but the process requires multiple steps. See the 'build_and_package.sh`
script for details.

**_Note:_** 32 bits builds are not supported

## Integration

When integrating OSD_Lite into other projects with cmake, the following dependencies
are provided for either static or dynamic linking:
- osd::osd_lite_static
- osd::osd_lite_shared

**_Note:_** OSD_Lite automatically exports cmake configuration files, so it should
never be necessary to write custom CMake `Find` modules.

**_Note:_** the `osd::` namespace is only pre-fixed when using the "exported" build
of the library through `find_package`. If OSD_Lite is imported via git submodule or
CMake's `FetchContent` mechanism, the `osd::` prefix should not be used to name 
OSD_Lite dependencies (this appears to be a somewhat inconsistent policy of CMake's
EXPORT feature that may change in the future).
