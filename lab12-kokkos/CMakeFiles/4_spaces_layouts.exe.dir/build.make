# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/cmake-3.21.4-mnpy2irg5xtqp547zxkyzm2mltu6kq6a/bin/cmake

# The command to remove a file.
RM = /hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/cmake-3.21.4-mnpy2irg5xtqp547zxkyzm2mltu6kq6a/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos

# Include any dependencies generated for this target.
include CMakeFiles/4_spaces_layouts.exe.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/4_spaces_layouts.exe.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/4_spaces_layouts.exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/4_spaces_layouts.exe.dir/flags.make

CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o: CMakeFiles/4_spaces_layouts.exe.dir/flags.make
CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o: 4_spaces_layouts.cpp
CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o: CMakeFiles/4_spaces_layouts.exe.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o"
	/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/kokkos-3.6.01-qu45u5vcjttvscfh76nijf2keukcj5z7/bin/kokkos_launch_compiler /hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/kokkos-3.6.01-qu45u5vcjttvscfh76nijf2keukcj5z7/bin/nvcc_wrapper g++ g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o -MF CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o.d -o CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o -c /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos/4_spaces_layouts.cpp

CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos/4_spaces_layouts.cpp > CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.i

CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos/4_spaces_layouts.cpp -o CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.s

# Object files for target 4_spaces_layouts.exe
4_spaces_layouts_exe_OBJECTS = \
"CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o"

# External object files for target 4_spaces_layouts.exe
4_spaces_layouts_exe_EXTERNAL_OBJECTS =

4_spaces_layouts.exe: CMakeFiles/4_spaces_layouts.exe.dir/4_spaces_layouts.cpp.o
4_spaces_layouts.exe: CMakeFiles/4_spaces_layouts.exe.dir/build.make
4_spaces_layouts.exe: /hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/kokkos-3.6.01-qu45u5vcjttvscfh76nijf2keukcj5z7/lib64/libkokkoscontainers.so.3.6.01
4_spaces_layouts.exe: /hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/kokkos-3.6.01-qu45u5vcjttvscfh76nijf2keukcj5z7/lib64/libkokkoscore.so.3.6.01
4_spaces_layouts.exe: /hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/cuda-11.1.0-t3z5kvfalh6n7rizd7scrtgskvaz7dxa/lib64/stubs/libcuda.so
4_spaces_layouts.exe: /hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/cuda-11.1.0-t3z5kvfalh6n7rizd7scrtgskvaz7dxa/lib64/libcudart.so
4_spaces_layouts.exe: /usr/lib64/libdl.so
4_spaces_layouts.exe: CMakeFiles/4_spaces_layouts.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 4_spaces_layouts.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/4_spaces_layouts.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/4_spaces_layouts.exe.dir/build: 4_spaces_layouts.exe
.PHONY : CMakeFiles/4_spaces_layouts.exe.dir/build

CMakeFiles/4_spaces_layouts.exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/4_spaces_layouts.exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/4_spaces_layouts.exe.dir/clean

CMakeFiles/4_spaces_layouts.exe.dir/depend:
	cd /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos /users/xiny/GitHubRepository/Math6370-Yang/lab12-kokkos/CMakeFiles/4_spaces_layouts.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/4_spaces_layouts.exe.dir/depend

