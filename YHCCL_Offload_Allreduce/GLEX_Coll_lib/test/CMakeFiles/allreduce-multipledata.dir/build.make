# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_COMMAND = /mnt/share/hyy/cmake/cmake-3.23.1-linux-x86_64/bin/cmake

# The command to remove a file.
RM = /mnt/share/hyy/cmake/cmake-3.23.1-linux-x86_64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/share/pjt/GLEX_Coll_lib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/share/pjt/GLEX_Coll_lib

# Include any dependencies generated for this target.
include test/CMakeFiles/allreduce-multipledata.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/allreduce-multipledata.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/allreduce-multipledata.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/allreduce-multipledata.dir/flags.make

test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o: test/CMakeFiles/allreduce-multipledata.dir/flags.make
test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o: test/allreduce-multipledata.cpp
test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o: test/CMakeFiles/allreduce-multipledata.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/share/pjt/GLEX_Coll_lib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o"
	cd /mnt/share/pjt/GLEX_Coll_lib/test && /usr/local/ompi-x/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o -MF CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o.d -o CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o -c /mnt/share/pjt/GLEX_Coll_lib/test/allreduce-multipledata.cpp

test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.i"
	cd /mnt/share/pjt/GLEX_Coll_lib/test && /usr/local/ompi-x/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/share/pjt/GLEX_Coll_lib/test/allreduce-multipledata.cpp > CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.i

test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.s"
	cd /mnt/share/pjt/GLEX_Coll_lib/test && /usr/local/ompi-x/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/share/pjt/GLEX_Coll_lib/test/allreduce-multipledata.cpp -o CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.s

# Object files for target allreduce-multipledata
allreduce__multipledata_OBJECTS = \
"CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o"

# External object files for target allreduce-multipledata
allreduce__multipledata_EXTERNAL_OBJECTS =

build/test/allreduce-multipledata: test/CMakeFiles/allreduce-multipledata.dir/allreduce-multipledata.cpp.o
build/test/allreduce-multipledata: test/CMakeFiles/allreduce-multipledata.dir/build.make
build/test/allreduce-multipledata: test/CMakeFiles/allreduce-multipledata.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/share/pjt/GLEX_Coll_lib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../build/test/allreduce-multipledata"
	cd /mnt/share/pjt/GLEX_Coll_lib/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/allreduce-multipledata.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/allreduce-multipledata.dir/build: build/test/allreduce-multipledata
.PHONY : test/CMakeFiles/allreduce-multipledata.dir/build

test/CMakeFiles/allreduce-multipledata.dir/clean:
	cd /mnt/share/pjt/GLEX_Coll_lib/test && $(CMAKE_COMMAND) -P CMakeFiles/allreduce-multipledata.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/allreduce-multipledata.dir/clean

test/CMakeFiles/allreduce-multipledata.dir/depend:
	cd /mnt/share/pjt/GLEX_Coll_lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/share/pjt/GLEX_Coll_lib /mnt/share/pjt/GLEX_Coll_lib/test /mnt/share/pjt/GLEX_Coll_lib /mnt/share/pjt/GLEX_Coll_lib/test /mnt/share/pjt/GLEX_Coll_lib/test/CMakeFiles/allreduce-multipledata.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/allreduce-multipledata.dir/depend

