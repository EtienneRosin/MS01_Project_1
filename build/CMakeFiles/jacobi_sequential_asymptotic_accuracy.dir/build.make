# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/etienne/Documents/Developer/MS01_Project_1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/etienne/Documents/Developer/MS01_Project_1/build

# Include any dependencies generated for this target.
include CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/flags.make

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/flags.make
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o: /Users/etienne/Documents/Developer/MS01_Project_1/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/etienne/Documents/Developer/MS01_Project_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o -MF CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o.d -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o -c /Users/etienne/Documents/Developer/MS01_Project_1/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/etienne/Documents/Developer/MS01_Project_1/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp > CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.i

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/etienne/Documents/Developer/MS01_Project_1/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.s

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/flags.make
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o: /Users/etienne/Documents/Developer/MS01_Project_1/utils/Vector.cpp
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/etienne/Documents/Developer/MS01_Project_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o -MF CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o.d -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o -c /Users/etienne/Documents/Developer/MS01_Project_1/utils/Vector.cpp

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/etienne/Documents/Developer/MS01_Project_1/utils/Vector.cpp > CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.i

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/etienne/Documents/Developer/MS01_Project_1/utils/Vector.cpp -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.s

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/flags.make
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o: /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_MPI.cpp
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/etienne/Documents/Developer/MS01_Project_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o -MF CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o.d -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o -c /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_MPI.cpp

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_MPI.cpp > CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.i

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_MPI.cpp -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.s

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/flags.make
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o: /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_seq.cpp
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/etienne/Documents/Developer/MS01_Project_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o -MF CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o.d -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o -c /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_seq.cpp

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_seq.cpp > CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.i

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_seq.cpp -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.s

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/flags.make
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o: /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_sequential.cpp
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/etienne/Documents/Developer/MS01_Project_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o -MF CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o.d -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o -c /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_sequential.cpp

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_sequential.cpp > CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.i

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/etienne/Documents/Developer/MS01_Project_1/solvers/jacobi_sequential.cpp -o CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.s

# Object files for target jacobi_sequential_asymptotic_accuracy
jacobi_sequential_asymptotic_accuracy_OBJECTS = \
"CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o" \
"CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o" \
"CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o" \
"CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o" \
"CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o"

# External object files for target jacobi_sequential_asymptotic_accuracy
jacobi_sequential_asymptotic_accuracy_EXTERNAL_OBJECTS =

jacobi_sequential_asymptotic_accuracy: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp.o
jacobi_sequential_asymptotic_accuracy: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/utils/Vector.cpp.o
jacobi_sequential_asymptotic_accuracy: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_MPI.cpp.o
jacobi_sequential_asymptotic_accuracy: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_seq.cpp.o
jacobi_sequential_asymptotic_accuracy: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/solvers/jacobi_sequential.cpp.o
jacobi_sequential_asymptotic_accuracy: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/build.make
jacobi_sequential_asymptotic_accuracy: /opt/homebrew/Cellar/open-mpi/5.0.3_1/lib/libmpi.dylib
jacobi_sequential_asymptotic_accuracy: CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/etienne/Documents/Developer/MS01_Project_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable jacobi_sequential_asymptotic_accuracy"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/build: jacobi_sequential_asymptotic_accuracy
.PHONY : CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/build

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/cmake_clean.cmake
.PHONY : CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/clean

CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/depend:
	cd /Users/etienne/Documents/Developer/MS01_Project_1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/etienne/Documents/Developer/MS01_Project_1 /Users/etienne/Documents/Developer/MS01_Project_1 /Users/etienne/Documents/Developer/MS01_Project_1/build /Users/etienne/Documents/Developer/MS01_Project_1/build /Users/etienne/Documents/Developer/MS01_Project_1/build/CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/jacobi_sequential_asymptotic_accuracy.dir/depend

