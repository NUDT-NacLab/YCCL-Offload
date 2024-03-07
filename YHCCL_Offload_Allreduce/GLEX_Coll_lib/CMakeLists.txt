cmake_minimum_required(VERSION 2.8.12)
project(GLEX_coll)

#set(MPI_DIR "/mnt/share/pjt/mpi3/")
#set(MPI_INCLUDE_DIRS "/mnt/share/pjt/mpi3/include/")
#set(INTERFACE_INCLUDE_DIRECTORIES "/mnt/share/pjt/mpi3/include/")



# set(MPI_DIR "/usr/local/ompi-x/")
# set(MPI_INCLUDE_DIRS "/usr/local/ompi-x/include/")
# set(INTERFACE_INCLUDE_DIRECTORIES "/usr/local/ompi-x/include/")

# find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

if(DEFINED ENV{PJT_GDB_DEBUG})
    set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -ggdb")
    set(CMAKE_C_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -ggdb")
endif()
set(CMAKE_C_FLAGS_DEBUG "$ENV{CFLAGS} -std=c99")

link_directories("/usr/local/glex/include")
link_directories("/vol8/home/nudtgcy/PJT-WorkDir/numactl-master/build/lib")
link_directories("/vol8/home/nudtgcy/PJT-WorkDir/numactl-master/build/include")
include_directories("${PROJECT_SOURCE_DIR}/IntCompresser")
link_directories("${PROJECT_SOURCE_DIR}/IntCompresser")

ADD_SUBDIRECTORY(./src)
ADD_SUBDIRECTORY(./test)
