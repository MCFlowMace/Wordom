##---------------------------------------------------------------------------
## Author:      Michele Seeber
## Copyright:   (c) Michele Seeber
## License:     GPL
##---------------------------------------------------------------------------

project(WORDOM)
cmake_minimum_required(VERSION 2.6.2)

 #CUDA for clustering with leader
find_package(CUDA)
if (CUDA_FOUND)

	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_20 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_52,code=compute_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_62,code=sm_62)
	cuda_add_executable(wordom src/main.c src/tsne.cpp src/sptree.cpp src/fileio.c src/tools.c src/datahandler.c src/analysis.c src/geometry.c src/moldiff.c src/cluster.c src/pca.c src/qentropy.c src/kga.c src/ssa.c src/enm.c src/psn.c src/corr.c src/surf.c src/xdrfile.c src/xdrfile_xtc.c src/qcprot.c src/volumes.c src/com.c src/ring.c src/tilt.c src/twist.c src/hole.c src/flux.c src/gCluster.cu)
	add_definitions( -D CUDA )
	set_source_files_properties(src/gCluster.cu PROPERTIES COMPILE_FLAGS --maxrregcount=32)
	
else (CUDA_FOUND)
	add_executable(wordom src/main.c src/tsne.cpp src/sptree.cpp src/fileio.c src/tools.c src/datahandler.c src/analysis.c src/geometry.c src/moldiff.c src/cluster.c src/pca.c src/qentropy.c src/kga.c src/ssa.c src/enm.c src/psn.c src/corr.c src/surf.c src/xdrfile.c src/xdrfile_xtc.c src/qcprot.c src/volumes.c src/com.c src/ring.c src/tilt.c src/twist.c src/hole.c src/flux.c )
endif (CUDA_FOUND)

set(CMAKE_C_FLAGS "-D_GNU_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D GMX -O2 -g -w")
set_source_files_properties(src/xdrfile.c src/xdrfile_xtc.c PROPERTIES COMPILE_FLAGS -w)

# lapack for fast matrix diagonalization
enable_language (Fortran OPTIONAL)
if (CMAKE_Fortran_COMPILER_WORKS)
  find_package(LAPACK)
  find_package(BLAS)
  if (LAPACK_FOUND AND BLAS_FOUND)
    add_definitions( -D LAPACK )
    set(HAVE_LAPACK_LIB 1)
    find_library(LAPACK_LIB lapack)
    target_link_libraries (wordom ${LAPACK_LIB})
    find_library(BLAS_LIB lapack)
    target_link_libraries (wordom ${BLAS_LIB})
  endif (LAPACK_FOUND AND BLAS_FOUND)
endif (CMAKE_Fortran_COMPILER_WORKS)

# threading
find_package (Threads)
if (THREADS_FOUND)
  target_link_libraries (wordom ${CMAKE_THREAD_LIBS_INIT})
  add_definitions( -D THREADED )
endif (THREADS_FOUND)

# platform checks
if(APPLE)
  add_definitions( -D MACOSX )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
endif(APPLE)
if(WIN32)
  add_definitions( -D FWIN32 )
   if(MSVC)
      add_definitions(-DWIN32_LEAN_AND_MEAN)
      add_definitions(-D_CRT_SECURE_NO_WARNINGS)
   endif(MSVC)
   if(CMAKE_C_COMPILER_ID STREQUAL GNU)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-subsystem,windows")
   endif(CMAKE_C_COMPILER_ID STREQUAL GNU)
   if(MINGW AND EXISTS /mingw)
      list(APPEND CMAKE_PREFIX_PATH /mingw)
	endif(MINGW AND EXISTS /mingw)
endif(WIN32)

# math libraries - no if
find_library(M_LIB m)
target_link_libraries (wordom ${M_LIB})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND UNIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local/" CACHE STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif()
#set (CMAKE_INSTALL_PREFIX /usr/local)
INSTALL(FILES build/wordom DESTINATION bin PERMISSIONS WORLD_EXECUTE)

# Add regression tests
enable_testing()
add_subdirectory(tests/)


