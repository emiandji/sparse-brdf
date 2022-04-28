


#FILE (GLOB_RECURSE project_SRCS *.cpp *.cxx *.cc *.C *.c *.h)


#
# Project Properties
#
CMAKE_MINIMUM_REQUIRED (VERSION 2.6.4)


PROJECT (ccslib)
SET (APPLICATION_NAME "ccslib")
SET (APPLICATION_CODENAME "${PROJECT_NAME}")
SET (APPLICATION_COPYRIGHT_YEARS "2021 2022")
SET (APPLICATION_VERSION_STRING "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}.${APPLICATION_VERSION_PATCH}-${APPLICATION_VERSION_TYPE}")
SET (APPLICATION_VENDOR_ID "MIT/LIU")
SET (APPLICATION_VENDOR_NAME "Computer Graphics and Image Processing Laboratory - Linkopings University")
SET (APPLICATION_VENDOR_URL "vcl.itn.liu.se")
SET (APPLICATION_ID "${APPLICATION_VENDOR_ID}.${PROJECT_NAME}")
 
#
# BUILD TYPE
#
set(CMAKE_BUILD_TYPE Release)
         
#
# Debugging Options
#
SET (CMAKE_VERBOSE_MAKEFILE 0) # Use 1 for debugging, 0 for release



SET(MAINFOLDER ${CMAKE_CURRENT_SOURCE_DIR})
FILE(GLOB libData_SRCS ${MAINFOLDER}/cssData/*.h)
FILE(GLOB libDict_SRCS ${MAINFOLDER}/cssDict/*.cpp)
FILE(GLOB libRecon_SRCS ${MAINFOLDER}/cssRecon/*.h)
FILE(GLOB libUtil_SRCS ${MAINFOLDER}/cssUtil/*.cpp ${MAINFOLDER}/cssUtil/*.c)

SET(CCSLIB_SRCS ${libData_SRCS} ${libDict_SRCS} ${libRecon_SRCS})
SET(proj_EXECS ${MAINFOLDER}/cssMain)




SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O2 -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS -DHAVE_MAT_INT64_T -DHAVE_MAT_UINT64_T -std=c++11 -fopenmp") #-DEIGEN_DONT_PARALLELIZE

SET(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${CMAKE_CXX_FLAGS})


find_package( HDF5 REQUIRED )
FIND_PACKAGE( Boost 1.53 COMPONENTS system serialization filesystem random REQUIRED )

SET(EIGEN3 "/usr/include/eigen3")


INCLUDE_DIRECTORIES(${MAINFOLDER})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${EIGEN3})
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})


link_directories(${Boost_LIBRARY_DIR})


SET(project_LIBS 
${Boost_LIBRARIES}
${HDF5_LIBS}
z
matio
pthread 
m 
dl
${LIBS})




SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin/)

ADD_LIBRARY(ccsutils STATIC ${libUtil_SRCS})
#set_target_properties(ccsutils PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(ccsutils ${project_LIBS})

ADD_LIBRARY(ccsdict STATIC ${libDict_SRCS})
set_target_properties(ccsdict PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(ccsdict ccsutils ${project_LIBS})



# BRDF ####################################################################
ADD_EXECUTABLE(brdf3Dcomp ${proj_EXECS}/brdf3Dcomp.cpp)
TARGET_LINK_LIBRARIES(brdf3Dcomp ccsutils ccsdict ${project_LIBS})
ADD_EXECUTABLE(brdf3Drecon ${proj_EXECS}/brdf3Drecon.cpp)
TARGET_LINK_LIBRARIES(brdf3Drecon ccsutils ccsdict ${project_LIBS})
ADD_EXECUTABLE(brdf3Dinterp ${proj_EXECS}/brdf3Dinterp.cpp)
TARGET_LINK_LIBRARIES(brdf3Dinterp ccsutils ccsdict ${project_LIBS})


# for debugging ###########################################################
ADD_EXECUTABLE(temptest ${proj_EXECS}/temptest.cpp )
TARGET_LINK_LIBRARIES(temptest ccsutils ccsdict ${project_LIBS})

