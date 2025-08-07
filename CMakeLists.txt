#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.17)
cmake_host_system_information(RESULT host QUERY HOSTNAME)
cmake_host_system_information(RESULT ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
project(HypTPC)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
#find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build HypTPC with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS EG RIO)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${PROJECT_SOURCE_DIR}/.build")
# add_compile_options(-D MemoryLeak)

#----------------------------------------------------------------------------
# Add the executables, and link it to the Geant4 libraries
add_executable(G4HypTPC main.cc ${sources} ${headers})
target_link_libraries(G4HypTPC ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

if(UNIX AND NOT APPLE)
  root_generate_dictionary(G4HypTPC_dict TParticle.h TVector3.h
		LINKDEF ${PROJECT_SOURCE_DIR}/linkdef/LinkDef.h
		MODULE G4HypTPC)
  add_library(G4HypTPC_lib SHARED G4HypTPC_dict)
  target_link_libraries(G4HypTPC G4HypTPC_lib)
endif()

#----------------------------------------------------------------------------
install(TARGETS G4HypTPC RUNTIME DESTINATION bin)
install(FILES
  ${PROJECT_BINARY_DIR}/libG4HypTPC_lib.so
  ${PROJECT_BINARY_DIR}/libG4HypTPC_rdict.pcm
  DESTINATION bin)
