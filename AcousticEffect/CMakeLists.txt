cmake_minimum_required(VERSION 3.9)

project(AcousticEffect)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(ROOT REQUIRED)
#include_directories(${ROOT_INCLUDE_DIRS})
find_package(OpenMP)


include_directories(
	${CMAKE_SOURCE_DIR}/src
	/home/delirium/root/include
	${CMAKE_CURRENT_SOURCE_DIR}/../Common
)

#set(LIBS ${LIBS}
#	${ROOT_LIBRARIES}
#)

link_directories(
	/home/delirium/root/lib
)

include("CMakeLists.Files.txt")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

source_group("Headers" FILES ${HEADERS})
set(SOURCES ${SOURCES} ${HEADERS})

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${LIBS} OpenMP::OpenMP_CXX)
#if(OpenMP_CXX_FOUND)
#	  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
#endif()
