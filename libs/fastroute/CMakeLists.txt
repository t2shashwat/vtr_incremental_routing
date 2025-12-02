# libs/fastroute/CMakeLists.txt

add_library(fastroute STATIC
    src/dist.cpp
    src/mst2.cpp
    src/neighbors.cpp
    src/heap.cpp
    src/maze.cpp
    src/RipUp.cpp
    src/route.cpp
    src/maze3D.cpp
    src/DataProc.cpp
    src/FastRoute.cpp
    src/EdgeShift.cpp
    src/dl.cpp
    src/RSMT.cpp
    src/utility.cpp
    src/rand-pts.cpp
    ${CMAKE_SOURCE_DIR}/vpr/src/util/flute.cpp

)

target_include_directories(fastroute PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/vpr/src/util
)
