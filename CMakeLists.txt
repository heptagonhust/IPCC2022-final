cmake_minimum_required(VERSION 3.1)
project(sparsify)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
add_compile_options(-Wall -Wno-unused-result -Wno-sign-compare -O3 -g -mavx2)
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)


add_executable(sparsify main.cpp)
add_executable(ng-sparsify new_main.cpp timer.hpp)

find_package(TBB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

# find_package(MPI REQUIRED)
# target_link_libraries(pivot MPI::MPI_CXX)
target_link_libraries(ng-sparsify Threads::Threads)
target_link_libraries(ng-sparsify OpenMP::OpenMP_CXX)
target_link_libraries(ng-sparsify TBB::tbb)
