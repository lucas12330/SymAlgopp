cmake_minimum_required(VERSION 3.16)
project(consommateur_symalgopp LANGUAGES CXX)

find_package(symalgopp REQUIRED)

add_executable(consommateur main.cpp)
target_link_libraries(consommateur PRIVATE symalgopp::symalgopp)
