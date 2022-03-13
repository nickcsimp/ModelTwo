cmake_minimum_required(VERSION 3.20)
project(ModelTwo)

set(CMAKE_CXX_STANDARD 14)

add_executable(ModelTwo main.cpp Conglomerate.cpp Conglomerate.h Connection.cpp Connection.h Polymer.cpp Polymer.h FreeSite.cpp FreeSite.h UnconnectedNeighbours.cpp UnconnectedNeighbours.h ConnectedNeighbours.cpp ConnectedNeighbours.h System.cpp System.h Tests.cpp Tests.h settings.cpp settings.h)