cmake_minimum_required(VERSION 3.6)
project(GraphicsEngine)

set(CMAKE_CXX_STANDARD 11)

add_executable(GraphicsEngine
        easy_image.cc
        easy_image.h
        engine.cc
        ini_configuration.cc
        ini_configuration.h
        l_parser.cc
        l_parser.h Intro.cpp EngineUtilities.cpp session1.cpp session2.cpp session4.cpp session5.cpp session3.cpp session7.cpp)
