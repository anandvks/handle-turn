# Copyright: (C) 2016 iCub Facility Istituto Italiano di Tecnologia
# Authors: Ugo Pattacini
# CopyPolicy: Released under the terms of the GNU GPL v2.0.

cmake_minimum_required(VERSION 2.8.9)
project(handle-turn)

find_package(YARP REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${YARP_MODULE_PATH})

find_package(ICUB REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${ICUB_MODULE_PATH})

find_package(IPOPT REQUIRED)

include(AddInstallRPATHSupport)
add_install_rpath_support(LIB_DIRS "${CMAKE_INSTALL_FULL_LIBDIR}"
                          BIN_DIRS "${CMAKE_INSTALL_FULL_BINDIR}"
                                   "${CMAKE_INSTALL_PREFIX}/${CER_DYNAMIC_PLUGINS_INSTALL_DIR}"
                          USE_LINK_PATH)

add_subdirectory(cer_kinematics)
add_subdirectory(modules)

