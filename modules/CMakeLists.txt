# Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
# Author: Ugo Pattacini
# email:  ugo.pattacini@iit.it

get_property(cer_kinematics_INCLUDE_DIRS
             TARGET cer_kinematics
             PROPERTY BUILD_INTERFACE_INCLUDE_DIRECTORIES)

include_directories(${YARP_INCLUDE_DIRS} ${cer_kinematics_INCLUDE_DIRS})

add_definitions(-D_USE_MATH_DEFINES)
add_executable(compute-trajectory compute-trajectory.cpp)
target_link_libraries(compute-trajectory ${YARP_LIBRARIES} cer_kinematics)

