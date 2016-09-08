# Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
# Author: Ugo Pattacini
# email:  ugo.pattacini@iit.it

get_property(cer_kinematics_INCLUDE_DIRS
             TARGET cer_kinematics
             PROPERTY BUILD_INTERFACE_INCLUDE_DIRECTORIES)

include_directories(${YARP_INCLUDE_DIRS} ${cer_kinematics_INCLUDE_DIRS})

add_definitions(-D_USE_MATH_DEFINES)
add_executable(handle-turn-compute-trajectory handle-turn-compute-trajectory.cpp)
target_link_libraries(handle-turn-compute-trajectory ${YARP_LIBRARIES} cer_kinematics)
install(TARGETS handle-turn-compute-trajectory DESTINATION bin)
