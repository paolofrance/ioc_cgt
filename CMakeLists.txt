cmake_minimum_required(VERSION 3.1)

project(ioc_cgt)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
  urdf
  cnr_controller_interface
  state_space_filters
  cnr_hardware_interface
  pluginlib
  sensor_msgs
  std_msgs
  name_sorting
  joint_limits_interface
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)




catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ioc_cgt
 CATKIN_DEPENDS urdf cnr_controller_interface state_space_filters cnr_hardware_interface pluginlib sensor_msgs std_msgs name_sorting joint_limits_interface
 DEPENDS EIGEN3
)

include_directories( include ${catkin_INCLUDE_DIRS} )




#add_library(${PROJECT_NAME} src/${PROJECT_NAME}/ioc_cgt.cpp )
#add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
#target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} Eigen3::Eigen)
#target_compile_options(${PROJECT_NAME} PUBLIC -Wall -faligned-new
#      $<$<CONFIG:Release>:-Ofast -funroll-loops -ffast-math >)

#if(${CMAKE_VERSION} VERSION_GREATER  "3.16.0")
#target_precompile_headers(${PROJECT_NAME} PUBLIC
#  <cmath>
#  <Eigen/Core>
#   <cnr_logger/cnr_logger.h>
#   <ros/time.h>
#   <sensor_msgs/JointState.h>
# )
#endif()




