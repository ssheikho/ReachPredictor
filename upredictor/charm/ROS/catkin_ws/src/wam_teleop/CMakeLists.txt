cmake_minimum_required(VERSION 2.8.3)
project(wam_teleop)

find_package(catkin
	REQUIRED COMPONENTS
	message_generation geometry_msgs wam_msgs wam_srvs roscpp)
# std_msgs sensor_msgs

#catkin_package(
#	INCLUDE_DIRS include
#)

include_directories(include
	${catkin_INCLUDE_DIRS}
)
	
# Declare the service files to be built
add_executable(wam_joystick_teleop
	src/wam_joystick_teleop.cpp
)

target_link_libraries(wam_joystick_teleop ${catkin_LIBRARIES})

add_dependencies(wam_joystick_teleop ${catkin_EXPORTED_TARGETS})
add_dependencies(wam_joystick_teleop ${catkin_EXPORTED_TARGETS})

# Declare that this catkin package's runtime dependencies
catkin_package(CATKIN_DEPENDS message_runtime)
# std_msgs sensor_msgs
