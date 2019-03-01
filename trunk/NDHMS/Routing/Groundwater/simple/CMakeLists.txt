cmake_minimum_required (VERSION 2.8)

add_library(nwm_routing_groundwater_simple STATIC
	module_simple_basin.F
	module_simple_data.F              
)

target_link_libraries(nwm_routing_groundwater_simple PUBLIC nwm_routing_groundwater)
target_link_libraries(nwm_routing_groundwater_simple PUBLIC nwm_routing_groundwater_bucket)
