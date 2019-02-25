cmake_minimum_required (VERSION 2.8)

add_library(nwm_routing_groundwater_simple STATIC
	module_simple_basin.F
	module_simple_data.F              
)
