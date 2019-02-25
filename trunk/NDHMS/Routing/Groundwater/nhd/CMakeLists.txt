cmake_minimum_required (VERSION 2.8)

add_library(nwm_routing_groundwater_nhd STATIC
	module_nhd_basin.F
	module_nhd_data.F
	module_nhd_tests.F              
)
