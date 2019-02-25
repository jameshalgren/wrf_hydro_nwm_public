cmake_minimum_required (VERSION 2.8)

add_library(nwm_routing_groundwater STATIC
	module_UDMAP.F
	module_groundwater_base.F       
	module_groundwater_static_data.F          
)

add_subdirectory("bucket")
add_subdirectory("nhd")
add_subdirectory("simple")

