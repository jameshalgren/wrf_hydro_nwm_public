cmake_minimum_required (VERSION 2.8)

add_library(nwm_routing_groundwater_bucket STATIC
	module_bucket_impl.F
	module_bucket_properties.F
	module_bucket_state.F              
)
