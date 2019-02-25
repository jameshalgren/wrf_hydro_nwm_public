cmake_minimum_required (VERSION 2.8)

add_library(nwm_routing STATIC
	module_channel_routing.F
	module_date_utilities_rt.F
	module_GW_baseflow.F
	module_gw_gw2d.F
	module_HYDRO_io.F
	module_HYDRO_utils.F
	module_lsm_forcing.F
	module_noah_chan_param_init_rt.F
	module_NWM_io_dict.F
	module_NWM_io.F
	module_reservoir_routing.F
	module_RT.F
	Noah_distr_routing.F
	Noah_distr_routing_overland.F
	Noah_distr_routing_subsurface.F        
)
