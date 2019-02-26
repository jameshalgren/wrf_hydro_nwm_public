cmake_minimum_required (VERSION 2.8)

# build the version static library
add_library(nwm_driver STATIC
	module_HYDRO_drv.F
)
