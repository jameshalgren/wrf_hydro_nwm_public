
module module_dynamic_wave_channel

    use module_hydrualic_channel_base

    !> @brief class that holds data about the geometry of channel objects
    type :: channel_geometry

    end type channel_geometry

    type, extends(hydraulic_channel_base_struct) :: dynamic_wave_channel_struct

    end type dynamic_wave_channel_struct

end module module_dynamic_wave_channel
