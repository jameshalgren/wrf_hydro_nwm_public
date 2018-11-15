
module module_hydrologic_base

    use module_channel_base

    type, abstract, extends(channel_base_properties) :: hydrologic_base_properties

    end type hydrologic_base_properties

    type, abstract, extends(channel_base_struct) :: hydrologic_channel_base_struct

    end type hydrologic_channel_base_struct

end module module_hydrologic_base
