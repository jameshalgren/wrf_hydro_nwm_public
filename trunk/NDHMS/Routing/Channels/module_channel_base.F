!> @brief this module holds the data structures used to define the base class for the hierarchy of channel objects

module module_channel_base

    !> @brief class that holds data about the geometry of channel objects
    type :: channel_geometry

        real :: channel_length   !< the length of the channel or waterbody along the primary direction of flow
        real :: slope            !< the slope of channel
        real :: side_slope       !< the side slope of the channel
        real :: top_width        !< the width of channel at top bank
        real :: bottom_width     !< the width at the bottom of the channel

    end type channel_geometry

    type channel_base_state

    end type

    !> &brief class that is root of the hierarchy of channel objects. It holds data that is useful to most channel types
    !  and waterbody types
    type, abstract :: channel_base_struct

        type (channel_geometry) :: geometry !< structure that holds variables that describe channel geometry

        real :: mannings_number  !< the Mannings Number that approximates friction ect along the channel


    end type channel_base_struct

end module module_channel_base
