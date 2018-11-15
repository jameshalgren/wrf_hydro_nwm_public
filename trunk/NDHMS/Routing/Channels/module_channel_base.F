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

    !> @breif class that holds state variables for a channel base object
    type channel_base_state
        real :: current_flow    !< the flow on the channel at the current timestep
        real :: prev_flow       !< the flow on the channel at the previous timestep
        real :: depth           !< current depth of the channel
        real :: retention_depth !< minimal depth of water for this channel
        real :: nundge          !< observed value to nudge flow closer to
    end type


    !> @breif class that holds properties for a channel base object
    type channel_base_properties
        real :: mannings_num    !< the Mannings Number that approximates friction ect along the channel

        integer :: typel        !< the type of this channel (depreciated!)
        integer :: stream_order !< the steam order of this channel

        integer :: to_node      !< the destination node for this channel link
        integer :: from_node    !< the source of this channel link
    end type

    !> @breif class that holds input data needed by the physics code that processes channel base objects
    type channel_base_input
        double precision :: accumulated_surface_lateral_runnoff  !< accumulated surface runnoff for this link
        double precision :: accumulated_bucket_flux              !< accumulated groundwater flux for this link
    end type

    !> @breif class that holds output data produced by the physics code that processes channel base objects
    type channel_base_output
        real :: surface_lateral_flux  !< lateral surface flux from this link
        real :: bucket_flux           !< flux to groundwater from this link
    end type

    !> &brief class that is root of the hierarchy of channel objects. It holds data that is useful to most channel types
    !  and waterbody types
    type, abstract :: channel_base_struct

        type (channel_geometry), pointer :: geometry          !< structure that holds variables that describe channel geometry
        type (channel_base_state), pointer :: state           !< structure that holds the state of channel objects
        type (channel_base_properties), pointer :: properties !< structure that holds properties of the channel object
        type (channel_base_input), pointer :: input           !< structure that holds the input data for a channel object
        type (channel_base_output), pointer :: output         !< structure that holds the output data for a channel object


    end type channel_base_struct

end module module_channel_base
