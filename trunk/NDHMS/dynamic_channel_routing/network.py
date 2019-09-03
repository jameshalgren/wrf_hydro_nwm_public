import constants
import helpers

class Network:
    '''Class definition for reaches related as part of a computational scheme for
       open channel routing '''
    def __init__(self):
        '''initialize a new Network of sections/reaches'''
        self.sections = []
        self.time_list = [] # TODO: this initialization could be for a datetime series to contain the timestamps
        self.upstream_flow_ts = []
        self.downstream_stage_ts = []

    def input_and_initialize_simple(self):
        pass

    def compute_initial_state(self):
        pass

    def compute_time_steps(self):
        '''This function can operate with
        1) Nt and dt (number of time steps and size of time step) and a pointer to boundary information
        2) List of times and a pointer to boundary information
        3) an initial time, a list of time deltas, and a corresponding list of boundary conditions
         but since they really all boil down to the last situation, we'll just
         make it work for #3 and then have other translator methods that create these.'''

        for j, t in enumerate(self.time_list):
            # print(j+1 , len(self.time_list), len(self.upstream_flow_ts), len(self.downstream_stage_ts))
            if j+1 < len(self.time_list):
                self.compute_next_time_step_state(j_current = j
                                                  , j_next = j + 1
                                                  , upstream_flow_current = self.upstream_flow_ts[j]
                                                  , upstream_flow_next = self.upstream_flow_ts[j+1]
                                                  , downstream_stage_current = self.downstream_stage_ts[j]
                                                  , downstream_stage_next = self.downstream_stage_ts[j+1])
                # self.compute_next_time_step_state(t, j)

    def compute_next_time_step_state(self, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next): # Will be defined in child classes
                                         #TODO: find out how important it is to have the variables defined in this dummy function
                                         #TODO: change the function to be dependednt on the time value instead of simply time-step index.

        pass

    class TimeStep:
    #TODO: QUESTION FOR Nick
        ''' When we are passing the time steps out to the Fortran module,
        we only want to pass one timestep at a time and receive another one back.
        How does that happen best? '''
        def __init__(self, time_step = None, new_flow = None, new_depth = None):
            # Per-time-step at-a-section properties
            self.time = time_step
            self.flow = new_flow
            self.depth = new_depth

            # Per-time-step downstream reach properties
            self.friction_slope_ds = 0

    def add_time_step(self, section, new_flow, new_depth):
        section.time_steps.append(self.TimeStep(new_flow = new_flow, new_depth=new_depth))

    def add_upstream_boundary_condition_time_step(self, section, upstream_flow):
        section.time_steps.append(self.TimeStep(new_flow = upstream_flow))

    def add_downstream_boundary_condition_time_step(self, section, downstream_depth):
        section.time_steps.append(self.TimeStep(new_depth = downstream_depth))

    def add_normal_depth_time_step(self, section, new_flow):
        new_depth = helpers.y_direct(section.bottom_width, section.manning_n_ds, section.bed_slope_ds, new_flow)
        section.time_steps.append(self.TimeStep(new_flow=new_flow, new_depth=new_depth))

    class Section:
        #TODO: The Section Class needs to be sub-classed with Different types,
        #e.g., SectionRectangle, SectionTrapezoid, SectionTrapFlood (for the type that
        #currently used in the National Water Model), SectionDepthArea, SectionDepthWidth, ...
        #def __init__(self, bottom_width, side_slope):
        def __init__(self, bottom_z, comid=None, station=None, dx_ds = 10):
            #Time independent at-a-station properties
            self.comid = comid
            self.station = station
            self.bottom_z = bottom_z
            self.time_steps = [] # array of values

            #Time independent downstream reach properties
            self.dx_ds = 0 # Distance to downstream section
            self.loss_coeff_ds = 0 # Contraction and other loss coefficients to downstream section
                                # C in the following equation
                                # hl = Sf * dx + C * (V1**2/2g - V2**2/2g)
            self.bed_slope_ds = 0 # Bed slope (S0) to downstream section
            #ADD NEIGHBOR Concept

    class RectangleSection(Section):
        def __init__(self, bottom_width, manning_n_ds = 0.015, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.bottom_width = bottom_width
            self.manning_n_ds = constants.MANNING_SI
            #self.sk = constants.MANNING_SI

        def get_area_depth(self, depth):
            return self.bottom_width * depth

        def get_area_j(self, j):
            return self.bottom_z * self.time_steps[j].depth

        def get_wetted_perimeter_depth(self, depth):
            return self.bottom_width + 2.0 * depth

        def get_wetted_perimeter_j(self, j):
            return self.bottom_z + 2.0 * self.time_steps[j].depth
