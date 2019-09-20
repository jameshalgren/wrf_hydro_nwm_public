import constants
import helpers
import numpy as np
import csv
from math import sqrt

class Network:
    '''Class definition for reaches related as part of a computational scheme for
       open channel routing '''
    #TODO: Somewhere, there will need to be a de-tangling of how we call and initialize a rectangular channel vs.
    #      vs. trapezoidal, vs. generalized, etc. Perhaps that can be handled within the input files...
    #      For now, assume rectangular
    def __init__(self, input_type = 'simple', input_vars = None):
        '''initialize a new Network of sections/reaches'''
        self.sections = []
        self.time_list = [] # TODO: this initialization could be for a datetime series to contain the timestamps
        self.upstream_flow_ts = []
        self.downstream_stage_ts = []

        self.I_UPSTREAM = 1
        self.I_DOWNSTREAM = 0
        self.gravity = constants.GRAVITY_SI
        self.manning_m = constants.MANNING_M

        if input_vars:
            if input_type is 'simple':
                self.input_and_initialize_simple(**input_vars)
                #helped by this SO post:
                #https://stackoverflow.com/questions/334655/passing-a-dictionary-to-a-function-as-keyword-parameters
                               #
                               # n_sections = input_vars['n_sections']
                               # , n_timesteps = input_vars['n_timesteps']
                               # , station_downstream = input_vars['station_downstream']
                               # , station_upstream = input_vars['station_upstream']
                               # , bottom_width_downstream = input_vars['bottom_width_downstream']
                               # , bottom_width_upstream = input_vars['bottom_width_upstream']
                               # , bottom_z_downstream = input_vars['bottom_z_downstream']
                               # , bottom_z_upstream = input_vars['bottom_z_upstream']
                               # , dx_ds_boundary = input_vars['dx_ds_boundary']
                               # , S0_ds_boundary = input_vars['S0_ds_boundary']
                               # , manning_n_ds_all = input_vars['manning_n_ds_all']
                               # , loss_coeff_all = input_vars['loss_coeff_all']
                               # , hydrograph_steady_time = input_vars['hydrograph_steady_time']
                               # , hydrograph_event_width = input_vars['hydrograph_event_width']
                               # , hydrograph_skewness = input_vars['hydrograph_skewness']
                               # , hydrograph_qpeak = input_vars['hydrograph_qpeak']
                               # )
            elif input_type == 'file':
                filetype = input_vars['filetype']
                if filetype == 'json': #TODO: Replace json psuedocode
                    file = read(input_path)
                    for v in file:
                        pass
                elif filetype == 'mesh.py':
                    self.input_and_initialize_meshpyfile(**input_vars)
                else: print('not-yet-implemented input_type')
            else: #If nothing was defined, prepare a simple network with some default parameters.
            #TODO: This should be properly error trapped/handled.
                self.input_and_initialize_simple(n_sections = 11
                                            , n_timesteps = 22
                                            , station_downstream = 0
                                            , station_upstream = 10000000
                                            , bottom_width_downstream = 100
                                            , bottom_width_upstream = 1000
                                            , bottom_z_downstream = 0
                                            , bottom_z_upstream = 100
                                            , dx_ds_boundary = 1000
                                            , S0_ds_boundary = 0.0001
                                            , manning_n_ds_all = 0.035
                                            , loss_coeff_all = 0.01
                                            , hydrograph_steady_time = 0
                                            , hydrograph_event_width = 7
                                            , hydrograph_skewness = 4
                                            , hydrograph_qpeak = 5000)


    def input_and_initialize_meshpyfile(self, filetype = None, input_path=None):
        pass

    def input_and_initialize_simple(self, n_sections = 11
                                        , n_timesteps = 22
                                        , station_downstream = 0
                                        , station_upstream = 10000000
                                        , bottom_width_downstream = 100
                                        , bottom_width_upstream = 1000
                                        , bottom_z_downstream = 0
                                        , bottom_z_upstream = 100
                                        , dx_ds_boundary = 1000
                                        , S0_ds_boundary = 0.0001
                                        , manning_n_ds_all = 0.035
                                        , loss_coeff_all = 0.01
                                        , hydrograph_steady_time = 0
                                        , hydrograph_event_width = 7
                                        , hydrograph_skewness = 4
                                        , hydrograph_qpeak = 5000):
        ''' This input option is intended to be an extremely simple channel for testing and plotting development'''
        self.I_UPSTREAM = n_sections - 1
        self.I_DOWNSTREAM = 0

        stations = np.linspace(station_downstream, station_upstream, n_sections, False)
        bottom_widths = np.linspace(bottom_width_upstream, bottom_width_downstream, len(stations), False)
        bottom_zs = np.linspace(bottom_z_downstream,bottom_z_upstream, len(stations), False)

        # print(NCOMP, len(stations))

        for i, bw in enumerate(bottom_widths):
            # continue
            self.sections.append(Network.RectangleSection(station=stations[i]
                                    , bottom_width=bottom_widths[i]
                                    , bottom_z = bottom_zs[i]
                                    , manning_n_ds = manning_n_ds_all))
            self.sections[i].loss_coeff_ds = loss_coeff_all
            # print(sections[i].bed_slope_ds, sections[i].dx_ds, sections[i].bottom_z)
            if i == 0:
                self.sections[i].dx_ds = dx_ds_boundary #Irrelevant with the slope defined
                self.sections[i].bed_slope_ds = S0_ds_boundary
            else:
                self.sections[i].dx_ds = self.sections[i].station - self.sections[i-1].station
                self.sections[i].bed_slope_ds = (self.sections[i].bottom_z \
                                            - self.sections[i-1].bottom_z) \
                                            / self.sections[i].dx_ds

        #TODO: clean up this code to generate intial upstream flow and downstream stage boundary time series
        self.upstream_flow_ts = helpers.Generate_Hydrograph(n_timesteps , hydrograph_steady_time
                                                                                , hydrograph_event_width
                                                                                , hydrograph_skewness
                                                                                , hydrograph_qpeak)
#        self.time_list, self.downstream_stage_ts = [zip(j, 5*helpers.y_direct(self.sections[self.I_DOWNSTREAM].bottom_width
#                                             , self.sections[self.I_DOWNSTREAM].manning_n_ds
#                                             , self.sections[self.I_DOWNSTREAM].bed_slope_ds
#                                             , q )) for j, q in enumerate(self.upstream_flow_ts)]

        self.time_list = [j for j, _ in enumerate(self.upstream_flow_ts)]
        #TODO: Convert all timesteps to time_stamps
        # import pandas as pd
        # pandas.date_range("11:00", "21:30", freq="30min")
        # datelist = pd.date_range(pd.datetime.today(), periods=100).tolist()


        self.downstream_stage_ts = [5*helpers.y_direct(self.sections[self.I_DOWNSTREAM].bottom_width
                                              , self.sections[self.I_DOWNSTREAM].manning_n_ds
                                              , self.sections[self.I_DOWNSTREAM].bed_slope_ds
                                              , q ) for q in self.upstream_flow_ts]

        # print(self.upstream_flow_ts)
        # print(self.downstream_stage_ts)

    def compute_initial_state(self):
        pass

    def compute_time_steps(self, verbose=False):
        '''This function can operate with
        1) Nt and dt (number of time steps and size of time step) and a pointer to boundary information
        2) List of times and a pointer to boundary information
        3) an initial time, a list of time deltas, and a corresponding list of boundary conditions
         but since they really all boil down to the last situation, we'll just
         make it work for #3 and then have other translator methods that create these.'''

        # print(self.time_list)
        for j, t in enumerate(self.time_list):
            if verbose: print(j+1 , len(self.time_list), len(self.upstream_flow_ts), len(self.downstream_stage_ts))
            if verbose: print(f'timestep {j} {t}')
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

    def write_state(self, type, file):
        output = self.output_state(type, file)
        with open(file, 'w') as output_file:
            output_file.write(output)

    def output_state(self, type):
        if type is 'pickle':
            state = self.pickle_output(self.network)
        else:
            print("only 'pickle' output is implemented")
        return state

    def pickle_output(self, network):
        return pickle.dumps(network.sections)

    def output_writer(time_step, output_opt = 1, output_path = None):
        ##WARNING THIS FUNCTION DOESN'T WORK at this point!
        #TODO: convert this to Section-based code
        if output_opt == 1:
            # Open files for output
            with open(f'{output_path}/flow_area.txt', 'w') as area_output: # unit 8
                with open(f'{output_path}/flow.txt', 'w') as flow_output: # unit 9
                    # Output initial condition
                    area_output.write(f'{t} {" ".join((y[0][i] - z[i]) * bo[i] for i in range(n_sections))}\n')
                    flow_output.write(f'{t}  {" ".join(q[0][i] for i in range(n_sections))}\n')
        else:
            print(f'{t} {" ".join((y[0][i] - z[i]) * bo[i] for i in range(n_sections))}\n')
            print(f'{t}  {" ".join(q[0][i] for i in range(n_sections))}\n')

    class TimeStep:
        # TODO: Juzer find how to pass the time steps out to the Fortran module,
        # we only want to pass one timestep at a time and receive another one back.
        # How does that happen best? '''
        def __init__(self, new_time = None, new_flow = None, new_depth = None):
            # Per-time-step at-a-section properties
            self.time = new_time
            self.flow = new_flow
            self.depth = new_depth

            # Per-time-step downstream reach properties
            self.friction_slope_ds = 0

    def add_time_step(self, section, new_flow, new_depth): #TODO: the Self and Section inputs are probably redundant
        section.time_steps.append(self.TimeStep(new_flow = new_flow, new_depth=new_depth))

    def add_upstream_boundary_condition_time_step(self, section, upstream_flow): #TODO: The Self and Section inputs are probably redundant
        section.time_steps.append(self.TimeStep(new_flow = upstream_flow))

    def add_downstream_boundary_condition_time_step(self, section, downstream_depth): #TODO: the Self and Section inputs are probably redundant
        section.time_steps.append(self.TimeStep(new_depth = downstream_depth))

    def add_normal_depth_time_step(self, section, new_flow): #TODO: the Self and Section inputs are probably redundant
        new_depth = helpers.y_direct(section.bottom_width, section.manning_n_ds, section.bed_slope_ds, new_flow)
        section.time_steps.append(self.TimeStep(new_flow=new_flow, new_depth=new_depth))

    class Section:
        #TODO: The Section Class needs to be sub-classed with Different types,
        #e.g., RectangleSection, TrapezoidSection, TrapFloodSection (for the type that
        #currently used in the National Water Model), DepthAreaSection, DepthWidthSection, ...
        #def __init__(self, bottom_width, side_slope):
        def __init__(self, bottom_z, comid=None, station=None, dx_ds = 10):
            #Time independent at-a-station properties
            self.comid = comid
            self.station = station
            self.bottom_z = bottom_z

            self.time_steps = [] # array of values

            #Time independent downstream reach properties
            self.dx_ds = dx_ds # Distance to downstream section
            self.loss_coeff_ds = 0 # Contraction and other loss coefficients to downstream section
                                # C in the following equation
                                # hl = Sf * dx + C * (V1**2/2g - V2**2/2g)
            self.bed_slope_ds = 0 # Bed slope (S0) to downstream section
            #ADD NEIGHBOR Concept

    class RectangleSection(Section):
        def __init__(self, bottom_width, manning_n_ds = 0.015, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.bottom_width = bottom_width
            self.manning_n_ds = manning_n_ds
            self.dbdx_ds = 0 # change in width to the next downstream section

        def get_celerity_area(self, area, gravity, debug = False):
            if debug: print(gravity, area, self.bottom_width)
            return sqrt(gravity * area / self.bottom_width)

        def get_depth_area(self, area):
            return area / self.bottom_width

        def get_area_depth(self, depth):
            return self.bottom_width * depth

        def get_area_j(self, j):
            return self.bottom_z * self.time_steps[j].depth

        def get_wetted_perimeter_area(self, area):
            return self.bottom_width + 2.0 * area / self.bottom_width

        def get_wetted_perimeter_depth(self, depth):
            return self.bottom_width + 2.0 * depth

        def get_wetted_perimeter_j(self, j):
            return self.bottom_z + 2.0 * self.time_steps[j].depth

    class IrregularSection(Section):

        def get_wetted_perimeter_area(self, area):
            return self.bottom_width + 2.0 * area / self.bottom_width
