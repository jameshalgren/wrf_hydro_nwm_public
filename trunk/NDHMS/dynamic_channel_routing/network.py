import constants
import helpers
import numpy as np
import csv

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
        with open(input_path, newline='') as f:
    
            # Put the first chunk of each line into a lsit
            data = list(map(lambda x:x.split(' ')[0] , f.read().split('\n')))

            #TODO: Get rid of this kludge to cast the input numbers into the right datatype
            for i, item in enumerate (data[0:23]):
                data[i] = float(item)

            data[3] = int(data[3])
            data[4] = int(data[4])

            # Assign all the input values into the variables
            dtini, dxini, tfin, n_sections, ntim, phi, theta, thetas, thesinv, alfa2,\
                alfa4, f, skk, yy, qq, cfl, time_step_optimization, yw, bw, w, option, yn, qn, igate,\
                bed_elevation_path, upstream_path, downstream_path, channel_width_path,\
                output_path, option_dsbc, null = data

        I_UPSTREAM = n_sections - 1
        I_DOWNSTREAM = 0

        # Read in bed elevation and bottom width input
        with open(bed_elevation_path, newline='') as file1:
            with open(channel_width_path, newline='') as file2:
                read_data1 = list(csv.reader(file1, delimiter=' '))
                read_data2 = list(csv.reader(file2, delimiter=' '))
                for i in range(n_sections): # use usual convention of i as spatial dimension
                    z = float(read_data1[i][1])
                    y0 = z + yy #TODO: This seems like a clunky definition of the initial water surface
                                #      Elevation and I think we can do better.
                    b0 = float(read_data2[i][1])
                    self.sections.append(Network.RectangleSection(i, b0, z, dxini)) # dx is a static value in the test cases

        # Read hydrograph input Upstream and Downstream
        with open(upstream_path, newline='') as file3:
            with open(downstream_path, newline='') as file4:
                read_data3 = list(csv.reader(file3, delimiter=' '))
                read_data4 = list(csv.reader(file4, delimiter=' '))
                for j in range(ntim): # use usual convention of j as time dimension
                    self.upstream_flow_ts.append(float(read_data3[j][1]))
                    self.downstream_stage_ts.append(float(read_data4[j][1]))
                    self.time_list.append(j)
                    #self.sections[I_UPSTREAM].time_steps.append(self.TimeStep(new_flow = q_upstream))
                    #self.sections[I_DOWNSTREAM].time_steps.append(self.TimeStep(new_depth = y_downstream))
#                else:
#                    #TODO: Work with Nick to get this data dictionary thing passing
#                    #into and out of this function properly
#                    #TODO: INSERT code to generate intial sections and boundary time series
#                    input_data.update({"dtini": 10.0})
#                    input_data.update({"dxini": 20.0})
#                    input_data.update({"tfin": 5000.})
#                    input_data.update({"ncomp": 501})
#                    input_data.update({"ntim": 5000})
#                    input_data.update({"phi": 1.0})
#                    input_data.update({"theta": 1.0})
#                    input_data.update({"thetas": 1.0})
#                    input_data.update({"thesinv": 1.0})
#                    input_data.update({"alfa2": 0.5})
#                    input_data.update({"alfa4": 0.1})
#                    input_data.update({"f": 1.0})
#                    input_data.update({"skk": 20.0 })
#                    input_data.update({"yy": 6.0})
#                    input_data.update({"qq": 100.0})
#                    input_data.update({"cfl": 1.0})
#                    input_data.update({"ots": 0.0})
#                    input_data.update({"yw": 0.0})
#                    input_data.update({"bw": 20.0})
#                    input_data.update({"w": 1.1})
#                    input_data.update({"option": 1.0})
#                    input_data.update({"yn": 0.1000})
#                    input_data.update({"qn": 0.0085})
#                    input_data.update({"igate": 700})
#
#
#                    time_steps = range(100)
#                    stations = range(10000,11001, 100)
#                    bottom_widths = range(100, 1001, 100)
#                    bottom_zs = range(0,100,10)
#
#                    upstream_flows = Generate_Hydrograph(100 , 20 , 2 , 4 , 5000)
#                    # for i, station, bottom_width, bottom_z in enumerate(stations):
#                    #     print(f'{i} {station} {bottom_width} {bottom_z}')
#                    # for i, flow in enumerate(upstream_flows):
#                    #     print(f'{i} {flow}')
#                    sections = []
#
#                    for i, bw in enumerate(bottom_widths):
#                        sections.append(Section(stations[i], bottom_widths[i], bottom_zs[i]))
#                        if i == 0:
#                            sections[i].dx_ds = 10
#                            sections[i].bed_slope_ds = .1
#                        else:
#                            sections[i].dx_ds = sections[i].station - sections[i-1].station
#                            sections[i].bed_slope_ds = (sections[i].bottom_z - \
#                                                        sections[i-1].bottom_z)/ \
#                                                        sections[i].dx_ds
#
#
#                    return section_arr, input_data
#
    #TODO: These Input and Initialize methods could be different methods within the Network class
    #TODO: Make GRAVITY and MANNING_SI constants consistent with anticipated units in the input step and
    #get them to be called/passsed consistently.
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
        self.time_list = range(n_timesteps)
        #TODO: Convert all timesteps to time_stamps
        # import pandas as pd
        # pandas.date_range("11:00", "21:30", freq="30min")
        # datelist = pd.date_range(pd.datetime.today(), periods=100).tolist()

        I_UPSTREAM = n_sections - 1
        I_DOWNSTREAM = 0

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
        self.upstream_flow_ts = helpers.Generate_Hydrograph(len(self.time_list) , hydrograph_steady_time
                                                                                , hydrograph_event_width
                                                                                , hydrograph_skewness
                                                                                , hydrograph_qpeak)
#        self.time_list, self.downstream_stage_ts = [zip(j, 5*helpers.y_direct(self.sections[I_DOWNSTREAM].bottom_width
#                                             , self.sections[I_DOWNSTREAM].manning_n_ds
#                                             , self.sections[I_DOWNSTREAM].bed_slope_ds
#                                             , q )) for j, q in enumerate[self.upstream_flow_ts]]
#
        self.time_list = [j for j, _ in enumerate(self.upstream_flow_ts)]

        self.downstream_stage_ts = [5*helpers.y_direct(self.sections[I_DOWNSTREAM].bottom_width
                                              , self.sections[I_DOWNSTREAM].manning_n_ds
                                              , self.sections[I_DOWNSTREAM].bed_slope_ds
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
