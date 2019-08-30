from DynamicNetwork import Network
from DynamicNetwork import Section

class RectangleSection(Section):
    def __init__(self, comid, bottom_width, bottom_z, dx_ds = 10, manning_n_ds = 0.015):
        #Time independent at-a-station properties
        self.station = None
        self.bottom_width = bottom_width
        self.bottom_z = bottom_z
        self.manning_n_ds = manning_n_ds
        self.time_steps = [] # array of values
        self.sk = 1.49
        
        #Time independent downstream reach properties
        self.dx_ds = 0 # Distance to downstream section
        self.dbdx_ds = 0 # Distance to downstream section
        self.bed_slope_ds = 0 # Bed slope (S0) to downstream section
#ADD NEIGHBOR Concept

class MESH_LookupNetwork(Network):
    #TODO: These Input and Initialize methods could be different methods within the Network class
    def input_and_initialize(self, input_opt=None, input_path=None, output_path=None, upstream_flow_ts=None, downstream_stage_ts=None):
        if input_opt == 1:
            '''Use File Input from the Mesh.py application'''

        elif input_opt == 2:
            ''' This input option is intended to be an extremely simple channel for testing and plotting development''' 
        
        else:
            input_vars = {}

            self.time_list = range(22)
            # import pandas as pd
            # pandas.date_range("11:00", "21:30", freq="30min")
            # datelist = pd.date_range(pd.datetime.today(), periods=100).tolist()

            n_sections = 11
            I_UPSTREAM = n_sections - 1
            I_DOWNSTREAM = 0

            station_downstream = 0
            station_upstream = 1000000
            stations = np.linspace(station_downstream, station_upstream, n_sections, False)
            bottom_widths = np.linspace(1000, 100, len(stations), False)
            bottom_zs = np.linspace(0,100, len(stations), False)

            # print(NCOMP, len(stations))

            input_vars.update({"dxini": 1000})
            input_vars.update({"manning_n_ds": 0.035})
            input_vars.update({"loss_coeff": 0.1})

            for i, bw in enumerate(bottom_widths):
                # continue
                self.sections.append(Section(station=stations[i]
                                        , bottom_width=bottom_widths[i]
                                        , bottom_z = bottom_zs[i]
                                        , manning_n_ds = input_vars['manning_n_ds']))
                self.sections[i].loss_coeff_ds = input_vars['loss_coeff']
                # print(sections[i].bed_slope_ds, sections[i].dx_ds, sections[i].bottom_z)
                if i == 0:
                    self.sections[i].dx_ds = input_vars['dxini'] #Irrelevant with the slope defined
                    self.sections[i].bed_slope_ds = 0.0001
                else:
                    self.sections[i].dx_ds = self.sections[i].station - self.sections[i-1].station
                    self.sections[i].bed_slope_ds = (self.sections[i].bottom_z \
                                                - self.sections[i-1].bottom_z) \
                                                / self.sections[i].dx_ds

            #TODO: clean up this code to generate intial upstream flow and downstream stage boundary time series
            self.upstream_flow_ts = Generate_Hydrograph(len(self.time_list) , 0 , 7 , 4 , 5000)
            self.downstream_stage_ts = [5*y_direct(self.sections[I_DOWNSTREAM].bottom_width
                                                 , self.sections[I_DOWNSTREAM].manning_n_ds
                                                 , self.sections[I_DOWNSTREAM].bed_slope_ds
                                                 , q ) for q in self.upstream_flow_ts]
            # print(self.upstream_flow_ts)
            # print(self.downstream_stage_ts)

        return input_vars

    def compute_initial_state(self):
        ''' Compute a steady initial state (this uses the same math as the next-
            time-step-state, only we simply assume we are using the first timestep
            of the boundary time-series.)
        '''
        # print(self.upstream_flow_ts)
        # print(self.downstream_stage_ts)
        self.add_steady_time_step(0, 0, self.sections, self.downstream_stage_ts[0], self.upstream_flow_ts[0])

    def compute_next_time_step_state(self, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next):
        ''' the Steady Network performs the Standard Step method to compute a steady profile
            for each flow and downstream stage in the series.
        '''

    def add_steady_time_step(self, j_current, j_next, sections, downstream_stage_next, upstream_flow_next):
        # TODO: Ask Nick -- should this be inside the Steady Timestep class?

        for i, section in enumerate(sections):
            ''' Step through using the standard step method
            '''
            # print(f'dssn {downstream_stage_next}')
            if i == 0: #Add known downstream boundary
                section.time_steps.append(self.TimeStep(new_flow = upstream_flow_next
                                                       ,new_depth = downstream_stage_next))
                # print(f'stsd_jcurr {section.time_steps[j_current].depth}')
                # print(f'sdstsd_jcurr {section_ds.time_steps[j_current].depth}')
                # print(f's0tsd_jcurr {sections[0].time_steps[j_current].depth}')
                # stage.append(downstream_stage_next)
                # WS.append(section.bottom_z + stage[0])
                # A.append(section.get_area_depth(stage[0]))
                continue

            section_ds = sections[i-1]
            #section
            # print(f'jnext {j_next}')
            # print(f'jcurr {j_current}')
            # print(f'stsd_jcurr {section.time_steps[j_current].depth}')
            # print(f'sdstsd_jcurr {section_ds.time_steps[j_current].depth}')
            # print(f's0tsd_jcurr {sections[0].time_steps[j_current].depth}')
            # print(f'stsd_jnext {section.time_steps[j_next].depth}')
            # print(f'sdstsd_jnext {section_ds.time_steps[j_next].depth}')
            # print(f's0tsd_jnext {sections[0].time_steps[j_next].depth}')
            #Use normal depth as an seed estimate convergence solution for standard depth
            y_guess = y_direct(section.bottom_width, section.manning_n_ds, section.bed_slope_ds, upstream_flow_next)
            y_standard = (self.y_standard_step(j_next, section_ds, section, upstream_flow_next, y_guess))

            section.time_steps.append(self.TimeStep(new_flow = upstream_flow_next, new_depth = y_standard))



    class TimeStep(MESH_LookupNetwork.TimeStep):
        def __init__(self, *args, **kwargs):
            # super(Network.TimeStep, self).__init__(*args, **kwargs)
            super().__init__(*args, **kwargs)

def main():
    
    # network = DummyNetwork()
    # network = SimpleFlowTrace() #DongHa's method.
    network = MESH_LookupNetwork()
    # network = MuskCNetwork()
    # network = MESHDNetwork()

    network.input_and_initialize()
    network.compute_initial_state()
    network.compute_time_steps()

if __name__ == "__main__":
    main()
