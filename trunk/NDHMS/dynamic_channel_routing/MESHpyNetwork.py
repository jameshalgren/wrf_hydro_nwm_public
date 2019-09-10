from __future__ import division
import helpers
import constants
from network import Network
import os

class MESHpyNetwork(Network):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def compute_initial_state(self):
        ''' Compute a steady initial state (this uses the same math as the next-
            time-step-state, only we simply assume we are using the first timestep
            of the boundary time-series.)
        '''
        #print(self.upstream_flow_ts)
        #print(self.downstream_stage_ts)
        self.add_normal_time_step(0, 0, self.sections, self.downstream_stage_ts[0], self.upstream_flow_ts[0])

    def compute_next_time_step_state(self, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next):
        ''' the Steady Network performs the Standard Step method to compute a steady profile
            for each flow and downstream stage in the series.
        '''
        self.add_normal_time_step(j_current, j_next, self.sections, downstream_stage_next, upstream_flow_next)

    def add_normal_time_step(self, j_current, j_next, sections, downstream_stage_next, upstream_flow_next):
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
            y_guess = helpers.y_direct(section.bottom_width, section.manning_n_ds, section.bed_slope_ds, upstream_flow_next)
            section.time_steps.append(self.TimeStep(new_flow = upstream_flow_next, new_depth = y_guess))



    class TimeStep(Network.TimeStep):
        def __init__(self, *args, **kwargs):
            # super(Network.TimeStep, self).__init__(*args, **kwargs)
            super().__init__(*args, **kwargs)

            # Per-time-step at-a-section properties
            self.flow = 0
            self.delta_flow_corrector = 0
            self.delta_flow_predictor = 0 
            self.delta_area_corrector = 0 
            self.delta_area_predictor = 0 
            self.depth = 0
            self.water_z = 0
            self.flow_area = 0
            self.ci1 = 0
            self.hy = 0 #(temporary variable for computing co)
            self.c0 = 0

            # Per-time-step downstream reach properties
            self.ci2_ds = 0
            # self.friction_slope_ds = 0 # Derived from parent Class, Network.TimeStep
            self.as0_ds = 0
            self.gs0_ds = 0
            self.sigma_ds = 0 # Sigma is approximately equivalent to the courant parameter, but not quite
            self.dbdx = 0

def main():
    
    input_type = 'file'
    input_vars = {}
    input_vars['filetype'] = 'mesh.py'
    root = os.path.abspath(r'c:/Users/james.halgren/Downloads/MESH_test/')
    Main_Example_Path = os.path.join(root , 'US')
    Sub_Example_Path = os.path.join(Main_Example_Path , 'BW')
    This_Example_Path = os.path.join(Sub_Example_Path, 'Q')

    #C:\Users\james.halgren\Downloads\MESH_test\US\BW\Q\Qvar_us_2YNorm\Qvar_us_0033_5.0-10000.0_0100_0000-0004-0200_2NormalDepth
    input_path = os.path.join(This_Example_Path,'Qvar_us_2YNorm','Qvar_us_0033_5.0-10000.0_0100_0000-0004-0200_2NormalDepth',"input.txt")
    input_vars['input_path'] = input_path
    # print(input_path)
    
    # if len(sys.argv) > 1:
    #     input_path = sys.argv[1]
    # else:
    #     input_path = os.path.join(This_Example_Path,"input.txt")
    #     print(input_path)
    #     #input_path = "./input.txt"

    # network = DummyNetwork()
    # network = SimpleFlowTrace() #DongHa's method.
    # network = SteadyNetwork(input_type = input_type, input_vars = input_vars)
    #input_and_initialize(sections, input_path, input_opt)
    network = MESHpyNetwork(input_type = input_type, input_vars = input_vars)
    # network = MuskCNetwork()
    # network = MESHDNetwork()

    network.compute_initial_state()
    network.compute_time_steps(verbose = True)

if __name__ == "__main__":
    main()
