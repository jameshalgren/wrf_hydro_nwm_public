from __future__ import division
import helpers
import constants
from network import Network
import os
import meshfunc as meshfunc

class MESHpyNetwork(Network):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def input_and_initialize_simple(section_arr, input_path=None, input_opt=1, output_path=None):
        '''Override the simple input_and_initialize function from network.py'''
        input_data = {}
        #TODO: Work with Nick to get this data dictionary thing passing
        #into and out of this function properly
        #TODO: INSERT code to generate intial sections and boundary time series
        input_data.update({"dtini": 10.0})
        input_data.update({"dxini": 20.0})
        input_data.update({"tfin": 5000.})
        input_data.update({"ncomp": 501})
        input_data.update({"ntim": 5000})
        input_data.update({"phi": 1.0})
        input_data.update({"theta": 1.0})
        input_data.update({"thetas": 1.0})
        input_data.update({"thesinv": 1.0})
        input_data.update({"alfa2": 0.5})
        input_data.update({"alfa4": 0.1})
        input_data.update({"f": 1.0})
        input_data.update({"skk": 20.0 })
        input_data.update({"yy": 6.0})
        input_data.update({"qq": 100.0})
        input_data.update({"cfl": 1.0})
        input_data.update({"ots": 0.0})
        input_data.update({"yw": 0.0})
        input_data.update({"bw": 20.0})
        input_data.update({"w": 1.1})
        input_data.update({"option": 1.0})
        input_data.update({"yn": 0.1000})
        input_data.update({"qn": 0.0085})
        input_data.update({"igate": 700})


        time_steps = range(100)
        stations = range(10000,11001, 100)
        bottom_widths = range(100, 1001, 100)
        bottom_zs = range(0,100,10)

        upstream_flows = Generate_Hydrograph(100 , 20 , 2 , 4 , 5000)
        # for i, station, bottom_width, bottom_z in enumerate(stations):
        #     print(f'{i} {station} {bottom_width} {bottom_z}')
        # for i, flow in enumerate(upstream_flows):
        #     print(f'{i} {flow}')
        sections = []

        for i, bw in enumerate(bottom_widths):
            sections.append(Section(stations[i], bottom_widths[i], bottom_zs[i]))
            if i == 0:
                sections[i].dx_ds = 10
                sections[i].bed_slope_ds = .1
            else:
                sections[i].dx_ds = sections[i].station - sections[i-1].station
                sections[i].bed_slope_ds = (sections[i].bottom_z - \
                                            sections[i-1].bottom_z)/ \
                                            sections[i].dx_ds

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

        # secpred, apply_corrector, and matrixc should be replaced by
        # section, apply_predictor, and matrixp, respectively.
        # There would need to be a flag of some sort to handle the direction.
        meshfunc.compute_sections(self.sections, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next)
        meshfunc.matrixp(self.sections, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next)
        meshfunc.apply_predictor(self.sections, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next)
        meshfunc.secpred(self.sections, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next)
        meshfunc.matrixc(self.sections, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next)
        meshfunc.apply_corrector(self.sections, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next)

        self.add_normal_time_step(j_current, j_next, self.sections, downstream_stage_next, upstream_flow_next)

    def add_normal_time_step(self, j_current, j_next, sections, downstream_stage_next, upstream_flow_next):
        # TODO: Ask Nick -- should this be inside the Steady Timestep class?

        for i, section in enumerate(sections):
            # print(f'downstream stage next {downstream_stage_next}')
            if i == 0: #Add known downstream boundary
                section.time_steps.append(self.TimeStep(new_flow = upstream_flow_next
                                                        , new_depth = downstream_stage_next))
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
            self.delta_flow_corrector = 0
            self.delta_flow_predictor = 0
            self.delta_area_corrector = 0
            self.delta_area_predictor = 0
            self.water_z = 0
            self.flow_area = 0
            self.ci1 = 0
            self.hy = 0 # Hydraulic Radius (used to compute co)
            self.conveyance = 0

            # Per-time-step downstream reach properties
            self.ci2_ds = 0
            # self.friction_slope_ds = 0 # Derived from parent Class, Network.TimeStep
            self.as0_ds = 0
            self.gs0_ds = 0
            self.sigma_ds = 0 # Sigma is approximately equivalent to the courant parameter, but not quite
            self.dbdx = 0
            self.velocity = 0
            self.celerity = 0
            self.cour = 0
            self.b11 = 0
            self.b12 = 0
            self.b21 = 0
            self.b22 = 0
            self.g11inv = 0
            self.g12inv = 0
            self.g21inv = 0
            self.g22inv = 0

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
