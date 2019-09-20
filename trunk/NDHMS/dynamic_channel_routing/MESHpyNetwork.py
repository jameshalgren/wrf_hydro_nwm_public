# import required modules
from __future__ import division
import helpers
import constants
import meshconstants
from network import Network
import csv
import os

class MESHpyNetwork(Network):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        '''USE A Global Declaration here to manage these values'''
        self.debug = False
        self.dx_tolerance = meshconstants.DX_TOLERANCE
        self.depth_tolerance = meshconstants.DEPTH_TOLERANCE
        self.celerity_epsilon = meshconstants.CELERITY_EPSILON
        self.crmax = 0.0
        self.crmin = 100.0
        self.predictor_step = True
        self.yy = 0.0
        self.qq = 0.0
        self.phi = meshconstants.PHI         # source term treatment (0:explicit, 1:implicit)
        self.theta = meshconstants.THETA     # ?
        self.thetas = meshconstants.THETAS   # ?
        self.thesinv = meshconstants.THESINV # ?
        self.alfa2 = meshconstants.ALFA2     # emp parameter for artificial diffusion (lit)
        self.alfa4 = meshconstants.ALFA4     # maximum value for artificial diffusion

    def input_and_initialize_meshpysimple(section_arr, input_path=None, input_opt=1, output_path=None):
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

        upstream_flows = helpers.Generate_Hydrograph(100 , 20 , 2 , 4 , 5000)
        # for i, station, bottom_width, bottom_z in enumerate(stations):
        #     print(f'{i} {station} {bottom_width} {bottom_z}')
        # for i, flow in enumerate(upstream_flows):
        #     print(f'{i} {flow}')

        for i, bw in enumerate(bottom_widths):
            self.sections.append(Section(stations[i], bottom_widths[i], bottom_zs[i]))
            if i == self.I_DOWNSTREAM:
                self.sections[i].dx_ds = 10
                self.sections[i].bed_slope_ds = .1
            else:
                self.sections[i].dx_ds = self.sections[i].station - self.sections[i-1].station
                self.sections[i].bed_slope_ds = (self.sections[i].bottom_z - \
                                            self.sections[i-1].bottom_z)/ \
                                            self.sections[i].dx_ds

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
            dtini, dxini, tfin, n_sections, ntim, self.phi, self.theta, self.thetas\
                , self.thesinv, self.alfa2, self.alfa4, f, skk, self.yy, self.qq, cfl\
                , time_step_optimization, yw, bw, w, option, yn, qn, igate\
                , bed_elevation_path, upstream_path, downstream_path, channel_width_path\
                , output_path, option_dsbc, null = data

        self.I_UPSTREAM = 0
        self.I_DOWNSTREAM = n_sections - 1
        self.manning_m = f

        # Read in bed elevation and bottom width input
        with open(bed_elevation_path, newline='') as file1:
            with open(channel_width_path, newline='') as file2:
                read_data1 = list(csv.reader(file1, delimiter=' '))
                read_data2 = list(csv.reader(file2, delimiter=' '))
                for i in range(n_sections): # use usual convention of i as spatial dimension
                    z = float(read_data1[i][1])
                    y0 = z + self.yy #TODO: This seems like a clunky definition of the initial water surface
                                #      Elevation and I think we can do better.
                    b0 = float(read_data2[i][1])
                    # print(f'b0 {b0}')
                    self.sections.append(self.RectangleSection(station = i
                                         , bottom_z = z
                                         , bottom_width = b0
                                         , dx_ds = dxini
                                         , manning_n_ds = 1/skk)) # dx is a static value in the test cases

        # Read hydrograph input Upstream and Downstream
        with open(upstream_path, newline='') as file3:
            with open(downstream_path, newline='') as file4:
                read_data3 = list(csv.reader(file3, delimiter=' '))
                read_data4 = list(csv.reader(file4, delimiter=' '))
                for j in range(ntim): # use usual convention of j as time dimension
                    self.upstream_flow_ts.append(float(read_data3[j][1]))
                    self.downstream_stage_ts.append(float(read_data4[j][1]))
                    self.time_list.append(dtini * j)

        for i, section in enumerate(self.sections):
            for j, new_time in enumerate(self.time_list):
                section.time_steps.append(self.TimeStep(new_time = new_time
                                    , new_flow = self.qq
                                    , new_depth = self.yy
                                    , new_water_z = section.bottom_z + self.yy))
                    #self.sections[self.I_UPSTREAM].time_steps.append(self.TimeStep(new_flow = q_upstream))
                    #self.sections[self.I_DOWNSTREAM].time_steps.append(self.TimeStep(new_depth = y_downstream))
    def compute_initial_state(self):
        ''' Initial state computed in initialization function.'''
        #print(self.upstream_flow_ts)
        #print(self.downstream_stage_ts)
        # self.add_normal_time_step(0, 0, self.sections, self.downstream_stage_ts[0], self.upstream_flow_ts[0])
        pass

    def compute_next_time_step_state(self, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next):

        # There is the possibility that the Theta for the Predictor and
        # Corrector steps could be numerically different -- but we don't need to
        # take advantage of that here, so theta and theta_s are set to be
        # equal.

        # secpred, apply_corrector, and matrixc should be replaced by
        # section, apply_predictor, and matrixp, respectively.
        # There would need to be a flag of some sort to handle the direction.
        self.compute_sections(self.sections, j_current , j_next
                                         , upstream_flow_current , upstream_flow_next
                                         , downstream_stage_current , downstream_stage_next
                                         , predictor_step = True)
        self.thes = self.thetas
        self.matrix_pc(self.sections, j_current , j_next
                                         , upstream_flow_current , upstream_flow_next
                                         , downstream_stage_current , downstream_stage_next
                                         , predictor_step = True)
        self.apply_predictor(self.sections, j_current , j_next
                                         , upstream_flow_current , upstream_flow_next
                                         , downstream_stage_current , downstream_stage_next)
        self.compute_sections(self.sections, j_current , j_next
                                         , upstream_flow_current , upstream_flow_next
                                         , downstream_stage_current , downstream_stage_next
                                         , predictor_step = False)
        self.thes = self.thesinv
        self.matrix_pc(self.sections, j_current , j_next
                                         , upstream_flow_current , upstream_flow_next
                                         , downstream_stage_current , downstream_stage_next
                                         , predictor_step = False)
        self.apply_corrector(self.sections, j_current , j_next
                                         , upstream_flow_current , upstream_flow_next
                                         , downstream_stage_current , downstream_stage_next)

        #self.add_normal_time_step(j_current, j_next, self.sections, downstream_stage_next, upstream_flow_next)

    # def add_normal_time_step(self, j_current, j_next, sections, downstream_stage_next, upstream_flow_next):
    #     # TODO: Determine if this should be called from the Timestep class
    #
    #     for i, section in enumerate(sections):
    #         # print(f'downstream stage next {downstream_stage_next}')
    #         if i == self.I_DOWNSTREAM: #Add known downstream boundary
    #             section.time_steps.append(self.TimeStep(new_flow = upstream_flow_next
    #                                                     , new_depth = downstream_stage_next))
    #             # print(f'stsd_jcurr {section.time_steps[j_current].depth}')
    #             # print(f'sdstsd_jcurr {section_ds.time_steps[j_current].depth}')
    #             # print(f's0tsd_jcurr {sections[0].time_steps[j_current].depth}')
    #             # stage.append(downstream_stage_next)
    #             # WS.append(section.bottom_z + stage[0])
    #             # A.append(section.get_area_depth(stage[0]))
    #             continue
    #
    #         section_ds = sections[i-1]
    #         #section
    #         # print(f'jnext {j_next}')
    #         # print(f'jcurr {j_current}')
    #         # print(f'stsd_jcurr {section.time_steps[j_current].depth}')
    #         # print(f'sdstsd_jcurr {section_ds.time_steps[j_current].depth}')
    #         # print(f's0tsd_jcurr {sections[0].time_steps[j_current].depth}')
    #         # print(f'stsd_jnext {section.time_steps[j_next].depth}')
    #         # print(f'sdstsd_jnext {section_ds.time_steps[j_next].depth}')
    #         # print(f's0tsd_jnext {sections[0].time_steps[j_next].depth}')
    #         #Use normal depth as an seed estimate convergence solution for standard depth
    #         y_guess = helpers.y_direct(section.bottom_width, section.manning_n_ds, section.bed_slope_ds, upstream_flow_next)
    #         section.time_steps.append(self.TimeStep(new_flow = upstream_flow_next, new_depth = y_guess))
    #         section.time_steps[0].water_z = section.time_steps[0].depth + section.bottom_z
    #         # print(section.time_steps[0].depth, section.time_steps[0].water_z)
    #
    def apply_predictor(self
            , section_arr
            , j_current
            , j_next
            , upstream_flow_current
            , upstream_flow_next
            , downstream_stage_current
            , downstream_stage_next):

        dt = self.time_list[j_next] - self.time_list[j_current]
        for i, section in enumerate(section_arr):
            # Use the flow time series for the upstream boundary
            section_j = section.time_steps[j_current]
            section_j.sigma_ds = dt / section.dx_ds
            if i == self.I_UPSTREAM:
                #Set the upstream boundary flow
                #This is set as a delta -- a correction factor -- to work in the
                #predictor sweep.
                #So we set the Delta Q predictor, which is set from the input
                #time series for the upstream-most point.
                # TODO: Confirm that DAP(1) is never used... Why not?
                # section_j.delta_area_predictor = 0.0
                section_j.delta_flow_predictor = upstream_flow_next - upstream_flow_current
            elif i == self.I_DOWNSTREAM:
                section_j.delta_flow_predictor = upstream_flow_next - upstream_flow_current
                section_j.delta_area_predictor = 0.0
                #TODO: Ask Ehab why these values are touched in the predictor step
                section_j.delta_area_corrector = 0.0
                section_j.delta_flow_corrector = section_j.delta_flow_predictor
            else:
                section_US = section_arr[i-1]
                section_US_j = section_US.time_steps[j_current]
                section_j.rhs1 = -1.0 * section_j.sigma_ds\
                        * (section_j.f1 - section_US_j.f1\
                            - section_j.d1 + section_US_j.d1)
                section_j.rhs2 = -1.0 * section_j.sigma_ds\
                        * (section_j.f2 - section_US_j.f2\
                            - section_j.d2 + section_US_j.d2)\
                        + dt * self.gravity\
                             * (section_j.ci2_ds + section_j.as0_ds)
                c11 = section_j.g11inv * section_US_j.b11 + section_j.g12inv * section_US_j.b21
                c12 = section_j.g11inv * section_US_j.b12 + section_j.g12inv * section_US_j.b22
                c21 = section_j.g21inv * section_US_j.b11 + section_j.g22inv * section_US_j.b21
                c22 = section_j.g21inv * section_US_j.b12 + section_j.g22inv * section_US_j.b22
                section_j.delta_area_predictor =\
                    section_j.g11inv * section_j.rhs1\
                    + section_j.g12inv * section_j.rhs2\
                    - c11 * section_US_j.delta_area_predictor\
                    - c12 * section_US_j.delta_flow_predictor
                section_j.delta_flow_predictor =\
                    section_j.g21inv * section_j.rhs1\
                    + section_j.g22inv * section_j.rhs2\
                    - c21 * section_US_j.delta_area_predictor\
                    - c22 * section_US_j.delta_flow_predictor

            #Update via predictor
            self.debug = True
            if self.debug and (i == self.I_UPSTREAM or i == self.I_DOWNSTREAM): print (f'area {section_j.flow_area}, dap {section_j.delta_area_predictor}')
            section_j.areap = section_j.flow_area + section_j.delta_area_predictor
            if self.debug and (i == self.I_UPSTREAM or i == self.I_DOWNSTREAM): print (f'flow {section_j.flow_area}, dqp {section_j.delta_flow_predictor}')
            section_j.qp = section_j.flow + section_j.delta_flow_predictor
            self.debug = False

    def dsbc():
        '''
          c
          c----------------------------------------------------------------------
          c     Downstream boundary condition
                subroutine dsbc(n)
          c----------------------------------------------------------------------
          c
                parameter(grav=9.81)
                common/arrays/ flow_area(1000),y(100000,1000),q(100000,1000),bo(1000),
              1    areap(1000),qp(1000),z(1000),delta_flow_predictor(1000),av11(1000),av12(1000),
              1                                       av21(1000),av22(1000),
              2              delta_flow_corrector(1000),delta_area_predictor(1000),delta_area_corrector(1000),ci1(1000),ci2(1000),
              3           aS0(1000),d(1000),f1(1000),g11inv(1000),g12inv(1000),
              4g21inv(1000),g22inv(1000),f2(1000),b11(1000),b12(1000),b21(1000),
              5b22(1000),eps2(1000),eps4(1000),d1(1000),d2(1000),u(1000),c(1000),
              6                        sk(1000),co(1000),gS0(1000),dbdx(1000),
              7                      dt(1000),ityp(1000),dx(1000)
                common/var/n_sections,cfl,f,qus,yds,rhs1,rhs2,c11,c12,c21,c22,us,time_step_optimization,
              1           thes,vv,dtini
                common/matrix/phi,theta,alfa2,alfa4,w
                common/sgate/ag,bg,option,yn,eps,dmeu,yw,bw,ydsn
          c
          c compute conjugate depth at dwon stream end
                ads=(y(j,I_UPSTREAM)-z(I_UPSTREAM))*bo(I_UPSTREAM)
                frds=q(j,I_UPSTREAM)/sqrt(grav*ads**3.0/bo(I_UPSTREAM))
                yconj=0.5*(y(j,I_UPSTREAM)-z(I_UPSTREAM))*(sqrt(1.0+8.0*frds*frds)-1.0)
                yconj=yconj+z(I_UPSTREAM)
                write(*,*)yn,y(j,I_UPSTREAM),yconj,z(I_UPSTREAM)
                if(yconj.lt.yn) then
          c      write(*,*)'no'
                if(option.eq.1.0) then
          c
          c downstream water level imposed (option 1)
                    delta_area_corrector(I_UPSTREAM)=(yn-y(j,I_UPSTREAM))*bo(I_UPSTREAM)
                    delta_area_predictor(I_UPSTREAM)=(yn-y(j,I_UPSTREAM))*bo(I_UPSTREAM)
                    delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
          c          write(*,*)'yes',delta_area_corrector(I_UPSTREAM),yn,y(j,I_UPSTREAM),delta_flow_corrector(I_UPSTREAM)
                  elseif(option.eq.2.0) then
          c
          c downstream flow imposed  (option 2)
                    delta_flow_corrector(I_UPSTREAM)=0.0
                    delta_flow_predictor(I_UPSTREAM)=0.0
                    delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
                  elseif(option.eq.3.0) then
          c
          c downstream rating curve imposed (option 3)
                    delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
                    yn=(flow_area(I_UPSTREAM)+delta_area_predictor(I_UPSTREAM))/bo(I_UPSTREAM)
                    qn=0.65*10*1.0*sqrt(2.0*grav*(yn-0.5))
                    delta_flow_predictor(I_UPSTREAM)=qn-q(j,I_UPSTREAM)
                    delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
                endif
                else
          c
          c super critical flow exist at downstream exit
                  delta_area_corrector(I_UPSTREAM)=delta_area_predictor(I_UPSTREAM)
                  delta_flow_corrector(I_UPSTREAM)=delta_flow_predictor(I_UPSTREAM)
                endif
                return
                end
          c----------------------------------------------------------------------
        '''
        pass

    def compute_sections(self
            , section_arr
            , j_current
            , j_next
            , upstream_flow_current
            , upstream_flow_next
            , downstream_stage_current
            , downstream_stage_next
            , predictor_step = False):

        for i, section in enumerate(section_arr):
            section_j = section.time_steps[j_current]
            if predictor_step: # If we are on the second step, applying the predictors, use the areap and qp
                if i == self.I_DOWNSTREAM:
                    section_j.depth = downstream_stage_current
                else:
                    section_j.depth = section_j.water_z - section.bottom_z
                section_j.flow_area = section.get_area_depth(section_j.depth)
                area = section_j.flow_area
                flow = section_j.flow
            else:
                area = section_j.areap # computed via 'apply_predictor' method
                flow = section_j.qp # computed via 'apply_predictor' method
                section_j.depthp = section.get_depth_area(area)

            # print (f'flow_area {section_j.flow_area}')
            # print (f'i, z: {i} {section_j.water_z}')
            # print(f'depth {section_j.depth}')
            # print(f'bw {section.bottom_width}')
            # print(f'current stage {downstream_stage_current}')
            # print(f'next stage {downstream_stage_next}')
            #TODO: make ci1 compatible with generalized sections.
            section_j.ci1 = section.bottom_width * (section_j.depth ** 2.0) / 2.0
            #TODO: Is there a reason for the difference between this Pw calculation and the matrixc calculations?
            section_j.hy = area / section.get_wetted_perimeter_depth(section_j.depth)
            section_j.conveyance_ds = self.manning_m * area * section_j.hy ** (2.0/3.0)
            # print (f'c0 {section_j.c0}')
            # print (f'manning_m {self.manning_m}')
            # print (f'flow_area {section_j.flow_area}')
            # print (f'Rw {section_j.hy}')
            if predictor_step:
                if i == self.I_UPSTREAM: continue
                section_US = section_arr[i-1] #Upstream is increasing in index for our scheme
                section_US_j = section_US.time_steps[j_current]
                    #TODO: make ci2 compatible with generalized sections.
                section_US_j.ci2_ds = ((section_j.depth ** 2.0) * (section_US_j.depth ** 2.0)) \
                              * (section.bottom_width - section_US.bottom_width) \
                              / (section_US.dx_ds * 4.0)
                # TODO: STRIKE THE BED SLOPE CALCULATION BECAUSE IT IS A STATIC VALUE
                # print('Upstream n, flow, conveyance_ds {} {} {}', \
                            # section_US.manning_n_ds\
                            # , section_US_j.flow\
                            # , section_US_j.conveyance_ds)
                # print('This section n, flow, conveyance_ds {} {} {}', \
                            # section_US.manning_n_ds\
                            # , section_j.flow\
                            # , section_j.conveyance_ds)
                #TODO: DongHa  Determine if we need to make the 0.5 to compute fs a variable
                section_US_j.friction_slope_ds = section_US.manning_n_ds \
                               * 0.5 * section_US_j.flow \
                               * abs(section_US_j.flow) / (section_US_j.conveyance_ds ** 2.0) \
                               + section_US.manning_n_ds * 0.5 * section_j.flow \
                               * abs(section_j.flow) / (section_j.conveyance_ds ** 2.0)
                section_US_j.as0_ds = (section_US_j.flow_area + section_j.flow_area) \
                               / 2.0 * (section_US.bed_slope_ds \
                                        - section_US_j.friction_slope_ds)
                section_US_j.gs0_ds = self.gravity * (section_US.bed_slope_ds \
                                        - section_US_j.friction_slope_ds)
                section_US.dbdx_ds = (section.bottom_width - section_US.bottom_width) \
                                      / section_US.dx_ds
            if not predictor_step:
                if i == self.I_DOWNSTREAM: continue
                section_DS = section_arr[i+1]
                section_DS_j = section_DS.time_steps[j_current]
                    #TODO: not compatible with generalized sections.
                section_j.ci2_ds = ((section_DS_j.depth ** 2.0) * (section_j.depth ** 2.0)) \
                              * (section_DS.bottom_width - section.bottom_width) \
                              / (section.dx_ds * 4.0)
                #TODO: remove bed_slope_ds calculation from the timestep loop: it does not change
                section.bed_slope_ds = (section.bottom_z - section_DS.bottom_z) \
                                        / section.dx_ds
                section_j.friction_slope_ds = section.manning_n_ds \
                               * 0.5 * section_j.flow \
                               * abs(section_j.flow) / (section_j.conveyance_ds ** 2.0) \
                               + section.manning_n_ds * 0.5 * section_DS_j.flow \
                               * abs(section_DS_j.flow) / (section_DS_j.conveyance_ds ** 2.0)
                section_j.as0_ds = (section_j.flow_area + section_DS_j.flow_area) \
                               / 2.0 * (section.bed_slope_ds \
                                        - section_j.friction_slope_ds)
                section_j.gs0_ds = self.gravity * (section.bed_slope_ds \
                                        - section_j.friction_slope_ds)
                #TODO: remove dbdx_ds calculation from the timestep loop: it does not change
                section.dbdx_ds = (section_DS.bottom_width - section.bottom_width) \
                                      / section.dx_ds

    def matrix_pc(self
            , section_arr
            , j_current
            , j_next
            , upstream_flow_current
            , upstream_flow_next
            , downstream_stage_current
            , downstream_stage_next
            , predictor_step = False):

        #crmax = 0.0
        #crmin = 100.0
        for i, section in enumerate(section_arr):
            section_j = section.time_steps[j_current]
            if self.predictor_step: # If we are on the second step, applying the predictors, use the areap and qp
                area = section_j.flow_area
                flow = section_j.flow
            else:
                area = section_j.areap
                flow = section_j.qp
            # print(f'area, flow {area}, {flow}')
            section_j.velocity = flow / area
            section_j.celerity = section.get_celerity_area(area, self.gravity, debug = self.debug)
            # print(section_j.celerity)
            if section_j.velocity == section_j.celerity:
                section_j.celerity = section_j.celerity + self.celerity_epsilon
            # print(section_j.celerity_epsilon)

            #c     This is the matrix L (left eigenvector matrix - eq 13)
            e11 = 1.0
            e12 = -1.0 / (section_j.velocity - section_j.celerity)
            e21 = 1.0
            e22 = -1.0 / (section_j.velocity + section_j.celerity)

            #c       L^{-1} (inverse of Left eigenvector matrix)
            f11 = -(section_j.velocity - section_j.celerity) / (2.0 * section_j.celerity)
            f12 = (section_j.velocity + section_j.celerity) / (2.0 * section_j.celerity)
            f21 = -(section_j.velocity ** 2.0 - section_j.celerity ** 2.0) / (2.0 * section_j.celerity)
            f22 = (section_j.velocity ** 2.0 - section_j.celerity ** 2.0) / (2.0 * section_j.celerity)

            #c       Diagonal wave matrix D (eq 12)
            section_j.d11 = abs(section_j.velocity + section_j.celerity)
            section_j.d22 = abs(section_j.velocity - section_j.celerity)

            #c       Equation 11 (L^{-1} D L)
            a11 = e11 * f11 * section_j.d11 + e21 * f12 * section_j.d22
            a12 = e12 * f11 * section_j.d11 + e22 * f12 * section_j.d22
            a21 = e11 * f21 * section_j.d11 + e21 * f22 * section_j.d22
            a22 = e12 * f21 * section_j.d11 + e22 * f22 * section_j.d22

            dt = self.time_list[j_next] - self.time_list[j_current]

            #c     Calculating dK/dA (eq 15)
            dkda = section.get_dkda_area(area)

            #c     Matrix S (eq 14)
            st11 = 0.0
            st12 = 0.0
            st21 = section.get_st21_area (area, flow, section_j.gs0_ds
                                   , section_j.conveyance_ds, dkda
                                   , self.gravity, self .manning_m)
            # TODO: Determine if st22 is a section-property-dependent value
            # and if so, move it to a method within the RectangleSection class
            st22 = -2.0 * self.manning_m * area \
                 * area \
                 / section_j.conveyance_ds

            section_US = section_arr[i-1]
            if predictor_step:
                # TODO: Confirm the order of these statements with Ehab
                # The following if statements run in reverse order in matrixp vs. matrixc
                # TODO: combine the two inner statements if the order remains important
                # TODO: Combine all three statments if the order is not important
                thessign = -1.0
                if section.dx_ds < self.dx_tolerance:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section.dx_ds
                    self.crmax = max(self.crmax , section_j.sigma_ds * max(section_j.d11, section_j.d22))
                    self.crmin = min(self.crmin , section_j.sigma_ds * max(section_j.d11, section_j.d22))

                # c     LHS of eq 7
                section_j.b11 = (0.5 - self.phi - self.theta  * section_j.sigma_ds * a11
                                    + thessign * 0.5 * self.thetas * st11 * dt)
                section_j.b12 = (-self.theta  * section_j.sigma_ds * a12
                                    + thessign * 0.5 * self.thetas * st12 * dt)
                section_j.b21 = (-self.theta * section_j.sigma_ds * a21
                                    + thessign * 0.5 * self.thetas * st21 * dt)
                section_j.b22 = (0.5 - self.phi - self.theta  * section_j.sigma_ds * a22
                                    + thessign * 0.5 * self.thetas * st22 * dt)

                if i == self.I_UPSTREAM:
                    section_j.sigma_ds = dt
                elif section_US.dx_ds < self.dx_tolerance:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section_US.dx_ds

            elif not predictor_step:
                thessign = 1.0
                # TODO: Confirm that these non-parallel if statements are congruent
                if i == self.I_UPSTREAM:
                    section_j.sigma_ds = dt
                elif section.dx_ds < self.dx_tolerance:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section_US.dx_ds

                # c     LHS of eq 7
                section_j.b11 = (0.5 - self.phi - self.theta  * section_j.sigma_ds * a11
                                    + thessign * 0.5 * self.thetas * st11 * dt)
                section_j.b12 = (-self.theta  * section_j.sigma_ds * a12
                                    + thessign * 0.5 * self.thetas * st12 * dt)
                section_j.b21 = (-self.theta * section_j.sigma_ds * a21
                                    + thessign * 0.5 * self.thetas * st21 * dt)
                section_j.b22 = (0.5 - self.phi - self.theta  * section_j.sigma_ds * a22
                                    + thessign * 0.5 * self.thetas * st22 * dt)

                if section.dx_ds == 0.0:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section.dx_ds
                    self.crmax = max(self.crmax , section_j.sigma_ds * max(section_j.d11, section_j.d22))
                    self.crmin = min(self.crmin , section_j.sigma_ds * max(section_j.d11, section_j.d22))

            if predictor_step:
                #TODO: ASK Ehab if 'dt' missing from matrixc computation is intentional
                #If it is an omission, these two equations can be combined...
                thessign = -1.0
                g11 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a11
                              + thessign * 0.5 * self.thetas * st11 * dt)
                g12 = (self.theta * section_j.sigma_ds * a12 + thessign * 0.5
                              * self.thetas * st12 * dt)
                g21 = (self.theta * section_j.sigma_ds * a21 + thessign * 0.5
                              * self.thetas * st21 * dt)
                g22 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a22
                              + thessign * 0.5 * self.thetas * st22 * dt)
            elif not predictor_step: # This reads more clearly that the equivalent simple 'else'
                thessign = 1.0
                g11 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a11
                              + thessign * 0.5 * self.thetas * st11)
                g12 = (self.theta * section_j.sigma_ds * a12 + thessign * 0.5
                              * self.thetas * st12)
                g21 = (self.theta * section_j.sigma_ds * a21 + thessign * 0.5
                              * self.thetas * st21)
                g22 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a22
                              + thessign * 0.5 * self.thetas * st22)

            section_j.g11inv =  g22 / (g11 * g22 - g12 * g21)
            section_j.g12inv = -g12 / (g11 * g22 - g12 * g21)
            section_j.g21inv = -g21 / (g11 * g22 - g12 * g21)
            section_j.g22inv =  g11 / (g11 * g22 - g12 * g21)

            section_j.f1 = flow
            section_j.f2 = flow ** 2.0 \
                        / area \
                        + self.gravity * section_j.ci1

          #      if(i.ge.2.and.i.lt.n_sections) then
            if i == self.I_DOWNSTREAM:
                section_US = section_arr[i-1]
                section_US_j = section_US.time_steps[j_current]
                section_j.eps2 = section_US_j.eps2
            elif i == self.I_UPSTREAM:
                section_DS = section_arr[i+1]
                section_DS_j = section_DS.time_steps[j_current]
                section_j.eps2 = section_DS_j.eps2
            else: # if i > self.I_DOWNSTREAM and i < self.I_UPSTREAM:
                section_DS = section_arr[i+1]
                section_DS_j = section_DS.time_steps[j_current]
                section_US = section_arr[i-1]
                section_US_j = section_US.time_steps[j_current]
                if self.predictor_step: # If we are on the second step, applying the predictors, use the areap and qp
                    area_DS = section_DS_j.flow_area
                    flow_DS = section_DS_j.flow
                    area_US = section_US_j.flow_area
                    flow_US = section_US_j.flow
                else:
                    area_DS = section_DS_j.areap
                    flow_DS = section_DS_j.qp
                    area_US = section_US_j.areap
                    flow_US = section_US_j.qp

                dip1 = section_DS.get_depth_area(area_DS)
                di = 2 * section.get_depth_area(area)
                dim1 = section_US.get_depth_area(area_US)
                if abs(dip1 + di + dim1) < self.depth_tolerance:
                    section_j.eps2 = self.depth_tolerance
                else:
                    section_j.eps2 = self.alfa2 * abs(dip1 - di + dim1) / (dip1 + di + dim1)

          #c
          #TODO: WHAT IS THIS STEP --- something to do with the boundary?
          #      do 20 i=2,n_sections-1
          #        if(ityp(i).ne.1) then
          #          eps2(i)=eps2(i-1)
          #          eps2(i+1)=eps2(i+2)
          #        endif
          #20    continue
          #c

        if predictor_step:
            for i, section in enumerate(section_arr):
                if i == self.I_UPSTREAM or i == self.I_DOWNSTREAM:
                    continue
                section_j = section.time_steps[j_current]
                section_DS = section_arr[i+1]
                section_DS_j = section_DS.time_steps[j_current]
                section_j.eps2 = max(section_DS_j.eps2, section_j.eps2)
                # print('0.0, alfa4, eps2, velocity, celerity {} {} {} {} {}'.format(0.0, self.alfa4, section_j.eps2
                #                 , section_j.velocity, section_j.celerity))
                section_j.eps4 = max(0.0, self.alfa4 - section_j.eps2
                                / (section_j.velocity + section_j.celerity))
              #      eps4(i+1)=max(0.,alfa4-eps2(i+1)/(u(i+1)+c(i+1)))
              #c      write(*,*)'corr',i,eps2(i)

        elif not predictor_step:
            for i, section in enumerate(reversed(section_arr)):
                if i == self.I_UPSTREAM or i == self.I_DOWNSTREAM:
                    continue
                section_j = section.time_steps[j_current]
                section_DS = section_arr[i+1]
                section_DS_j = section_DS.time_steps[j_current]
                section_DS_j.eps2 = max(section_DS_j.eps2, section_j.eps2)
                # print('0.0, alfa, eps2, velocity, celerity {} {} {} {} {}'.format(0.0, self.alfa4, section_j.eps2
                #                 , section_j.velocity, section_j.celerity))
                section_DS_j.eps4 = max(0.0, self.alfa4 - section_DS_j.eps2
                                / (section_DS_j.velocity + section_DS_j.celerity))
          #c      write(*,*)'corr',i,eps2(i)

          #      d1(1)=0.0
          #      d2(1)=0.0
          #      d1(I_UPSTREAM)=0.0
          #      d2(I_UPSTREAM)=0.0

        for i, section in enumerate(section_arr):
            section_j = section.time_steps[j_current] #TODO: Implement direct iterator for timestep
            if i == self.I_UPSTREAM:
                section_j.d1 = 0.0
                section_j.d2 = 0.0
            elif i == self.I_DOWNSTREAM:
                section_j.d1 = 0.0
                section_j.d2 = 0.0
            else:
                if predictor_step: # If we are on the second step, applying the predictors, use the areap and qp
                    section_DS = section_arr[i+1]
                    section_DS_j = section_DS.time_steps[j_current]
                    area_DS = section_DS_j.flow_area
                    flow_DS = section_DS_j.flow
                else:
                    section_US = section_arr[i-1]
                    section_US_j = section_US.time_steps[j_current]
                    area_US = section_US_j.areap
                    flow_US = section_US_j.qp
                ei = max(abs(section_j.velocity + section_j.celerity),
                                    (section_j.velocity - section_j.celerity))
                if predictor_step:
                    ei1 = max(abs(section_DS_j.velocity + section_DS_j.celerity),
                                    (section_DS_j.velocity - section_DS_j.celerity))
                elif not predictor_step:
                    ei1 = max(abs(section_US_j.velocity + section_US_j.celerity),
                                    (section_US_j.velocity - section_US_j.celerity))
                eia = (ei + ei1) / 2.0
              #      if(ityp(i-1).ne.1) then
              #        d1(i)=0.0
              #        d2(i)=0.0
              #      elseif(i.eq.2.or.i.eq.(n_sections-1)) then
                if predictor_step: # This if could technically be combined with the above if statement
                                   # but it reads more clearly to have it here.
                    section_US = section_arr[i-1]
                    section_US_j = section_US.time_steps[j_current]
                    area_US = section_DS_j.flow_area
                    flow_US = section_DS_j.flow
                    area_difference = area_DS - area
                    flow_difference = flow_DS - flow
                elif not predictor_step:
                    section_DS = section_arr[i+1]
                    section_DS_j = section_DS.time_steps[j_current]
                    area_DS = section_US_j.areap
                    flow_DS = section_US_j.qp
                    area_difference = area - area_US
                    flow_difference = flow - flow_US
                section_j.d1 = section_j.eps2 * eia * area_difference
                section_j.d2 = section_j.eps2 * eia * flow_difference
                if (i > self.I_UPSTREAM + 1 and i < self.I_DOWNSTREAM - 1):
                    # print(i, self.I_DOWNSTREAM, self.I_UPSTREAM)
              #        d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))
              #        d2(i)=eps2(i)*eia*(qp(i)-qp(i-1)) #TODO: Determine if this reduction is unnecessary (i.e., Could this be replaced by simply using the next lines and letting the equation reduce automatically based on the value in eps4?
                    if predictor_step:
                        section_2DS = section_arr[i+2]
                        section_2DS_j = section_DS.time_steps[j_current]
                        area_2DS = section_2DS_j.flow_area
                        flow_2DS = section_2DS_j.flow
                        d1_eps4_diff = - 1.0 * section_j.eps4 * (area_2DS - 3.0 * area_DS + 3.0 * area - area_US)
                        d2_eps4_diff = - 1.0 * section_j.eps4 * (flow_2DS - 3.0 * flow_DS + 3.0 * flow - flow_US)
                    elif not predictor_step:
                        section_2US = section_arr[i-2]
                        section_2US_j = section_2US.time_steps[j_current]
                        area_2US = section_2US_j.areap
                        flow_2US = section_2US_j.qp
                        d1_eps4_diff = - 1.0 * section_j.eps4 * (area_DS - 3.0 * area + 3.0 * area_US - area_2US)
                        d2_eps4_diff = - 1.0 * section_j.eps4 * (flow_DS - 3.0 * flow + 3.0 * flow_US - flow_2US)
                    section_j.d1 = section_j.d1 + d1_eps4_diff
                    section_j.d2 = section_j.d1 + d2_eps4_diff
              #      else
              #      d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))-eps4(i)*
              #    1        (areap(i+1)-3*areap(i)+3*areap(i-1)-areap(i-2))
              #      d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))-eps4(i)*(qp(i+1)-3*qp(i)+
              #    1      3*qp(i-1)-qp(i-2))
              #      endif
              #c      write(*,*)'corr',i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)

    def apply_corrector(self
            , section_arr
            , j_current
            , j_next
            , upstream_flow_current
            , upstream_flow_next
            , downstream_stage_current
            , downstream_stage_next):
        '''docstring here'''
        dt = self.time_list[j_next] - self.time_list[j_current]
        for i, section in enumerate(section_arr):
        # do i=ncomp-1,1,-1
        #TODO: CONFIRM That the reversed array is necessary in this case (I don't think it is...)
        #     cour=dt(i)/dx(i)
            # Use the flow time series for the upstream boundary
            section_j = section.time_steps[j_current]
            if i == self.I_UPSTREAM or i == self.I_DOWNSTREAM:
                #Set the upstream boundary flow
                #This is set as a delta -- a correction factor -- to work in the
                #predictor sweep.
                #So we set the Delta Q predictor, which is set from the input
                #time series for the upstream-most point.
                section_j.delta_flow_corrector = section_j.delta_flow_predictor
                section_j.delta_area_corrector = 0.0
            else:
                section_DS = section_arr[i+1]
                section_DS_j = section_DS.time_steps[j_current]
                section_j.sigma_ds = dt / section.dx_ds
                section_j.rhs1 = -1.0 * section_j.sigma_ds\
                        * (section_DS_j.f1 - section_j.f1\
                            - section_DS_j.d1 + section_j.d1)
            #     rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))
                section_j.rhs2 = -1.0 * section_j.sigma_ds\
                        * (section_DS_j.f2 - section_j.f2\
                            - section_DS_j.d2 + section_j.d2)\
                        + dt * self.gravity\
                             * (section_j.ci2_ds + section_j.as0_ds)
            #     rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i))
                c11 = section_j.g11inv * section_DS_j.b11 + section_j.g12inv * section_DS_j.b21
                c12 = section_j.g11inv * section_DS_j.b12 + section_j.g12inv * section_DS_j.b22
                c21 = section_j.g21inv * section_DS_j.b11 + section_j.g22inv * section_DS_j.b21
                c22 = section_j.g21inv * section_DS_j.b12 + section_j.g22inv * section_DS_j.b22
            #     c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
            #     c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
            #     c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
            #     c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
                section_j.delta_area_predictor =\
                    section_j.g11inv * section_j.rhs1\
                    + section_j.g12inv * section_j.rhs2\
                    - c11 * section_DS_j.delta_area_corrector\
                    - c12 * section_DS_j.delta_flow_corrector
            #     dac(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dac(i+1)-c12*dqc(i+1)
                section_j.delta_flow_predictor =\
                    section_j.g21inv * section_j.rhs1\
                    + section_j.g22inv * section_j.rhs2\
                    - c21 * section_DS_j.delta_area_corrector\
                    - c22 * section_DS_j.delta_flow_corrector
            #     dqc(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dac(i+1)-c22*dqc(i+1)

            # end do

    class RectangleSection(Network.RectangleSection):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def get_dkda_area(self, area):
            wetted_perimeter = self.get_wetted_perimeter_area(area)
            return (1 / self.manning_n_ds *
                  (( 5.0 / 3.0 * area ** (2.0/3.0) * wetted_perimeter) -
                  (area ** (5.0/3.0) * 2.0 / self.bottom_width)) /
                  (wetted_perimeter ** 2.0))

        def get_st21_area(self, area, flow, gs0, conveyance, dkda, gravity, manning_m):
            return (gravity * area / self.bottom_width \
                / self.bottom_width * self.dbdx_ds + gs0 \
                + manning_m * 2.0 * gravity \
                * area * flow \
                * abs (flow) / conveyance ** 3.0 \
                * dkda)

    class TimeStep(Network.TimeStep):
        '''MESH-specific time-step values'''
        def __init__(self, new_water_z = 0.0, *args, **kwargs):
            super().__init__(*args, **kwargs)

            # Per-time-step at-a-section properties
            self.delta_flow_corrector = 0
            self.delta_flow_predictor = 0
            self.delta_area_corrector = 0
            self.delta_area_predictor = 0
            self.water_z = new_water_z
            self.flow_area = 0
            self.areap = 10
            self.qp = 10
            self.ci1 = 0
            self.hy = 0 # Hydraulic Radius (used to compute co)

            # Per-time-step downstream reach properties
            self.conveyance_ds = 0
            self.ci2_ds = 0
            # self.friction_slope_ds = 0 # Derived from parent Class, Network.TimeStep
            self.as0_ds = 0
            self.gs0_ds = 0
            self.sigma_ds = 0 # Sigma is related to the courant parameter: CFL = celerity * sigma
            #self.cour = 0
            self.dbdx = 0
            self.velocity = 0
            self.celerity = 0
            self.f1 = 0
            self.f2 = 0
            self.d1 = 0
            self.d2 = 0
            self.d11 = 0
            self.d22 = 0
            self.b11 = 0
            self.b12 = 0
            self.b21 = 0
            self.b22 = 0
            self.g11inv = 0
            self.g12inv = 0
            self.g21inv = 0
            self.g22inv = 0
            self.eps2 = 0
            self.eps4 = 0
            self.rhs1 = 0
            self.rhs2 = 0

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
