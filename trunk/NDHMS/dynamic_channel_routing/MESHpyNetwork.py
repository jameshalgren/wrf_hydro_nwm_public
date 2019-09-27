# import required modules
from __future__ import division
import helpers
import constants
import meshconstants
from network import Network
import csv
import os

'''This python module implements  (and to a degree, extends) the
   MESH code first developed by Dr. Ehab Meselhe.
   A number of helpful references are available to explain the
   derivation of the method and may be requested by contacting
   the authors or Dr. Meselhe.
   We gratefully acknowlege his contribution to this work.'''
class MESHpyNetwork(Network):
    '''USE Global Declarations here to manage these values'''
    # TODO Determine why these values do not persist when declared in __init__
    debug = False
    dx_tolerance = meshconstants.DX_TOLERANCE
    depth_tolerance = meshconstants.DEPTH_TOLERANCE
    celerity_epsilon = meshconstants.CELERITY_EPSILON
    area_tolerance = meshconstants.AREA_TOLERANCE
    crmax = 0.0
    crmin = 100.0
    predictor_step = True
    yy = 0.0
    qq = 0.0
    phi = meshconstants.PHI         # source term treatment (0:explicit, 1:implicit)
    theta = meshconstants.THETA     # ?
    thetas = meshconstants.THETAS   # ?
    thesinv = meshconstants.THESINV # ?
    alfa2 = meshconstants.ALFA2     # emp parameter for artificial diffusion (lit)
    alfa4 = meshconstants.ALFA4     # maximum value for artificial diffusion

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def input_and_initialize_meshpyfile(self, filetype = None, input_path=None):
        with open(input_path, newline='') as f:

            # Put the first chunk of each line into a lsit
            data = list(map(lambda x:x.split(' ')[0] , f.read().split('\n')))
            # print(data)

            #TODO: Get rid of this kludge to cast the input numbers into the right datatype
            for i, item in enumerate (data[0:23]):
                data[i] = float(item)

            data[3] = int(data[3])
            data[4] = int(data[4])

            # Assign all the input values into the variables
            dtini, dxini, tfin, n_sections, ntim, self.phi, self.theta, self.thetas\
                , self.thetasinv, self.alfa2, self.alfa4, f, skk, self.yy, self.qq, cfl\
                , time_step_optimization, yw, bw, w, option, yn, qn, igate\
                , bed_elevation_path, upstream_path, downstream_path, channel_width_path\
                , output_path, option_dsbc, null = data

            print(data)

        self.I_UPSTREAM = 0
        self.I_DOWNSTREAM = n_sections - 1
        self.manning_m = f

        # print(f'yy {self.yy}')

        root = os.path.dirname(input_path)
        # print(root)
        # Read in bed elevation and bottom width input
        with open(os.path.join(root, bed_elevation_path), newline='') as file1:
            with open(os.path.join(root,channel_width_path), newline='') as file2:
                read_data1 = list(csv.reader(file1, delimiter=' '))
                read_data2 = list(csv.reader(file2, delimiter=' '))
                for i in range(n_sections): # use usual convention of i as spatial dimension
                    z = float(read_data1[i][1])
                    b0 = float(read_data2[i][1])
                    # print(f'i: {i} -- b0 {b0}, z {z}')
                    self.sections.append(self.RectangleSection(station = i
                                         , bottom_z = z
                                         , bottom_width = b0
                                         , dx_ds = dxini
                                         , manning_n_ds = 1/skk)) # dx is a static value in the test cases
                us_section = None
                for section in self.sections:
                    if us_section:
                        us_section.ds_section = section
                    section.us_section = us_section
                    us_section = section
        # Read hydrograph input Upstream and Downstream
        with open(os.path.join(root, upstream_path), newline='') as file3:
            with open(os.path.join(root, downstream_path), newline='') as file4:
                read_data3 = list(csv.reader(file3, delimiter=' '))
                read_data4 = list(csv.reader(file4, delimiter=' '))
                for j in range(ntim): # use usual convention of j as time dimension
                    self.upstream_flow_ts.append(float(read_data3[j][1]))
                    self.downstream_stage_ts.append(float(read_data4[j][1]))
                    self.time_list.append(dtini * j)

    def compute_initial_state(self, verbose = False
                                    , write_output = False
                                    , output_path = None):
        ''' Initial state computed in initialization function.'''
        q_init = self.qq
        y_init = self.yy
        j_init = 0
        t_init = self.time_list[j_init]
        for i, section in enumerate(self.sections):
            # print(f'yy, qq {self.yy} {self.qq}')
            if i == self.I_UPSTREAM:
                q_init = self.upstream_flow_ts[j_init]
            elif i == self.I_DOWNSTREAM:
                y_init = self.downstream_stage_ts[j_init]
            section.time_steps.append(self.TimeStep(new_time = t_init
                                , new_flow = q_init
                                , new_depth = y_init
                                , new_water_z = section.bottom_z + y_init
                                , new_area = section.get_area_depth(y_init)))
            if section.ds_section:
                section.bed_slope_ds = (section.bottom_z - \
                                section.ds_section.bottom_z) / section.dx_ds
            else:
                section.bed_slope_ds = .00001
        #print(self.upstream_flow_ts)
        #print(self.downstream_stage_ts)

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

        # print(f'yy {self.yy}')
        self.compute_sections(section_arr = self.sections
                        , j_current = j_current
                        , j_next = j_next
                        , upstream_flow_current = upstream_flow_current
                        , upstream_flow_next = upstream_flow_next
                        , downstream_stage_current = downstream_stage_current
                        , downstream_stage_next = downstream_stage_next
                        , predictor_step = True)
        #TODO: Determine the need for thes AND thetas
        # print(f'self.thes; self.thetas: {self.thes} {self.thetas}')
        # self.thes = self.thetas
        self.matrix_pc(section_arr = self.sections
                        , j_current = j_current
                        , j_next = j_next
                        , upstream_flow_current = upstream_flow_current
                        , upstream_flow_next = upstream_flow_next
                        , downstream_stage_current = downstream_stage_current
                        , downstream_stage_next = downstream_stage_next
                        , predictor_step = True)
        self.apply_predictor(section_arr = self.sections
                        , j_current = j_current
                        , j_next = j_next
                        , upstream_flow_current = upstream_flow_current
                        , upstream_flow_next = upstream_flow_next
                        , downstream_stage_current = downstream_stage_current
                        , downstream_stage_next = downstream_stage_next)
        self.thetas = self.thetasinv
        self.compute_sections(section_arr = self.sections
                        , j_current = j_current
                        , j_next = j_next
                        , upstream_flow_current = upstream_flow_current
                        , upstream_flow_next = upstream_flow_next
                        , downstream_stage_current = downstream_stage_current
                        , downstream_stage_next = downstream_stage_next
                        , predictor_step = False)
        self.matrix_pc(section_arr = self.sections
                        , j_current = j_current
                        , j_next = j_next
                        , upstream_flow_current = upstream_flow_current
                        , upstream_flow_next = upstream_flow_next
                        , downstream_stage_current = downstream_stage_current
                        , downstream_stage_next = downstream_stage_next
                        , predictor_step = False)
        self.apply_corrector(section_arr = self.sections
                        , j_current = j_current
                        , j_next = j_next
                        , upstream_flow_current = upstream_flow_current
                        , upstream_flow_next = upstream_flow_next
                        , downstream_stage_current = downstream_stage_current
                        , downstream_stage_next = downstream_stage_next)

    def dsbc():
          # c----------------------------------------------------------------------
          # c     Downstream boundary condition
          #       subroutine dsbc(n)
          # c----------------------------------------------------------------------
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
                    #TODO: Make sure we are handling the DS Boundary correctly
                    section_j.depth = downstream_stage_current
                self.debug = False
                if self.debug and i == self.I_DOWNSTREAM:
                    print(f'predictor step {predictor_step} i {i}')
                    print(f'depth {section_j.depth} depthp {section_j.depthp}')
                self.debug = False
                area = section_j.flow_area
                flow = section_j.flow
                depth = section_j.depth
            elif not predictor_step:
                self.debug = False
                if self.debug and i == self.I_DOWNSTREAM:
                    print(f'corrector step i {i}')
                    print(f'depth {section_j.depth} depthp {section_j.depthp}')
                self.debug = False
                area = section_j.areap # computed via 'apply_predictor' method
                flow = section_j.qp # computed via 'apply_predictor' method
                depth = section_j.depthp

            if self.debug:
                print (f'flow_area {section_j.flow_area}')
                print (f'i, z: {i} {section_j.water_z}')
                print(f'depth {section_j.depth}')
                print(f'depthp {section_j.depthp}')
                print(f'yy {self.yy}')
                print(f'bw {section.bottom_width}')
                print(f'current stage {downstream_stage_current}')
                print(f'next stage {downstream_stage_next}')
            section_j.ci1 = section.get_ci1_depth(depth)
            section_j.hy = area / section.get_wetted_perimeter_depth(depth)
            section_j.conveyance_ds = section.manning_n_ds * area * \
                                        section_j.hy ** (2.0/3.0)
            self.debug = True
            if self.debug:
                if i == self.I_UPSTREAM or i == self.I_DOWNSTREAM:
                    print (f'ci1 {section_j.ci1: 9g} Rw {section_j.hy: 9g} conveyance {section_j.conveyance_ds: 9g}')
                self.debug = False
            if predictor_step:
                if i == self.I_UPSTREAM:
                    pass
                else:
                    section_US = section.us_section
                    section_US_j = section_US.time_steps[j_current]
                    #TODO: USE DBDX in Ci2 computation
                    section_j.dbdx_ds = section_US.get_dbdx_ds_depth(section_US_j.depth, depth)                                  / section.dx_ds
                    section_j.ci2_ds = \
                        section_US.get_ci2_depth_depth_ds(section_US_j.depth, depth)
                    if self.debug:
                        print('Upstream n, flow, conveyance_ds {} {} {}', \
                                    section_US.manning_n_ds\
                                    , section_US_j.flow\
                                    , section_US_j.conveyance_ds)
                        print('This section n, flow, conveyance_ds {} {} {}', \
                                    section.manning_n_ds\
                                    , section_j.flow\
                                    , section_j.conveyance_ds)
                    #TODO: DongHa  Determine if we need to make the 0.5 to compute fs a variable
                    #TODO: Move manning_m into conveyance calculation
                    #TODO: consider Nazmul suggestion to combine Manning N and Manning M
                    section_j.friction_slope_ds = 0.5 * self.manning_m * (
                                     section_US_j.flow * abs(section_US_j.flow)
                                   / (section_US_j.conveyance_ds ** 2.0)
                                   + section_j.flow * abs(section_j.flow)
                                   / (section_j.conveyance_ds ** 2.0))
                    section_j.as0_ds = (section_US_j.flow_area + section_j.flow_area) \
                                   / 2.0 * (section_US.bed_slope_ds \
                                            - section_j.friction_slope_ds)
                    section_j.gs0_ds = self.gravity * (section_US.bed_slope_ds \
                                            - section_j.friction_slope_ds)
            if not predictor_step:
                if i == self.I_DOWNSTREAM:
                    pass
                else:
                    section_DS = section.ds_section
                    section_DS_j = section_DS.time_steps[j_current]
                    section_j.dbdx_ds = section.get_dbdx_ds_depth(depth, section_DS_j.depthp)                                  / section.dx_ds
                    section_j.ci2_ds = \
                        section.get_ci2_depth_depth_ds(depth, section_DS_j.depthp)
                        #TODO: Verify that depthp is calculated in time.
                    self.debug = True
                    if self.debug:
                        if  i == self.I_DOWNSTREAM or i == self.I_UPSTREAM:
                            print ('depthp{: 9g} depthp_DS{: 9g} {: 9g}'\
                                 .format(section_j.ci2_ds, section_j.dbdx_ds
                                 , section_j.friction_slope_ds))
                        self.debug = False
                    section_j.friction_slope_ds = 0.5 * self.manning_m * (
                                     section_j.flow * abs(section_j.flow)
                                   / (section_j.conveyance_ds ** 2.0)
                                   + section_DS_j.flow * abs(section_DS_j.flow)
                                   / (section_DS_j.conveyance_ds ** 2.0))
                    section_j.as0_ds = (section_j.areap + section_DS_j.areap) \
                           / 2.0 * (section.bed_slope_ds \
                                    - section_j.friction_slope_ds)
                    section_j.gs0_ds = self.gravity * (section.bed_slope_ds \
                                            - section_j.friction_slope_ds)
                self.debug = True
                if self.debug:
                    if  i == self.I_DOWNSTREAM or i == self.I_UPSTREAM:
                        print ('ci2 {: 9g} dbdx {: 9g} fs{: 9g}'\
                             .format(section_j.ci2_ds, section_j.dbdx_ds
                             , section_j.friction_slope_ds))
                    self.debug = False
    def matrix_pc(self
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
            if self.predictor_step:
                area = section_j.flow_area
                flow = section_j.flow
            else:
                # If we are on the second step, applying the predictors,
                # use the areap and qp
                area = section_j.areap
                flow = section_j.qp
            # print(f'area, flow {area}, {flow}')
            section_j.velocity = flow / area
            section_j.celerity = section.get_celerity_area(area, self.gravity
                                                        , debug = self.debug)
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
            f11 = -(section_j.velocity - section_j.celerity) \
                    / (2.0 * section_j.celerity)
            f12 = (section_j.velocity + section_j.celerity) \
                    / (2.0 * section_j.celerity)
            f21 = -(section_j.velocity ** 2.0 - section_j.celerity ** 2.0) \
                    / (2.0 * section_j.celerity)
            f22 = (section_j.velocity ** 2.0 - section_j.celerity ** 2.0) \
                    / (2.0 * section_j.celerity)

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
                                   , section_j.dbdx_ds
                                   , self.gravity, self .manning_m)
            # TODO: Determine if st22 is a section-property-dependent value
            # and if so, move it to a method within the RectangleSection class
            st22 = -2.0 * self.manning_m * area \
                 * area \
                 / section_j.conveyance_ds

            section_US = section.us_section
            if predictor_step:
                # TODO: Confirm the order of these statements with Ehab
                # The following if statements run in reverse order in matrixp
                # vs. matrixc
                # TODO: combine the two inner statements if order is important
                # TODO: Combine all three statments if order is not important
                thetassign = -1.0
                if section.dx_ds < self.dx_tolerance:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section.dx_ds
                    self.crmax = max(self.crmax , section_j.sigma_ds \
                            * max(section_j.d11, section_j.d22))
                    self.crmin = min(self.crmin , section_j.sigma_ds \
                            * max(section_j.d11, section_j.d22))

                # c     LHS of eq 7
                section_j.b11 = (0.5 - self.phi
                                    - self.theta  * section_j.sigma_ds * a11
                                    + thetassign * 0.5 * self.thetas * st11 * dt)
                section_j.b12 = (-self.theta  * section_j.sigma_ds * a12
                                    + thetassign * 0.5 * self.thetas * st12 * dt)
                section_j.b21 = (-self.theta * section_j.sigma_ds * a21
                                    + thetassign * 0.5 * self.thetas * st21 * dt)
                section_j.b22 = (0.5 - self.phi
                                    - self.theta  * section_j.sigma_ds * a22
                                    + thetassign * 0.5 * self.thetas * st22 * dt)

                if i == self.I_UPSTREAM:
                    section_j.sigma_ds = dt
                elif section_US.dx_ds < self.dx_tolerance:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section_US.dx_ds

            elif not predictor_step:
                thetassign = 1.0
                # TODO: Confirm that these non-parallel if statements are congruent
                if i == self.I_UPSTREAM:
                    section_j.sigma_ds = dt
                elif section.dx_ds < self.dx_tolerance:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section_US.dx_ds

                # c     LHS of eq 7
                section_j.b11 = (0.5 - self.phi - self.theta  * section_j.sigma_ds * a11
                                    + thetassign * 0.5 * self.thetas * st11 * dt)
                section_j.b12 = (-self.theta  * section_j.sigma_ds * a12
                                    + thetassign * 0.5 * self.thetas * st12 * dt)
                section_j.b21 = (-self.theta * section_j.sigma_ds * a21
                                    + thetassign * 0.5 * self.thetas * st21 * dt)
                section_j.b22 = (0.5 - self.phi - self.theta  * section_j.sigma_ds * a22
                                    + thetassign * 0.5 * self.thetas * st22 * dt)

                if section.dx_ds == 0.0:
                    section_j.sigma_ds = dt
                else:
                    section_j.sigma_ds = dt / section.dx_ds
                    self.crmax = max(self.crmax , section_j.sigma_ds * max(section_j.d11, section_j.d22))
                    self.crmin = min(self.crmin , section_j.sigma_ds * max(section_j.d11, section_j.d22))

            if predictor_step:
                #TODO: ASK Ehab if 'dt' missing from matrixc computation is intentional
                #If it is an omission, these two equations can be combined...
                thetassign = -1.0
                g11 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a11
                              + thetassign * 0.5 * self.thetas * st11 * dt)
                g12 = (self.theta * section_j.sigma_ds * a12 + thetassign * 0.5
                              * self.thetas * st12 * dt)
                g21 = (self.theta * section_j.sigma_ds * a21 + thetassign * 0.5
                              * self.thetas * st21 * dt)
                g22 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a22
                              + thetassign * 0.5 * self.thetas * st22 * dt)
            elif not predictor_step: # This reads more clearly that the equivalent simple 'else'
                thetassign = 1.0
                g11 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a11
                              + thetassign * 0.5 * self.thetas * st11)
                g12 = (self.theta * section_j.sigma_ds * a12 + thetassign * 0.5
                              * self.thetas * st12)
                g21 = (self.theta * section_j.sigma_ds * a21 + thetassign * 0.5
                              * self.thetas * st21)
                g22 = (0.5 + self.phi + self.theta * section_j.sigma_ds * a22
                              + thetassign * 0.5 * self.thetas * st22)

            section_j.g11inv =  g22 / (g11 * g22 - g12 * g21)
            section_j.g12inv = -g12 / (g11 * g22 - g12 * g21)
            section_j.g21inv = -g21 / (g11 * g22 - g12 * g21)
            section_j.g22inv =  g11 / (g11 * g22 - g12 * g21)

            section_j.f1 = flow
            section_j.f2 = flow ** 2.0 \
                        / area \
                        + self.gravity * section_j.ci1

          #      if(i.ge.2.and.i.lt.n_sections) then
            if i > self.I_DOWNSTREAM and i < self.I_UPSTREAM:
                section_DS = section.ds_section
                section_DS_j = section_DS.time_steps[j_current]
                section_US = section.us_section
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
                # Assign values for eps2 at boundaries
                if i == self.I_DOWNSTREAM - 1:
                    section_DS_j.eps2 = section_j.eps2
                elif i == self.I_UPSTREAM + 1:
                    section_US_j.eps2 = section_j.eps2

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
                section_DS = section.ds_section
                section_DS_j = section_DS.time_steps[j_current]
                section_j.eps2 = max(section_DS_j.eps2, section_j.eps2)
                # print('0.0, alfa4, eps2, velocity, celerity {} {} {} {} {}'.format(0.0, self.alfa4, section_j.eps2
                #                 , section_j.velocity, section_j.celerity))
                section_j.eps4 = max(0.0, self.alfa4 - section_j.eps2
                                / (section_j.velocity + section_j.celerity))
              #      eps4(i+1)=max(0.,alfa4-eps2(i+1)/(u(i+1)+c(i+1)))
              #c      write(*,*)'corr',i,eps2(i)

        elif not predictor_step:
            for k, section in enumerate(reversed(section_arr)):
                i = self.I_DOWNSTREAM - k
                # print(k)
                if i == self.I_DOWNSTREAM:
                    continue
                section_j = section.time_steps[j_current]
                section_DS = section.ds_section
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
                    section_DS = section.ds_section
                    section_DS_j = section_DS.time_steps[j_current]
                    area_DS = section_DS_j.flow_area
                    flow_DS = section_DS_j.flow
                else:
                    section_US = section.us_section
                    section_US_j = section_US.time_steps[j_current]
                    area_US = section_US_j.areap
                    flow_US = section_US_j.qp
                ei = max(abs(section_j.velocity + section_j.celerity),
                                    abs(section_j.velocity - section_j.celerity))
                if predictor_step:
                    ei1 = max(abs(section_DS_j.velocity + section_DS_j.celerity),
                                    abs(section_DS_j.velocity - section_DS_j.celerity))
                elif not predictor_step:
                    ei1 = max(abs(section_US_j.velocity + section_US_j.celerity),
                                    abs(section_US_j.velocity - section_US_j.celerity))
                eia = (ei + ei1) / 2.0
              #      if(ityp(i-1).ne.1) then
              #        d1(i)=0.0
              #        d2(i)=0.0
              #      elseif(i.eq.2.or.i.eq.(n_sections-1)) then
                if predictor_step:
                    area_difference = area_DS - area
                    flow_difference = flow_DS - flow
                elif not predictor_step: # This if could technically be combined with the above if statement
                                   # but it reads more clearly to have it here.
                    area_difference = area - area_US
                    flow_difference = flow - flow_US
                section_j.d1 = section_j.eps2 * eia * area_difference
                section_j.d2 = section_j.eps2 * eia * flow_difference
                if (i > self.I_UPSTREAM + 1 and i < self.I_DOWNSTREAM - 1):
                    # print(i, self.I_DOWNSTREAM, self.I_UPSTREAM)
              #        d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))
              #        d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))
              #TODO: Determine if the reduction for the next-to-boundary sections
              # is unnecessary (i.e., Could the if statement be replaced by simply
              # using the next lines and letting the equation reduce automatically
              # based on the value in eps4?
                    if predictor_step:
                        section_US = section.us_section
                        section_US_j = section_US.time_steps[j_current]
                        area_US = section_US_j.flow_area
                        flow_US = section_US_j.flow
                        section_2DS = section.ds_section.ds_section
                        section_2DS_j = section_DS.time_steps[j_current]
                        area_2DS = section_2DS_j.flow_area
                        flow_2DS = section_2DS_j.flow
                        d1_eps4_diff = - 1.0 * section_j.eps4 * (area_2DS
                                    - 3.0 * area_DS + 3.0 * area - area_US)
                        d2_eps4_diff = - 1.0 * section_j.eps4 * (flow_2DS
                                    - 3.0 * flow_DS + 3.0 * flow - flow_US)
                    elif not predictor_step:
                        section_DS = section.ds_section
                        section_DS_j = section_DS.time_steps[j_current]
                        area_DS = section_DS_j.areap
                        flow_DS = section_DS_j.qp
                        section_2US = section.us_section.us_section
                        section_2US_j = section_2US.time_steps[j_current]
                        area_2US = section_2US_j.areap
                        flow_2US = section_2US_j.qp
                        d1_eps4_diff = - 1.0 * section_j.eps4 * (area_DS
                                - 3.0 * area + 3.0 * area_US - area_2US)
                        d2_eps4_diff = - 1.0 * section_j.eps4 * (flow_DS
                                - 3.0 * flow + 3.0 * flow_US - flow_2US)
                    section_j.d1 = section_j.d1 + d1_eps4_diff
                    section_j.d2 = section_j.d2 + d2_eps4_diff
              #      else
              #      d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))-eps4(i)*
              #    1        (areap(i+1)-3*areap(i)+3*areap(i-1)-areap(i-2))
              #      d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))-eps4(i)*(qp(i+1)-3*qp(i)+
              #    1      3*qp(i-1)-qp(i-2))
              #      endif
              #c      write(*,*)'corr',i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)

    def apply_predictor(self
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
            # Use the flow time series for the upstream boundary
            section_j = section.time_steps[j_current]
            section_j.sigma_ds = dt / section.dx_ds
            if i == self.I_UPSTREAM:
                #Set the upstream boundary flow
                #This is set as a delta -- a correction factor -- to work in the
                #predictor sweep.
                #So we set the Delta Q predictor, which is set from the input
                #time series for the upstream-most point.
                # TODO: Confirm that DAP(1) is never used... Why not? ASK EHAB
                # section_j.delta_area_predictor = 0.0
                section_j.delta_flow_predictor = upstream_flow_next - upstream_flow_current
            elif i == self.I_DOWNSTREAM:
                # section_j.delta_area_predictor = 0.0
                #TODO: Ask Ehab why these values are touched in the predictor step
                section_j.delta_area_corrector = 0.0
                section_j.delta_flow_corrector = section_j.delta_flow_predictor
            else:
                section_US = section.us_section
                section_US_j = section_US.time_steps[j_current]
                section_j.rhs1 = -1.0 * section_j.sigma_ds\
                        * (section_j.f1 - section_US_j.f1\
                            - section_j.d1 + section_US_j.d1)
                section_j.rhs2 = -1.0 * section_j.sigma_ds\
                        * (section_j.f2 - section_US_j.f2\
                            - section_j.d2 + section_US_j.d2)\
                        + dt * self.gravity\
                             * (section_j.ci2_ds + section_j.as0_ds)
                c11 = section_j.g11inv * section_US_j.b11 + section_j.g12inv \
                                                            * section_US_j.b21
                c12 = section_j.g11inv * section_US_j.b12 + section_j.g12inv \
                                                            * section_US_j.b22
                c21 = section_j.g21inv * section_US_j.b11 + section_j.g22inv \
                                                            * section_US_j.b21
                c22 = section_j.g21inv * section_US_j.b12 + section_j.g22inv \
                                                            * section_US_j.b22
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
            if self.debug and (i == self.I_UPSTREAM or i == self.I_DOWNSTREAM): print (f'area {section_j.flow_area}, dap {section_j.delta_area_predictor}')
            section_j.areap = section_j.flow_area + section_j.delta_area_predictor
            section_j.depthp = section.get_depth_area(section_j.areap)
            if self.debug and (i == self.I_UPSTREAM or i == self.I_DOWNSTREAM): print (f'flow {section_j.flow_area}, dqp {section_j.delta_flow_predictor}')
            section_j.qp = section_j.flow + section_j.delta_flow_predictor

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
        for k, section in enumerate(reversed(self.sections)):
        # for i, section in enumerate((self.sections)):
            i = self.I_DOWNSTREAM - k
        #TODO: CONFIRM That the reversed array is necessary in this case (I don't think it is...)
            # Use the flow time series for the upstream boundary
            section_j = section.time_steps[j_current]
            if i == self.I_DOWNSTREAM:
                #Set the upstream boundary flow
                #This is set as a delta -- a correction factor -- to work in the
                #predictor sweep.
                #So we set the Delta Q predictor, which is set from the input
                #time series for the upstream-most point.
                section_j.delta_area_corrector = 0.0
                section_j.delta_flow_corrector = section_j.delta_flow_predictor
            else:
                #### WHY DOES THIS NOT WORK WITH THE APPROPRIATE DS SECTION????
                section_DS = section.ds_section
                #### WHY DOES THIS NOT WORK WITH THE APPROPRIATE DS SECTION????

                section_DS_j = section_DS.time_steps[j_current]
                section_j.sigma_ds = dt / section.dx_ds
                section_j.rhs1 = -1.0 * section_j.sigma_ds\
                        * (section_DS_j.f1 - section_j.f1\
                            - section_DS_j.d1 + section_j.d1)
                section_j.rhs2 = -1.0 * section_j.sigma_ds\
                        * (section_DS_j.f2 - section_j.f2\
                            - section_DS_j.d2 + section_j.d2)\
                        + dt * self.gravity\
                             * (section_j.ci2_ds + section_j.as0_ds)
                c11 = section_j.g11inv * section_DS_j.b11 + section_j.g12inv \
                                                            * section_DS_j.b21
                c12 = section_j.g11inv * section_DS_j.b12 + section_j.g12inv \
                                                            * section_DS_j.b22
                c21 = section_j.g21inv * section_DS_j.b11 + section_j.g22inv \
                                                            * section_DS_j.b21
                c22 = section_j.g21inv * section_DS_j.b12 + section_j.g22inv \
                                                            * section_DS_j.b22
                section_j.delta_area_corrector =\
                    section_j.g11inv * section_j.rhs1\
                    + section_j.g12inv * section_j.rhs2\
                    - c11 * section_DS_j.delta_area_corrector\
                    - c12 * section_DS_j.delta_flow_corrector
                section_j.delta_flow_corrector =\
                    section_j.g21inv * section_j.rhs1\
                    + section_j.g22inv * section_j.rhs2\
                    - c21 * section_DS_j.delta_area_corrector\
                    - c22 * section_DS_j.delta_flow_corrector

            da = (section_j.delta_area_predictor + section_j.delta_area_corrector) \
                    / 2.0
            dq = (section_j.delta_flow_predictor + section_j.delta_flow_corrector) \
                    / 2.0
            if (da + section_j.flow_area) > self.area_tolerance:
                next_flow_area = section_j.flow_area + da
            else:
                next_flow_area = self.area_tolerance
            next_depth = section.get_depth_area(next_flow_area)
            next_flow = section_j.flow + dq
            self.debug = True
            section.time_steps.append(self.TimeStep(new_time = self.time_list[j_next]
                                , new_flow = next_flow
                                , new_depth = next_depth
                                , new_water_z = section.bottom_z + next_depth
                                , new_area = next_flow_area))
            section_jnext = section.time_steps[j_next]
            # if self.debug:
            if self.debug and (i == self.I_UPSTREAM or i == self.I_UPSTREAM + 1):
                print('current depth area flow {: 10g} {: 9g} {: 9g}'.format\
                        (section_j.depth, section_j.flow_area
                            , section_j.flow))
                print(f'              da    dq             {da: 9g} {dq: 9g}')
                print('next    depth area flow {: 10g} {: 9g} {: 9g}'.format\
                        (next_depth, next_flow_area, next_flow))
                print('next    depth area flow {: 10g} {: 9g} {: 9g}\n(in new array)'.format\
                        (section_jnext.depth, section_jnext.flow_area
                            , section_jnext.flow))
            self.debug = False

    class RectangleSection(Network.RectangleSection):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def get_ci1_depth(self, depth):
            return self.bottom_width * (depth ** 2.0) / 2.0

        def get_dbdx_ds_depth(self, depth, depth_ds):
            #FOR RectangleSection -- the depth is trivial
            return (self.ds_section.bottom_width - self.bottom_width) \
                      / (self.dx_ds)

        def get_ci2_depth_depth_ds(self, depth, depth_ds):
            # print((depth_ds ** 2.0) + (depth ** 2.0))
            # print(self.ds_section.bottom_width - self.bottom_width)
            # print(self.dx_ds)
            return (((depth_ds ** 2.0) + (depth ** 2.0))
                      * (self.ds_section.bottom_width - self.bottom_width) \
                      / (self.dx_ds * 4.0))

        def get_dkda_area(self, area):
            wetted_perimeter = self.get_wetted_perimeter_area(area)
            return (1 / self.manning_n_ds *
                  (( 5.0 / 3.0 * area ** (2.0/3.0) * wetted_perimeter) -
                  (area ** (5.0/3.0) * 2.0 / self.bottom_width)) /
                  (wetted_perimeter ** 2.0))

        def get_st21_area(self, area, flow, gs0, conveyance,
                                dkda, dbdx, gravity, manning_m):
            return (gravity * area / self.bottom_width \
                / self.bottom_width * dbdx + gs0 \
                + manning_m * 2.0 * gravity \
                * area * flow \
                * abs (flow) / conveyance ** 3.0 \
                * dkda)

    class TimeStep(Network.TimeStep):
        '''MESH-specific time-step values'''
        def __init__(self, new_water_z = 0.0, *args, **kwargs):
            super().__init__(*args, **kwargs)

            # Per-time-step at-a-section properties
            self.delta_flow_corrector = 0.0
            self.delta_flow_predictor = 0.0
            self.delta_area_corrector = 0.0
            self.delta_area_predictor = 0.0
            self.water_z = new_water_z
            self.areap = 0.0
            self.qp = 0.0
            self.depthp = 0.0
            self.ci1 = 0.0
            self.hy = 0.0 # Hydraulic Radius (used to compute co)

            # Per-time-step downstream reach properties
            self.conveyance_ds = 0
            self.ci2_ds = 0
            # self.friction_slope_ds = 0 # Derived from parent Class, Network.TimeStep
            self.as0_ds = 0
            self.gs0_ds = 0
            self.sigma_ds = 0 # Sigma is related to the courant parameter: CFL = celerity * sigma
            #self.cour = 0
            self.dbdx_ds = 0
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
    # root = os.path.abspath(r'c:/Users/james.halgren/Downloads/MESH_test/')
    # Main_Example_Path = os.path.join(root , 'US')
    # Sub_Example_Path = os.path.join(Main_Example_Path , 'BW')
    # This_Example_Path = os.path.join(Sub_Example_Path, 'Q')

    #C:\Users\james.halgren\Downloads\MESH_test\US\BW\Q\Qvar_us_2YNorm\Qvar_us_0033_5.0-10000.0_0100_0000-0004-0200_2NormalDepth
    # input_path = os.path.join(This_Example_Path,'Qvar_us_2YNorm','Qvar_us_0033_5.0-10000.0_0100_0000-0004-0200_2NormalDepth',"input.txt")
    root = os.path.abspath(os.path.dirname(__file__))
    test_folder = os.path.join(root, r'test')
    output_folder = os.path.join(test_folder, r'output')
    input_path = os.path.join(test_folder, r'input.txt')
    output_path = os.path.join(output_folder, r'out.txt')

    input_vars[r'input_path'] = input_path
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

    network.compute_initial_state(write_output = True
                                                    , output_path = output_path)
    network.debug = False
    network.compute_time_steps(verbose = True, write_output = True
                                                    , output_path = output_path)
    network.output_dump(output_path = output_path, verbose = True)

if __name__ == "__main__":
    main()
