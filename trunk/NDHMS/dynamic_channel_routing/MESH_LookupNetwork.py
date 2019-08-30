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


class TimeStep:
    def __init__(self, new_flow = None, new_depth = None):
        # Per-time-step at-a-section properties
        self.flow = 0
        self.delta_flow_corrector = 0
        self.delta_flow_predictor = 0 
        self.delta_area_corrector = 0 
        self.delta_area_corrector = 0 
        self.depth = 0
        self.water_z = 0
        self.flow_area = 0
        self.ci1 = 0
        self.hy = 0 #(temporary variable for computing co)
        self.c0 = 0

        # Per-time-step downstream reach properties
        self.ci2_ds = 0
        self.friction_slope_ds = 0
        self.as0_ds = 0
        self.gs0_ds = 0
        self.sigma_ds = 0 # Sigma is approximately equivalent to the courant parameter, but not quite
        self.dbdx = 0

class MESH_LookupNetwork(Network):
    #TODO: These Input and Initialize methods could be different methods within the Network class
    def input_and_initialize(self, input_opt=None, input_path=None, output_path=None, upstream_flow_ts=None, downstream_stage_ts=None):
        input_data = {}
        if input_opt == 1:
            '''Use File Input from the Mesh.py application'''

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
                        sections.append(Section(i, b0, z, dxini)) # dx is a static value in the test cases
        
            # Read hydrograph input Upstream and Downstream
            with open(upstream_path, newline='') as file3:
                with open(downstream_path, newline='') as file4:
                    read_data3 = list(csv.reader(file3, delimiter=' '))
                    read_data4 = list(csv.reader(file4, delimiter=' '))
                    for j in range(ntim): # use usual convention of j as time dimension
                        q_upstream = float(read_data3[j][1])
                        y_downstream = float(read_data4[j][1])
                        sections[I_UPSTREAM].add_upstream_boundary_condition_time_step(q_upstream)
                        sections[I_DOWNSTREAM].add_downstream_boundary_condition_time_step(y_downstream) 
    
        elif input_opt == 2:
            ''' This input option is intended to be an extremely simple channel for testing and plotting development''' 
    else:
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
    

    return section_arr, input_data

        
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



    def secpred(section_arr):
        '''
                subroutine secpred 
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
          c     
                do 10 i=1,n_sections
                d(i)=areap(i)/bo(i)      
                ci1(i)=bo(i)*d(i)*d(i)/2.0
                hy=areap(i)/(2.0*d(i)+bo(i))
                co(i)=sk(i)*areap(i)*hy**(2.0/3.0)
          10    continue    
                do 20 i=1,n_sections-1
                if(ityp(i).eq.1) then      
                ci2(i)=(d(i)*d(i)+d(i+1)*d(i+1))*(bo(i+1)-bo(i))/(dx(i)*4.)
                beds=(z(i)-z(i+1))/dx(i)
                fs=f*0.5*qp(i)*abs(qp(i))/(co(i)*co(i))+
              1   f*0.5*qp(i+1)*abs(qp(i+1))/(co(i+1)*co(i+1))
                aS0(i)=(areap(i)+areap(i+1))/2.0*(beds-fs)
                gS0(i)=grav*(beds-fs)
                dbdx(i)=(bo(i+1)-bo(i))/dx(i) 
                endif
          20    continue                      
                return
                end
        '''
        return


    def compute_sections(section_arr, time_step):
        for i, section in enumerate(section_arr):
            section_US = sections[i-1]
            section_j = section.time_steps[time_step]
            section_US_j = section_US.time_steps[time_step]
            section_j.depth = sect_j.water_z - section.bottom_z
            section_j.flow_area = section.bottom_width * section_j.depth
            #TODO: All of these 
            #Calculations based on the geometry should be 
            #in class methods
            
            section_j.ci1 = section.bottom_width * (section_j.depth ** 2.0) / 2.0
            section_j.hy = section_j.flow_area / (2.0 * section_j.depth + section.bottom_width)
            section_j.c0 = section.sk * section_j.flow_area * section_j.hy ** (2.0/3.0)
            #C       write(*,*)d(i),flow_area(i),ci1(i),y(j,i),z(i)
            #TODO: Where Does this WRITE TO???
            if i == 0: #treat 
                continue
            else: 
                section_j.ci2_ds = ((section_j.depth ** 2.0) * (section_US_j.depth ** 2.0)) \
                              * (section.bottom_width - section_US.bottom_width) \
                              / (section_US.dx_ds * 4.0)
                section_US.bed_slope_ds = (section_US.bottom_z - section.bottom_z) \
                                        / section_US.dx_ds
                section_US_j.friction_slope_ds = section_US.manning_n_ds \
                               * 0.5 * section_US_j.flow \
                               * abs(section_US_j.flow) / (section_US_j.c0 ** 2.0) \
                               + section_US.manning_n_ds * 0.5 * section_j.flow \
                               * abs(section_j.flow) / (section_j.c0 ** 2.0)
                section_US_j.as0_ds = (section_US_j.area + section_j.area) \
                               / 2.0 * (section_US.bed_slope_ds \
                                        - section_US_j.bed_slope_ds)
                section_US_j.gs0_ds = GRAVITY * (section_US.bed_slope_ds \
                                        - section_US_j.bed_slope_ds)
                section_US.dbdx_ds = (section.bottom_width - section_US.bottom_width) \
                                      / section_US.dx_ds
        '''
          C in this version, we are computing the section attributes as they were
          C in the previous time step -- essentially as the initial condition for
          C the predictor step to advance to the corrector step and then to the
          C next time step
                subroutine section(n)
                parameter(grav=9.81)
                common/arrays/ area(1000),y(100000,1000),q(100000,1000),bo(1000),
              1    areap(1000),qp(1000),z(1000),dqp(1000),av11(1000),av12(1000),
              1                                       av21(1000),av22(1000),
              2              dqc(1000),dap(1000),dac(1000),ci1(1000),ci2(1000),
              3           aso(1000),d(1000),f1(1000),g11inv(1000),g12inv(1000),
              4g21inv(1000),g22inv(1000),f2(1000),b11(1000),b12(1000),b21(1000),
              5b22(1000),eps2(1000),eps4(1000),d1(1000),d2(1000),u(1000),c(1000),
              6                        sk(1000),co(1000),gso(1000),dbdx(1000),
              7                      dt(1000),ityp(1000),dx(1000)
                common/var/ncomp,cfl,f,qus,yds,rhs1,rhs2,c11,c12,c21,c22,us,ots,
              1           thes,vv,dtini
                common/matrix/phi,theta,alfa2,alfa4,w
                common/sgate/ag,bg,option,yn,eps,dmeu,yw,bw,ydsn
          c                                
                
                do 10 i=1,ncomp
                d(i)=y(n,i)-z(i)
                area(i)=bo(i)*d(i)
                ci1(i)=bo(i)*d(i)*d(i)/2.0
                hy=area(i)/(2.0*d(i)+bo(i))
                co(i)=sk(i)*area(i)*hy**(2.0/3.0)
          C       write(*,*)d(i),area(i),ci1(i),y(n,i),z(i)
          10    continue                            
          C       stop
                do 20 i=2,ncomp   
                if(ityp(i-1).eq.1) then               
                ci2(i)=(d(i)*d(i)+d(i-1)*d(i-1))*(bo(i)-bo(i-1))/(dx(i-1)*4.)
                beds=(z(i-1)-z(i))/dx(i-1)
                fs=f*0.5*q(n,i-1)*abs(q(n,i-1))/(co(i-1)*co(i-1))+
              1   f*0.5*q(n,i)*abs(q(n,i))/(co(i)*co(i))
                aso(i)=(area(i)+area(i-1))/2.0*(beds-fs)
                gso(i)=grav*(beds-fs)
                dbdx(i)=(bo(i)-bo(i-1))/dx(i-1) 
                endif
          20    continue                      
                return
                end

        '''
        return

    def matrixp():
        '''
                subroutine matrixp(j)
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
          c                                
                do 10 i=1,n_sections
                u(i)=q(j,i)/flow_area(i)
                c(i)=sqrt(grav*flow_area(i)/bo(i))
          c
          c     This is the matrix L (left eigenvector matrix - eq 13)
                e11=1.0 
                if(u(i).eq.c(i)) c(i)=c(i)+0.00001
                e12=-1.0/(u(i)-c(i))
                e21=1.0
                e22=-1.0/(u(i)+c(i))
          c
          c     L^{-1} (inverse of Left eigenvector matrix)
                f11=-(u(i)-c(i))/(2.0*c(i))
                f12=(u(i)+c(i))/(2.0*c(i))
                f21=-(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))   
                f22=(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))   
          c
          c     Diagonal wave matrix D (eq 12)
                d11=abs(u(i)+c(i))
                d22=abs(u(i)-c(i))     
          c
          c     Equation 11 (L^{-1} D L)
                a11=e11*f11*d11+e21*f12*d22         
                a12=e12*f11*d11+e22*f12*d22         
                a21=e11*f21*d11+e21*f22*d22         
                a22=e12*f21*d11+e22*f22*d22
          c      
          c
                dt(i)=dtini
          c
                ter=bo(i)+2.0*flow_area(i)/bo(i)
                dkda=sk(i)*((5.0/3.0*flow_area(i)**(2.0/3.0)*ter)-
              1     (flow_area(i)**(5.0/3.0)*2.0/bo(i)))/ter/ter
                st11=0.0
                st12=0.0
                st21=grav*flow_area(i)/bo(i)/bo(i)*dbdx(i)+gS0(i)+f*
              1     2.0*grav*flow_area(i)*q(j,i)*abs(q(j,i))/co(i)**3.0*dkda
                st22=-2*f*q(j,i)*flow_area(i)/co(i)/co(i)                
          c                                          
                if(dx(i).eq.0.0) then
                  cour=dt(i)
                else
                  cour=dt(i)/dx(i)
                  crmax=max(crmax,cour*max(d11,d22))
                  crmin=min(crmin,cour*max(d11,d22))
                endif

          c     LHS of eq 7
                b11(i)=0.5-phi-theta*cour*a11-0.5*thes*st11*dt(i)
                b12(i)=-theta*cour*a12-0.5*thes*st12*dt(i)
                b21(i)=-theta*cour*a21-0.5*thes*st21*dt(i)
                b22(i)=0.5-phi-theta*cour*a22-0.5*thes*st22*dt(i)
          c                          
                if(i == 1) then
                    cour=dt(i)
                else if (dx(i-1) == 0.0) then
                    cour=dt(i)
                else 
                  cour=dt(i) / dx(i-1)
                endif
          c     
                g11=0.5+phi+theta*cour*a11-0.5*thes*st11*dt(i)
                g12=theta*cour*a12-0.5*thes*st12*dt(i)
                g21=theta*cour*a21-0.5*thes*st21*dt(i)
                g22=0.5+phi+theta*cour*a22-0.5*thes*st22*dt(i)
          c     
                g11inv(i)= g22/(g11*g22-g12*g21)
                g12inv(i)=-g12/(g11*g22-g12*g21)
                g21inv(i)=-g21/(g11*g22-g12*g21)
                g22inv(i)= g11/(g11*g22-g12*g21) 
          c             
                f1(i)=q(j,i)
                f2(i)=q(j,i)*q(j,i)/flow_area(i)+grav*ci1(i)
          c                  
                if(i.ge.2.and.i.lt.I_UPSTREAM) then
                dip1=flow_area(i+1)/bo(i+1)
                di=2*flow_area(i)/bo(i)
                dim1=flow_area(i-1)/bo(i-1)
                eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1) 
          c      write(*,*)i,eps2(i),flow_area(i),flow_area(i+1),flow_area(i-1)
                endif              
          10    continue           
          c
                eps2(1)=eps2(2)
                eps2(I_UPSTREAM)=eps2(I_UPSTREAM-1)
          c
                do 20 i=2,n_sections-1
                  if(ityp(i).ne.1) then 
                    eps2(i)=eps2(i-1)
                    eps2(i+1)=eps2(i+2)
                  endif
          20    continue
                do 40 i=2,n_sections-1             
                eps2(i)=max(eps2(i+1),eps2(i))
          c      u(i)=(u(i+1)+u(i))/2.0
          c      c(i)=(c(i+1)+c(i))/2.0
                eps4(i)=max(0.,alfa4-eps2(i)/(u(i)+c(i))) 
          c      write(*,*)i,eps2(i)
          40    continue           
                d1(1)=0.0
                d2(1)=0.0
                d1(I_UPSTREAM)=0.0
                d2(I_UPSTREAM)=0.0
          c
                do 50 i=2,n_sections-1
                d11=abs(u(i)+c(i))
                d22=abs(u(i)-c(i))
                ei=max(d11,d22)
                d11=abs(u(i+1)+c(i+1))
                d22=abs(u(i+1)-c(i+1))
                ei1=max(d11,d22)
                eia=(ei+ei1)/2.0                              
                if(ityp(i).ne.1) then
                  d1(i)=0.0
                  d2(i)=0.0
                elseif(i.eq.2.or.i.eq.(I_UPSTREAM-1)) then ##TODO: Figure out if this is being treated correctly -- Would expect 'n_sections-1' to be more meaningful index, but the translation is confusing
                  d1(i)=eps2(i)*eia*(flow_area(i+1)-flow_area(i))
                  d2(i)=eps2(i)*eia*(q(j,i+1)-q(j,i))  
          c      write(*,*)i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)
                else
                d1(i)=eps2(i)*eia*(flow_area(i+1)-flow_area(i))-eps4(i)*(flow_area(i+2)-
              1       3*flow_area(i+1)+3*flow_area(i)-flow_area(i-1))                                       
                d2(i)=eps2(i)*eia*(q(j,i+1)-q(j,i))-eps4(i)*(q(j,i+2)-
              1       3*q(j,i+1)+3*q(j,i)-q(j,i-1))
                endif                                
          c      write(*,*)i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)
          50    continue 
                return
                end
        '''
        return

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
        return

    def matrixc():
        '''
                subroutine matrixc
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
                do 10 i=1,n_sections
                u(i)=qp(i)/areap(i)
                c(i)=sqrt(grav*areap(i)/bo(i))   
          c
                e11=1.0
                if(u(i).eq.c(i)) c(i)=c(i)+0.00001
                e12=-1.0/(u(i)-c(i))
                e21=1.0
                e22=-1.0/(u(i)+c(i))    
          c
                f11=-(u(i)-c(i))/(2.0*c(i))
                f12=(u(i)+c(i))/(2.0*c(i))
                f21=-(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))   
                f22=(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))   
          c
                d11=abs(u(i)+c(i))
                d22=abs(u(i)-c(i))
          c
                a11=e11*f11*d11+e21*f12*d22         
                a12=e12*f11*d11+e22*f12*d22         
                a21=e11*f21*d11+e21*f22*d22         
                a22=e12*f21*d11+e22*f22*d22         
          c         
                dt(i)=dtini
          c
          c     Calculating dK/dA (eq 15)                     
                ter=bo(i)+2.0*areap(i)/bo(i)
                dkda=sk(i)*((5.0/3.0*areap(i)**(2.0/3.0)*ter)-
              1     (areap(i)**(5.0/3.0)*2.0/bo(i)))/ter/ter
          c
          c     Matrix S (eq 14)
                st11=0.0
                st12=0.0
                st21=grav*areap(i)/bo(i)/bo(i)*dbdx(i)+gS0(i)+f*
              1     2.0*grav*areap(i)*qp(i)*abs(qp(i))/co(i)**3.0*dkda
                st22=-2*f*qp(i)*areap(i)/co(i)/co(i)                
          c       
          c     cour == sigma
                if(i == 1) then
                  cour=dt(i)
                else if(dx(i-1) == 0.0) then
                  cour=dt(i)
                else
                  cour=dt(i)/dx(i-1)
                endif
                b11(i)=0.5-phi-theta*cour*a11+0.5*thes*st11*dt(i)
                b12(i)=-theta*cour*a12+0.5*thes*st12*dt(i)
                b21(i)=-theta*cour*a21+0.5*thes*st21*dt(i)
                b22(i)=0.5-phi-theta*cour*a22+0.5*thes*st22*dt(i)
          c                
                if(dx(i).eq.0.0) then
                  cour=dt(i)
                else
                  cour=dt(i)/dx(i)
                  crmax=max(crmax,cour*max(d11,d22))
                  crmin=min(crmin,cour*max(d11,d22))
                endif
                g11=0.5+phi+theta*cour*a11+0.5*thes*st11
                g12=theta*cour*a12+0.5*thes*st12
                g21=theta*cour*a21+0.5*thes*st21
                g22=0.5+phi+theta*cour*a22+0.5*thes*st22
          c     
                g11inv(i)= g22/(g11*g22-g12*g21)
                g12inv(i)=-g12/(g11*g22-g12*g21)
                g21inv(i)=-g21/(g11*g22-g12*g21)
                g22inv(i)= g11/(g11*g22-g12*g21)
          c                                        
                f1(i)=qp(i)
                f2(i)=qp(i)*qp(i)/areap(i)+grav*ci1(i)        
          c
                if(i.ge.2.and.i.lt.n_sections) then                  
                dip1=areap(i+1)/bo(i+1)
                di=2*areap(i)/bo(i)
                dim1=areap(i-1)/bo(i-1)
                eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
                endif
          10    continue
                eps2(1)=eps2(2)
                eps2(I_UPSTREAM)=eps2(I_UPSTREAM-1)
          c
                do 20 i=2,n_sections-1
                  if(ityp(i).ne.1) then 
                    eps2(i)=eps2(i-1)
                    eps2(i+1)=eps2(i+2)
                  endif
          20    continue
          c
                do 40 i=n_sections-1,1,-1               
                eps2(i+1)=max(eps2(i+1),eps2(i))
          c      u(i+1)=(u(i+1)+u(i))/2.0
          c      c(i+1)=(c(i+1)+c(i))/2.0
                eps4(i+1)=max(0.,alfa4-eps2(i+1)/(u(i+1)+c(i+1)))
          c      write(*,*)'corr',i,eps2(i)
          40    continue
                d1(1)=0.0
                d2(1)=0.0
                d1(I_UPSTREAM)=0.0
                d2(I_UPSTREAM)=0.0
          c
                do 50 i=2,n_sections-1
                d11=abs(u(i)+c(i))
                d22=abs(u(i)-c(i))
                ei=max(d11,d22)
                d11=abs(u(i-1)+c(i-1))
                d22=abs(u(i-1)-c(i-1))
                ei1=max(d11,d22)
                eia=(ei+ei1)/2.0
                if(ityp(i-1).ne.1) then
                  d1(i)=0.0
                  d2(i)=0.0
                elseif(i.eq.2.or.i.eq.(n_sections-1)) then
                  d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))                                      
                  d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))                                            
                else
                d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))-eps4(i)*
              1        (areap(i+1)-3*areap(i)+3*areap(i-1)-areap(i-2))
                d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))-eps4(i)*(qp(i+1)-3*qp(i)+
              1      3*qp(i-1)-qp(i-2))
                endif
          c      write(*,*)'corr',i,d1(i),d2(i),eps2(i),flow_area(i+1),flow_area(i)
          50    continue  
                return
                end
        '''
        return


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
