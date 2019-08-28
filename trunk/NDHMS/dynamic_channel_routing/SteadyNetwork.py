# import required modules
from __future__ import division
import sys
import numpy as np
import pandas as pd
from scipy.optimize import fmin
from scipy import stats
from math import ceil
#from plotly.offline import plot
#import plotly.graph_objs as go
import matplotlib.pyplot as plt
import csv
import os
GRAVITY = 9.81
MANNING_SI = 1.00
MANNING_UK = 1.49
'''Actual value: 1.4859185775 = (1 / 0.3048) ** (1/3)'''

def Bernoulli_Energy(WSE, V, hl = 0, gravity = GRAVITY):
    ''' compute flow using Manning's equation '''
    return (WSE) + (V ** 2.0) / (2 * gravity) - hl

def y_direct(B, n, S0, Q, k=MANNING_SI):
    ''' function to compute error in depth calculation for a guessed depth compared to a calculated depth for a given flow.
        Uses scipy.optimize fmin '''
    y_optimum = fmin(flow_min, Q/B/3, args=(n, S0, Q, B, k), full_output=True, disp=False)
    return float(y_optimum[0])

def flow_min(y, n, S0, Q, B, k = MANNING_SI):
    ''' computes the error term comparing the Manning's computed depth with the given Q '''
    epsilon = np.abs(Manning_Q(y, n, S0, B, k) - Q)
    return epsilon

def Manning_Slope(n, Q, A, Rw, k = MANNING_SI):
    #print(f'n * Q / (k * A * (Rw ** (2/3)))) ** 2.0 {n} * {Q} / ({k} * {A} * ({Rw} ** (2/3)))) ** 2.0')
    return (n * Q / (k * A * (Rw ** (2.0/3.0)))) ** 2.0

def Manning_Q(y, n, S0, B, k = MANNING_SI):
    ''' compute flow using Manning's equation '''
    return (k/n)*(S0**0.5)*((B*y)**(5/3))/((B+y)**(2/3))

# Generate_Hydrograph from Maryam Asgari Lamjiri and Kelly Flint
# NWC Summer institute coordinators 2019
def Generate_Hydrograph (Num_TimeStep,Steady_time,width,Skewness,QMax):
    # Steady_time += 359
    XX = np.linspace(-5,5,Num_TimeStep)
    Steady_t = -5 + Steady_time/Num_TimeStep*10
    t = (XX - Steady_t)/width
    P = 2.0 / width * stats.norm.pdf(t) * stats.norm.cdf(Skewness*t)
    Q = (QMax-100) * P/np.max(P)+100
    # plt.plot(np.arange(0,Num_TimeStep),Q)
    return Q


#########################################################
#            FINITE DIFFERENCE METHOD                   #
#                                                       #
#  A program for one dimensional flow in open channels  #
#                                                       #
#########################################################

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
        new_depth = y_direct(section.bottom_width, section.manning_n_ds, section.bed_slope_ds, new_flow)
        section.time_steps.append(self.TimeStep(new_flow=new_flow, new_depth=new_depth))

class SteadyNetwork(Network):
    #TODO: These Input and Initialize methods could be different methods within the Network class
    def input_and_initialize(self, input_opt=1, input_path=None, output_path=None, upstream_flow_ts=None, downstream_stage_ts=None):
        ''' This input option is intended to be an extremely simple channel for testing and plotting development'''
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
        self.add_steady_time_step(j_current, j_next, self.sections, downstream_stage_next, upstream_flow_next)

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

    def y_standard_step(self, j, section_ds, section_us, Q, y_guess = 1.0, k=MANNING_SI, gravity=GRAVITY):
    #TODO: Make GRAVITY and MANNING_K constants consistent with anticipated units and
    #get them to be called/passsed consistently.
        ''' function to compute error in depth calculation for a guessed depth compared to a calculated depth for a given flow.
            Uses scipy.optimize fmin '''
        y_ds = section_ds.time_steps[j].depth
        z_us = section_us.bottom_z
        z_ds = section_ds.bottom_z
        n_us = section_us.manning_n_ds
        n_ds = section_ds.manning_n_ds

        loss_coeff = section_us.loss_coeff_ds # The loss coefficient is a reach property so we use only one of them from the upstream reach
        dx = section_us.dx_ds # The downstream distance is a reach property so we use only one of them from the upstream reach

        # y_standard = fmin(self.stage_standard_min, y_guess, args=(y_ds, z_us, z_ds, V_us, V_ds, hl_us2ds, gravity), full_output=True, disp=True)
        y_standard = fmin(self.stage_standard_min, y_guess, args=(section_us, section_ds, loss_coeff, dx, Q, y_ds, z_us, z_ds, n_us, n_ds, gravity, k), full_output=True, disp=False)

        y_us = float(y_standard[0])

        print(f'{dx}')
        print(f'y_us: {y_us}   y_ds: {y_ds}')
        print(f'z_us: {z_us}   z_ds: {z_ds}')
        print(f'n_us: {n_us}   n_ds: {n_ds}')
        print(f'{loss_coeff}')

        Area_us = section_us.get_area_depth(y_us)
        Area_ds = section_ds.get_area_depth(y_ds)
        print(f'Area_us: {Area_us}   Area_ds: {Area_ds}')
        V_us = Q / Area_us
        V_ds = Q / Area_ds
        print(f'V_us: {V_us}   V_ds: {V_ds}')
        V_head_us = V_us ** 2.0 / (2 * gravity) 
        V_head_ds = V_ds ** 2.0 / (2 * gravity)
        print(f'V_head_us: {V_head_us}   V_head_ds: {V_head_ds}')
        Pw_us = section_us.get_wetted_perimeter_depth(y_us)
        Pw_ds = section_ds.get_wetted_perimeter_depth(y_ds)
        print(f'Pw_us: {Pw_us}   Pw_ds: {Pw_ds}')
        Rw_us = Area_us / Pw_us
        Rw_ds = Area_ds / Pw_ds
        print(f'Rw_us: {Rw_us}   Rw_ds: {Rw_ds}')
        Area_mean = (Area_ds + Area_us) / 2.0
        n_mean = (n_ds + n_us) / 2.0
        Rw_mean = (Rw_ds + Rw_us) / 2.0
        Sf_ds = Manning_Slope(Q = Q, A = Area_ds, Rw = Rw_ds, n = n_ds, k = k)
        Sf_us = Manning_Slope(Q = Q, A = Area_us, Rw = Rw_us, n = n_us, k = k)
        Sf_mean = Manning_Slope(Q = Q, A = Area_mean, Rw = Rw_mean, n = n_mean, k = k)
        #Sf_mean = (Sf_ds + Sf_us) / 2.0
        print (f'Q {Q}')
        print(f'Sf_ds = {Sf_ds} Sf_us = {Sf_us} Sf_mean = {Sf_mean}')
        hl_us2ds = Sf_mean * dx + loss_coeff * abs(V_head_us -  V_head_ds)
        print('hl_us2ds = {hl_us2ds}')
        print('y_us, y_guess, y_ds, z_us, z_ds, Q, V_us, V_ds, hl_us2ds, MANNING_SI')
        print(y_us, y_guess, y_ds, z_us, z_ds, Q, V_us, V_ds, hl_us2ds, MANNING_SI)
        WSE_ds = y_ds + z_ds
        WSE_us = y_us + z_us
        E_ds = Bernoulli_Energy(WSE_ds, V_ds, 0, gravity) 
        E_us = Bernoulli_Energy(WSE_us, V_us, hl_us2ds, gravity)
        print(f'E_ds: {E_ds}')
        print(f'E_us: {E_us}')

        return float(y_standard[0])

    def stage_standard_min(self, y_us, section_us, section_ds, loss_coeff, dx, Q, y_ds, z_us, z_ds, n_us, n_ds, gravity = GRAVITY, k=MANNING_SI):
        ''' computes the error term comparing the Manning's computed depth with the given Q '''
        Area_us = section_us.get_area_depth(y_us)
        Area_ds = section_ds.get_area_depth(y_ds)
        # print(f'{Area_us} {Area_ds}')
        V_us = Q / Area_us
        V_ds = Q / Area_ds
        # print(f'{V_us} {V_ds}')
        V_head_us = V_us ** 2.0 / (2 * gravity) 
        V_head_ds = V_ds ** 2.0 / (2 * gravity)
        # print(f'{V_head_us} {V_head_ds}')
        Pw_us = section_us.get_wetted_perimeter_depth(y_us)
        Pw_ds = section_ds.get_wetted_perimeter_depth(y_ds)
        # print(f'{Pw_us} {Pw_ds}')
        Rw_us = Area_us / Pw_us
        Rw_ds = Area_ds / Pw_ds
        # print(f'{Rw_us} {Rw_ds}')

        # print (f'Q {Q}')

        Area_mean = (Area_ds + Area_us) / 2.0
        n_mean = (n_ds + n_us) / 2.0
        Rw_mean = (Rw_ds + Rw_us) / 2.0
        Sf_ds = Manning_Slope(Q = Q, A = Area_ds, Rw = Rw_ds, n = n_ds, k = k)
        Sf_us = Manning_Slope(Q = Q, A = Area_us, Rw = Rw_us, n = n_us, k = k)
        #Sf_mean = Manning_Slope(Q = Q, A = Area_mean, Rw = Rw_mean, n = n_mean, k = k)
        Sf_mean = (Sf_ds + Sf_us) / 2.0
        # print(f'Sf_ds = {Sf_ds} Sf_us = {Sf_us} Sf_mean = {Sf_mean}')

        hl_us2ds = Sf_mean * dx + loss_coeff * abs(V_head_us -  V_head_ds)
        # print('hl_us2ds')
        # print(hl_us2ds)

        # print(y_ds, z_ds)
        # print(y_us, z_us)
        WSE_ds = y_ds + z_ds
        WSE_us = y_us + z_us
        E_ds = Bernoulli_Energy(WSE_ds, V_ds, 0, gravity) 
        E_us = Bernoulli_Energy(WSE_us, V_us, hl_us2ds, gravity)
        # print(f'E_ds: {E_ds}')
        # print(f'E_us: {E_us}')
        epsilon = np.abs(E_us - E_ds)
        #epsilon = np.abs(Bernoulli_Energy(WSE_ds, V_ds, 0, gravity) - Bernoulli_Energy(WSE_us, V_us, hl_us2ds, gravity))
        return epsilon


    class TimeStep(Network.TimeStep):
        def __init__(self, *args, **kwargs):
            # super(Network.TimeStep, self).__init__(*args, **kwargs)
            super().__init__(*args, **kwargs)

class Section:
    #TODO: The Section Class needs to be sub-classed with Different types,
    #e.g., SectionRectangle, SectionTrapezoid, SectionTrapFlood (for the type that
    #currently used in the National Water Model), SectionDepthArea, SectionDepthWidth, ...
    #def __init__(self, bottom_width, side_slope):
    def __init__(self, bottom_width, bottom_z, comid=None, station=None, dx_ds = 10, manning_n_ds = 0.015):
        #Time independent at-a-station properties
        self.comid = comid
        self.station = station
        self.bottom_width = bottom_width
        self.bottom_z = bottom_z
        self.manning_n_ds = manning_n_ds
        self.time_steps = [] # array of values
        self.sk = MANNING_SI

        #Time independent downstream reach properties
        self.dx_ds = 0 # Distance to downstream section
        self.loss_coeff_ds = 0 # Contraction and other loss coefficients to downstream section
                            # C in the following equation
                            # hl = Sf * dx + C * (V1**2/2g - V2**2/2g)
        self.dbdx_ds = 0 # Distance to downstream section
        self.bed_slope_ds = 0 # Bed slope (S0) to downstream section
        #ADD NEIGHBOR Concept

    def get_area_depth(self, depth):
        return self.bottom_width * depth

    def get_area_j(self, j):
        return self.bottom_z * self.time_steps[j].depth

    def get_wetted_perimeter_depth(self, depth):
        return self.bottom_width + 2.0 * depth

    def get_wetted_perimeter_j(self, j):
        return self.bottom_z + 2.0 * self.time_steps[j].depth

def main():
    # network = DummyNetwork()
    # network = SimpleFlowTrace() #DongHa's method.
    network = SteadyNetwork()
    # network = MuskCNetwork()
    # network = MESHDNetwork()

    network.input_and_initialize()
    network.compute_initial_state()
    network.compute_time_steps()

if __name__ == "__main__":
    main()
