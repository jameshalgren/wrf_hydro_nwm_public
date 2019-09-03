# import required modules
from __future__ import division
import helpers
import constants
from network import Network
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

class DummyNetwork(Network):
    #TODO: These Input and Initialize methods could be different methods within the Network class
    def input_and_initialize(self, input_opt=1, input_path=None, output_path=None, upstream_flow_ts=None, downstream_stage_ts=None):
        ''' This input option is intended to be an extremely simple channel for testing and plotting development'''
        input_vars = {}

        self.time_list = range(101)
#        import pandas as pd
#        pandas.date_range("11:00", "21:30", freq="30min")
#        datelist = pd.date_range(pd.datetime.today(), periods=100).tolist()

        n_sections = 40
        I_UPSTREAM = n_sections - 1
        I_DOWNSTREAM = 0

        station_downstream = 10000
        station_upstream = 11000
        stations = np.linspace(station_downstream, station_upstream, n_sections, False)
        bottom_widths = np.linspace(100, 1000, len(stations), False)
        bottom_zs = np.linspace(0,100, len(stations), False)

        #print(NCOMP, len(stations))

        #/input_vars.update({"dxini": 1000})
        input_vars.update({"manning_n_ds": 0.035})

        for i, bw in enumerate(bottom_widths):
            # continue
            self.sections.append(Section(station=stations[i]
                                    , bottom_width=bottom_widths[i]
                                    , bottom_z = bottom_zs[i]
                                    , manning_n_ds = input_vars['manning_n_ds']))
            #print(sections[i].bed_slope_ds, sections[i].dx_ds, sections[i].bottom_z)
            if i == 0:
                # self.sections[i].dx_ds = input_vars['dxini'] #Irrelevant with the slope defined
                self.sections[i].bed_slope_ds = .0001
            else:
                self.sections[i].dx_ds = self.sections[i].station - self.sections[i-1].station
                self.sections[i].bed_slope_ds = (self.sections[i].bottom_z \
                                            - self.sections[i-1].bottom_z) \
                                            / self.sections[i].dx_ds

        #TODO: clean up this code to generate intial upstream flow and downstream stage boundary time series
        self.upstream_flow_ts = helpers.Generate_Hydrograph(len(self.time_list) , 20 , 2 , 4 , 5000)
        self.downstream_stage_ts = [helpers.y_direct(self.sections[I_DOWNSTREAM].bottom_width
                                             , self.sections[I_DOWNSTREAM].manning_n_ds
                                             , self.sections[I_DOWNSTREAM].bed_slope_ds
                                             , q ) for q in self.upstream_flow_ts]
        # print(self.upstream_flow_ts)
        # print(self.downstream_stage_ts)

        return input_vars

    def compute_initial_state(self):
        ''' Compute a dummy initial state
        '''
        #print(self.upstream_flow_ts)
        #print(self.downstream_stage_ts)
        for section in self.sections:
            self.add_normal_depth_time_step(section, self.upstream_flow_ts[0])

    def compute_next_time_step_state(self, j_current
                                         , j_next
                                         , upstream_flow_current
                                         , upstream_flow_next
                                         , downstream_stage_current
                                         , downstream_stage_next):
        ''' the Dummy Network simply copies the current channel state to the next time step
            flow
        '''
        for section in self.sections:
            #print(j)
            self.add_normal_depth_time_step(section, upstream_flow_next)

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
        self.sk = constants.MANNING_SI

        #Time independent downstream reach properties
        self.dx_ds = 0 # Distance to downstream section
        self.dbdx_ds = 0 # Distance to downstream section
        self.bed_slope_ds = 0 # Bed slope (S0) to downstream section
        #ADD NEIGHBOR Concept

def main():
    network = DummyNetwork()
    # network = SimpleFlowTrace() #DongHa's method.
    # network = SteadyNetwork()
    # network = MuskCNetwork()
    # network = MESHDNetwork()

    network.input_and_initialize()

    network.compute_initial_state()
    network.compute_time_steps()

    cols_for_subplots = 4
    fig, axes = plt.subplots(nrows=ceil(len(network.sections)/cols_for_subplots), ncols=cols_for_subplots, squeeze=False)

    #TODO: make plotting another method within the network
    for i, section in enumerate(network.sections):
        a = pd.Series(time_step.depth for i, time_step in enumerate(section.time_steps))
        # print(m, n)
        m = i//cols_for_subplots
        n = i%cols_for_subplots
        a.plot(ax = axes[m,n])

    return

if __name__ == "__main__":
    main()
