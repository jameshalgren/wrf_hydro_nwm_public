# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
import sys
import os
import time
import multiprocessing
from functools import partial
connections = None
networks = None

ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL: 
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')
    sys.path.append(r'../fortran_routing/mc_pylink_v00/MC_singleCH_singleTS')
    sys.setrecursionlimit(4000)

import pickle
import networkbuilder
import nhd_network_traversal as nnt
import mc_sc_stime as mc
import numpy as np

def set_networks(
    supernetwork = ''
    , geo_input_folder = None
    , verbose = True
    , debuglevel = 0
    ):

    supernetwork_data = nnt.set_supernetwork_data(
      supernetwork = supernetwork
      , geo_input_folder = geo_input_folder
      ) 
    supernetwork_values = nnt.get_nhd_connections(
      supernetwork_data = supernetwork_data
      , verbose = verbose
      , debuglevel = debuglevel
      )
    return supernetwork_data, supernetwork_values

def recursive_junction_read (
                             segments
                             , order_iter
                             , connections
                             , network
                             , terminal_code = 0
                             , verbose = False
                             , debuglevel = 0
                            ):
    #global connections, network

    #Greatest to least ordering of the river system for computation.
    #IDs are in order of magnitude so largest order IDs will be computed first from the list reading left to right
    #Can easily change to great sublists for each magnitude of order if necessary
    #temporary function to gather unique order listings

    for segment in segments:
        csegment = segment
        if 1 == 1:
        #try:
            reach = {}
            reach.update({'reach_tail':csegment})
            reach.update({'downstream_reach':connections[csegment]['downstream']})
            reachset = set()
            reachset.add(csegment)
            usegments = connections[segment]['upstreams']
            while True: 
                if usegments == {terminal_code}: # HEADWATERS
                    if debuglevel <= -3: print(f"headwater found at {csegment}")
                    network['total_segment_count'] += 1
                    if debuglevel <= -3: print(f"segs at csegment {csegment}: {network['total_segment_count']}")
                    reachset.add(csegment)
                    reach.update({'reach_head':csegment})
                    reach.update({'order':order_iter})
                    if order_iter == 0: network.update({'terminal_reach':csegment})#; import pdb; pdb.set_trace() #TODO: FIX THIS; SEEMS FRAGILE
                    network.update({'maximum_order':max(network['maximum_order'],order_iter)})
                    reach.update({'segments':reachset})
                    network['reaches'].update({csegment:reach})
                    network['headwaters'].add(csegment)
                    break
                elif len(usegments) >= 2: # JUNCTIONS
                    if debuglevel <= -3: print(f"junction found at {csegment} with upstreams {usegments}")
                    network['total_segment_count'] += 1
                    if debuglevel <= -3: print(f"segs at csegment {csegment}: {network['total_segment_count']}")
                    reachset.add(csegment)
                    reach.update({'reach_head':csegment})
                    reach.update({'order':order_iter})
                    if order_iter == 0: network.update({'terminal_reach':csegment})#; import pdb; pdb.set_trace() #TODO: FIX THIS; SEEMS FRAGILE
                    network.update({'maximum_order':max(network['maximum_order'],order_iter)})
                    reach.update({'segments':reachset})
                    network['reaches'].update({csegment:reach})
                    network['total_junction_count'] += 1 #the Terminal Segment
                    network['junctions'].add(csegment)
                    recursive_junction_read (
                            usegments
                            , order_iter + 1
                            , connections
                            , network
                            , terminal_code = terminal_code
                            , verbose = verbose
                            , debuglevel = debuglevel) 
                    break
                if debuglevel <= -3: print(f"segs at csegment {csegment}: {network['total_segment_count']}")
                # the terminal code will indicate a headwater
                if debuglevel <= -4: print(usegments)
                (csegment,) = usegments
                usegments = connections[csegment]['upstreams']
                network['total_segment_count'] += 1
                reachset.add(csegment)

                # print(usegments)
        #except:
            #if debuglevel <= -2: 
                #print(f'There is a problem with connection: {segment}: {connections[segment]}')

def network_trace(
        terminal_segment
        , order_iter
        , connections
        , terminal_code = 0
        , verbose= False
        , debuglevel = 0
    ):

    network = {}
    us_length_total = 0
    
    if debuglevel <= -1: print(f'\ntraversing upstream on network {terminal_segment}:')
    # try:
    if 1 == 1:
        network.update({'total_segment_count': 0}) 
        network.update({'total_junction_count': 0})
        network.update({'maximum_order':0})
        network.update({'junctions':set()})
        network.update({'headwaters':set()})
        network.update({'reaches':{}}) 
        recursive_junction_read(
                  [terminal_segment]
                  , order_iter
                  , connections
                  , network
                  , verbose = verbose
                  , terminal_code = terminal_code
                  , debuglevel = debuglevel)
        if verbose: print(f"junctions: {network['total_junction_count']}")
        if verbose: print(f"segments: {network['total_segment_count']}")
    # except Exception as exc:
    #     print(exc)
    #TODO: compute upstream length as a surrogate for the routing computation
    return {terminal_segment: network, 'upstream_length': us_length_total}

def compose_reaches(
        supernetwork_values = None
        , terminal_code = 0
        , debuglevel = 0
        , verbose = False
        , showtiming = False
    ):

    terminal_segments = supernetwork_values[4] 
    circular_segments = supernetwork_values[6]
    terminal_segments_super = terminal_segments - circular_segments
    connections = supernetwork_values[0]
        
    networks = {terminal_segment:{}
                      for terminal_segment in terminal_segments_super 
                     }  

    if verbose: print('verbose output')
    if verbose: print(f'number of Independent Networks to be analyzed is {len(networks)}')
    if verbose: print(f'Multi-processing will use {multiprocessing.cpu_count()} CPUs')
    if verbose: print(f'debuglevel is {debuglevel}')

    results_serial = {}
    init_order = 0
    for terminal_segment, network in networks.items():
        network.update(network_trace(terminal_segment, init_order, connections, terminal_code = terminal_code, verbose = verbose, debuglevel = debuglevel)[terminal_segment])
        up_reaches = networkbuilder.get_up_connections(
            network['reaches']
            , terminal_code
            , network['headwaters']
            , {network['terminal_reach']}
            , r'upstream_reaches'
            , r'downstream_reach'
            , verbose = verbose
            , debuglevel = debuglevel
            )

        if debuglevel <= -1:
            if debuglevel <=-1: print(f'terminal_segment: {terminal_segment}')
            if debuglevel <=-2: 
                for k, v in network.items():
                    if type(v) is dict:
                        print (f'{k}:')
                        for k1, v1 in v.items():
                            print(f'\t{k1}: {v1}')
                    else: print(f'{k}: {v}')
    if debuglevel <= -1: print(f'Number of networks in the Supernetwork: {len(networks.items())}')

    return networks

def compute_network(
        terminal_segment = None
        , network = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    global connections

    if verbose: print(f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_order']}")

    ordered_reaches = {}
    for head_segment, reach in network['reaches'].items():
        if reach['order'] not in ordered_reaches:
            ordered_reaches.update({reach['order']:[]}) #TODO: Should this be a set/dictionary?
        ordered_reaches[reach['order']].append([head_segment
                  , reach
                  ])
    for x in range(network['maximum_order'],-1,-1):
        for head_segment, reach in ordered_reaches[x]:
            compute_reach_up2down(
                head_segment = head_segment
                , reach = reach
                # , connections = connections
                , supernetwork_data = supernetwork_data
                , verbose = verbose
                , debuglevel = debuglevel
                )

#TODO: generalize with a direction flag
def compute_reach_up2down(
        head_segment = None
        , reach = None
        # , connections = None
        , supernetwork_data = None
        , verbose = False
        , debuglevel = 0
        ):
    global connections
    if debuglevel <= -1: print(f"\nreach: {head_segment} (order: {reach['order']} n_segs: {len(reach['segments'])})")
    current_segment = reach['reach_head']
    next_segment = connections[current_segment]['downstream'] 
    while True:
        # Thanks to SO post for a reminder of this "Loop-and-a-half" construct
        # https://stackoverflow.com/questions/1662161/is-there-a-do-until-in-python
        compute_segment(
            current_segment = current_segment
            , supernetwork_data = supernetwork_data
            , verbose = verbose
            , debuglevel = debuglevel
            )
        if current_segment == reach['reach_tail']:
            if debuglevel <= -1: print(f'{current_segment} (tail)')
            break
        if debuglevel <= -1: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream'] 

def compute_segment(
    current_segment = None
    , supernetwork_data = None
    # , connections = None
    , verbose = False
    , debuglevel = 0
    ):
    global connections
    
    if debuglevel <= -2:
        print(f'data for segment: {current_segment}')
        for k, v in connections[current_segment].items():
            if type(v) is list: 
                print (f'{k}:')
                print (f"manningn: {v[supernetwork_data['manningn_col']]}")
                print (f"slope: {v[supernetwork_data['slope_col']]}")
                print (f"bottom width: {v[supernetwork_data['bottomwidth_col']]}")
                print (f"Muskingum K: {v[supernetwork_data['MusK_col']]}")
                print (f"Muskingum X: {v[supernetwork_data['MusX_col']]}")
                print (f"Channel Side Slope: {v[supernetwork_data['ChSlp_col']]}")
            else: print(f'{k}: {v}')

    #"compute MC"
    computeMC = True
    #"compute Dummy"
    computeDummy = False
    #"compute MSH"
    computeMSH= False

     
    if computeMC:
       compute_mc(connections[current_segment]
            , supernetwork_data = supernetwork_data
            , verbose = verbose
            , debuglevel = debuglevel)
        
    elif computeDummy:
        print ('computeDummy')
                   
    elif computeMSH:
        print ('computeMSH')                   
                 


# ### Psuedocode
# 
# ```
# Call Compute Network
#     for each reach in the network
#         Call compute reach
#             For each segment in the reach
#                 Import the mc object/module
#                 Call prepare segment
#                     Populate the Mc segment array with the individual segment properties
#                     Obtain and set any lateral inflows
#                     obtain and set the upstream and downstream (requrires logic to deal with junctions or boundaries at headwaters)
#             With the populated arrays, execute MC for the reach
# ```     
#         

def get_upstream_inflow():
    pass

def get_lateral_inflow():
    pass

def compute_junction_downstream():
    pass

def compute_mc(
        current_segment = None
        , supernetwork_data = None
        , verbose = False
        , debuglevel = 0
        ):
    #global connections

    
 
    #print (current_segment['data'][supernetwork_data['manningn_col']])
    #print (current_segment['data'][supernetwork_data['slope_col']])
    # print (current_segment['data'][supernetwork_data['bottomwidth_col']])
    #print (current_segment['data'][supernetwork_data['MusK_col']])
    #print (current_segment['data'][supernetwork_data['MusX_col']])
    #print (current_segment['data'][supernetwork_data['ChSlp_col']])

    #if verbose: 
              

    #print (f"manningn: {connections[current_segment]}")
    #print ('computeMC')

    # Thanks to SO post for a reminder of this "Loop-and-a-half" construct
    # https://stackoverflow.com/questions/1662161/is-there-a-do-until-in-python

    ntim=2;       #the number of time steps necessary for variables passed to mc module to compute successfully
    nlinks=2;     #the number of links needed to define varialbe qd. ** nlinks is not used in fortran source code.

    ncomp0=3; mc.var.ncomp0=ncomp0  #the number of segments of a reach upstream of the current reach
    ncomp=3; mc.var.ncomp=ncomp  #the number of segments of the current reach 
    mxseg=max(ncomp0,ncomp)
    mc.var.uslinkid=1
    mc.var.linkid=2

    #MC model outputs
    mc.var.qd=np.zeros((ntim,mxseg,nlinks))  #will store MC output qdc (flow downstream current timestep) 
    mc.var.vela=np.zeros((ntim,ncomp)) 
    mc.var.deptha=np.zeros((ntim,ncomp))

    #print
    #import pdb; pdb.set_trace()
    
    #lateral flow
    mc.var.qlat=np.zeros((ncomp))

    dt=30.0 ;      mc.var.dt= dt  #60.0;
    dx=current_segment['data'][supernetwork_data['length_col']] ;     mc.var.dx= dx  #20.0
    bw=current_segment['data'][supernetwork_data['bottomwidth_col']];       mc.var.bw= bw #50
    tw= 0.01*bw; mc.var.tw= tw
    twcc=tw;     mc.var.twcc=twcc
    n=current_segment['data'][supernetwork_data['manningn_col']] ;      mc.var.n=n #0.03
    ncc=n;       mc.var.ncc=ncc
    cs=current_segment['data'][supernetwork_data['ChSlp_col']] ;    mc.var.cs=cs #1.0e6
    so=current_segment['data'][supernetwork_data['slope_col']];    mc.var.so=so #0.002
    ck= current_segment['data'][supernetwork_data['MusK_col']];   mc.var.ck = ck 
    cx= current_segment['data'][supernetwork_data['MusX_col']];   mc.var.cx = cx

    #import pdb; pdb.set_trace()
    #run M-C model
    nts=1  #the number of timestep in simulation
    wnlinks=20 #the number of all links in simulation
    wmxseg=5  #max number of segments among all links

    #variable storing all outputs in time
    wqd= np.zeros((nts,wmxseg,wnlinks))   
    wvela= np.zeros((nts,wmxseg,wnlinks)) 
    wdeptha= np.zeros((nts,wmxseg,wnlinks))

    for k in range (0,nts):        
        #input lateral flow for current reach; input flow to upstream reach of current reach
        for i in range(0,ncomp):
            mc.var.qlat[i]= (k+1)*2.0
            mc.var.qd[1,i,0]= (k+1)*10.0

        mc.mc.main()

        for i in range(0,ncomp):
            #current link(=reach)
            #qd[k,i,j]: k=0/1: previous/current timestep; i: node ID; j=0/1: upstream/current reach
            mc.var.qd[0,i,1]= mc.var.qd[1,i,1]
            mc.var.vela[0,i]= mc.var.vela[1,i]
            mc.var.deptha[0,i]= mc.var.deptha[1,i]
            #upstream link(=reach)        
            mc.var.qd[0,i,0]=  mc.var.qd[1,i,0] 

        #output keeping
        j=1 #temporarily assigned link ID for the current reach
        for i in range(0,ncomp):
            wqd[k,i,j]= mc.var.qd[1,i,1]
            wvela[k,i,j]= mc.var.vela[1,i]
            wdeptha[k,i,j]= mc.var.deptha[1,i]

    #test output    
    j=1
    #print('MC_TEST_OUTPUT')
#     for k in range (0,nts):
#         for i in range(0,ncomp):
#             print(wqd[k,i,j])

#     for k in range (0,nts):
#         for i in range(0,ncomp):
#             print(wvela[k,i,j])

#     for k in range (0,nts):
#         for i in range(0,ncomp):
#             print(wdeptha[k,i,j])

def main():

    global connections
    global networks

    verbose = True
    debuglevel = 0
    showtiming = True

    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo', r'Channels')

    #TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_ge5'
    """##NHD CONUS order 5 and greater"""
    supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    # supernetwork = 'Mainstems_CONUS'
    # supernetwork = 'CONUS_FULL_RES_v20'
    # supernetwork = 'CONUS_Named_Streams' #create a subset of the full resolution by reading the GNIS field
    # supernetwork = 'CONUS_Named_combined' #process the Named streams through the Full-Res paths to join the many hanging reaches


    if verbose: print('creating supernetwork connections set')
    if showtiming: start_time = time.time()
    #STEP 1
    supernetwork_data, supernetwork_values = set_networks(
        supernetwork = supernetwork
        , geo_input_folder = geo_input_folder
        , verbose = False
        # , verbose = verbose
        , debuglevel = debuglevel
        )
    if verbose: print('supernetwork connections set complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 2
    if showtiming: start_time = time.time()
    if verbose: print('organizing connections into reaches ...')
    networks = compose_reaches(
        supernetwork_values
        , verbose = False
        # , verbose = verbose
        , debuglevel = debuglevel
        , showtiming = showtiming
        )
    if verbose: print('reach organization complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 3
    if showtiming: start_time = time.time()
    if verbose: print('executing computation on ordered reaches ...')
    connections = supernetwork_values[0]
    for terminal_segment, network in networks.items():
        compute_network(
            terminal_segment = terminal_segment
            , network = network
            , supernetwork_data = supernetwork_data
            # , connections = connections
            , verbose = False
            # , verbose = verbose
            , debuglevel = debuglevel
        )
    if verbose: print('ordered reach computation complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

if __name__ == '__main__':
    main()
