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

def set_network():
    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo', r'Channels')

    """##NHD Subset (Brazos/Lower Colorado)"""
    Brazos_LowerColorado_ge5 = True
    """##NHD CONUS order 5 and greater"""
    CONUS_ge5 = False
    """These are large -- be careful"""
    CONUS_FULL_RES_v20 = False
    CONUS_Named_Streams = False #create a subset of the full resolution by reading the GNIS field
    CONUS_Named_combined = False #process the Named streams through the Full-Res paths to join the many hanging reaches

    debuglevel = -1
    verbose = True

    if Brazos_LowerColorado_ge5:
        Brazos_LowerColorado_ge5_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'Brazos_LowerColorado_ge5'
        , geo_input_folder = geo_input_folder)
        return Brazos_LowerColorado_ge5_supernetwork, nnt.get_nhd_connections(
        supernetwork = Brazos_LowerColorado_ge5_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    elif CONUS_ge5:
        CONUS_ge5_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_ge5'
        , geo_input_folder = geo_input_folder)
        return CONUS_ge5_supernetwork, nnt.get_nhd_connections(
        supernetwork = CONUS_ge5_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    elif CONUS_FULL_RES_v20:
        CONUS_FULL_RES_v20_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_FULL_RES_v20'
        , geo_input_folder = geo_input_folder)
        return CONUS_FULL_RES_v20_supernetwork, nnt.get_nhd_connections(
        supernetwork = CONUS_FULL_RES_v20_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    elif CONUS_Named_Streams:
        CONUS_Named_Streams_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_Named_Streams'
        , geo_input_folder = geo_input_folder)

        if not CONUS_Named_combined:
            return CONUS_Named_Streams_supernetwork, nnt.get_nhd_connections(
            supernetwork = CONUS_Named_Streams_supernetwork
            , verbose = verbose
            , debuglevel = debuglevel)

        else:
            ''' NOW Combine the two CONUS analyses by starting with the Named Headwaters
                but trace the network down the Full Resolution NHD. It should only work
                if the other two datasets have been computed.
                ANY OTHER Set of Headerwaters could be substituted'''
            
            if not (CONUS_Named_Streams and CONUS_FULL_RES_v20):
                print('\n\nWARNING: If this works, you may be using old data...')

            
            # Use only headwater segments that are in the full dataset.
            headwater_segments_combined = CONUS_FULL_RES_v20_values[3] & \
                                        CONUS_Named_Streams_values[3]
            # Need to make sure that these objects are independent -- we will modify them a bit.
            connections_combined = pickle.loads(pickle.dumps(CONUS_FULL_RES_v20_values[0]))
            terminal_segments_combined = pickle.loads(pickle.dumps(CONUS_FULL_RES_v20_values[4]))
            terminal_code_combined = CONUS_Named_Streams_supernetwork['terminal_code']
            
            for segment in connections_combined: #Clear the upstreams and rebuild it with just named streams
                connections_combined[segment].pop('upstreams',None)


           #TODO: THIS IS INCOMPLETE FOR THE COMBINED

            (junction_segments_combined
            , visited_segments_combined
            , visited_terminal_segments_combined
            , junction_count_combined) = networkbuilder.get_up_connections(
                        connections = connections_combined
                        , terminal_code = terminal_code_combined
                        , headwater_segments = headwater_segments_combined
                        , terminal_segments = terminal_segments_combined
                        , verbose = verbose
                        , debuglevel = debuglevel)

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
                    if order_iter == 0: network.update({'terminal_reach':csegment}) #TODO: FIX THIS; SEEMS FRAGILE
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
    
    if verbose: print(f'\ntraversing upstream on network {terminal_segment}:')
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

    start_time = time.time()
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
    print("--- %s seconds: serial compute ---" % (time.time() - start_time))
    if debuglevel <= -1: print(f'Number of networks in the Supernetwork: {len(networks.items())}')

    return networks

def compute_network(
        terminal_segment = None
        , network = None
        , supernetwork = None
        , connections = None
        , verbose = False
        , debuglevel = 0
        ):
    if verbose: print(f"\n\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_order']}")

    for x in range(network['maximum_order'],-1,-1):
        for head_segment, reach in network['reaches'].items():
            if x == reach['order']:

                compute_reach_up2down(
                    head_segment = head_segment
                    , reach = reach
                    , connections = connections
                    , supernetwork = supernetwork
                    , verbose = verbose
                    , debuglevel = debuglevel
                    )

#TODO: generalize with a direction flag
def compute_reach_up2down(
        head_segment = None
        , reach = None
        , connections = None
        , supernetwork = None
        , verbose = False
        , debuglevel = 0
        ):
    if debuglevel <= -1: print(f"\nreach: {head_segment} (order: {reach['order']} n_segs: {len(reach['segments'])})")
    current_segment = reach['reach_head']
    next_segment = connections[current_segment]['downstream'] 
    while True:
        # Thanks to SO post for a reminder of this "Loop-and-a-half" construct
        # https://stackoverflow.com/questions/1662161/is-there-a-do-until-in-python
        compute_segment(
            current_segment = current_segment
            , supernetwork = supernetwork
            , verbose = verbose
            , debuglevel = debuglevel
            )
        if current_segment == reach['reach_tail']:
            if debuglevel <= -1: print(f'{current_segment} (tail)')
            break
        if debuglevel <= -1: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream'] 

def compute_reach_up2down_mc(
        head_segment = None
        , reach = None
        , connections = None
        , supernetwork = None
        , verbose = False
        , debuglevel = 0
        ):

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


    #lateral flow
    mc.var.qlat=np.zeros((ncomp))

    dt=60.0;     mc.var.dt= dt
    dx=20.0;     mc.var.dx= dx
    bw=50;       mc.var.bw= bw
    tw= 0.01*bw; mc.var.tw= tw
    twcc=tw;     mc.var.twcc=twcc
    n=0.03;      mc.var.n=n
    ncc=n;       mc.var.ncc=ncc
    cs=1.0e6;    mc.var.cs=cs
    so=0.002;    mc.var.so=so

    #run M-C model
    nts=10  #the number of timestep in simulation
    wnlinks=20 #the number of all links in simulation
    wmxseg=5   #max number of segments among all links

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
    print('MC_TEST_OUTPUT')
    for k in range (0,nts):
      for i in range(0,ncomp):
          print(wqd[k,i,j])

    for k in range (0,nts):
      for i in range(0,ncomp):
          print(wvela[k,i,j])

    for k in range (0,nts):
      for i in range(0,ncomp):
          print(wdeptha[k,i,j])
                                              
def compute_segment(
    current_segment = None
    , supernetwork = None
    , verbose = False
    , debuglevel = 0
    ):
    global connections
    
    if debuglevel <= -1: print(f'data for segment: {current_segment}')
    if debuglevel <= -2:
        for k, v in connections[current_segment].items():
            if type(v) is list: # print the 
                print (f'{k}:')
                print (f"manningn: {v[supernetwork['manningn_col']]}")
                print (f"slope: {v[supernetwork['slope_col']]}")
                print (f"bottom width: {v[supernetwork['bottomwidth_col']]}")
                print (f"Muskingum K: {v[supernetwork['MusK_col']]}")
                print (f"Muskingum X: {v[supernetwork['MusX_col']]}")
                print (f"Channel Side Slope: {v[supernetwork['ChSlp_col']]}")
            else: print(f'{k}: {v}')

def get_upstream_inflow():
    pass

def get_lateral_inflow():
    pass

def compute_junction_downstream():
    pass

def main():

    global connections
    global networks

    debuglevel = -2
    verbose = True


    if verbose: print('creating supernetwork connections set')
    start_time = time.time()
    #STEP 1
    supernetwork, supernetwork_values = set_network()
    if verbose: print('supernetwork connections set complete')
    if debuglevel <= -1: print("--- in %s seconds: ---" % (time.time() - start_time))
    if verbose: print('ordering reaches ...')

    start_time = time.time()
    #STEP 2
    networks = compose_reaches(
        supernetwork_values
        , verbose = verbose
        , debuglevel = debuglevel
        )
    if verbose: print('ordered reaches complete')
    if debuglevel <= -1: print("--- in %s seconds: ---" % (time.time() - start_time))

    start_time = time.time()
    if verbose: print('executing computation on ordered reaches ...')
    #STEP 3
    connections = supernetwork_values[0]
    for terminal_segment, network in networks.items():
        compute_network(
            terminal_segment = terminal_segment
            , network = network
            , supernetwork = supernetwork
            , connections = connections
            , verbose = verbose
            , debuglevel = debuglevel
        )
    if verbose: print('ordered reach computation complete')
    if debuglevel <= -1: print("--- in %s seconds: ---" % (time.time() - start_time))

    #start_time = time.time()
    #if verbose: print('executing computation on ordered reaches ...')
    #connections = supernetwork_values[0]
    #for terminal_segment, network in networks.items():
        #compute_network(terminal_segment, network, connections, verbose, debuglevel)
    #if verbose: print('ordered reach computation complete')
    #if debuglevel <= -1: print("--- in %s seconds: ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()
