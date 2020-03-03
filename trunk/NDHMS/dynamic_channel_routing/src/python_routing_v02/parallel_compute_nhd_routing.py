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

# WARNING: These global declarations cause the parallel implementation to 
# crash when executed on Windows
connections = None
networks = None
from sys import platform
if platform == "linux" or platform == "linux2":
    pass
elif platform == "darwin":
    pass
elif platform == "win32":
    print('The parallel version of compute_nhd_routing.py will not execute as currently')
    print('written due to the lack of a fork() capability in the windows OS.')
    print('For parallel execution, please us a *nix OS.')
    print('\nexiting...')
    sys.exit()
    # Some useful references:
    # https://stackoverflow.com/questions/985281/what-is-the-closest-thing-windows-has-to-fork/985525#985525
    # https://stackoverflow.com/questions/8220108/how-do-i-check-the-operating-system-in-python
    # https://stackoverflow.com/questions/6596617/python-multiprocess-diff-between-windows-and-linux

ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL: 
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')
    sys.setrecursionlimit(4000)

import pickle
import networkbuilder
import nhd_network_traversal as nnt

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

def compute_network_parallel(
        large_networks = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    
    global connections

    overall_max = -1
    for terminal_segment, network in large_networks.items():
        overall_max = max(network['maximum_order'], overall_max)

    if verbose: print(f"Executing simulation for all large network beginning with maximum order {overall_max}")

    if 1 == 1:
        ordered_reaches = {}
        for terminal_segment, network in large_networks.items():
            for head_segment, reach in network['reaches'].items():
                if reach['order'] not in ordered_reaches:
                    ordered_reaches.update({reach['order']:[]}) #TODO: Should this be a set/dictionary?
                ordered_reaches[reach['order']].append([head_segment
                          , reach
                          , supernetwork_data
                          , verbose
                          , debuglevel])
        with multiprocessing.Pool() as netpool:
            for x in range(overall_max,-1,-1):
                nslist3 = ordered_reaches[x]
                if verbose: print(f"Executing simulation for {len(nslist3)} large network reaches of order {x}")
                results = netpool.starmap(compute_reach_up2down, nslist3)

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

    MC_CALL = False
    if MC_CALL:
        print (connections[current_segment]['data'][supernetwork_data['manningn_col']])
        print (connections[current_segment]['data'][supernetwork_data['slope_col']])
        print (connections[current_segment]['data'][supernetwork_data['bottomwidth_col']])
        print (connections[current_segment]['data'][supernetwork_data['MusK_col']])
        print (connections[current_segment]['data'][supernetwork_data['MusX_col']])
        print (connections[current_segment]['data'][supernetwork_data['ChSlp_col']])

def get_upstream_inflow():
    pass

def get_lateral_inflow():
    pass

def compute_junction_downstream():
    pass

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

    parallel_split = 10 # -1 turns off the splitting and runs everything through the lumped execution
    #STEP 3a -- Large Networks
    if showtiming: start_time = time.time()
    if verbose: print(f'executing computation on ordered reaches for networks of order greater than {parallel_split} ...')
    connections = supernetwork_values[0]
    large_networks = {terminal_segment: network \
                      for terminal_segment, network in networks.items() \
                      if network['maximum_order'] > parallel_split}
    # print(large_networks)
    compute_network_parallel(
        large_networks
        , supernetwork_data = supernetwork_data
        # , connections = connections
        , verbose = False
        # , verbose = verbose
        , debuglevel = debuglevel
    )
    if verbose: print(f'ordered reach computation complete for networks of order greater than {parallel_split}')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 3b -- Small Networks
    if showtiming: start_time = time.time()
    if verbose: print(f'executing parallel computation on ordered reaches of order less than {parallel_split} ...')
    nslist = ([terminal_segment
                  , network 
                  , supernetwork_data #TODO: This should probably be global...
                  , False
                  , debuglevel] 
                  for terminal_segment, network in networks.items() 
                  if network['maximum_order'] <= parallel_split)
    with multiprocessing.Pool() as pool:
        results = pool.starmap(compute_network, nslist)
    if verbose: print(f'parallel ordered reach computation complete of order less than {parallel_split}')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

if __name__ == '__main__':
    main()
