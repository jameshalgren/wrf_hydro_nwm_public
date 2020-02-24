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

ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL: 
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')
    sys.setrecursionlimit(4000)

import pickle
import networkbuilder
import nhd_network_traversal as nnt

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
        return nnt.get_nhd_connections(
        supernetwork = Brazos_LowerColorado_ge5_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)
        return Brazos_LowerColorado_ge5_values

    elif CONUS_ge5:
        CONUS_ge5_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_ge5'
        , geo_input_folder = geo_input_folder)
        return nnt.get_nhd_connections(
        supernetwork = CONUS_ge5_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    elif CONUS_FULL_RES_v20:
        CONUS_FULL_RES_v20_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_FULL_RES_v20'
        , geo_input_folder = geo_input_folder)
        return nnt.get_nhd_connections(
        supernetwork = CONUS_FULL_RES_v20_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    elif CONUS_Named_Streams:
        CONUS_Named_Streams_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_Named_Streams'
        , geo_input_folder = geo_input_folder)

        if not CONUS_Named_combined:
            return nnt.get_nhd_connections(
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

            
            # Use only headwater keys that are in the full dataset.
            headwater_keys_combined = CONUS_FULL_RES_v20_values[3] & \
                                        CONUS_Named_Streams_values[3]
            # Need to make sure that these objects are independent -- we will modify them a bit.
            connections_combined = pickle.loads(pickle.dumps(CONUS_FULL_RES_v20_values[0]))
            terminal_keys_combined = pickle.loads(pickle.dumps(CONUS_FULL_RES_v20_values[4]))
            terminal_code_combined = CONUS_Named_Streams_supernetwork['terminal_code']
            
            for key in connections_combined: #Clear the upstreams and rebuild it with just named streams
                connections_combined[key].pop('upstreams',None)

            (junction_keys_combined
            , visited_keys_combined
            , visited_terminal_keys_combined
            , junction_count_combined) = networkbuilder.get_up_connections(
                        connections = connections_combined
                        , terminal_code = terminal_code_combined
                        , headwater_keys = headwater_keys_combined
                        , terminal_keys = terminal_keys_combined
                        , verbose = verbose
                        , debuglevel = debuglevel)

def recursive_junction_read (
                             keys
                             , con
                             , network
                             , terminal_code = 0
                             , verbose = False
                             , debuglevel = 0
                            ):
    #global con, network
    for key in keys:
        ckey = key
        if 1 == 1:
        #try:
            ukeys = con[key]['upstreams']
            while not len(ukeys) >= 2 and not (ukeys == {terminal_code}):
                if debuglevel <= -3: print(f"segs at ckey {ckey}: {network['total_segment_count']}")
                # the terminal code will indicate a headwater
                if debuglevel <= -4: print(ukeys)
                (ckey,) = ukeys
                ukeys = con[ckey]['upstreams']
                network['total_segment_count'] += 1
            if ukeys == {terminal_code}:
                if debuglevel <= -3: print(f"headwater found at {ckey}")
                network['total_segment_count'] += 1
                if debuglevel <= -3: print(f"segs at ckey {ckey}: {network['total_segment_count']}")
            elif len(ukeys) >= 2:
                network['total_segment_count'] += 1
                if debuglevel <= -3: print(f"junction found at {ckey} with upstreams {ukeys}")
                network['total_segment_count'] += 1
                if debuglevel <= -3: print(f"segs at ckey {ckey}: {network['total_segment_count']}")
                network['total_junction_count'] += 1 #the Terminal Segment
                recursive_junction_read (
                        ukeys
                        , con
                        , network
                        , terminal_code = terminal_code
                        , verbose = verbose
                        , debuglevel = debuglevel) 
                # print(ukeys)
                ukeys = con[ckey]['upstreams']
                ckey = ukeys
        #except:
            #if debuglevel <= -2: 
                #print(f'There is a problem with connection: {key}: {con[key]}')

def network_trace(
                        nid
                        , con
                        , terminal_code = 0
                        , verbose= False
                        , debuglevel = 0
                        ):

    network = {}
    us_length_total = 0
    
    if verbose: print(f'\ntraversing upstream on network {nid}:')
    # try:
    if 1 == 1:
        network.update({'total_segment_count': 0}) 
        network.update({'total_junction_count': 0})
        network.update({'junctions':{}})
        network.update({'reaches': {}}) #the Terminal Segment
        recursive_junction_read(
                  [nid]
                  , con
                  , network
                  , verbose = verbose
                  , terminal_code = terminal_code
                  , debuglevel = debuglevel)
        if verbose: print(f"junctions: {network['total_junction_count']}")
        if verbose: print(f"segments: {network['total_segment_count']}")
    # except Exception as exc:
    #     print(exc)
    #TODO: compute upstream length as a surrogate for the routing computation
    return {nid: network, 'upstream_length': us_length_total}

def compose_reaches(
        network_values = None
        , terminal_code = 0
        , debug_level = 0
        , verbose = False
        ):


    terminal_keys = network_values[4] 
    circular_keys = network_values[6]
    terminal_keys_super = terminal_keys - circular_keys
    con = network_values[0]
        
    networks = {terminal_key:{}
                      for terminal_key in terminal_keys_super 
                     }  
    debuglevel = -2
    verbose = False

    if verbose: print('verbose output')
    if verbose: print(f'number of Independent Networks to be analyzed is {len(super_networks)}')
    if verbose: print(f'Multi-processing will use {multiprocessing.cpu_count()} CPUs')
    if verbose: print(f'debuglevel is {debuglevel}')

    start_time = time.time()
    results_serial = {}
    for nid, network in networks.items():
        network.update(network_trace(nid, con, terminal_code = terminal_code, verbose = verbose, debuglevel = debuglevel)[nid])
    print("--- %s seconds: serial compute ---" % (time.time() - start_time))
    if debuglevel <= -1: print(len(networks.items()))
    if debuglevel <= -2: print(networks)

def main():
    network_values = set_network()
    reach_values = compose_reaches(network_values)


if __name__ == '__main__':
    main()

    # # nhd_conus_networks = []
    # # conus_supernetwork = SuperNetwork()
    #
    # print([network.networkID for network in nhd_conus_networks])
    # print([network.reachCollection for network in nhd_conus_networks])
    # print((network.reachCollection for network in nhd_conus_networks))
    # print(network.reachCollection for network in nhd_conus_networks)
    # for network in nhd_conus_networks:
    #     print([reach.reachID for reach in network.reachCollection])

    # print([reach.order for reach in network.reachCollection
                                    # for network in nhd_conus_networks])
