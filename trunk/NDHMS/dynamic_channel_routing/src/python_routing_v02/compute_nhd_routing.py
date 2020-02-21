# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
import sys
import os
ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL: 
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')

import pickle
import networkbuilder
import nhd_network_traversal as nnt

def main():

    
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
        Brazos_LowerColorado_ge5_values = nnt.get_nhd_connections(
        supernetwork = Brazos_LowerColorado_ge5_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    if CONUS_ge5:
        CONUS_ge5_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_ge5'
        , geo_input_folder = geo_input_folder)
        CONUS_ge5_values = nnt.get_nhd_connections(
        supernetwork = CONUS_ge5_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    if CONUS_FULL_RES_v20:
        CONUS_FULL_RES_v20_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_FULL_RES_v20'
        , geo_input_folder = geo_input_folder)
        CONUS_FULL_RES_v20_values = nnt.get_nhd_connections(
        supernetwork = CONUS_FULL_RES_v20_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    if CONUS_Named_Streams:
        CONUS_Named_Streams_supernetwork = \
        nnt.set_supernetwork_data(supernetwork = 'CONUS_Named_Streams'
        , geo_input_folder = geo_input_folder)
        #print(CONUS_Named_Streams_supernetwork)
        CONUS_Named_Streams_values = nnt.get_nhd_connections(
        supernetwork = CONUS_Named_Streams_supernetwork
        , verbose = verbose
        , debuglevel = debuglevel)

    if CONUS_Named_combined:
        ''' NOW Combine the two CONUS analyses by starting with the Named Headwaters
            but trace the network down the Full Resolution NHD. It should only work
            if the other two datasets have been computed.
            ANY OTHER Set of Headerwaters could be substituted'''
        
        if not (CONUS_Named_Streams and CONUS_FULL_RES_v20):
            print('\n\nWARNING: If this works, you are using old data...')

        
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
