# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
import networkbuilder
import sys
import os
import geopandas as gpd

# NOTE: these methods can lose the "connections" and "rows" arguments when
# implemented as class methods where those arrays are members of the class.

if 1 == 1:
    """##NHD Subset (Brazos/Lower Colorado)"""
    # import numpy as np

    # root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # test_folder = os.path.join(root, r'test')
    # input_folder = os.path.join(test_folder, r'input', r'geo')
    # input_path = os.path.join(input_folder, r'input.txt')
    # nhd_conus = gpd.read_file(os.path.join('..','test','input','geo','Channels','NHD_BrazosLowerColorado_Channels.shp'))
    # nhd_conus = gpd.read_file('/content/small/trunk/NDHMS/dynamic_channel_routing/test/input/text/Channels/NHD_BrazosLowerColorado_Channels.shp')
    # nhd_conus = gpd.read_file('/content/small/trunk/NDHMS/dynamic_channel_routing/test/input/text/Channels/Channels.shp')

    CONUS = False
    if CONUS:
        nhd_conus = gpd.read_file('../test/input/geo/Channels/NHD_Conus_Channels.shp')
        key_col_NHD = 1
        downstream_col_NHD = 6
        length_col_NHD = 5
        terminal_code_NHD = 0

    else:
        nhd_conus = gpd.read_file('../test/input/geo/Channels/NHD_BrazosLowerColorado_Channels.shp')
        key_col_NHD = 2
        downstream_col_NHD = 7
        length_col_NHD = 6
        terminal_code_NHD = 0

    nhd_conus_rows = nhd_conus.to_numpy()

    # Kick off recursive call for all connections and keys
    (all_keys_NHD, ref_keys_NHD, headwater_keys_NHD
     , terminal_keys_NHD) = networkbuilder.determine_keys(
                 rows = nhd_conus_rows
                 , key_col = key_col_NHD
                 , downstream_col = downstream_col_NHD
                 , terminal_code = terminal_code_NHD
                 , verbose = True
                 , debuglevel = 0)

    (down_connections_NHD
     , up_connections_NHD) = networkbuilder.traverse_segments(
                 rows = nhd_conus_rows
                 , key_col = key_col_NHD
                 , downstream_col = downstream_col_NHD
                 , length_col = length_col_NHD
                 , terminal_code = terminal_code_NHD
                 , headwater_keys = headwater_keys_NHD
                 , terminal_keys = terminal_keys_NHD
                 , verbose = True
                 , debuglevel = 0)

# print_connections(headwater_keys = test_headwater_keys
#                 , terminal_keys = test_terminal_keys
#                 , down_connections = test_down_connections
#                 , up_connections = test_up_connections
#                 , terminal_code = test_terminal_code)

    for key, value in up_connections_NHD.items():
        if len(value['upstreams']) > 1:
        # if 1 == 1:
            # print (link, link['upstreams'][:])
            print (key, value['upstreams'])


    print(len(down_connections_NHD))
    print(len(up_connections_NHD))

    nhd_conus_networks = []
    # conus_supernetwork = SuperNetwork()

    networkbuilder.buildNetwork(
                rows = nhd_conus_rows
                , networks = nhd_conus_networks
                , down_connections = down_connections_NHD
                , up_connections = up_connections_NHD
                , key_col = key_col_NHD
                , downstream_col = downstream_col_NHD
                , length_col = length_col_NHD
                , terminal_code = terminal_code_NHD
                , headwater_keys = headwater_keys_NHD
                , terminal_keys = terminal_keys_NHD
                , verbose = True
                , debuglevel = -2
                )

    print([network.networkID for network in nhd_conus_networks])
    print([network.reachCollection for network in nhd_conus_networks])
    print((network.reachCollection for network in nhd_conus_networks))
    print(network.reachCollection for network in nhd_conus_networks)
    for network in nhd_conus_networks:
        print([reach.reachID for reach in network.reachCollection])

    print([reach.order for reach in network.reachCollection
                                    for network in nhd_conus_networks])
