# -*- coding: utf-8 -*-
"""NHD Network traversal

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
import network
import reach
import nexus

# class Reach():
#     def __init__(self, reachID = None, order = -1):
#         self.reachID = reachID
#         self.segmentCollection = []
#         self.order = order
#         self.upstreamNexus = None
#         self.downstreamNexus = None
#
# class Segment():
#     def __init__(self, segmentID = None, length = -1):
#         self.segmentID = segmentID
#         self.length = length
#         self.upstreamSegment = None
#         self.downstreamSegment = None

def get_leaves_from_trunk():
    '''
    Define a function to assign leaves of a given nodeset
    If the leaves are already present this should return self.leaves of the trunk
    If the leaves have not yet been calculated, this should call the recursive
    search to define the leaves.
    '''
    pass

# NOTE: these methods can lose the "connections" and "rows" arguments when
# implemented as class methods where those arrays are members of the class.
def recurse_downstream(key, rows, key_col, downstream_col, length_col, down_connections, terminal_code):
    for row in rows:
        # print(f'{key_col}  {key}')
        if row[key_col] == key:
            down_connections[key] = {'downstream': row[downstream_col], 'length': row[length_col] }
            if not key == terminal_code:
                recurse_downstream(row[downstream_col], rows, key_col, downstream_col, length_col, down_connections, terminal_code)

def determine_keys(rows, key_col, downstream_col, terminal_code
            , verbose = False, debuglevel = 0):
    # Get the upstream nodes to start with
    if verbose: print('starting build')
    if verbose: print('all_keys')
    all_keys = [row[key_col] for row in rows]
    if debuglevel == -1: print(f'length = {len(all_keys)}')
    if verbose: print('all_keys complete')
    if debuglevel == -2: print(all_keys)
    if verbose: print('ref_keys')
    ref_keys = [row[downstream_col] for row in rows]
    if verbose: print('ref_keys complete')
    if debuglevel == -2: print(ref_keys)
    if verbose: print('headwater_keys')
    headwater_keys = [x for x in all_keys if x not in ref_keys]
    if debuglevel == -2: print(headwater_keys)
    if verbose: print('headwater_keys complete')

    # Get the downstream terminating nodes
    if verbose: print('terminal_keys')
    terminal_keys = []
    for row in rows:
        if row[downstream_col] == terminal_code:
            terminal_keys.append(row[key_col])
    if debuglevel == -2: print(terminal_keys)
    if verbose: print('terminal_keys complete')

    return all_keys, ref_keys, headwater_keys, terminal_keys

def build_connections(rows, key_col, downstream_col, length_col, terminal_code
                    , headwater_keys, terminal_keys
                    , verbose = False, debuglevel = 0):
    ''' Recursive call to go all the way down the relationship '''
    up_connections = {}
    down_connections = {}

    if verbose: print('\nStart Recursion')
    if verbose: print('down_connections')
    for key in headwater_keys:
        recurse_downstream(key = key
             , rows = rows
             , key_col = key_col
             , downstream_col = downstream_col
             , length_col = length_col
             , down_connections = down_connections
             , terminal_code = terminal_code)
    if debuglevel == -2: print(down_connections)
    if verbose: print('down_connections complete')

    # Create inverse of connections looking upstream
    if verbose: print('up_connections')
    for i in headwater_keys:
        up_connections[i] = { 'upstreams': [terminal_code], 'length': down_connections[i]['length'] }

    for i in down_connections.keys():
        if not down_connections[i]['downstream'] == terminal_code:
            up_connections[down_connections[i]['downstream']] = { 'upstreams': [], 'length': down_connections[i]['length'] }
                                    #  `down_connections[down_connections[i]['downstream']]['length']` to just `down_connections[i]['length']`
            for j in down_connections.keys():
                if down_connections[i]['downstream'] == down_connections[j]['downstream']:
                    up_connections[down_connections[i]['downstream']]['upstreams'].append(j)
    if debuglevel == -2: print(up_connections)
    if verbose: print('up_connections complete')

    return down_connections, up_connections

def rec_print_down(key, down_connections, terminal_code, debuglevel = 0):
    if key == terminal_code: return

    print(f"{key} with length {down_connections[key]['length']}")
    rec_print_down(down_connections[key]['downstream'], down_connections, terminal_code)

def rec_print_up(keys, tab_count, up_connections, down_connections, terminal_code, debuglevel = 0):
    if not isinstance(keys, list): keys = [keys]
    tab_count += 1
    for key in keys:
        if not key == terminal_code:
            print(f"{'.' * (tab_count)}\\{key} with length {down_connections[key]['length']}\\")
            rec_print_up(up_connections[key]['upstreams'], tab_count, up_connections, down_connections, terminal_code)

def print_connections(headwater_keys = None, terminal_keys = None
                    , down_connections = None, up_connections = None
                    , terminal_code = None):
    try:
        if headwater_keys:
            print("########################")
            print("Downstream Connections")
            print("########################")
            for key in headwater_keys:
                rec_print_down(key, down_connections, terminal_code)
                print("########################")

        if terminal_keys:
            print("########################")
            print("Upstream Connections")
            print("########################")
            for key in terminal_keys:
                rec_print_up([key], -1, up_connections, down_connections, terminal_code)
                print("########################")
    except:
        if verbose: print('''provide headwater_keys, down_connections, and a terminal code
to print Downstream Connections.

Provide terminal_keys, up_connections, down_connections, and a terminal code
to print Upstream Connections.''')
import uuid
def buildSuperNetwork(
                networks = None
                ):
    for iter, network in enumerate(networks):
        pass

def buildNetwork(
                rows = None
                , networks = None
                , down_connections = None
                , up_connections = None
                , key_col = None
                , downstream_col = None
                , length_col = None
                , terminal_code = None
                , headwater_keys = None
                , terminal_keys = None
                , verbose = False
                , debuglevel = 0
                ):
    for iter, trunk_key in enumerate(terminal_keys):
        # networks.head_segments = get_leaves_from_trunk(trunk)
        networks.append(network.Network(uuid.uuid4()))
        print(f'network = {networks[iter].networkID}')

        reach_iter = 0
        print(reach_iter)
        networks[iter].reachCollection.append(reach.Reach(reachID = uuid.uuid4()))
        print(f'reach = {networks[iter].reachCollection[reach_iter].reachID}')
        rec_Network_up(networks[iter].reachCollection
                            , reach_iter
                            , [trunk_key]
                            , -1
                            , up_connections
                            , down_connections
                            , terminal_code
                            , debuglevel)
        print("########################")

def rec_Network_up(reaches
                    , reach_iter
                    , keys
                    , tab_count
                    , up_connections
                    , down_connections
                    , terminal_code
                    , debuglevel = 0):
    tab_count += 1
    for key_iter, key in enumerate(keys):
        if len(keys) > 1:
            reaches.append(reach.Reach(reachID = uuid.uuid4(), order = reach_iter - key_iter))
            reach_iter += 1
            if debuglevel < 0: print(reach_iter)
            if debuglevel < -1: print(f"{'.' * int(tab_count/10)}reach = {reaches[reach_iter].reachID}")
        if not key == terminal_code:
            reaches[reach_iter].segmentCollection.append(reach.Reach.Segment(key, down_connections[key]['length']))
            if debuglevel < -2: print(f"{'.' * (tab_count)}\\{key} with length {down_connections[key]['length']}\\")
            rec_Network_up(
                        reaches
                        , reach_iter
                        , up_connections[key]['upstreams']
                        , tab_count
                        , up_connections
                        , down_connections
                        , terminal_code
                        , debuglevel)

if 1 == 1:
    """##TEST"""
    print("")
    print ('Executing Test')
    # Test data
    test_rows = [
        [0,456,-999,0],
        [1,178,4,0],
        [2,394,0,0],
        [3,301,2,0],
        [4,798,0,0],
        [5,679,4,0],
        [6,523,0,0],
        [7,815,2,0],
        [8,841,-999,0],
        [9,514,8,0],
        [10,458,9,0],
        [11,832,10,0],
        [12,543,11,0],
        [13,240,12,0],
        [14,548,13,0],
        [15,920,14,0],
        [16,920,15,0],
        [17,514,16,0],
        [18,458,17,0],
        [19,832,18,0],
        [20,543,19,0],
        [21,240,16,0],
        [22,548,21,0],
        [23,920,22,0],
        [24,240,23,0],
        [25,548,12,0],
        [26,920,25,0],
        [27,920,26,0],
        [28,920,27,0],
    ]

    test_key_col = 0
    test_downstream_col = 2
    test_length_col = 1
    test_terminal_code = -999

    (test_all_keys, test_ref_keys, test_headwater_keys
     , test_terminal_keys) = determine_keys(
                 rows = test_rows
                 , key_col = test_key_col
                 , downstream_col = test_downstream_col
                 , terminal_code = test_terminal_code
                 , verbose = True
                 , debuglevel = -2)

    (test_down_connections
     , test_up_connections) = build_connections(
                 rows = test_rows
                 , key_col = test_key_col
                 , downstream_col = test_downstream_col
                 , length_col = test_length_col
                 , terminal_code = test_terminal_code
                 , headwater_keys = test_headwater_keys
                 , terminal_keys = test_terminal_keys
                 , verbose = True
                 , debuglevel = -2)

    print_connections(headwater_keys = test_headwater_keys
                    , terminal_keys = test_terminal_keys
                    , down_connections = test_down_connections
                    , up_connections = test_up_connections
                    , terminal_code = test_terminal_code)

    test_networks = []
    test_supernetwork = network.SuperNetwork()

    buildNetwork(
                rows = test_rows
                , networks = test_networks
                , down_connections = test_down_connections
                , up_connections = test_up_connections
                , key_col = test_key_col
                , downstream_col = test_downstream_col
                , length_col = test_length_col
                , terminal_code = test_terminal_code
                , headwater_keys = test_headwater_keys
                , terminal_keys = test_terminal_keys
                , verbose = True
                , debuglevel = -3
                )

if 1 == 1:
    """##NHD Subset (Brazos/Lower Colorado)"""
    import sys
    import os
    import geopandas as gpd
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
     , terminal_keys_NHD) = determine_keys(
                 rows = nhd_conus_rows
                 , key_col = key_col_NHD
                 , downstream_col = downstream_col_NHD
                 , terminal_code = terminal_code_NHD
                 , verbose = True
                 , debuglevel = 0)

    (down_connections_NHD
     , up_connections_NHD) = build_connections(
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

    buildNetwork(
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
