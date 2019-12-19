# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
import networkbuilder as networkbuilder
import recursive_print
import os
import geopandas as gpd

# NOTE: these methods can lose the "connections" and "rows" arguments when
# implemented as class methods where those arrays are members of the class.

if 1 == 1:

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo', r'Channels')

    Brazos_FULL_RES = False
    LowerColorado_FULL_RES = False
    LowerColorado_CONCHOS_FULL_RES = False
    Mainstems_CONUS = False
    """##NHD CONUS order 5 and greater"""
    CONUS_ge5 = False
    """##NHD Subset (Brazos/Lower Colorado)"""
    Brazos_LowerColorado_ge5 = True
    CONUS_FULL_RES = False

    debuglevel = -1
    verbose = True

    # The following datasets are extracts from the feature datasets available
    # from https://www.nohrsc.noaa.gov/pub/staff/keicher/NWM_live/web/data_tools/
    # the CONUS_ge5 and Brazos_LowerColorado_ge5 datasets are included
    # in the github test folder
    if Brazos_FULL_RES:
        nhd_conus_file_path = os.path.join(geo_input_folder
                , r'Export_For_Test_ONLYBrazos_ALLORDERS.shp')
        key_col_NHD = 1
        downstream_col_NHD = 6
        length_col_NHD = 5
        terminal_code_NHD = 0
        title_string = 'Brazos \nFull Res'
        driver_string = 'ESRI Shapefile'
        layer_string = 0

    elif LowerColorado_FULL_RES:
        nhd_conus_file_path = os.path.join(geo_input_folder
                , r'Export_For_Test_ONLYLowerColorado_ALLORDERS.shp')
        key_col_NHD = 1
        downstream_col_NHD = 6
        length_col_NHD = 5
        terminal_code_NHD = 0
        title_string = 'Lower Colorado\nFull Res'
        driver_string = 'ESRI Shapefile'
        layer_string = 0

    elif LowerColorado_CONCHOS_FULL_RES:
        nhd_conus_file_path = os.path.join(geo_input_folder
                , r'Export_For_Test_ONLYLowerColorado_CONCHOS_ALLORDERS.shp')
        key_col_NHD = 1
        downstream_col_NHD = 6
        length_col_NHD = 5
        terminal_code_NHD = 0
        title_string = 'Conchos sub-basin\nLower Colorado Full Res '
        driver_string = 'ESRI Shapefile'
        layer_string = 0

    elif Mainstems_CONUS:
        nhd_conus_file_path = os.path.join(geo_input_folder
                , r'downstream_reaches_v1_GCS.shp')
        key_col_NHD = 0
        downstream_col_NHD = 2
        length_col_NHD = 10
        terminal_code_NHD = 0
        title_string = 'CONUS "Mainstem"'
        driver_string = 'ESRI Shapefile'
        layer_string = 0

    elif CONUS_ge5:
        nhd_conus_file_path = os.path.join(geo_input_folder
                , r'NHD_Conus_Channels.shp')
        key_col_NHD = 1
        downstream_col_NHD = 6
        length_col_NHD = 5
        terminal_code_NHD = 0
        title_string = 'CONUS Order 5 and Greater '
        driver_string = 'ESRI Shapefile'
        layer_string = 0

    elif CONUS_FULL_RES:
        nhd_conus_file_path = '../../../../../../GISTemp/nwm_v12.gdb'
        key_col_NHD = 0
        downstream_col_NHD = 5
        length_col_NHD = 4
        terminal_code_NHD = 0
        title_string = 'CONUS Full Resolution NWM v1.2'
        driver_string = 'FileGDB'
        layer_string = 'channels_nwm_v12_routeLink'

    elif Brazos_LowerColorado_ge5:
        nhd_conus_file_path = os.path.join(geo_input_folder
                , r'NHD_BrazosLowerColorado_Channels.shp')
        key_col_NHD = 2
        downstream_col_NHD = 7
        length_col_NHD = 6
        terminal_code_NHD = 0
        title_string = 'Brazos + Lower Colorado\nNHD stream orders 5 and greater\n'
        title_string = 'CONUS Order 5 and Greater '
        driver_string = 'ESRI Shapefile'
        layer_string = 0

    if verbose: print(title_string)
    if debuglevel <= -1: print(f'reading -- dataset: {nhd_conus_file_path}; layer: {layer_string}; fiona driver: {driver_string}')
    nhd_conus = gpd.read_file(nhd_conus_file_path, driver=driver_string, layer=layer_string)
    if debuglevel <= -1: print(nhd_conus.head())
    nhd_conus_rows = nhd_conus.to_numpy()

    # Kick off recursive call for all connections and keys
    (connections_NHD) = networkbuilder.get_down_connections(
                 rows = nhd_conus_rows
                 , key_col = key_col_NHD
                 , downstream_col = downstream_col_NHD
                 , length_col = length_col_NHD
                 , verbose = verbose
                 , debuglevel = debuglevel)

    (connections_NHD, all_keys_NHD, ref_keys_NHD, headwater_keys_NHD
     , terminal_keys_NHD
     , terminal_ref_keys_NHD
     , circular_keys_NHD) = networkbuilder.determine_keys(
                 connections = connections_NHD
                 , rows = nhd_conus_rows
                 , key_col = key_col_NHD
                 , downstream_col = downstream_col_NHD
                 , terminal_code = terminal_code_NHD
                 , verbose = verbose
                 , debuglevel = debuglevel)

    (junction_keys_NHD) = networkbuilder.get_up_connections(
                 connections = connections_NHD
                 , terminal_code = terminal_code_NHD
                 , headwater_keys = headwater_keys_NHD
                 , terminal_keys = terminal_keys_NHD
                 , verbose = verbose
                 , debuglevel = debuglevel)

    recursive_print.print_basic_network_info(
                 connections = connections_NHD
                 , headwater_keys = headwater_keys_NHD
                 , junction_keys = junction_keys_NHD
                 , terminal_keys = terminal_keys_NHD
                 , terminal_code = terminal_code_NHD
                 , verbose = True
                 )

    if 1 == 0: #THE RECURSIVE PRINT IS NOT A GOOD IDEA WITH LARGE NETWORKS!!!
        recursive_print.print_connections(
                    headwater_keys = headwater_keys_NHD
                    , down_connections = connections_NHD
                    , up_connections = connections_NHD
                    , terminal_code = terminal_code_NHD
                    , terminal_keys = terminal_keys_NHD
                    , terminal_ref_keys = terminal_ref_keys_NHD
                    , debuglevel = -2
                    )

    # # nhd_conus_networks = []
    # # conus_supernetwork = SuperNetwork()
    #
    # networkbuilder.buildNetwork(
    #             rows = nhd_conus_rows
    #             , networks = nhd_conus_networks
    #             , down_connections = connections_NHD
    #             , up_connections = connections_NHD
    #             , key_col = key_col_NHD
    #             , downstream_col = downstream_col_NHD
    #             , length_col = length_col_NHD
    #             , terminal_code = terminal_code_NHD
    #             , headwater_keys = headwater_keys_NHD
    #             , terminal_keys = terminal_keys_NHD
    #             , verbose = True
    #             , debuglevel = -1
    #             )

    # print([network.networkID for network in nhd_conus_networks])
    # print([network.reachCollection for network in nhd_conus_networks])
    # print((network.reachCollection for network in nhd_conus_networks))
    # print(network.reachCollection for network in nhd_conus_networks)
    # for network in nhd_conus_networks:
    #     print([reach.reachID for reach in network.reachCollection])

    # print([reach.order for reach in network.reachCollection
                                    # for network in nhd_conus_networks])
