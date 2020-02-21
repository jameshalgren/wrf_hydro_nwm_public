# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
import nhd_network_utilities as nnu
import recursive_print
import os

# NOTE: these methods can lose the "connections" and "rows" arguments when
# implemented as class methods where those arrays are members of the class.

def set_supernetwork_data(
    supernetwork = ''
    # Note: supernetworks may contain:
    # Brazos_FULL_RES
    # LowerColorado_FULL_RES
    # LowerColorado_CONCHOS_FULL_RES
    # Mainstems_CONUS
    # CONUS_ge5
    # Brazos_LowerColorado_ge5
    # CONUS_Named_Streams
    # CONUS_FULL_RES_v12
    # CONUS_FULL_RES_v20
    , geo_input_folder = None
    , verbose = True
    , debuglevel = 0
    ):

    # The following datasets are extracts from the feature datasets available
    # from https://www.nohrsc.noaa.gov/pub/staff/keicher/NWM_live/web/data_tools/
    # the CONUS_ge5 and Brazos_LowerColorado_ge5 datasets are included
    # in the github test folder
    if supernetwork == 'Brazos_FULL_RES':
        return {
            'geofile_path': os.path.join(geo_input_folder
                    , r'Export_For_Test_ONLYBrazos_ALLORDERS.shp')
            , 'key_col' : 1
            , 'downstream_col' : 6
            , 'length_col' : 5
            , 'terminal_code' : 0
            , 'title_string' : 'Brazos \nFull Res'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'LowerColorado_FULL_RES':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'Export_For_Test_ONLYLowerColorado_ALLORDERS.shp')
            , 'key_col' : 1
            , 'downstream_col' : 6
            , 'length_col' : 5
            , 'terminal_code' : 0
            , 'title_string' : 'Lower Colorado\nFull Res'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'LowerColorado_CONCHOS_FULL_RES':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'Export_For_Test_ONLYLowerColorado_CONCHOS_ALLORDERS.shp')
            , 'key_col' : 1
            , 'downstream_col' : 6
            , 'length_col' : 5
            , 'terminal_code' : 0
            , 'title_string' : 'Conchos sub-basin\nLower Colorado Full Res '
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'Brazos_LowerColorado_ge5':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'NHD_BrazosLowerColorado_Channels.shp')
            , 'key_col' : 2
            , 'downstream_col' : 7
            , 'length_col' : 6
            , 'terminal_code' : 0
            , 'title_string' : 'NHD Subset including Brazos + Lower Colorado\nNHD stream orders 5 and greater\n'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'Mainstems_CONUS':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'downstream_reaches_v1_GCS.shp')
            , 'key_col' : 0
            , 'downstream_col' : 2
            , 'length_col' : 10
            , 'terminal_code' : 0
            , 'title_string' : 'CONUS "Mainstem"'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'CONUS_ge5':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'NHD_Conus_Channels.shp')
            , 'key_col' : 1
            , 'downstream_col' : 6
            , 'length_col' : 5
            , 'terminal_code' : 0
            , 'title_string' : 'NHD CONUS Order 5 and Greater'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'CONUS_Named_Streams':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'channels_nwm_v12_routeLink_NamedOnly.shp')
            , 'key_col' : 0
            , 'downstream_col' : 5
            , 'length_col' : 4
            , 'terminal_code' : 0
            , 'title_string' : 'NHD v1.2 segments corresponding to NHD 2.0 GNIS labeled streams\n'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'CONUS_FULL_RES_v20':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'RouteLink_NWMv2.0_20190517_cheyenne_pull.nc')
            , 'key_col' : 0
            , 'downstream_col' : 2
            , 'length_col' : 8
            , 'terminal_code' : 0
            , 'title_string' : 'CONUS Full Resolution NWM v2.0'
            , 'driver_string' : 'NetCDF'
            , 'layer_string' : 0
          }

    elif supernetwork == 'CONUS_FULL_RES_v12':
        return {
            'geofile_path' : os.path.join(geo_input_folder
                    , r'channels_nwm_v12_routeLink_all.shp')
            , 'key_col' : 0
            , 'downstream_col' : 5
            , 'length_col' : 4
            , 'terminal_code' : 0
            , 'title_string' : 'CONUS Full Resolution NWM v1.2'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

def get_nhd_connections(
    supernetwork = {}
    , debuglevel = 0
    , verbose = False
    ):
    return nnu.do_network(
        geofile_path = supernetwork['geofile_path']
          , key_col = supernetwork['key_col']
          , downstream_col = supernetwork['downstream_col']
          , length_col = supernetwork['length_col']
          , terminal_code = supernetwork['terminal_code']
          , title_string = supernetwork['title_string']
          , driver_string = supernetwork['driver_string']
          , layer_string = supernetwork['layer_string']
          , debuglevel = debuglevel
          , verbose = verbose
        )

def main():
    # find the path of the test scripts, several levels above the script path
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    test_folder = os.path.join(root, r'test')
    
    supernetworks = {}
    # NOT IN GIT REPO # supernetworks.update({'Brazos_FULL_RES':{}}) 
    # NOT IN GIT REPO # supernetworks.update({'LowerColorado_FULL_RES':{}}) 
    # NOT IN GIT REPO # supernetworks.update({'LowerColorado_CONCHOS_FULL_RES':{}}) 
    # NOT IN GIT REPO # supernetworks.update({'Mainstems_CONUS':{}})
    supernetworks.update({'CONUS_ge5':{}}) ##NHD CONUS order 5 and greater"""
    supernetworks.update({'Brazos_LowerColorado_ge5':{}}) ##NHD Subset (Brazos/Lower Colorado)"""
    # NOT IN GIT REPO # supernetworks.update({'CONUS_Named_Streams':{}})
    # NOT IN GIT REPO # supernetworks.update({'CONUS_FULL_RES_v12':{}}) 
    supernetworks.update({'CONUS_FULL_RES_v20':{}}) # = False


    for supernetwork in supernetworks:
        supernetworks[supernetwork] = set_supernetwork_data(
          supernetwork = supernetwork
          , geo_input_folder = os.path.join(test_folder, r'input', r'geo', r'Channels')
        )
        network_out_values = \
            get_nhd_connections(
              supernetworks[supernetwork]
              , debuglevel = -1
              , verbose = True
        )

        recursive_print.print_basic_network_info(
          connections = network_out_values[0]
            , headwater_keys = network_out_values[3]
            , junction_keys = network_out_values[7]
            , terminal_keys = network_out_values[4]
            , terminal_code = supernetworks[supernetwork]['terminal_code']
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
