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

def main():
    # find the path of the test scripts, several levels above the script path
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    test_folder = os.path.join(root, r'test')
    
    supernetworks = {}
    # NOT IN GIT REPO # supernetworks.update({'Brazos_FULL_RES':{}}) 
    # NOT IN GIT REPO # supernetworks.update({'LowerColorado_FULL_RES':{}}) 
    # NOT IN GIT REPO # supernetworks.update({'LowerColorado_CONCHOS_FULL_RES':{}}) 
    supernetworks.update({'Pocono_TEST1':{}})
    # supernetworks.update({'Mainstems_CONUS':{}})
    # REMOVED FROM GIT REPO USE Mainstems Instead # supernetworks.update({'CONUS_ge5':{}}) ##NHD CONUS order 5 and greater"""
    # supernetworks.update({'Brazos_LowerColorado_ge5':{}}) ##NHD Subset (Brazos/Lower Colorado)"""
    # supernetworks.update({'CONUS_Named_Streams':{}})
    # # NOT IN GIT REPO # supernetworks.update({'CONUS_FULL_RES_v12':{}}) 
    # supernetworks.update({'CONUS_FULL_RES_v20':{}}) # = False

    debuglevel = -3
    verbose = True

    for supernetwork in supernetworks:
        supernetworks[supernetwork] = nnu.set_supernetwork_data(
          supernetwork = supernetwork
            , geo_input_folder = os.path.join(test_folder, r'input', r'geo', r'Channels')
            , debuglevel = debuglevel
            , verbose = verbose
        )
        if debuglevel <= -1: 
            if verbose: print(f'\n\nSupernetwork:')
            print(f'{supernetwork}')
        if debuglevel <= -2: 
            if verbose: print(r'All items in the above supernetwork:')
            for k,v in supernetworks[supernetwork].items():
                 print(f"{{'{k}' : {v}}}")

        network_out_values = \
          nnu.get_nhd_connections(
            supernetworks[supernetwork]
            , debuglevel = debuglevel
            , verbose = verbose
        )

        recursive_print.print_basic_network_info(
          connections = network_out_values[0]
            , headwater_keys = network_out_values[3]
            , junction_keys = network_out_values[7]
            , terminal_keys = network_out_values[4]
            , terminal_code = supernetworks[supernetwork]['terminal_code']
            , verbose = verbose
        )

        if 1 == 1: #THE RECURSIVE PRINT IS NOT A GOOD IDEA WITH LARGE NETWORKS!!!
            recursive_print.print_connections(
                        headwater_keys = network_out_values[3]
                        , down_connections = network_out_values[0]
                        , up_connections = network_out_values[0]
                        , terminal_code = supernetworks[supernetwork]['terminal_code']
                        , terminal_keys = network_out_values[4]
                        , terminal_ref_keys = network_out_values[5]
                        , debuglevel = debuglevel
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
