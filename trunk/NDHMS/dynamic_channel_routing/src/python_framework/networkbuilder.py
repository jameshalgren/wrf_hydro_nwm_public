import recursive_print

def get_down_connections(
    rows
    , key_col
    , downstream_col
    , length_col
    , mask_set = None
    , verbose = False
    , debuglevel = 0
    ):
    # TODO: Consider moving debug and verbose prints to the calling function
    if debuglevel <= -100: breakpoint()
    if verbose: print('down connections ...')
    # connections = {row[key_col]: { 'data': list(row)
    #                     , 'visited': False
    #                     , 'downstream': row[downstream_col]
    #                     , 'length': row[length_col]}
    #                             for row in rows}
    if mask_set: connections = {row[key_col]: { 'downstream': row[downstream_col]
                        , 'length': row[length_col]}
                                for row in rows 
                                if row[key_col] in mask_set}
    else: connections = {row[key_col]: { 'downstream': row[downstream_col]
                        , 'length': row[length_col]}
                                for row in rows}
    if debuglevel <= -1: print(f'found {len(connections.keys())} segments')
    if debuglevel <= -2: print(connections)
    if verbose: print('down_connections complete')

    return connections

def determine_keys(connections
            , key_col, downstream_col
            , terminal_code
            , verbose = False, debuglevel = 0):

    if verbose: print('ref_keys ...')
    # ref_keys = {row[downstream_col] for row in rows}
    ref_keys = {con['downstream'] for key, con in connections.items()}
    if debuglevel <= -1: print(f'found {len(ref_keys)} ref_keys')
    if debuglevel <= -2: print(ref_keys)
    if verbose: print('ref_keys complete')
    if debuglevel <= -2: print(ref_keys)
    if verbose: print('headwater_keys ...')
    # headwater_keys = {x for x in connections.keys() if x not in ref_keys}
    headwater_keys = {x for x in connections.keys() if x not in ref_keys}
    if debuglevel <= -1: print(f'found {len(headwater_keys)} headwater segments')
    if debuglevel <= -2: print(headwater_keys)
    if verbose: print('headwater_keys complete')

    # Get the downstream terminating nodes
    if verbose: print('terminal_keys ...')
    # Find the pointing-to keys not found in the key dataset.
    terminal_ref_keys = {x for x in ref_keys if x not in connections.keys()} 
    # import pdb; pdb.set_trace()

    # Then collect the keys associated with those 'pointing-tos'
    terminal_keys = set()
    # for row in rows:
    for key, con in connections.items():
        # curr_term_ref_key = row[downstream_col]
        curr_term_ref_key = con['downstream']
        if curr_term_ref_key in terminal_ref_keys:
            if curr_term_ref_key != terminal_code:
                if debuglevel <= -2:
                    # print(f'Non-standard terminal key {row[downstream_col]} found in segment {row[key_col]}')
                    print(f"Non-standard terminal key {con['downstream']} found in segment {key}")
            elif curr_term_ref_key == terminal_code:
                if debuglevel <= -2:
                    # print(f'Standard terminal key {row[downstream_col]} found in segment {row[key_col]}')
                    print(f"Standard terminal key {con['downstream']} found in segment {key}")
            # terminal_keys.add(row[key_col])
            terminal_keys.add(key)
    if debuglevel <= -1: print(f'found {len(terminal_keys)} terminal segments')
    if debuglevel <= -1: print(f'of those, {len([x for x in terminal_ref_keys if x != terminal_code])} had non-standard terminal keys')
    if debuglevel <= -2: print(terminal_keys)
    if verbose: print('terminal_keys complete')

    if verbose: print('circular_keys ...')
    circular_keys = set()
    for key, value in connections.items():
        try:
        #TODO: benchmark try/except vs. nested if statment on 'in' to handle terminal keys
        # e.g., "if key not in terminal_keys: ... etc.
            if connections[connections[key]['downstream']]['downstream'] == key:
                circular_keys.add(key)
            elif connections[connections[connections[key]['downstream']]['downstream']]['downstream'] == key:
                circular_keys.add(key)
            elif connections[connections[connections[connections[key]['downstream']]['downstream']]['downstream']]['downstream'] == key:
                circular_keys.add(key)
            elif connections[connections[connections[connections[connections[key]['downstream']]['downstream']]['downstream']]['downstream']]['downstream'] == key:
                circular_keys.add(key)
        except: pass

    if debuglevel <= -1: print(f'identified at least {len(circular_keys)} segments with circular references testing to the fourth level')
    if debuglevel <= -2: print(circular_keys)
    if verbose: print('circular_keys complete')


    return connections.keys(), ref_keys, headwater_keys \
                , terminal_keys, terminal_ref_keys, circular_keys

def get_up_connections(connections
                    , terminal_code
                    , headwater_keys
                    , terminal_keys
                    , verbose = False, debuglevel = 0):

    # Create inverse of connections looking upstream
    if verbose: print('identifying upstream connections and junctions ...')

    # Using Sets for Junction and Visited keys is REALLY, REALLY, REALLY, FAST!!!
    junction_keys = set()
    visited_keys = set()
    visited_terminal_keys = set()
    junction_count = 0
    for hkey in headwater_keys:
        # TODO: Create a dictionary key identifying relationship to the terminal segment. 

        # Start with the headwater keys and label the upstream connections
        # with the terminal_code...
        # connections[hkey].update({'upstreams' : [terminal_code]})
        connections[hkey].update({'upstreams' : {terminal_code}})
        visited_keys.add(hkey)
        # Then iterate through the list and search for the other values
        ukey = hkey
        # print(ukey, hkey)
        # print(visited_keys)
        # print(ukey not in terminal_keys)
        # print(ukey not in junction_keys)
        while (ukey not in terminal_keys) and (ukey not in junction_keys):
            dkey = connections[ukey]['downstream']
            if 'upstreams' not in connections[dkey]: # Check for key in dictionary https://able.bio/rhett/check-if-a-key-exists-in-a-python-dictionary--73iajoz
                # connections[dkey].update({'upstreams': []}) #TODO: Consider making this a set/hash for the unusual possibility of many upstream segments
                connections[dkey].update({'upstreams': set()}) #TODO: alternate Upstreams-as-a-set method
                # connections[dkey]['upstreams'].append(ukey)
                connections[dkey]['upstreams'].add(ukey)
                visited_keys.add(dkey)
            else:
                if terminal_code in connections[dkey]['upstreams']:
                    # If the downstream node here is labeled as a headwater (because it
                    # has an upstream set with the terminal code), it means
                    # that the network had a break and that the traversal has
                    # spanned the gap and the headwater is not actually not a terminating node.
                    # In that case, reset the node to be a blank list (or set, if using 
                    # that method), then proceed downstream.
                    # THIS IS A DANGEROUS STEP AND DESERVES ADDITIONAL REVIEW
                    # TO MAKE SURE IT IS DOING WHAT WE INTEND AS DESCRIBED ABOVE
                    # connections[dkey].update({'upstreams': []}) 
                    connections[dkey].update({'upstreams' : set()})

                # connections[dkey]['upstreams'].append(ukey)
                connections[dkey]['upstreams'].add(ukey)
                visited_keys.add(dkey)
                # print(dkey, connections[dkey]['upstreams'], visited_keys)
                if len(connections[dkey]['upstreams'])  == 2:
                    # breakpoint()
                    if dkey not in junction_keys:
                        junction_keys.add(dkey)
                        # print(len(junction_keys), dkey)
                    if debuglevel <= -2: print (f"Junction found above/into Segment {dkey} with upstream Segments {connections[dkey]['upstreams']}")
                    junction_count += 1
                elif len(connections[dkey]['upstreams']) > 2:
                    if dkey not in junction_keys:
                        #At this point, the logic does not allow for this to be a non-junction
                        #TODO: raise/handle error/warning
                        print('key error -- junction analysis has an undetermined anomaly!')
#                         print(dkey in visited_keys)
#                         for temp_ukey in connections[dkey]['upstreams']:
#                             print(temp_ukey, temp_ukey in visited_keys)
                    if debuglevel <= -2: print (f"revisited Junction above/into Segment {dkey} now with upstream Segments {connections[dkey]['upstreams']}")
                    junction_count += 1
            ukey = dkey
        if ukey in terminal_keys:
            visited_terminal_keys.add(ukey)

        #len(visited_keys) > 14266: breakpoint()

    if debuglevel <= -1: print(f'visited {len(visited_keys)} segments')
    if debuglevel <= -1: print(f'found {junction_count} junctions in {len(junction_keys)} junction nodes')
    if debuglevel <= -2: print(junction_keys)
    if debuglevel <= -2: print(connections)
    if verbose: print('up_connections complete')
    if verbose: print('')

    return junction_keys, visited_keys, visited_terminal_keys, junction_count

def main():
    """##TEST"""
    print("")
    print ('Executing Test')
    # Test data
    test_rows = [
        [50,178,51,0],
        [51,178,50,0],
        [60,178,61,0],
        [61,178,62,0],
        [62,178,60,0],
        [70,178,71,0],
        [71,178,72,0],
        [72,178,73,0],
        [73,178,70,0],
        [80,178,81,0],
        [81,178,82,0],
        [82,178,83,0],
        [83,178,84,0],
        [84,178,80,0],
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

    (test_connections) = get_down_connections(
                rows = test_rows
                , key_col = test_key_col
                , downstream_col = test_downstream_col
                , length_col = test_length_col
                , verbose = True
                , debuglevel = -2
                )

    (test_all_keys, test_ref_keys, test_headwater_keys
     , test_terminal_keys
     , test_terminal_ref_keys
     , test_circular_keys) = determine_keys(
                connections = test_connections
                , key_col = test_key_col
                , downstream_col = test_downstream_col
                , terminal_code = test_terminal_code
                , verbose = True
                , debuglevel = -2
                )

    test_junction_keys = get_up_connections(
                connections = test_connections
                , terminal_code = test_terminal_code
                , headwater_keys = test_headwater_keys
                , terminal_keys = test_terminal_keys
                , verbose = True
                , debuglevel = -2
                )

    recursive_print.print_connections(
                headwater_keys = test_headwater_keys
                , down_connections = test_connections
                , up_connections = test_connections
                , terminal_code = test_terminal_code
                , terminal_keys = test_terminal_keys
                , terminal_ref_keys = test_terminal_ref_keys
                , debuglevel = -2
                )

    recursive_print.print_basic_network_info(
                connections = test_connections
                , headwater_keys = test_headwater_keys
                , junction_keys = test_junction_keys
                , terminal_keys = test_terminal_keys
                , terminal_code = test_terminal_code
                , verbose = True
                , debuglevel = -2
                )

if __name__ == '__main__':
    main()
