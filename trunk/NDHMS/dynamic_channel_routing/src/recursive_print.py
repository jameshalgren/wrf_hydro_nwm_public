
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
        if verbose: print('''Provide headwater_keys, down_connections, and a terminal code
to print Downstream Connections.

Provide terminal_keys, up_connections, down_connections, and a terminal code
to print Upstream Connections.''')
