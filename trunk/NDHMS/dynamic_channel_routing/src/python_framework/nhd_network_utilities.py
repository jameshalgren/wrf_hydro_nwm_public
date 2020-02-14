import networkbuilder as networkbuilder
import recursive_print
import geopandas as gpd
import xarray as xr

def do_network(
        geofile_path = None
        , title_string = None
        , layer_string = None
        , driver_string = None
        , key_col = None
        , downstream_col = None
        , length_col = None
        , terminal_code = None
        , verbose = False
        , debuglevel = 0
        ):

    if verbose: print(title_string)
    geofile_rows = get_geofile_table_rows(
        geofile_path = geofile_path
        , layer_string = layer_string
        , driver_string = driver_string
        , verbose = verbose
        , debuglevel = debuglevel
    )

    return build_connections_object(
        geofile_rows = geofile_rows
        , key_col = key_col
        , downstream_col = downstream_col
        , length_col = length_col
        , terminal_code = terminal_code
        , verbose = verbose
        , debuglevel = debuglevel
        )

    # return connections, all_keys, ref_keys, headwater_keys \
    #     , terminal_keys, terminal_ref_keys \
    #     , circular_keys, junction_keys


def get_geofile_table_rows(
        geofile_path = None
        , layer_string = None
        , driver_string = None
        , verbose = False
        , debuglevel = 0
        ):

    # NOTE: these methods can lose the "connections" and "rows" arguments when
    # implemented as class methods where those arrays are members of the class.
    if driver_string == 'NetCDF': # Use Xarray to read a netcdf table
        geofile = xr.open_dataset(geofile_path)
        geofile_rows = (geofile.to_dataframe()).values
        # The xarray method for NetCDFs was implemented after the geopandas method for 
        # GIS source files. It's possible (probable?) that we are doing something 
        # inefficient by converting away from the Pandas dataframe.
        # TODO: Check the optimal use of the Pandas dataframe
        if debuglevel <= -1: print(f'reading -- dataset: {geofile_path}; layer: {layer_string}; driver: {driver_string}')
    else: # Read Shapefiles, Geodatabases with Geopandas/fiona
        if debuglevel <= -1: print(f'reading -- dataset: {geofile_path}; layer: {layer_string}; fiona driver: {driver_string}')
        geofile = gpd.read_file(geofile_path, driver=driver_string, layer=layer_string)
        geofile_rows = geofile.to_numpy()
        if debuglevel <= -2: geofile.plot() 
    if debuglevel <= -1: print(geofile.head()) 

    return geofile_rows

def build_connections_object(
        geofile_rows = None
        , key_col = None
        , downstream_col = None
        , length_col = None
        , terminal_code = None
        , verbose = False
        , debuglevel = 0
        ):
    (connections) = networkbuilder.get_down_connections(
                    rows = geofile_rows
                    , key_col = key_col
                    , downstream_col = downstream_col
                    , length_col = length_col
                    , verbose = verbose
                    , debuglevel = debuglevel)
    
    (all_keys, ref_keys, headwater_keys
        , terminal_keys
        , terminal_ref_keys
        , circular_keys) = networkbuilder.determine_keys(
                    connections = connections
                    , key_col = key_col
                    , downstream_col = downstream_col
                    , terminal_code = terminal_code
                    , verbose = verbose
                    , debuglevel = debuglevel)
    
    (junction_keys, visited_keys
     , visited_terminal_keys
     , junction_count) = networkbuilder.get_up_connections(
                    connections = connections
                    , terminal_code = terminal_code
                    , headwater_keys = headwater_keys
                    , terminal_keys = terminal_keys
                    , verbose = verbose
                    , debuglevel = debuglevel)
    return connections, all_keys, ref_keys, headwater_keys \
        , terminal_keys, terminal_ref_keys \
        , circular_keys, junction_keys \
        , visited_keys, visited_terminal_keys \
        , junction_count

