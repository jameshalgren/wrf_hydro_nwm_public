import networkbuilder as networkbuilder
import recursive_print
import geopandas as gpd
import pandas as pd
import zipfile
import xarray as xr

def do_network(
        geo_file_path = None
        , title_string = None
        , layer_string = None
        , driver_string = None
        , key_col = None
        , downstream_col = None
        , length_col = None
        , terminal_code = None
        , mask_file_path = None
        , mask_driver_string = None
        , mask_layer_string = None
        , mask_key_col = None
        , verbose = False
        , debuglevel = 0
        ):

    if verbose: print(title_string)
    geo_file_rows = get_geo_file_table_rows(
        geo_file_path = geo_file_path
        , layer_string = layer_string
        , driver_string = driver_string
        , verbose = verbose
        , debuglevel = debuglevel
    )

    if debuglevel <= -1: print(f'MASK: {mask_file_path}')
    if mask_file_path:
        mask_file_rows = get_geo_file_table_rows(
            geo_file_path = mask_file_path
            , layer_string = mask_layer_string
            , driver_string = mask_driver_string
            , verbose = verbose
            , debuglevel = debuglevel
            )
        #TODO: make mask dict with additional attributes, e.g., names
        mask_set = {row[mask_key_col] for row in mask_file_rows}

        return build_connections_object(
            geo_file_rows = geo_file_rows
            , mask_set = mask_set
            , key_col = key_col
            , downstream_col = downstream_col
            , length_col = length_col
            , terminal_code = terminal_code
            , verbose = verbose
            , debuglevel = debuglevel
            )

    else:
        return build_connections_object(
            geo_file_rows = geo_file_rows
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


def get_geo_file_table_rows(
        geo_file_path = None
        , layer_string = None
        , driver_string = None
        , verbose = False
        , debuglevel = 0
        ):

    # NOTE: these methods can lose the "connections" and "rows" arguments when
    # implemented as class methods where those arrays are members of the class.
    if driver_string == 'NetCDF': # Use Xarray to read a netcdf table
        if debuglevel <= -1: print(f'reading -- dataset: {geo_file_path}; layer: {layer_string}; driver: {driver_string}')
        geo_file = xr.open_dataset(geo_file_path)
        geo_file_rows = (geo_file.to_dataframe()).values
        # The xarray method for NetCDFs was implemented after the geopandas method for 
        # GIS source files. It's possible (probable?) that we are doing something 
        # inefficient by converting away from the Pandas dataframe.
        # TODO: Check the optimal use of the Pandas dataframe
    elif driver_string == 'zip': # Use Pandas to read zipped csv
        if debuglevel <= -1: print(f'reading -- dataset: {geo_file_path}; layer: {layer_string}; driver: {driver_string}')
        with zipfile.ZipFile(geo_file_path, 'r') as zcsv:
            with zcsv.open(layer_string) as csv:
                geo_file = pd.read_csv(csv)
        geo_file_rows = geo_file.to_numpy()
    else: # Read Shapefiles, Geodatabases with Geopandas/fiona
        if debuglevel <= -1: print(f'reading -- dataset: {geo_file_path}; layer: {layer_string}; fiona driver: {driver_string}')
        geo_file = gpd.read_file(geo_file_path, driver=driver_string, layer=layer_string)
        geo_file_rows = geo_file.to_numpy()
        if debuglevel <= -2: geo_file.plot() 
    if debuglevel <= -1: print(geo_file.head()) # Preview the first 5 lines of the loaded data

    return geo_file_rows

def build_connections_object(
        geo_file_rows = None
        , mask_set = None
        , key_col = None
        , downstream_col = None
        , length_col = None
        , terminal_code = None
        , verbose = False
        , debuglevel = 0
        ):
    (connections) = networkbuilder.get_down_connections(
                    rows = geo_file_rows
                    , mask_set = mask_set
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

