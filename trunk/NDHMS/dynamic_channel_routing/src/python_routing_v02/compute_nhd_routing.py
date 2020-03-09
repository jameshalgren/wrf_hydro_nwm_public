## Basic imports
import sys
import os
import time

# WARNING: These global declarations cause the parallel implementation to 
# crash when executed on Windows
connections = None
networks = None

from sys import platform
if platform == "linux" or platform == "linux2":
    pass
elif platform == "darwin":
    pass
elif platform == "win32":
    print('The parallel version of compute_nhd_routing.py will not execute as currently')
    print('written due to the lack of a fork() capability in the windows OS.')
    print('For parallel execution, please us a *nix OS.')
    print('\nexiting...')
    sys.exit()
    # Some useful references:
    # https://stackoverflow.com/questions/985281/what-is-the-closest-thing-windows-has-to-fork/985525#985525
    # https://stackoverflow.com/questions/8220108/how-do-i-check-the-operating-system-in-python
    # https://stackoverflow.com/questions/6596617/python-multiprocess-diff-between-windows-and-linux

ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL: 
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')
    sys.path.append(r'../fortran_routing/mc_pylink_v00/MC_singleCH_singleTS')
    sys.setrecursionlimit(4000)

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru

## Muskingum Cunge
import mc_sc_stime as mc
import numpy as np

def compute_network(
        terminal_segment = None
        , network = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    global connections

    if verbose: print(f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_reach_seqorder']}")

    ordered_reaches = {}
    for head_segment, reach in network['reaches'].items():
        if reach['seqorder'] not in ordered_reaches:
            ordered_reaches.update({reach['seqorder']:[]}) #TODO: Should this be a set/dictionary?
        ordered_reaches[reach['seqorder']].append([head_segment
                  , reach
                  ])
    for x in range(network['maximum_reach_seqorder'],-1,-1):
        for head_segment, reach in ordered_reaches[x]:
            compute_reach_up2down(
                head_segment = head_segment
                , reach = reach
                # , connections = connections
                , supernetwork_data = supernetwork_data
                , verbose = verbose
                , debuglevel = debuglevel
                )

#TODO: generalize with a direction flag
def compute_reach_up2down(
        head_segment = None
        , reach = None
        # , connections = None
        , supernetwork_data = None
        , verbose = False
        , debuglevel = 0
        ):
    global connections

    if debuglevel <= -1: print(f"\nreach: {head_segment} (seqorder: {reach['seqorder']} n_segs: {len(reach['segments'])})")
    current_segment = reach['reach_head']
    next_segment = connections[current_segment]['downstream'] 
    while True:
        # Thanks to SO post for a reminder of this "Loop-and-a-half" construct
        # https://stackoverflow.com/questions/1662161/is-there-a-do-until-in-python
        compute_segment(
            current_segment = current_segment
            , supernetwork_data = supernetwork_data
            , verbose = verbose
            , debuglevel = debuglevel
            )
        if current_segment == reach['reach_tail']:
            if debuglevel <= -1: print(f'{current_segment} (tail)')
            break
        if debuglevel <= -1: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream'] 

def compute_segment(
    current_segment = None
    , supernetwork_data = None
    # , connections = None
    , verbose = False
    , debuglevel = 0
    ):
    global connections
    
    if debuglevel <= -2:
        print(f'data for segment: {current_segment}')
        for k, v in connections[current_segment].items():
            if type(v) is list: 
                print (f'{k}:')
                print (f"manningn: {v[supernetwork_data['manningn_col']]}")
                print (f"slope: {v[supernetwork_data['slope_col']]}")
                print (f"bottom width: {v[supernetwork_data['bottomwidth_col']]}")
                print (f"Muskingum K: {v[supernetwork_data['MusK_col']]}")
                print (f"Muskingum X: {v[supernetwork_data['MusX_col']]}")
                print (f"Channel Side Slope: {v[supernetwork_data['ChSlp_col']]}")
            else: print(f'{k}: {v}')

    #"compute MC"
    computeMC = True
    #"compute Dummy"
    computeDummy = False
    #"compute MSH"
    computeMSH= False

     
    if computeMC:
       compute_mc(connections[current_segment]
            , supernetwork_data = supernetwork_data
            , verbose = verbose
            , debuglevel = debuglevel)
        
    elif computeDummy:
        print ('computeDummy')
                   
    elif computeMSH:
        print ('computeMSH')                   
                 


# ### Psuedocode
# 
# ```
# Call Compute Network
#     for each reach in the network
#         Call compute reach
#             For each segment in the reach
#                 Import the mc object/module
#                 Call prepare segment
#                     Populate the Mc segment array with the individual segment properties
#                     Obtain and set any lateral inflows
#                     obtain and set the upstream and downstream (requrires logic to deal with junctions or boundaries at headwaters)
#             With the populated arrays, execute MC for the reach
# ```     
#         

def get_upstream_inflow():
    pass

def get_lateral_inflow():
    pass

def compute_junction_downstream():
    pass

def compute_mc(
        current_segment = None
        , supernetwork_data = None
        , verbose = False
        , debuglevel = 0
        ):
    #global connections

    
 
    #print (current_segment['data'][supernetwork_data['manningn_col']])
    #print (current_segment['data'][supernetwork_data['slope_col']])
    # print (current_segment['data'][supernetwork_data['bottomwidth_col']])
    #print (current_segment['data'][supernetwork_data['MusK_col']])
    #print (current_segment['data'][supernetwork_data['MusX_col']])
    #print (current_segment['data'][supernetwork_data['ChSlp_col']])

    #if verbose: 
              

    #print (f"manningn: {connections[current_segment]}")
    #print ('computeMC')

    # Thanks to SO post for a reminder of this "Loop-and-a-half" construct
    # https://stackoverflow.com/questions/1662161/is-there-a-do-until-in-python

    ntim=2;       #the number of time steps necessary for variables passed to mc module to compute successfully
    nlinks=2;     #the number of links needed to define varialbe qd. ** nlinks is not used in fortran source code.

    ncomp0=3; mc.var.ncomp0=ncomp0  #the number of segments of a reach upstream of the current reach
    ncomp=3; mc.var.ncomp=ncomp  #the number of segments of the current reach 
    mxseg=max(ncomp0,ncomp)
    mc.var.uslinkid=1
    mc.var.linkid=2

    #MC model outputs
    mc.var.qd=np.zeros((ntim,mxseg,nlinks))  #will store MC output qdc (flow downstream current timestep) 
    mc.var.vela=np.zeros((ntim,ncomp)) 
    mc.var.deptha=np.zeros((ntim,ncomp))

    #print
    #import pdb; pdb.set_trace()
    
    #lateral flow
    mc.var.qlat=np.zeros((ncomp))

    dt=30.0 ;      mc.var.dt= dt  #60.0;
    dx=current_segment['data'][supernetwork_data['length_col']] ;     mc.var.dx= dx  #20.0
    bw=current_segment['data'][supernetwork_data['bottomwidth_col']];       mc.var.bw= bw #50
    tw= 0.01*bw; mc.var.tw= tw
    twcc=tw;     mc.var.twcc=twcc
    n=current_segment['data'][supernetwork_data['manningn_col']] ;      mc.var.n=n #0.03
    ncc=n;       mc.var.ncc=ncc
    cs=current_segment['data'][supernetwork_data['ChSlp_col']] ;    mc.var.cs=cs #1.0e6
    so=current_segment['data'][supernetwork_data['slope_col']];    mc.var.so=so #0.002
    ck= current_segment['data'][supernetwork_data['MusK_col']];   mc.var.ck = ck 
    cx= current_segment['data'][supernetwork_data['MusX_col']];   mc.var.cx = cx

    #import pdb; pdb.set_trace()
    #run M-C model
    nts=10  #the number of timestep in simulation
    wnlinks=20 #the number of all links in simulation
    wmxseg=5  #max number of segments among all links

    #variable storing all outputs in time
    wqd= np.zeros((nts,wmxseg,wnlinks))   
    wvela= np.zeros((nts,wmxseg,wnlinks)) 
    wdeptha= np.zeros((nts,wmxseg,wnlinks))

    for k in range (0,nts):        
        #input lateral flow for current reach; input flow to upstream reach of current reach
        for i in range(0,ncomp):
            mc.var.qlat[i]= (k+1)*2.0
            mc.var.qd[1,i,0]= (k+1)*10.0

        mc.mc.main()

        for i in range(0,ncomp):
            #current link(=reach)
            #qd[k,i,j]: k=0/1: previous/current timestep; i: node ID; j=0/1: upstream/current reach
            mc.var.qd[0,i,1]= mc.var.qd[1,i,1]
            mc.var.vela[0,i]= mc.var.vela[1,i]
            mc.var.deptha[0,i]= mc.var.deptha[1,i]
            #upstream link(=reach)        
            mc.var.qd[0,i,0]=  mc.var.qd[1,i,0] 

        #output keeping
        j=1 #temporarily assigned link ID for the current reach
        for i in range(0,ncomp):
            wqd[k,i,j]= mc.var.qd[1,i,1]
            wvela[k,i,j]= mc.var.vela[1,i]
            wdeptha[k,i,j]= mc.var.deptha[1,i]

    #test output    
    j=1
    #print('MC_TEST_OUTPUT')
#     for k in range (0,nts):
#         for i in range(0,ncomp):
#             print(wqd[k,i,j])

#     for k in range (0,nts):
#         for i in range(0,ncomp):
#             print(wvela[k,i,j])

#     for k in range (0,nts):
#         for i in range(0,ncomp):
#             print(wdeptha[k,i,j])

def main():

    global connections
    global networks

    verbose = True
    debuglevel = 0
    showtiming = True

    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo', r'Channels')

    #TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_ge5'
    """##NHD CONUS order 5 and greater"""
    supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    # supernetwork = 'Mainstems_CONUS'
    # supernetwork = 'CONUS_FULL_RES_v20'
    # supernetwork = 'CONUS_Named_Streams' #create a subset of the full resolution by reading the GNIS field
    # supernetwork = 'CONUS_Named_combined' #process the Named streams through the Full-Res paths to join the many hanging reaches


    if verbose: print('creating supernetwork connections set')
    if showtiming: start_time = time.time()
    #STEP 1
    supernetwork_data, supernetwork_values = nnu.set_networks(
        supernetwork = supernetwork
        , geo_input_folder = geo_input_folder
        , verbose = False
        # , verbose = verbose
        , debuglevel = debuglevel
        )
    if verbose: print('supernetwork connections set complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 2
    if showtiming: start_time = time.time()
    if verbose: print('organizing connections into networks and reaches ...')
    networks = nru.compose_reaches(
        supernetwork_values
        , verbose = False
        # , verbose = verbose
        , debuglevel = debuglevel
        , showtiming = showtiming
        )
    if verbose: print('reach organization complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 3
    if showtiming: start_time = time.time()
    if verbose: print('executing computation on ordered reaches ...')
    connections = supernetwork_values[0]
    from itertools import islice
    def take(iterable, n):
        return list(islice(iterable, n))
    import pdb; pdb.set_trace()

    for terminal_segment, network in networks.items():
        compute_network(
            terminal_segment = terminal_segment
            , network = network
            , supernetwork_data = supernetwork_data
            # , connections = connections
            , verbose = False
            # , verbose = verbose
            , debuglevel = debuglevel
        )
    if verbose: print('ordered reach computation complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

if __name__ == '__main__':
    main()
