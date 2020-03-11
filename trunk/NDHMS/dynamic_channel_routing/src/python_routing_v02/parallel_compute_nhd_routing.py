## Basic imports
import sys
import os
import time

## Parallel execution
import multiprocessing

# WARNING: These global declarations cause the parallel implementation to 
# crash when executed on Windows
connections = None
networks = None
flowdepthvel = None

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
    sys.setrecursionlimit(4000)
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.path.append(os.path.join(root, r'src', r'python_framework'))
    fortran_source_dir = os.path.join(root, r'src', r'fortran_routing', r'mc_pylink_v00', r'MC_singleCH_singleTS')
    sys.path.append(fortran_source_dir)
    try:
        import mc_sc_stime as mc
    except:
        import subprocess
        fortran_compile_call = []
        fortran_compile_call.append(r'f2py3')
        fortran_compile_call.append(r'-c')
        fortran_compile_call.append(r'varSingleChStime_f2py.f90')
        fortran_compile_call.append(r'MCsingleChStime_f2py_clean.f90')
        fortran_compile_call.append(r'-m')
        fortran_compile_call.append(r'mc_sc_stime')
        subprocess.run(fortran_compile_call, cwd=fortran_source_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            import mc_sc_stime as mc
        except Exception as e:
            print (e)

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru

## Muskingum Cunge
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
    global flowdepthvel 

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

def compute_network_parallel(
        large_networks = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    global connections
    global flowdepthvel 

    overall_max = -1
    for terminal_segment, network in large_networks.items():
        overall_max = max(network['maximum_reach_seqorder'], overall_max)

    if verbose: print(f"Executing simulation for all large network beginning with maximum order {overall_max}")

    nts = 1440 # number fof timestep
    
    for ts in range (0,nts):
        #print(f'timestep: {ts}\n')

    #if 1 == 1:
        ordered_reaches = {}
        for terminal_segment, network in large_networks.items():
            for head_segment, reach in network['reaches'].items():
                if reach['seqorder'] not in ordered_reaches:
                    ordered_reaches.update({reach['seqorder']:[]}) #TODO: Should this be a set/dictionary?
                ordered_reaches[reach['seqorder']].append([head_segment
                          , reach
                          , supernetwork_data
                          , verbose
                          , debuglevel])
        with multiprocessing.Pool() as netpool:
            for x in range(overall_max,-1,-1):
                nslist3 = ordered_reaches[x]
                if verbose: print(f"Executing simulation for {len(nslist3)} large network reaches of order {x}")
                results = netpool.starmap(compute_mc_reach_up2down, nslist3)

    if verbose: print(f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_reach_seqorder']}")



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
#                 With the populated arrays, execute MC for the reach
# ```     
#         

def read_segments():
    pass

def prepare_segments():
    pass

def handle_junctions():
    pass

def get_upstream_inflow():
    pass

def get_lateral_inflow():
    pass

def compute_junction_downstream():
    pass

#TODO: generalize with a direction flag
def compute_mc_reach_up2down(
        head_segment = None
        , reach = None
        #, connections = None
        , supernetwork_data = None
        , ts = 0
        , verbose = False
        , debuglevel = 0
        ):
    
    global connections
    global flowdepthvel
    # global network
    
    if verbose: print(f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})")
    
    ntim=2;       #the number of time steps necessary for variables passed to mc module to compute successfully
    nlinks=2;     #the number of links needed to define varialbe qd. ** nlinks is not used in fortran source code.

    mc.var.uslinkid=1
    mc.var.linkid=2
    ncomp0=1; mc.var.ncomp0=ncomp0  #the number of segments of a reach upstream of the current reach
    ncomp = len(reach['segments']) ;  mc.var.ncomp= ncomp  #the number of segments of the current reach 
    #mxseg=max(ncomp0,ncomp)    
    #MC model outputs
    mc.var.qd=np.zeros((ntim,ncomp,nlinks))  #will store MC output qdc (flow downstream current timestep) 
    mc.var.vela=np.zeros((ntim,ncomp)) 
    mc.var.deptha=np.zeros((ntim,ncomp))
    #lateral flow
    mc.var.qlat=np.zeros((ncomp))
     
    
    
    # upstream flow per reach
    qup_tmp = 0
    #import pdb; pdb.set_trace()
    if reach['upstream_reaches'] == {supernetwork_data['terminal_code']}: # Headwaters
        qup_tmp = (ts+1)*10.0 
    else: # Loop over upstream reaches
        #for us in reach['upstream_reaches']:
        for us in connections[reach['reach_head']]['upstreams']:
            #if us == 5507050 :
            #    import pdb; pdb.set_trace()
            #qup_tmp += flowdepthvel[network['reaches'][us]['reach_tail']]['flow']['curr']
            qup_tmp += flowdepthvel[us]['flow']['curr']
    
    flowdepthvel[reach['reach_head']]['flow']['curr'] = qup_tmp
    #print(qup_tmp)
            
    current_segment = reach['reach_head']
    next_segment = connections[current_segment]['downstream'] 
    #print(f'{current_segment}==>{next_segment} conections:{ncomp} timestep:{ts}')
    i = 0
    #input flow to upstream reach of current reach   
    mc.var.qd[1,i,0] = qup_tmp 
    
    while True:

        # for now treating as constant per reach 
        dt=300.0 ;      mc.var.dt= dt  #60.0;

        dx=connections[current_segment]['data'][supernetwork_data['length_col']] ;     mc.var.dx= dx  #20.0
        bw=connections[current_segment]['data'][supernetwork_data['bottomwidth_col']];       mc.var.bw= bw #50
        tw= 0.01*bw; mc.var.tw= tw
        twcc=tw;     mc.var.twcc=twcc
        n=connections[current_segment]['data'][supernetwork_data['manningn_col']] ;      mc.var.n=n #0.03
        ncc=n;       mc.var.ncc=ncc
        cs=connections[current_segment]['data'][supernetwork_data['ChSlp_col']] ;    mc.var.cs=cs #1.0e6
        so=connections[current_segment]['data'][supernetwork_data['slope_col']];    mc.var.so=so #0.002
        #ck= current_segment['data'][supernetwork['MusK_col']];   mc.var.ck = ck 
        #cx= current_segment['data'][supernetwork['MusX_col']];   mc.var.cx = cx
        #print (f'{current_segment}')
               
        flowdepthvel[current_segment]['qlat']['curr'] = (ts+1)*2.0      # lateral flow per segment
       
                      
        flowdepthvel[current_segment]['flow']['prev'] = flowdepthvel[current_segment]['flow']['curr']
        flowdepthvel[current_segment]['depth']['prev'] = flowdepthvel[current_segment]['depth']['curr']
        flowdepthvel[current_segment]['vel']['prev'] = flowdepthvel[current_segment]['vel']['curr']
        flowdepthvel[current_segment]['qlat']['prev'] = flowdepthvel[current_segment]['qlat']['curr']

        #print (f'counter = {i}')
        #if current_segment == 5559368 or i == 100:
        #    import pdb; pdb.set_trace()

        mc.var.qlat[i]= flowdepthvel[current_segment]['qlat']['curr']  # temporary assigned qlat 
        mc.var.qd[0,i,1]= flowdepthvel[current_segment]['flow']['prev']  # temporary assigned qd
        mc.var.vela[0,i] = flowdepthvel[current_segment]['vel']['prev']
        mc.var.deptha[0,i] = flowdepthvel[current_segment]['depth']['prev']
        
        i += 1
        
        if current_segment == reach['reach_tail']:
            if verbose: print(f'{current_segment} (tail)')
            break
        if verbose: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream'] 
        
    mc.mc.main()

    #print (f'{ts} end mc')
    
    current_segment = reach['reach_head']
    next_segment = connections[current_segment]['downstream'] 
    i = 0
    while True:
        flowdepthvel[current_segment]['flow']['curr'] = mc.var.qd[1,i,1] 
        flowdepthvel[current_segment]['depth']['curr'] = mc.var.deptha[1,i]
        flowdepthvel[current_segment]['vel']['curr'] = mc.var.vela[1,i]
        d = flowdepthvel[current_segment]['depth']['curr'] 
        q = flowdepthvel[current_segment]['flow']['curr']
        v = flowdepthvel[current_segment]['vel']['curr']
        ql = flowdepthvel[current_segment]['qlat']['curr']
        #print ( f'timestep: {ts} cur : {current_segment}  {q} {d} {v} {ql}')
        i += 1
        
        #print(f'timestep: {ts} {flowdepthvel[current_segment]}')
        #import pdb; pdb.set_trace()
        
            
        if current_segment == reach['reach_tail']:
            if verbose: print(f'{current_segment} (tail)')
            break
        if verbose: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream'] 

def main():

    global connections
    global networks
    global flowdepthvel

    verbose = True
    debuglevel = 0
    showtiming = True

    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo', r'Channels')

    #TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_ge5'
    """##NHD CONUS order 5 and greater"""
    # supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    supernetwork = 'Mainstems_CONUS'
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
    connections = supernetwork_values[0]
    #initialize flowdepthvel dict
    flowdepthvel = {connection:{'flow':{'prev':0, 'curr':0}
                                , 'depth':{'prev':-999, 'curr':0}
                                , 'vel':{'prev':0, 'curr':0}
                                , 'qlat':{'prev':0, 'curr':0}} for connection in connections} 
    parallel_split = -1 # -1 turns off the splitting and runs everything through the lumped execution

    #STEP 3a -- Large Networks
    if showtiming: start_time = time.time()
    if verbose: print(f'executing computation on ordered reaches for networks of order greater than {parallel_split} ...')
    large_networks = {terminal_segment: network \
                      for terminal_segment, network in networks.items() \
                      if network['maximum_reach_seqorder'] > parallel_split}
    # print(large_networks)
    compute_network_parallel(
        large_networks
        , supernetwork_data = supernetwork_data
        # , connections = connections
        , verbose = False
        # , verbose = verbose
        , debuglevel = debuglevel
    )
    if verbose: print(f'ordered reach computation complete for networks of order greater than {parallel_split}')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    ##STEP 3b -- Small Networks
    if parallel_split >= 0: print(r'DO NOT RUN WITH `parallel_split` >= 0')

if __name__ == '__main__':
    main()
