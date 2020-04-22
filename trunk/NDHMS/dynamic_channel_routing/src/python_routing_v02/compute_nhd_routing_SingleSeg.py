#!/usr/bin/env python
# coding: utf-8


# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
## Parallel execution
import multiprocessing
import os
import sys
import time
import csv
import numpy as np

ENV_IS_CL = False
if ENV_IS_CL:
    root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL:
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')
    sys.path.append(r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS')
    sys.setrecursionlimit(4000)

## Muskingum Cunge
COMPILE = True
if COMPILE:
    try:
        import subprocess

        fortran_compile_call = []
        fortran_compile_call.append(r'f2py3')
        fortran_compile_call.append(r'-c')
        fortran_compile_call.append(r'MCsingleSegStime_f2py_NOLOOP.f90')
        fortran_compile_call.append(r'-m')
        fortran_compile_call.append(r'mc_sseg_stime_NOLOOP')
        subprocess.run(fortran_compile_call, cwd=r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS',
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        import mc_sseg_stime_NOLOOP as mc
    except Exception as e:
        print(e)
else:
    import mc_sseg_stime_NOLOOP as mc

connections = None
networks = None
flowdepthvel = None
WRITE_OUTPUT = False  # True

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru

if WRITE_OUTPUT:
    # write to file 
    def writetoFile(file, writeString):
        file.write(writeString + '\n')
        file.flush()
        os.fsync(file.fileno())


def compute_network(
        terminal_segment=None
        , network=None
        , supernetwork_data=None
        # , connections = None
        , verbose=False
        , debuglevel=0
):
    global connections
    global flowdepthvel
    # global flowarray

    # = {connection:{'flow':{'prev':-999, 'curr':-999}
    #                            , 'depth':{'prev':-999, 'curr':-999}
    #                            , 'vel':{'prev':-999, 'curr':-999}} for connection in connections}

    # print(tuple(([x for x in network.keys()][i], [x for x in network.values()][i]) for i in range(len(network))))

    # if verbose: print(f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_order']}")

    ordered_reaches = {}
    for head_segment, reach in network['reaches'].items():
        if reach['seqorder'] not in ordered_reaches:
            ordered_reaches.update({reach['seqorder']: []})  # TODO: Should this be a set/dictionary?
        ordered_reaches[reach['seqorder']].append([head_segment
                                                      , reach
                                                   ])

    # initialize flowdepthvel dict
    nts = 50  # one timestep
    # nts = 1440 # number fof timestep = 1140 * 60(model timestep) = 86400 = day

    for ts in range(0, nts):
        # print(f'timestep: {ts}\n')

        for x in range(network['maximum_reach_seqorder'], -1, -1):
            for head_segment, reach in ordered_reaches[x]:
                # print(f'{{{head_segment}}}:{reach}')          

                compute_mc_reach_up2down(
                    head_segment=head_segment
                    , reach=reach
                    # , connections = connections
                    , supernetwork_data=supernetwork_data
                    , ts=ts
                    , verbose=verbose
                    , debuglevel=debuglevel
                )
                # print(f'{head_segment} {flowdepthvel[head_segment]}')          

    for x in range(network['maximum_reach_seqorder'], -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            # print(f'{{{head_segment}}}:{reach}')

            printarray(
                head_segment=head_segment
                , reach=reach
                , supernetwork_data=supernetwork_data
                , verbose=verbose
                , debuglevel=debuglevel
            )


# ### Psuedocode
# 


# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
        head_segment=None
        , reach=None
        # , connections = None
        , supernetwork_data=None
        , ts=0
        , verbose=False
        , debuglevel=0
):
    global connections
    global flowdepthvel
    # global network

    # if verbose: print(f"\nreach: {head_segment}")
    # if verbose: print(f"(reach: {reach})")
    # if verbose: print(f"(n_segs: {len(reach['segments'])})")
    if verbose: print(f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})")

    if WRITE_OUTPUT:
        filename = f'../../test/output/text/{head_segment}_{ts}.csv'
        file = open(filename, 'w+')
        writeString = f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])}  isterminal: {reach['upstream_reaches'] == {supernetwork_data['terminal_code']}} )  reach tail: {reach['reach_tail']}  upstream seg : "

    # upstream flow per reach
    qup_tmp = 0
    # import pdb; pdb.set_trace()
    if reach['upstream_reaches'] == {supernetwork_data['terminal_code']}:  # Headwaters
        qup_tmp = 0.0  # no flows
    else:  # Loop over upstream reaches
        # for us in reach['upstream_reaches']:
        for us in connections[reach['reach_head']]['upstreams']:
            if WRITE_OUTPUT: writeString = writeString + f" {us} "
            #if ts > 0:
            qup_tmp += flowdepthvel[us]['flowval'][-1]
            #else:
            #    qup_tmp += 0
    # print(qup_tmp)
    # writetoFile(file, writeString)

    qup = qup_tmp
    quc = qup

    current_segment = reach['reach_head']

    if WRITE_OUTPUT:
        writeString = writeString + f' timestep: {ts} cur : {current_segment}  upstream flow: {qup_tmp}'
        writetoFile(file, writeString)
        writeString = f"  , , , , , , "
        writetoFile(file, writeString)

    next_segment = connections[current_segment]['downstream']
    # print(f'{current_segment}==>{next_segment} conections:{ncomp} timestep:{ts}')
    # input flow to upstream reach of current reach     
    while True:

        # for now treating as constant per reach 
        dt = 60.0
        # import pdb; pdb.set_trace()
        bw = connections[current_segment]['data'][supernetwork_data['bottomwidth_col']]
        tw = 0.01 * bw
        twcc = tw
        dx = connections[current_segment]['data'][supernetwork_data['length_col']]
        bw = connections[current_segment]['data'][supernetwork_data['bottomwidth_col']]
        n_manning = connections[current_segment]['data'][supernetwork_data['manningn_col']]
        n_manning_cc = connections[current_segment]['data'][supernetwork_data['manningn_col']]
        cs = connections[current_segment]['data'][supernetwork_data['ChSlp_col']]
        s0 = connections[current_segment]['data'][supernetwork_data['slope_col']]

        # add some flow
        multiplier = 1
        if ts > 25:
            multiplier = 0
        qlat = 10.0 * (multiplier) # (ts + 1) * 10.0  # lateral flow per segment
        flowdepthvel[current_segment]['qlatval'].append(qlat)

        # print (f'counter = {i}')
        # if current_segment == 5559368 or i == 100:
        #    import pdb; pdb.set_trace()

        if ts > 0:
            qdp = flowdepthvel[current_segment]['flowval'][-1]
            velp = flowdepthvel[current_segment]['velval'][-1]
            depthp = flowdepthvel[current_segment]['depthval'][-1]
        else:
            qdp = 0
            velp = 0
            depthp = 0


        # run M-C model
        qdc, velc, depthc = singlesegment(
            dt=dt
            , qup=qup
            , quc=quc
            , qdp=qdp
            , qlat=qlat
            , dx=connections[current_segment]['data'][supernetwork_data['length_col']]
            , bw=connections[current_segment]['data'][supernetwork_data['bottomwidth_col']]
            , tw=tw
            , twcc=twcc
            , n_manning=connections[current_segment]['data'][supernetwork_data['manningn_col']]
            , n_manning_cc=connections[current_segment]['data'][supernetwork_data['manningn_col']]
            , cs=connections[current_segment]['data'][supernetwork_data['ChSlp_col']]
            , s0=connections[current_segment]['data'][supernetwork_data['slope_col']]
            , velp=velp
            , depthp=depthp
        )
        # print(qdc, velc, depthc)
        # print(qdc_expected, velc_expected, depthc_expected)

        if WRITE_OUTPUT:
            writeString = writeString + f',  {qdc},  {depthc},  {velc} '
            writetoFile(file, writeString)
        # for next segment qup / quc use the previous flow values
        if ts > 0:
            qup = flowdepthvel[current_segment]['flowval'][-1]  # input for next segment
        else:
            qup = 0
        quc = qup  # input for next segment
        flowdepthvel[current_segment]['flowval'].append(qdc)
        flowdepthvel[current_segment]['depthval'].append(depthc)
        flowdepthvel[current_segment]['velval'].append(velc)
        flowdepthvel[current_segment]['time'].append(ts)
        if current_segment == reach['reach_tail']:
            if verbose: print(f'{current_segment} (tail)')
            break
        if verbose: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream']
        # end loop initialized the MC vars 
        # print (flowdepthvel[current_segment]['flowval'])


def printarray(head_segment=None
               , reach=None
               , supernetwork_data=None
               , verbose=False
               , debuglevel=0
               ):
    global connections
    global flowdepthvel

    WRITE_OUTPUT = True
    header = [['time', 'qlat', 'q', 'd', 'v']]
    current_segment = reach['reach_head']
    next_segment = connections[current_segment]['downstream']
    # print(f'{current_segment}==>{next_segment} conections:{ncomp} timestep:{ts}')
    # input flow to upstream reach of current reach
    while True:

        if WRITE_OUTPUT:
            filename = f'../../test/output/text/{current_segment}.csv'
            with open(filename, 'w+') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_ALL)
                csvwriter.writerows(header)
                csvwriter.writerows(zip(flowdepthvel[current_segment]['time'],
                                        flowdepthvel[current_segment]['qlatval'],
                                        flowdepthvel[current_segment]['flowval'],
                                        flowdepthvel[current_segment]['depthval'],
                                        flowdepthvel[current_segment]['velval']))

        if current_segment == reach['reach_tail']:
            if verbose: print(f'{current_segment} (tail)')
            break
        if verbose: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream']


def singlesegment(
        dt  # dt
        , qup=None  # qup
        , quc=None  # quc
        , qdp=None  # qdp
        , qlat=None  # ql
        , dx=None  # dx
        , bw=None  # bw
        , tw=None  # tw
        , twcc=None  # twcc
        , n_manning=None  #
        , n_manning_cc=None  # ncc
        , cs=None  # cs
        , s0=None  # s0
        , velp=None  # velocity at previous time step
        , depthp=None  # depth at previous time step
):
    # call Fortran routine
    return mc.muskingcungenwm(
        dt, qup, quc, qdp, qlat, dx, bw, tw, twcc
        , n_manning, n_manning_cc, cs, s0, velp, depthp
    )
    # return qdc, vel, depth


def main():
    global connections
    global networks
    global flowdepthvel

    verbose = True
    debuglevel = 0
    showtiming = True

    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo')

    # TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_ge5'
    supernetwork = 'Pocono_TEST1'
    """##NHD CONUS order 5 and greater"""
    # supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    # supernetwork = 'Mainstems_CONUS'
    # supernetwork = 'CONUS_FULL_RES_v20'
    # supernetwork = 'CONUS_Named_Streams' #create a subset of the full resolution by reading the GNIS field
    # supernetwork = 'CONUS_Named_combined' #process the Named streams through the Full-Res paths to join the many hanging reaches

    if verbose: print('creating supernetwork connections set')
    if showtiming: start_time = time.time()
    # STEP 1
    supernetwork_data, supernetwork_values = nnu.set_networks(
        supernetwork=supernetwork
        , geo_input_folder=geo_input_folder
        , verbose=False
        # , verbose = verbose
        , debuglevel=debuglevel
    )
    if verbose: print('supernetwork connections set complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    # STEP 2
    if showtiming: start_time = time.time()
    if verbose: print('organizing connections into reaches ...')
    networks = nru.compose_reaches(
        supernetwork_values
        , verbose=False
        # , verbose = verbose
        , debuglevel=debuglevel
        , showtiming=showtiming
    )
    if verbose: print('reach organization complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    if showtiming: start_time = time.time()
    connections = supernetwork_values[0]

    flowdepthvel = {connection: {'qlatval': []
        , 'time': []
        , 'flowval': []
        , 'depthval': []
        , 'velval': []} for connection in connections}

    parallelcompute = True
    if not parallelcompute:
        if verbose: print('executing computation on ordered reaches ...')

        for terminal_segment, network in networks.items():
            compute_network(
                terminal_segment=terminal_segment
                , network=network
                , supernetwork_data=supernetwork_data
                # , connections = connections
                , verbose=False
                # , verbose = verbose
                , debuglevel=debuglevel
            )
            print(f'{terminal_segment}')
            if showtiming: print("... in %s seconds." % (time.time() - start_time))

    else:
        if verbose: print(f'executing parallel computation on ordered reaches .... ')
        # for terminal_segment, network in networks.items():
        #    print(terminal_segment, network)
        # print(tuple(([x for x in networks.keys()][i], [x for x in networks.values()][i]) for i in range(len(networks))))
        nslist = ([terminal_segment
                      , network
                      , supernetwork_data  # TODO: This should probably be global...
                      , False
                      , debuglevel]
                  for terminal_segment, network in networks.items())
        with multiprocessing.Pool() as pool:
            results = pool.starmap(compute_network, nslist)

    if verbose: print('ordered reach computation complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))


if __name__ == '__main__':
    main()
