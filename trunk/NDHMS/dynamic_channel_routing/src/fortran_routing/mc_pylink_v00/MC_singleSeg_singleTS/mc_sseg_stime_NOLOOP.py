import numpy as np
debuglevel = 0
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
        if debuglevel <= -2:
            subprocess.run(fortran_compile_call)
        else:
            subprocess.run(fortran_compile_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        import mc_sseg_stime_NOLOOP as mc
    except Exception as e:
        print (e)
else:
    import mc_sseg_stime_NOLOOP as mc

# # Method 1
# Python: time loop; segment loop; constant channel variables are passed to Fortran
# Fortran: Take constant variable values and then run MC for a single segment

def compute_mc_up2down_ReachbySegment():
    '''HANDLE LOOPING, 
        Then call single segment routine for each segment'''
    pass

def singlesegment(
        dt # dt
        , qup = None # qup
        , quc = None # quc
        , qdp = None # qdp
        , qlat = None # ql
        , dx = None # dx
        , bw = None # bw
        , tw = None # tw
        , twcc = None # twcc
        , n_manning = None #
        , n_manning_cc = None # ncc
        , cs = None # cs
        , s0 = None # s0
        , velp = None # velocity at previous time step
        , depthp = None # depth at previous time step
    ):

    # call Fortran routine
    return mc.muskingcungenwm(
        dt, qup, quc, qdp, qlat, dx, bw, tw, twcc
        ,n_manning, n_manning_cc, cs, s0, velp, depthp
    )
    #return qdc, vel, depth
        
def main ():

    dt = 60.0
    dx = 1800.0
    bw = 112.0
    tw = 448.0
    twcc = 623.5999755859375
    n_manning = 0.02800000086426735
    n_manning_cc = 0.03136000037193298
    cs = 1.399999976158142
    s0 = 0.0017999999690800905
    qlat = 40.0
    qup = 0.04598825052380562
    quc = 0.04598825052380562
    qdp = 0.21487340331077576
    velp = 0.070480190217494964
    depthp = 0.010033470578491688

    qdc_expected = 0.7570106983184814
    velc_expected = 0.12373604625463486
    depthc_expected = 0.02334451675415039

    #run M-C model
    qdc, velc, depthc = singlesegment(
        dt = dt
        , qup = qup
        , quc = quc
        , qdp = qdp
        , qlat = qlat
        , dx = dx
        , bw = bw
        , tw = tw
        , twcc = twcc
        , n_manning = n_manning
        , n_manning_cc = n_manning_cc
        , cs = cs
        , s0 = s0
        , velp = velp
        , depthp = depthp
    )
    print(qdc, velc, depthc)
    print(qdc_expected, velc_expected, depthc_expected)
            
if __name__ == '__main__':
    main()
