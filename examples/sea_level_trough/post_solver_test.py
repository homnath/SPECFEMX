from tester import Tester, parse_args, gen_proc_list


def run_test(tester, verb, vars, nsteps):
    '''
    Runs test to check variable values held within ensight files
    :param tester: Instance of Tester class containing details on stable and trial data
    :type tester: Tester class
    :param verb: Level of verbosity (0,1,2) - increased value = more verbose
    :type verb: integer
    :param vars: list of variables to be tested. Full names or abbreviations (given in brackets) may be used.
                 Available options: 'displacement' (dis), 'strain' (eps), 'gravity_potential' (gpot), 'oceanf' (of), 'iceload_phi' (lp),
                'iceload_sl' (lsl), 'iceload_u' (lu), 'ice' (ice), 'icerate' (irate), 'sea_level' (sl), 'gravity_acceleration' (grav)
    :type vars: list of strings
    :param nsteps: Number of timesteps to be tested
    :type nsteps: integer
    :return: None
    '''


    # Loop through timesteps and test each timestep:
    for ts in range(nsteps):
        if verb==0:
            print(f"\n------------------------ TIMESTEP {ts} ------------------------")

        for iproc in tester.proc_list:

            if verb>0:
                print(f"\n------------------------ TIMESTEP {ts} - PROCESSOR {iproc} ------------------------")

            # Set processor:
            tester.set_iproc(iproc)

            tester.read_stable()
            tester.read_trial()


            # Get the data for each variable (e.g. strain, displacement)
            tester.stable.get_vars(timestep=ts, vars=vars, verbose=verb)
            tester.trial.get_vars(timestep=ts,  vars=vars, verbose=verb)

            tester.compare_vars(verbose=verb)



if __name__ == "__main__":
    # Called from command line with cmd line args
    # e.g. python3 run_tests.py -test_id uptrough -nprocs 10 -verbosity=2 --pstable ./example_output/uptrough/ --ptrial ./example_output/uptrough/ --BV dis eps grav gpot --TS 2

    # Get args - holds user inputs
    a = parse_args()

    # TEST BODY VARIABLES
    if a.test_body:
        # Generate procs for body :
        PROC_LIST = gen_proc_list(typ='range', p0=a.Nmin, pmax=a.Nmax)

        # Create tester:
        TESTER = Tester(test_id=a.ID,
                        dpath_trial=a.ptrial,
                        dpath_stable=a.pstable,
                        proc_list=PROC_LIST,
                        casetype='body',
                        supress_warnings =a.SW, 
                        precision=a.PREC)

        run_test(tester=TESTER,  verb=a.V, vars=a.VAR, nsteps=a.TS)


    # TEST FREE SURFACE VARIABLES
    if a.test_fs:
        # Generate procs for free surface :
        PROC_LIST = gen_proc_list(typ='from_dir', label='free_surface', dir=a.pstable)

        # Create tester:
        TESTER = Tester(test_id=a.ID,
                        dpath_trial=a.ptrial,
                        dpath_stable=a.pstable,
                        proc_list=PROC_LIST,
                        casetype='fs',
                        supress_warnings =a.SW, 
                        precision=a.PREC)

        run_test(tester=TESTER,  verb=a.V, vars=a.VAR, nsteps=a.TS)