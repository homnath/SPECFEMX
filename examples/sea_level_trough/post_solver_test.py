from tester import Tester, parse_args, gen_proc_list
import sys

def run_test(tester, verb, vars, nsteps):
    # Loop through timesteps:
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
    # Get args:
    a = parse_args()

    if a.test_body:
        # Generate procs for body :
        PROC_LIST = gen_proc_list(typ='range', p0=0, pmax=a.N)

        # Create tester:
        TESTER = Tester(test_id=a.ID,
                        dpath_trial=a.ptrial,
                        dpath_stable=a.pstable,
                        proc_list=PROC_LIST,
                        casetype='body',
                        supress_warnings =True)

        run_test(tester=TESTER,  verb=a.V, vars=a.VAR, nsteps=a.TS)



    if a.test_fs:
        # Generate procs for free surface :
        PROC_LIST = gen_proc_list(typ='from_dir', label='free_surface', dir=a.pstable)

        # Create tester:
        TESTER = Tester(test_id=a.ID,
                        dpath_trial=a.ptrial,
                        dpath_stable=a.pstable,
                        proc_list=PROC_LIST,
                        casetype='fs',
                        supress_warnings =True)

        run_test(tester=TESTER,  verb=a.V, vars=a.VAR, nsteps=a.TS)