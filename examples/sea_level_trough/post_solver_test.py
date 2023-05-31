from tester import Tester, parse_args
import sys

def run_test(testid, nprocs, verb, body_vars, nsteps, dpath_trial, dpath_stable):
    # Create tester:
    tester = Tester(testid, nprocs, dpath_trial=dpath_trial, dpath_stable=dpath_stable)

    # Loop through timesteps:
    for ts in range(nsteps):
        for iproc in range(tester.nprocs):

            print(f"\n------------------------ TIMESTEP {ts} - PROCESSOR {iproc} ------------------------")

            # Set processor:
            tester.set_iproc(iproc)

            tester.read_stable()
            tester.read_trial()


            # Get the data for each variable (e.g. strain, displacement)
            # that is stored throughout the entire body (not just free surface):
            tester.stable.get_body_vars(timestep=ts, bodyvars=body_vars, verbose=verb)
            tester.trial.get_body_vars(timestep=ts,  bodyvars=body_vars, verbose=verb)
            tester.compare_body_vars(verbose=verb)


if __name__ == "__main__":


    # Get args:
    a = parse_args()
    run_test(testid=a.ID, nprocs=a.N, verb=a.V, body_vars=a.BV, nsteps=a.TS, dpath_trial=a.ptrial, dpath_stable=a.pstable)


