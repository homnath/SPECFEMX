from tester import Tester

TEST_ID   = 'uptrough'
NPROCS    = 10
verb      = 0
body_vars = ['dis', 'eps', 'gpot', 'grav']
nsteps    = 2


# Create tester:
tester = Tester(TEST_ID, NPROCS, dpath_trial='./output_uptrough/', dpath_stable='./stable_output')
tester.print_test_stats(nsteps=nsteps, verb=verb, bodyvars=body_vars)

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
