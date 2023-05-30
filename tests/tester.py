import numpy as np
from en_data import EnData
import ensightreader as er


class Tester():
    def __init__(self, test_id, nprocs, dpath_trial, dpath_stable, precision=1e-10):

        self.id             = test_id
        self.trial_fpath    = f'{dpath_trial}'
        self.stable_fpath   = f'{dpath_stable}'

        self.nprocs = nprocs
        self.lnp = len(str(self.nprocs))

        self.iproc = None
        self.precision = precision

        self.print_testing_precision()
        self.print_dirs()



    def read_stable(self):
        self.stable = self.read_ensight(self.stable_fpath, id='stable')


    def read_trial(self):
        self.trial = self.read_ensight(self.trial_fpath, id='trial')


    def set_iproc(self, ip):
        # Sets iproc number...note that needs a number of prefixed 0s depending on nproc:
        sip = str(ip)
        while len(sip)<self.lnp:
            sip = '0'+ sip

        self.iproc = sip



    def read_ensight(self, path, id):
        d = EnData(id)
        d.case       = er.read_case(path + f"/{self.id}_proc{self.iproc}.case")
        d.geofile    = d.case.get_geometry_model()
        d.part_names = d.geofile.get_part_names()
        d.part       = d.geofile.get_part_by_name(d.part_names[0])

        # Create dictionary of testable variables:
        # List defaults values to None
        d.testvarnames  = d.case.get_variables()
        d.testvars      = dict.fromkeys(d.testvarnames)
        d.ntestvars     = len(d.testvars)
        d.iproc = self.iproc

        return d


    def check_eq_num_vars(self):
        # Check the number of LOADED body variables is same for stable and test cases:
        if self.trial.nbody_vars_loaded != self.stable.nbody_vars_loaded:
            raise ValueError(f"Trial num body vars loaded  : {self.trial.nbody_vars_loaded} \n            Stable num body vars loaded : {self.stable.nbody_vars_loaded}")

        assert(len(self.trial.loaded_bodyvars)  == self.trial.nbody_vars_loaded)
        assert(len(self.stable.loaded_bodyvars) == self.stable.nbody_vars_loaded)



    def compare_body_vars(self, verbose=0):
        # First check there are equal numbers of loaded vars to test:
        self.check_eq_num_vars()

        nvars = self.trial.nbody_vars_loaded

        for i in range(nvars):
            # Get variable name, e.g. displacement (str)
            vs = self.stable.loaded_bodyvars[i]
            vt = self.trial.loaded_bodyvars[i]

            if vs!=vt:
                raise ValueError("self.stable.loaded_bodyvars[i] is not equal to self.trial.loaded_bodyvars[i]")


            if verbose==2:
                print(vt + ':      TRIAL         STABLE')
                print("     Min val:   ", np.min(self.trial.testvars[vt]), '    ', np.min(self.stable.testvars[vt]))
                print("     Max val:   ", np.max(self.trial.testvars[vt]), '    ', np.max(self.stable.testvars[vt]))


            try:
                assert ((np.abs(self.trial.testvars[vt] - self.stable.testvars[vt]) < self.precision).all())

                if verbose>0:
                    print(f'  ✔  Variable {vt} matches')
                    if verbose==2:
                        print()
            except:
                errstr = f'  ✖ Variable {vt} has different values'
                raise ValueError(errstr)


        print(f'-- Tested {nvars} variables for processor {self.iproc}')


    def print_dirs(self):
        print('Loading data from:')
        print(f'• TEST  : {self.trial_fpath}')
        print(f'• STABLE: {self.stable_fpath}')
        print()

    def print_testing_precision(self):

        print('*****************************************************************')
        print("""           _   _   _   _   _     _   _   _   _   _   _   _  
          / \ / \ / \ / \ / \   / \ / \ / \ / \ / \ / \ / \ 
         ( B | E | G | I | N ) ( T | E | S | T | I | N | G )
          \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/ \_/ 
          """)
        print('*****************************************************************')
        print('             TESTING WITH PRECISION: ', self.precision)
        print('*****************************************************************')
        print()

    def print_test_stats(self, nsteps=None, verb=None, bodyvars=None):
        print(f" Test ID         :     {self.id}")
        print(f" Num. procs      :     {self.nprocs}")
        print(f" Verbosity level :     {verb}")
        print(f" Num. timesteps  :     {nsteps}")
        print(f" Body variables  :     {bodyvars}")
