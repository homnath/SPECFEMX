import numpy as np
from en_data import EnData
import ensightreader as er
import argparse
import os

all_vars = ['displacement', 'strain', 'gravity_potential', 'oceanf', 'iceload_phi',
            'iceload_sl', 'iceload_u', 'ice', 'icerate', 'sea_level', 'gravity_acceleration']


class Tester():
    def __init__(self, test_id, casetype, proc_list, dpath_trial, dpath_stable, precision=1e-10, supress_warnings=False):

        self.id             = test_id
        self.id0            = test_id

        self.trial_fpath    = f'{dpath_trial}'
        self.stable_fpath   = f'{dpath_stable}'


        self.iproc = None
        self.precision = precision

        self.print_testing_precision()
        self.print_dirs()

        self.update_casetype(casetype)

        self.proc_list = proc_list
        self.nprocs    = len(self.proc_list)
        self.lnp       = len(str(self.nprocs))

        self.suppress_warnings = supress_warnings


    def update_casetype(self, cs):
        if np.logical_and(cs!='fs', cs!='body'):
            raise ValueError("Casetype must be 'fs' or 'body' ")
        # can be fs or body
        self.casetype = cs
        if self.casetype == 'fs':
            self.id = self.id0 + '_free_surface'

        print('Updated casetype to: ', cs)

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
        d = EnData(id, self.suppress_warnings)
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
        # Check the number of LOADED variables is same for stable and test cases:
        if self.trial.n_vars_loaded != self.stable.n_vars_loaded:
            raise ValueError(f"Trial num vars loaded  : {self.trial.n_vars_loaded} \n            Stable num vars loaded : {self.stable.n_vars_loaded}")

        assert(len(self.trial.loaded_vars)  == self.trial.n_vars_loaded)
        assert(len(self.stable.loaded_vars) == self.stable.n_vars_loaded)



    def compare_vars(self, verbose=0):
        # First check there are equal numbers of loaded vars to test:
        self.check_eq_num_vars()

        nvars = self.trial.n_vars_loaded

        for i in range(nvars):
            # Get variable name, e.g. displacement (str)
            vs = self.stable.loaded_vars[i]
            vt = self.trial.loaded_vars[i]

            if vs!=vt:
                raise ValueError("self.stable.loaded_vars[i] is not equal to self.trial.loaded_vars[i]")


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
                # Calculate avg error %
                err = np.abs(self.trial.testvars[vt] - self.stable.testvars[vt]).flatten()

                means = [np.abs(np.mean(self.trial.testvars[vt])), np.abs(np.mean(self.stable.testvars[vt]))]
                minmean = np.min(means)
                if minmean==0:
                    minmean = np.max(means)

                approx_err = np.mean(err)/minmean

                errstr = f'  ✖ Variable {vt} has different values -  mean error %: {approx_err}'
                raise ValueError(errstr)


        print(f'  ✔ Tested {nvars} variables for processor {self.iproc}')


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








def parse_args():

    parser = argparse.ArgumentParser(prog='SPECFEMX Tester', description='Run test case for SPECFEMX')
    parser.add_argument('-test_id', '--ID', type=str, nargs='?',
                        help='Test ID e.g. test1 or uptrough (str)')
    parser.add_argument('-nprocs', '--N', action='store', type=int,
                        help='Number of processors (int)')
    parser.add_argument('-ntsteps', '--TS', action='store', type=int,
                        help='Number of timesteps (int)')
    parser.add_argument('-verbosity', '--V', action='store', type=int, choices=range(3), default=0,
                        help='Verbosity (0-2: default=0) - higher val = more verbose (int)')
    parser.add_argument('-path_trial', '--ptrial', type=str, nargs='?',
                        help='File path to trial results directory (str)')
    parser.add_argument('-path_stable', '--pstable', type=str, nargs='?',
                        help='File path to stable results directory (str)')
    parser.add_argument('-vars', '--VAR', default=[], nargs='+', help='List of variables to test. Use "all" to test all available')
    parser.add_argument('-test_body', '--test_body', action='store_true', help='Test variable values in entire body')
    parser.add_argument('-test_fs', '--test_fs', action='store_true', help='Test variable values on free surface body')

    a = parser.parse_args()

    # Check if asking for all available variables:
    if np.logical_and(len(a.VAR) == 1, a.VAR[0] == 'all'):
        a.VAR = all_vars


    # Print test parameters:
    print("----------------- TESTING PARAMS -----------------")
    print(f" Test body vars. :     {a.test_body}")
    print(f" Test fs vars.   :     {a.test_fs}\n")
    print(f" Test ID          :     {a.ID}")
    print(f" Num. procs       :     {a.N}")
    print(f" Verbosity level  :     {a.V}")
    print(f" Num. timesteps   :     {a.TS}")
    print(f" Variables   :     {a.VAR}")
    print(f" Trial file path  :     {a.ptrial}")
    print(f" Stable file path :     {a.pstable}")
    print("--------------------------------------------------")

    return a





def gen_proc_list(typ, p0=0, pmax=None, dir=None, label=None):
    if typ=='range':

        if type(pmax)==type(None):
            raise ValueError("Max processor number not specified.")

        list = np.arange(p0, pmax)
        return list
    elif typ=='from_dir':

        assert(label!=None)

        # Initialise list:
        list = []

        # Looks for any processors that have a case file 'label' in the name
        for f in os.listdir(dir):
            if np.logical_and(f.find('case')!=-1 , f.find(label)!=-1):

                # Get processor value:
                ftmp = f[f.find('proc')+4:]
                proc = ftmp[:ftmp.find('.case')]

                list.append(int(proc))

        # Convert to numpy and sort:
        list = np.array(list)
        list.sort()

        return list