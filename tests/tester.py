import numpy as np
from en_data import EnData
import ensightreader as er
import argparse
import os

all_vars = ['displacement', 'strain', 'gravity_potential', 'oceanf', 'iceload_phi',
            'iceload_sl', 'iceload_u', 'ice', 'icerate', 'sea_level', 'gravity_acceleration']


class Tester():
    '''
    Class which performs the majority of testing.
    '''
    def __init__(self, test_id, casetype, proc_list, dpath_trial, dpath_stable, precision=1e-10, supress_warnings=False):
        '''
        Instantiates an instance of the Tester class
        :param test_id: unique ID for the test
        :type test_id: string
        :param casetype: signals whether testing free surface or entire body. Must be one of 'fs' or 'body'.
        :type casetype: string
        :param proc_list: list of processor IDs to be tested
        :type proc_list: list of integers
        :param dpath_trial: file path to the directory holding trial data outputted from SPECFEMX
        :type dpath_trial: string
        :param dpath_stable: file path to the directory holding stable data outputted from SPECFEMX to be tested against
        :type dpath_stable: string
        :param precision: Precision to which variables will be tested against one another. Default 1e-10.
        :type precision: int or float
        :param supress_warnings: If True, will not print warnings about requested test variables that couldnt be loaded. Default False.
        :type supress_warnings: Bool
        :returns None
        '''

        self.id             = test_id                   # Unique test ID
        self.id0            = test_id                   # Store of unique test ID

        self.trial_fpath    = f'{dpath_trial}'          # File path to trial data
        self.stable_fpath   = f'{dpath_stable}'         # File path to stable data

        self.iproc = None                               # ID of current processor
        self.precision = precision                      # Precision

        # Init. functions:
        # Print details on precision, directories to load from
        self.print_testing_precision()
        self.print_dirs()
        # Update the type of case file being searched for based on user input - edit value of self.casetype and self.id
        self.update_casetype(casetype)


        self.proc_list = np.array(proc_list)                      # List of processor IDs
        self.nprocs    = len(self.proc_list)                      # Number of processors to test
        self.lnp       = len(str(np.amax(self.proc_list)))        # Number of characters in highest proc ID string
                                                                  # --> needed since procs are stored as 02/002, not just '2'

        self.suppress_warnings = supress_warnings                 # Suppresses warnings from being printed



    def update_casetype(self, cs):
        '''
        Updates the value of self.casetype. If looking for free surface Case files then will edit the self.id tag
        :param cs: Case type - can be free surface ('fs') or for entire body ('body').
        :type cs: string
        :return:
        '''

        # Check input from user - can be fs or body
        if np.logical_and(cs!='fs', cs!='body'):
            raise ValueError("Casetype must be 'fs' or 'body' ")

        # Set new value
        self.casetype = cs

        # Edit id to search for free surface case files if desired
        if self.casetype == 'fs':
            self.id = self.id0 + '_free_surface'
        # Confirm to user
        print('Updated casetype to: ', cs)


    def read_stable(self):
        '''
        Reads ensight files from the file path stored in self.stable_fpath
        :return: none
        '''
        self.stable = self.read_ensight(self.stable_fpath, id='stable')


    def read_trial(self):
        '''
        Reads ensight files from the file path stored in self.trial_fpath
        :return: none
        '''
        self.trial = self.read_ensight(self.trial_fpath, id='trial')


    def set_iproc(self, ip):
        '''
        Appends the required number of 0s to the front of the processor ID value to be consistent with SPECFEMX files
        :param ip: Processor ID
        :return: None. updates self.iproc value
        '''

        # Sets iproc number...note that needs a number of prefixed 0s depending on nproc:
        sip = str(ip)
        while len(sip)<self.lnp:
            sip = '0'+ sip
        self.iproc = sip



    def read_ensight(self, path, id):
        '''
        Reads ensight data for a given processor from the specified file path
        :param path: File path to ensight data.
        :type path: string
        :param id: either 'stable' or 'trial' indicating the type of data being read
        :return: EnData object holding the relevant data
        '''

        # Instantiate EnData object with the given ID and warnings flag
        d = EnData(id, self.suppress_warnings)
        # Read the case file for this processor
        d.case       = er.read_case(path + f"/{self.id}_proc{self.iproc}.case")

        # Import other details
        d.geofile    = d.case.get_geometry_model()
        d.part_names = d.geofile.get_part_names()
        d.part       = d.geofile.get_part_by_name(d.part_names[0])

        # Create dictionary of testable variables:
        # List defaults values to None
        d.testvarnames  = d.case.get_variables()
        d.testvars      = dict.fromkeys(d.testvarnames) # Dictionary of testable variable names
        d.ntestvars     = len(d.testvars)               # Number of variables in file
        d.iproc = self.iproc                            # Assign the processor ID to the EnData obj.

        return d


    def check_eq_num_vars(self):
        '''
        Check the actual number of variables that could be loaded is same for stable and test cases:
        :return: none
        '''

        if self.trial.n_vars_loaded != self.stable.n_vars_loaded:
            raise ValueError(f"Trial num vars loaded  : {self.trial.n_vars_loaded} \n            Stable num vars loaded : {self.stable.n_vars_loaded}")

        assert(len(self.trial.loaded_vars)  == self.trial.n_vars_loaded)
        assert(len(self.stable.loaded_vars) == self.stable.n_vars_loaded)



    def compare_vars(self, verbose=0):
        '''
        Compares each loaded variable to see if stable and trial arrays are equal to within specified precision
        :param verbose: Specifies level of verbosity in output to user. Can be 0, 1, 2. Higher value = more verbose.
        :type verbose: integer (0, 1, 2)
        :return:
        '''

        # First check there are equal numbers of loaded vars to test and store that val
        self.check_eq_num_vars()
        nvars = self.trial.n_vars_loaded

        # Look through each variable
        for i in range(nvars):
            # Get variable name, e.g. displacement (str)
            vs = self.stable.loaded_vars[i]
            vt = self.trial.loaded_vars[i]

            # Ensure comparing same variables:
            if vs!=vt:
                raise ValueError("self.stable.loaded_vars[i] is not equal to self.trial.loaded_vars[i]")

            # For highly verbose - print min, max values of each array
            if verbose==2:
                print(vt + ':      TRIAL         STABLE')
                print("     Min val:   ", np.min(self.trial.testvars[vt]), '    ', np.min(self.stable.testvars[vt]))
                print("     Max val:   ", np.max(self.trial.testvars[vt]), '    ', np.max(self.stable.testvars[vt]))


            try:
                # Assert whether at each point in the array, the compared values are equal within specified precision
                # If assertion fails, will go to except block
                assert ((np.abs(self.trial.testvars[vt] - self.stable.testvars[vt]) < self.precision).all())

                # If quite verbose, print match for each variable
                if verbose>0:
                    print(f'  ✔  Variable {vt} matches')
                    if verbose==2:
                        print()
            except:
                # Block when values to not match

                # Calculate avg error %
                err = np.abs(self.trial.testvars[vt] - self.stable.testvars[vt]).flatten()
                means = [np.abs(np.mean(self.trial.testvars[vt])), np.abs(np.mean(self.stable.testvars[vt]))]

                # Use minimum mean value to divide by when calculating error - if val is 0 then use the other mean val.
                minmean = np.min(means)
                if minmean==0:
                    minmean = np.max(means)

                approx_err = np.mean(err)/minmean

                # Raise error that there are differences in values and state approximate error. %
                errstr = f'  ✖ Variable {vt} has different values -  mean error %: {approx_err}'
                raise ValueError(errstr)


        # If successful for ALL variables, print success of testing for this processor
        print(f'  ✔ Tested {nvars} variables for processor {self.iproc}')


    def print_dirs(self):
        '''
        Prints the directories being used for testing
        :return: none
        '''
        print('Loading data from:')
        print(f'• TEST  : {self.trial_fpath}')
        print(f'• STABLE: {self.stable_fpath}')
        print()

    def print_testing_precision(self):
        '''
        Prints the startup screen for testing with the precision
        :return: none
        '''
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
    '''
    Parses command line inputs from the user.
    :return: ArgumentParser object (from module argpass) holding inputs from user
    '''

    # Instantiate object
    parser = argparse.ArgumentParser(prog='SPECFEMX Tester', description='Run test case for SPECFEMX')

    # Assign arguments user can input
    parser.add_argument('-test_id', '--ID', type=str, nargs='?',
                        help='Test ID e.g. test1 or uptrough (str)')
    parser.add_argument('-Nmin', '--Nmin', action='store', type=int,
                        help='Number of processors (int)', default=0)
    parser.add_argument('-Nmax', '--Nmax', action='store', type=int,
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
    parser.add_argument('-suppress_warnings', '--SW', action='store_true', help='Supress printed warnings if vars not loaded')

    # Check the args that have been parsed by user
    a = parser.parse_args()

    # Check if asking for all available variables:
    if np.logical_and(len(a.VAR) == 1, a.VAR[0] == 'all'):
        a.VAR = all_vars


    # Print test parameters:
    print("----------------- TESTING PARAMS -----------------")
    print(f" Test body vars.  :     {a.test_body}")
    print(f" Test fs vars.    :     {a.test_fs}\n")
    print(f" Test ID          :     {a.ID}")
    print(f" Min proc in range:     {a.Nmin}")
    print(f" Max proc in range:     {a.Nmax}")
    print(f" Verbosity level  :     {a.V}")
    print(f" Num. timesteps   :     {a.TS}")
    print(f" Variables        :     {a.VAR}")
    print(f" Trial file path  :     {a.ptrial}")
    print(f" Stable file path :     {a.pstable}")
    print(f" Suppress warning :     {a.SW}")
    print("--------------------------------------------------")

    return a





def gen_proc_list(typ, p0=0, pmax=None, dir=None, label=None):
    '''
    Generates list of processors to be tested. For entire body, generally want all processors. For free surface case
    files, can only test processors on the free surface. Can either generate a list from a range of processor IDs, or
    search for all the processor IDs in a given file.

    If using typ='range' then pmax must be specified. This is the maximum processor ID. By default the minimum ID is 0,
    but can be specified with p0.

    If using typ='from_dir', must specify dir and label arguments. Function will seach in the specified directory 'dir'
    for any case files that include the string 'label', and append the processor from that file to the list.

    :param typ: Either 'range' or 'from_dir'
    :type typ: string
    :param p0: Minimum processor ID in range, if typ='range'. Default = 0
    :type p0: int
    :param pmax: Maximum processor ID in range, if typ='range'. Default=None so function will fail if not specified when using 'range'.
    :type pmax: int
    :param dir: file path of directory to search for case files
    :type dir: string
    :param label: string that specifies the case files to select e.g. 'free_surface' to find free
                  surface case files of the form 'uptrough_free_surface_proc00'
    :type label: string
    :return:
    '''


    # If using range type:
    if typ=='range':
        # Assert that max proc ID is specified
        if type(pmax)==type(None):
            raise ValueError("Max processor number not specified.")
        return np.arange(p0, pmax)

    # If using from_dir:
    elif typ=='from_dir':
        # Check label is specified
        if label==None:
            raise ValueError('label has not been specified')

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
        # Remove any duplicates and then sort
        list = np.unique(list)
        list.sort()

        return list
    else:
        raise ValueError("typ must have value 'range' or 'from_dir' ")