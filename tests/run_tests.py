import ensightreader as er
from ensightreader import read_case, ElementType
import matplotlib.pyplot as plt
import numpy as np

class Tester():

    def __init__(self, fpath): 
        self.fpath = fpath


        # Load parts: 
        self.read_ensight()

        # Create dictionary of testable variables: 
        # List defaults values to None
        self.testvarnames = self.case.get_variables()
        self.testvars = dict.fromkeys(self.testvarnames)
        self.ntestvars = len(self.testvars)


    def read_ensight(self): 
        self.case       = er.read_case(self.fpath)
        self.geofile    = self.case.get_geometry_model()
        self.part_names = self.geofile.get_part_names()     
        self.part       = self.geofile.get_part_by_name(self.part_names[0])




    def load_ensight_var(self, var_name, timestep):
        # Load variable: 
        var = self.case.get_variable(var_name, timestep=timestep)
        with open(var.file_path, "rb") as fp_var:
            var_out = var.read_node_data(fp_var, self.part.part_id) 

        return var_out

    def get_all_vars(self,timestep):
        # Function that loads and sets all available variables: 
        # For now let us load all of the variables present: 
        for v in trial.case.get_variables():
            # set each dictionary variable to the loaded np array of 
            # data for that  variable
            self.testvars[v] = trial.load_ensight_var(var_name=v, 
                                                       timestep=timestep
                                                      )


precision = 1e-16


print('*****************************************************************')
print("""   _   _   _   _   _     _   _   _   _   _   _   _  
  / \ / \ / \ / \ / \   / \ / \ / \ / \ / \ / \ / \ 
 ( B | E | G | I | N ) ( T | E | S | T | I | N | G )
  \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/ \_/ 
  """)
print('*****************************************************************')
print('             TESTING WITH PRECISION: ', precision)
print('*****************************************************************')

# ----------------------------------------------------------------------
# Test 1 

test_no = '1' 
fname = 'test1_trough'
ts = 0 # timestep

# File path for SPECFEMX_SL results we are testing
trial_fpath = f'./outputs/test{test_no}/{fname}_proc0.case'
# File path for results from SPECFEMX (stable) from HNG repo
stable_fpath = f'./stable_outputs/test{test_no}/{fname}_proc0.case'

print('Loading data from:')
print(f'• TEST  : {trial_fpath}')
print(f'• STABLE: {stable_fpath}')


# Get the data for each variable (e.g. strain, displacement)
trial = Tester(trial_fpath)
trial.get_all_vars(timestep=ts)

# Repeat for the stable simulation
stable = Tester(stable_fpath)
stable.get_all_vars(timestep=ts)



# Now check each datapoint is identical:
print()
# First check that the number of variables being tested is the same: 
try: 
    assert(trial.ntestvars==stable.ntestvars)
    confstr = f'Number of vars is equal: {trial.ntestvars}' 
    print(confstr)
except: 
    errstr = f'Number of variables outputted/loaded here is different: trial = {trial.ntestvars} stable = {stable.ntestvars} \n trial  variables: {trial.testvarnames} \n stable variables: {stable.testvarnames} \n ' 
    raise ValueError(errstr)


for i in range(trial.ntestvars):
    # Get variable name, e.g. displacement (str)
    v = trial.testvarnames[i]

    try: 
        assert((trial.testvars[v] - stable.testvars[v] < precision).all())
        confstr = f'✔  Variable {v} matches' 
        print(confstr)
    except: 
        errstr = f'✖ Variable {v} has different values' 
        raise ValueError(errstr)

# ----------------------------------------------------------------------
