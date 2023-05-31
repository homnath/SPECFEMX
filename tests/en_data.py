import numpy as np

VAR_ABBV = {'dis'   : 'displacement'   ,
            'eps'   : 'strain'         ,
            'gpot'  : 'gravity_potential',
            'of'    : 'oceanf',
            'lp'    : 'iceload_phi',
            'lsl'   : 'iceload_sl',
            'iu'    : 'iceload_u',
            'ice'   : 'ice',
            'irate' : 'icerate',
            'sl'    : 'sea_level',
            'grav'  : 'gravity_acceleration'
            }

class EnData():
    def __init__(self, id, warnings):
        self.id             = id
        self.case           = None
        self.geofile        = None
        self.part_names     = None
        self.part           = None
        self.testvarnames   = None
        self.testvars       = None
        self.ntestvars      = None

        self.n_vars_loaded = None
        self.iproc = None
        self.suppress_warnings = warnings

    def get_vars(self,timestep, vars, verbose):
        if verbose>0:
            print(f'Loading variables for {self.id}')

        # User to signal which to load:
        self.vars = vars
        self.loaded_vars = []
        # Function that loads and sets all available variables:
        # For now let us load all of the variables present:
        # set each dictionary variable to the loaded np array of
        # data for that  variable

        # Tracks if some desired variables are not available (e.g. free surface variables like ice are part of the
        # case file but not available inside the mesh (only on surface)
        unavail_vars = []
        unavail_switch = 0

        loaded_ctr = 0
        # Loops through variables user wants to test
        for bv in self.vars:
            # Will turn to 1 if this var has been set

            searchctr = 0
            ctr = 0
            # Loop through available:
            for v in self.case.get_variables():
                # Keep track of how many variables tried to match:
                ctr +=1

                if searchctr==1:
                    break;

                # Using abbreviated variable name?:
                try:
                    abbvvarname = VAR_ABBV[bv]
                except:
                    abbvvarname = bv


                if v == abbvvarname:
                    searchctr = 1   # Found so set to 1 to escape ```for v in self.case.get_variables()``` loop
                    try:
                        # Try to load:
                        var = self.case.get_variable(v, timestep=timestep)
                        with open(var.file_path, "rb") as fp_var:
                            self.testvars[v] = var.read_node_data(fp_var, self.part.part_id)

                            # Print conf
                            if verbose>1:
                                print(f'Loaded: {v} from {var.file_path}')
                            # Increase number of loaded vars:
                            loaded_ctr +=1
                            # Store loaded var:
                            self.loaded_vars.append(v)
                    except:
                        # Cant load it:
                        unavail_vars.append(v)
                        unavail_switch = 1


                # If it hasnt found what user is requesting:
                if np.logical_and(ctr == self.ntestvars, searchctr==0):
                    raise ValueError(f"You asked for variable {bv} but the available vars are: {self.case.get_variables()}")


        if self.suppress_warnings==False:
            if unavail_switch==1:
                print(f'  * WARNING: The following variables have no data in desired section - NOT LOADED:')
                for k in range(len(unavail_vars)):
                    print(f"    - {self.id} :  {unavail_vars[k]}")
                print()

        self.n_vars_loaded = loaded_ctr



