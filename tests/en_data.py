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
    def __init__(self, id):
        self.id             = id
        self.case           = None
        self.geofile        = None
        self.part_names     = None
        self.part           = None
        self.testvarnames   = None
        self.testvars       = None
        self.ntestvars      = None

        self.nbody_vars_loaded = None
        self.iproc = None

    def get_body_vars(self,timestep, bodyvars, verbose):
        if verbose>0:
            print(f'Loading body variables for {self.id}')

        # User to signal which to load:
        self.bodyvars = bodyvars
        self.loaded_bodyvars = []
        # Function that loads and sets all available variables:
        # For now let us load all of the variables present:
        # set each dictionary variable to the loaded np array of
        # data for that  variable


        loaded_ctr = 0
        # Loops through variables user wants to test
        for bv in self.bodyvars:
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
                            self.loaded_bodyvars.append(v)
                    except:
                        print(f'  * WARNING: {self.id}:{v}:proc{self.iproc} is a variable but file path doesnt exist - NOT LOADED')



                # If it hasnt found what user is requesting:
                if np.logical_and(ctr == self.ntestvars, searchctr==0):
                    raise ValueError(f"You asked for variable {bv} but the available vars are: {self.case.get_variables()}")


        self.nbody_vars_loaded = loaded_ctr
        if verbose>0:
            print()


