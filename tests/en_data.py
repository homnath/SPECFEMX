import numpy as np

# Dict of abbreviations for variables
VAR_ABBV = {'dis'   : 'displacement'   ,
            'eps'   : 'strain'         ,
            'gpot'  : 'gravity_potential',
            'of'    : 'oceanf',
            'lp'    : 'iceload_phi',
            'lsl'   : 'iceload_sl',
            'lu'    : 'iceload_u',
            'ice'   : 'ice',
            'irate' : 'icerate',
            'sl'    : 'sea_level',
            'grav'  : 'gravity_acceleration'
            }




class EnData():
    '''
    Class to store and load ensight data
    '''


    def __init__(self, id, warnings):
        '''
        :param id: Reference ID/name to indicate what data is stored. Typically 'stable' or 'trial'
        :type id: string
        :param warnings: Flag of whether to print warnings when variables cant be loaded
        :type warnings: bool
        '''

        # User prescribed
        self.id                  = id        # Ref ID specified by user
        self.suppress_warnings   = warnings

        # Initialise
        self.case           = None
        self.geofile        = None
        self.part_names     = None
        self.part           = None
        self.testvarnames   = None
        self.testvars       = None
        self.ntestvars      = None
        self.n_vars_loaded  = None
        self.iproc          = None



    def get_vars(self, timestep, vars, verbose):
        '''
        Load variable data from the .varname (e.g. .dis or .sl) files
        :param timestep: Timestep of data to be loaded
        :type timestep: Timestep of data to be loaded
        :param vars: List of variables to be loaded
        :type vars: List of strings
        :param verbose: Defines level of verbosity (0,1,2) - higher number = more verbose
        :type verbose: int
        :return: none
        '''

        # Signal the reference id for the data being loaded
        if verbose>0:
            print(f'Loading variables for {self.id}')

        # Loaded variables :
        self.vars = vars
        self.loaded_vars = []
        loaded_ctr = 0

        # Tracks if some desired variables are not available (e.g. free surface variables like ice are part of the
        # case file but not available inside the mesh (only on surface)
        unavail_vars = []
        unavail_switch = 0





        # Loops through variables user wants to test
        for bv in self.vars:

            # Will turn to 1 if this var has been loaded/set
            searchctr = 0

            # Keep track of how many variables tried to match so far:
            ctr = 0

            # Loop through available variables to try and match with user prescribed variable:
            for v in self.case.get_variables():
                ctr +=1

                # Has matched with the user prescribed variable name, escape loop
                if searchctr==1:
                    break;

                # Check if user input is using abbreviated variable name?:
                try:
                    abbvvarname = VAR_ABBV[bv]
                except:
                    abbvvarname = bv


                if v == abbvvarname:
                    # Found so set to 1 to escape ```for v in self.case.get_variables()``` loop
                    searchctr = 1
                    try:
                        # Try to load variable data with user specified name:
                        var = self.case.get_variable(v, timestep=timestep)

                        # set each dictionary variable to the loaded np array of data for that variable
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
                        # Cant load it - add to list of unavailable variables the user wanted
                        unavail_vars.append(v)
                        # Switch on to print warning that variables couldnt be printed
                        if self.suppress_warnings == False:
                            unavail_switch = 1


                # If it hasnt found what user is requesting:
                if np.logical_and(ctr == self.ntestvars, searchctr==0):
                    raise ValueError(f"You asked for variable {bv} but the available vars are: {self.case.get_variables()}")


            # If the variable does exist in the case file, but there is no data to load: (e.g. if trying to load a
            # var like ice_load (which is only defined on the free_surface) for the entire domain
            if unavail_switch==1:
                print(f'  * WARNING: The following variables have no data in desired section - NOT LOADED:')
                for k in range(len(unavail_vars)):
                    print(f"    - {self.id} :  {unavail_vars[k]}")
                print()

        # Set number of variables loaded
        self.n_vars_loaded = loaded_ctr



