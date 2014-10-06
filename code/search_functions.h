
void reflect(double *pvar);
void output_monolayer(char *configname, monolayer config, double energy, int frame);
void copymonolayer(monolayer new, monolayer *pchanging, int Nnodetypes);
void initialize_periodic_monolayer(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config);
void copybilayer(bilayer new, bilayer *pchanging, int Nnodetypes);
void initialize_periodic_bilayer(int Nmonomers, int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry);
void output_bilayer(char *configname, bilayer config, double energy, int frame);
