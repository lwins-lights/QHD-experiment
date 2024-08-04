extern const double compress_coef;
extern const double slope;
extern const int dim;
extern const double lb[];
extern const double ub[];

void get_potential_params(int &dim);
double get_potential(const double *x, double L);
void get_subgradient(const double *x, double *ret, double L);
double get_obj(const double *x);
void get_obj_subg(const double *x, double *ret);