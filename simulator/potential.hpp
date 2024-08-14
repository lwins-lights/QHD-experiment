extern const double compress_coef;
extern const double slope;
extern const int dim;
extern const double lb[];
extern const double ub[];
extern const double pinned[];

void get_potential_params(int &dim);
void get_pinned_point(double *ret, double L);
double get_potential(const double *x, double L);
void get_subgradient(const double *x, double *ret, double L);
double get_obj(const double *x);
void get_obj_subg(const double *x, double *ret);