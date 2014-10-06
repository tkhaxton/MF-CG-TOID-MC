# define rand_double (1.0*rand()/(1.0*RAND_MAX))
# define recenter(a, b) a-=floor(a/b+0.5)*b
# define fmod(a, b) a-=floor(a/b)*b
# define mod(a, b) ((a)>=0)?((a)%(b)):(((b)+(a)%(b))%(b))					//	% doesn't mod negative numbers correctly!

typedef struct int_triple_tag{
	int x;
	int y;
	int z;
} int_triple;

typedef struct int_double_tag{
	int x;
	int y;
} int_double;

typedef struct double_triple_tag{
	double x;
	double y;
	double z;
} double_triple;

double_triple add_double_triple(double_triple a, double_triple b);
double_triple scalar_multiply_double_triple(double_triple a, double s);
double dot_product(double_triple a, double_triple b);
double_triple rand_unit_cube();
double_triple rand_unit_ball();
double_triple rand_unit_sphere();
double_triple subtract_double_triple(double_triple a, double_triple b);
void recenter_double_triple(double_triple *pa, double_triple b);
void fmod_double_triple(double_triple *pa, double_triple b);
void normalize(double_triple *pvec);
double_triple cross_product(double_triple a, double_triple b);
double norm(double_triple a);
double largest_cubic_root(double a3, double a2, double a1, double a0);
double quartic_global_minimum(double p4, double p3, double p2, double p1);
double_triple normed(double_triple vec);
void forward_and_backward_matrix(double angle, double_triple axis_vector, double **forward, double **backward);
double_triple rotate_by_matrix(double_triple start, double **matrix);
void forward_matrix(double angle, double_triple axis_vector, double **forward);
void matrix_multiply(double **a, double **b, double **c);


