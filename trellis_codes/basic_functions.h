void conv_to_trellis(const int a, const int b, const int memory_order, 
		     const bool conv_array[], Node state_diagram[], 
		     const int nb_states);

void print_array(bool * input_array, int n);

void print_trellis(vector < vector <DecodingNode> > trellis, int n);

int hamming(const bool * str1, const bool * str2, int length);

void generate_rand_string(int * random_string, const int str_length, const int k,
			  gsl_rng * r);

void encoder(const int b, const int * input_str, const int input_length, 
	     const Node * state_diagram, bool * output);

int geom(const double prob);

int mismatched_blocks(const bool * str1, const bool * str2, int length, 
		      int block_size);

double entropy_vector(double * prob_vector, int vector_length);

void qc_interleaver_bool(const int length_vector, const bool * vector1, 
			 bool * vector2, const int condition); 

void qc_interleaver_double(const int length_vector, const double * vector1, 
			   double * vector2, const int condition); 

void qc_interleaver_vec(const int length_vector, 
			const double * vector1, 
			vector <double> vector2, const int condition);

double sum_log(const double x, const double y);

bool kd(const int b);

// unsigned int reverse_bits(unsigned int v, int length);

int bit_reverse(unsigned int v, unsigned int length);

void gen_rect(const int k, const int vpc, const int s);

void gen_markov_seq_init(const int length,  bool codeword[], const int order, 
			 const long double p[], const int init_state,
			 gsl_rng * r);

Node * rect_code(const int option, const int k, const int vpc, 
		 const int input_length_symb, const int b, 
		 const bool conv_array[], const long double markov_tran[], 
		 gsl_rng * r);
