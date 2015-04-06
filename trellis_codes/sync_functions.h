int gaps(const bool str1[], const bool str2[], int length_str1, int length_str2);

void delete_bits_table(bool codeword[], int initial_length, int positions[], 
		       int nb_deletions);

int binary_NSC(const bool codeword[], bool corrupted_codeword[], 
	       const int codeword_length, const double p_del, const double p_ins, 
	       const double p_dup, const double p_sub, gsl_rng * r);

int binary_NSC_drift(const bool codeword[], bool corrupted_codeword[], 
		     const int codeword_length, const double p_del, 
		     const double p_ins, const double p_dup, const double p_sub, 
		     int * final_drift, gsl_rng * r);
