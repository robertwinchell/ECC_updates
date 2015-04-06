
void MAP_NSC_rect(const int nb_branches, const int input_length_symb, 
		  const int half_corridor_width, const int b,
		  const double p_del, const double p_ins, const double p_dup, 
		  const double p_sub, const Node state_diagram[], const int vpc, 
		  const bool * corrupted_codeword, const int corrupted_length,
		  const vector < vector <double> > prior_prob,
		  vector < vector <double> > &decoded_app); 
