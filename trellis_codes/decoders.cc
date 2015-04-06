// ******************************************************************************

#include "definitions.h"
#include "basic_functions.h"
#include "sync_functions.h"

/* ******************************************************************************
  
MAP decoder for the NSC and rectangular graphs.

INPUT ARGUMENTS:
- nb_branches: the alphabet size, or the number of outgoing branches per state  
- input_length: the number of SYMBOLS in the input sequence
- half_corridor_width: the width of the half-corridor 
- b: the number of output bits per input symbols (code of rate a:b)
- p_del, p_ins, p_dup, p_sub: channel parameters
- state_diagram: the complete state diagram (rectangular graph)
- vpc: the number of vertices per column in the rectangular graph 
- corrupted_codeword: the received sequence
- corrupted_length: the length of the corrupted received sequence
- prior_prob: the prior probabilities of the input symbols
- decoded_app: the posterior probabilities of the input symbols after decoding

OPTIMIZATION:
- A matrix on the heap is used instead of the vector class for the large 
  multidimensional arrays required by the decoding algorithm

- The algorithm works on a narrow corridor in the lattices.

- The trellis of lattices use double instead of long double entries. This 
  reduces the memory requirements for large trellises and lattices. The 
  normalization has not caused problems during testing.

NOTES: 
- The decoding corridor used during the decoding process is based on 
  the length of the input and received sequences, whereas the half-corridor 
  width provided as an input argument might come from the parameters of the 
  NSC. The user has to make sure that the half-corridor width provided in 
  input is large enough to include the true synchronization path 
  Taking the corridor straight from the channel without expanding it might 
  exclude the true synchronization path from the decoder corridor, 
  thus decreasing performance.

- The columns of the lattices are only defined around the corridor width to 
  reduce memory requirements. The transitions from left to right in the lattices
  must therefore be dealt according to the drift of the corridors.

ASSUMPTIONS: 
- The codewords are terminated (in state 0)
- The encoders are rectangular graphs.

INEFFICIENCIES AND POSSIBLE IMPROVEMENTS
- Everything is done in float-point arithmetic; log domain is not used
- Max-Log-Map could be implemented

OUTPUT:
The algorithm returns the posterior probabilities of the input symbols 
after decoding in the argument decoded_app */

void MAP_NSC_rect(const int nb_branches, const int input_length_symb,
		  const int half_corridor_width, const int b, const double p_del,
		  const double p_ins, const double p_dup, const double p_sub,
		  const Node state_diagram[], const int vpc,  
		  const bool * corrupted_codeword, const int corrupted_length, 
		  const vector < vector <double> > prior_prob,
		  vector < vector <double> > &decoded_app) 
{

	const int codeword_length = input_length_symb * b;

  // **************************
  // *** Decoding Corridors ***
  // **************************
  
  // Corridor width
	int corridor_width = 2 * half_corridor_width + 1;
  
  // The center of the decoding corridor. 
	int corridor_center[codeword_length + 1];

  // Top and bottom of the decoding corridor
	int corridor_top[codeword_length + 1];
	int corridor_bottom[codeword_length + 1];
  
  // Corridor boundaries
	for(int i = 0; i < codeword_length + 1; i++)
	{   
    	    corridor_center[i] = int(double(corrupted_length) / double(codeword_length) * double(i));
    	    corridor_top[i] = corridor_center[i] - half_corridor_width;
    	    corridor_bottom[i] = corridor_center[i] + half_corridor_width;
      
      // Check if the corridor boundaries are in the lattice
    	    if(corridor_top[i] < 0) corridor_top[i] = 0;
    	    if(corridor_bottom[i] > corrupted_length) corridor_bottom[i] = corrupted_length;
	}
  
  // Corridor boundaries for deletions (forward and backward recursions)
	int del_top_f[codeword_length + 1];
	int del_bottom_f[codeword_length + 1];
	int del_top_b[codeword_length + 1];
	int del_bottom_b[codeword_length + 1];

	for(int i = 1; i < codeword_length + 1; i++)
	{   
    	    del_top_f[i] = corridor_top[i];
    	    del_bottom_f[i] = corridor_bottom[i - 1];
	}
 
	for(int i = 0; i < codeword_length; i++)
	{   
    	    del_top_b[i] = del_top_f[i + 1];
    	    del_bottom_b[i] = del_bottom_f[i + 1];
	}

  // Corridor boundaries for substitutions and correct transitions
  // (forward and backward transitions)
	int sub_top_f[codeword_length + 1];
	int sub_bottom_f[codeword_length + 1];
	int sub_top_b[codeword_length + 1];
	int sub_bottom_b[codeword_length + 1];

	for(int i = 1; i < codeword_length + 1; i++)
	{   
    	    if(corridor_top[i] == corridor_top[i - 1]) sub_top_f[i] = corridor_top[i] + 1;
    	    else sub_top_f[i] = corridor_top[i];
      
    	    if(corridor_bottom[i] == corridor_bottom[i - 1]) sub_bottom_f[i] = corridor_bottom[i];
    	    else sub_bottom_f[i] = corridor_bottom[i - 1] + 1;
	}

	for(int i = 0; i < codeword_length; i++)
	{   
    	    sub_top_b[i] = sub_top_f[i + 1] - 1;
    	    sub_bottom_b[i] = sub_bottom_f[i + 1] - 1;
	}

  // The corridor boundaries may be reduced when insertions / duplications 
  // and deletions cannot occur simultaneously

  // Only deletions
	if(p_del > 0 && p_ins + p_dup == 0)
	{
    	    corridor_top[0] = 0;
    	    corridor_bottom[0] = 0;

    	    for(int i = 1; i < codeword_length + 1; i++)
	    {
		if(corridor_top[i] < i - codeword_length + corrupted_length) corridor_top[i] = i - codeword_length + corrupted_length;
		if(corridor_bottom[i] > i) corridor_bottom[i] = i;
		if(del_top_f[i] < i - codeword_length + corrupted_length) del_top_f[i] = i - codeword_length + corrupted_length;
		if(del_bottom_f[i] > i - 1) del_bottom_f[i] = i - 1;
		if(sub_top_f[i] < i - codeword_length + corrupted_length) sub_top_f[i] = i - codeword_length + corrupted_length; 
		if(sub_bottom_f[i] > i) sub_bottom_f[i] = i;
	    }

    	    for(int i = 0; i < codeword_length; i++)
	    {   
		del_top_b[i] = del_top_f[i + 1];
		del_bottom_b[i] = del_bottom_f[i + 1];
	    }
      
    	    for(int i = 0; i < codeword_length; i++)
	    {   
		sub_top_b[i] = sub_top_f[i + 1] - 1;
		sub_bottom_b[i] = sub_bottom_f[i + 1] - 1;
	    }
	}
  
  // Only insertions and duplications
	if(p_del == 0 && p_ins + p_dup > 0)
	{
    	    for(int i = 0; i < codeword_length + 1; i++)
	    {
		if(corridor_top[i] < i) corridor_top[i] = i;
		if(corridor_bottom[i] > i + corrupted_length - codeword_length) corridor_bottom[i] = i + corrupted_length - codeword_length;
	    }

    	    for(int i = 1; i < codeword_length + 1; i++)
	    {
		sub_top_f[i] = corridor_top[i];
		sub_bottom_f[i] = corridor_bottom[i];
	    }
      
    	    for(int i = 0; i < codeword_length; i++)
	    {   
		sub_top_b[i] = sub_top_f[i + 1] - 1;
		sub_bottom_b[i] = sub_bottom_f[i + 1] - 1;
	    }
	}

  // Adjusted boundaries of all the decoding corridors 
	for(int i = 0; i < codeword_length + 1; i++)
	{
    	    corridor_top[i] = corridor_top[i] - corridor_center[i] + half_corridor_width;
    	    corridor_bottom[i] = corridor_bottom[i] - corridor_center[i] + half_corridor_width;
    	    del_top_f[i] = del_top_f[i] - corridor_center[i] + half_corridor_width;
    	    del_bottom_f[i] = del_bottom_f[i] - corridor_center[i] + half_corridor_width;
	    del_top_b[i] = del_top_b[i] - corridor_center[i] + half_corridor_width;
    	    del_bottom_b[i] = del_bottom_b[i] - corridor_center[i] + half_corridor_width;
    	    sub_top_f[i] = sub_top_f[i] - corridor_center[i] + half_corridor_width;
    	    sub_bottom_f[i] = sub_bottom_f[i] - corridor_center[i] + half_corridor_width;
    	    sub_top_b[i] = sub_top_b[i] - corridor_center[i] + half_corridor_width;
    	    sub_bottom_b[i] = sub_bottom_b[i] - corridor_center[i] + half_corridor_width;
	}

  // ****************************
  // *** Recursion parameters ***
  // ****************************  

	int temp_position; // Temp position in the codewords
	int diff_columns;  // Difference between the corridor center of 
                     // consecutive columns 

  // Probability that a bit is transmitted without synchronization errors
	const double p_trans = 1 - p_del - p_ins - p_dup; 
  // Probability that a transmitted bit is corrupted
	const double p_trans_error = (1 - p_del - p_ins - p_dup) * p_sub; 
  // Probability that a transmitted bit is transmitted without substitution error
	const double p_trans_ok = (1 - p_del - p_ins - p_dup) * (1 - p_sub); 

  // Trellis of lattices
	Matrix5 <double> trellis_of_lattices_f(input_length_symb, vpc, nb_branches, corridor_width, b+1);
	Matrix5 <double> trellis_of_lattices_b(input_length_symb, vpc, nb_branches, corridor_width, b+1);

  // *****************
  // Forward Recursion
  // *****************

  // Loop on time steps
	for(int i = 0; i < input_length_symb; i++)
	{  
      // Temp position in the codewords
    	    temp_position = i * b; 

      // First column of the lattices
      // Time 0
    	    if(i == 0)
	    {
	  // Only from state 0. Prior probabilities of the lattices. 
		for(int k = 0; k < nb_branches; k++)
		    trellis_of_lattices_f[0][0][k][half_corridor_width][0] = prior_prob[0][k];
	    }

      // First column of the lattices
      // Time > 0
    	    else
	    {
	  // Must use previous column!!!
		for(int j = (i - 1) * vpc; j < i * vpc; j++)
	    // Loop on the branches of the previous time step
		    for(int k = 0; k < nb_branches; k++)
	      // Loop on the branches of the next time step
	    		for(int k2 = 0; k2 < nb_branches; k2++)
			    for(int l = corridor_top[temp_position]; l <= corridor_bottom[temp_position]; l++)
				trellis_of_lattices_f[i][(state_diagram[j].transitions[k]) % vpc][k2][l][0] += trellis_of_lattices_f[i - 1][j % vpc][k][l][b] * prior_prob[i][k2];
	    };

      // Loop on the columns of the lattices
        for(int m = 0; m <= b; m++)
	{
	    temp_position = i * b + m; 
	    diff_columns = corridor_center[temp_position] - corridor_center[temp_position - 1];
	  
	  // Deletions (not in the first columm)
	    if(p_del > 0 && m > 0)
	    {
	        for(int l = del_top_f[temp_position]; l <= del_bottom_f[temp_position]; l++)
		    for(int j = i * vpc; j < (i + 1) * vpc; j++)
			for(int k = 0; k < nb_branches; k++)
			    trellis_of_lattices_f[i][j % vpc][k][l][m] += p_del * trellis_of_lattices_f[i][j % vpc][k][l + diff_columns][m - 1];
	    }
	  
	  // Correct transitions (not in the first column)(p_sub == 0)
	    if(p_sub == 0 && m > 0)
	    {
	        for(int l = sub_top_f[temp_position]; l <= sub_bottom_f[temp_position]; l++)
		    for(int j = i * vpc; j < (i + 1) * vpc; j++)	
			for(int k = 0; k < nb_branches; k++)
			    if(state_diagram[j].bits[k][m - 1] == corrupted_codeword[l + corridor_center[temp_position] - half_corridor_width - 1])
		    		trellis_of_lattices_f[i][j % vpc][k][l][m] += p_trans * trellis_of_lattices_f[i][j % vpc][k][l + diff_columns - 1][m - 1];
	    }
	
	  // Correct transitions and substitution errors 
	  // (not in the first column)(p_sub > 0)
	    else if (p_sub > 0 && m > 0)
	    {
	        for(int l = sub_top_f[temp_position]; l <= sub_bottom_f[temp_position]; l++)
		    for(int j = i * vpc; j < (i + 1) * vpc; j++)
			for(int k = 0; k < nb_branches; k++)
			    if(state_diagram[j].bits[k][m - 1] == corrupted_codeword[l + corridor_center[temp_position] - half_corridor_width - 1])
		    		trellis_of_lattices_f[i][j % vpc][k][l][m] += p_trans_ok * trellis_of_lattices_f[i][j % vpc][k][l + diff_columns - 1][m - 1];
			    else
		    		trellis_of_lattices_f[i][j % vpc][k][l][m] += p_trans_error * trellis_of_lattices_f[i][j % vpc][k][l + diff_columns - 1][m - 1];
	    }
  
	  // Insertions and duplications (not in the last column)
	    if(p_ins + p_dup > 0 && m < b)
	    {
	      // Loop on the rows of the corridor except the first one
	        for(int l = corridor_top[temp_position] + 1; l <= corridor_bottom[temp_position]; l++)
		{
		  // Insertions
		    if(p_ins > 0)
		    {
		        for(int j = i * vpc; j < (i + 1) * vpc; j++)
			    for(int k = 0; k < nb_branches; k++)
				trellis_of_lattices_f[i][j % vpc][k][l][m] += trellis_of_lattices_f[i][j % vpc][k][l - 1][m] * p_ins / 2;
		    }
		  
		  // Duplications
		    if(p_dup > 0)
		    {
		        for(int j = i * vpc; j < (i + 1) * vpc; j++)
			    for(int k = 0; k < nb_branches; k++)
				if(state_diagram[j].bits[k][m] == corrupted_codeword[l + corridor_center[temp_position] - half_corridor_width - 1])
				    trellis_of_lattices_f[i][j % vpc][k][l][m] += trellis_of_lattices_f[i][j % vpc][k][l - 1][m] * p_dup;
		    }
		}
	    }
	} // Loop on the columns 1 to m of the lattices

      // Normalization
        if(trellis_of_lattices_f[i][0][0][half_corridor_width][b] != 0 && trellis_of_lattices_f[i][0][0][half_corridor_width][b] < 1)
      	{
      	    for(int j = 0; j < vpc; j++)
      		for(int k = 0; k < nb_branches; k++)
      	    	    for(int l = 0; l < corridor_width; l++)
      			trellis_of_lattices_f[i][j][k][l][b] *= 10e100;
       	}

    } // Loop on time steps
    
  // ******************
  // Backward recursion
  // ******************

  // Loop on time steps
     for(int i = input_length_symb - 1; i >= 0; i--)
     {
        temp_position = (i + 1) * b; 
       
       // Last column of the lattices
       // Time t - 1
        if(i == input_length_symb - 1)
       	{
	    for(int j = i * vpc; j < (i + 1) * vpc; j++)
       		for(int k = 0; k < nb_branches; k++)
	      // Terminated conv. codes!!! Ending in state 0!!!
       		    if(state_diagram[j].transitions[k] == 0)
			trellis_of_lattices_b[i][j % vpc][k][half_corridor_width][b] = prior_prob[i][k];
	}
       
       // Last column of the lattices
       // Time < t - 1
    	 else
	 {
	    for(int j = i * vpc; j < (i + 1) * vpc; j++)
	     // Loop on the branches of the previous time step
	        for(int k = 0; k < nb_branches; k++)
	       // Loop on the branches of the next time step
	    	    for(int k2 = 0; k2 < nb_branches; k2++)
			for(int l = corridor_top[temp_position]; l <= corridor_bottom[temp_position]; l++)
			    trellis_of_lattices_b[i][j % vpc][k][l][b] += trellis_of_lattices_b[i + 1][(state_diagram[j].transitions[k]) % vpc][k2][l][0] * prior_prob[i][k];
	 }
       
       // Loop on the columns m - 1 to 0 in the lattice
    	 for(int m = b - 1; m >= 0; m--)
	 {
	     temp_position = i * b + m; 
	     diff_columns = corridor_center[temp_position + 1] - corridor_center[temp_position];

	   // Deletions
	     if(p_del > 0)
	     {
	        for(int l = del_top_b[temp_position]; l <= del_bottom_b[temp_position]; l++)
		    for(int j = i * vpc; j < (i + 1) * vpc; j++)
			for(int k = 0; k < nb_branches; k++)
		    	    trellis_of_lattices_b[i][j % vpc][k][l][m] += p_del * trellis_of_lattices_b[i][j % vpc][k][l - diff_columns][m + 1];
	     }
      
	   // Correct transitions (p_sub == 0)
	     if(p_sub == 0)
	     {
	        for(int l = sub_top_b[temp_position]; l <= sub_bottom_b[temp_position]; l++)
		    for(int j = i * vpc; j < (i + 1) * vpc; j++)
			for(int k = 0; k < nb_branches; k++)
		    	    if(state_diagram[j].bits[k][m] == corrupted_codeword[l + corridor_center[temp_position]- half_corridor_width])
		    		trellis_of_lattices_b[i][j % vpc][k][l][m] += p_trans * trellis_of_lattices_b[i][j % vpc][k][l - diff_columns + 1][m + 1];
	     }
	
	   // Correct transitions and substitution errors (p_sub > 0)
	     else
	     {
	        for(int l = sub_top_b[temp_position]; l <= sub_bottom_b[temp_position]; l++)
		    for(int j = i * vpc; j < (i + 1) * vpc; j++)
			for(int k = 0; k < nb_branches; k++)
		    	    if(state_diagram[j].bits[k][m] == corrupted_codeword[l + corridor_center[temp_position] - half_corridor_width])
		    		trellis_of_lattices_b[i][j % vpc][k][l][m] += p_trans_ok * trellis_of_lattices_b[i][j % vpc][k][l - diff_columns + 1][m + 1];
		    	    else
		    		trellis_of_lattices_b[i][j % vpc][k][l][m] += p_trans_error * trellis_of_lattices_b[i][j % vpc][k][l - diff_columns + 1][m + 1];
	     }

	   // Insertions and duplications
	     if(p_ins + p_dup > 0)
	     {
	       // Loop on the rows of the corridor except the first one
	    	 for(int l = corridor_bottom[temp_position] - 1; l >= corridor_top[temp_position]; l--)
		 {
		   // Insertions
		     if(p_ins > 0)
		     {
		        for(int j = i * vpc; j < (i + 1) * vpc; j++)
			    for(int k = 0; k < nb_branches; k++)
				trellis_of_lattices_b[i][j % vpc][k][l][m] += p_ins / 2 * trellis_of_lattices_b[i][j % vpc][k][l + 1][m];
		     }
	      
		   // Duplications
		     if(p_dup > 0)
		     {
		    	for(int j = i * vpc; j < (i + 1) * vpc; j++)
			    for(int k = 0; k < nb_branches; k++)
				if(state_diagram[j].bits[k][m] == corrupted_codeword[l + corridor_center[temp_position] - half_corridor_width])
			    	    trellis_of_lattices_b[i][j % vpc][k][l][m] += trellis_of_lattices_b[i][j % vpc][k][l + 1][m] * p_dup;
		     }
		 }
	     }

	 }  // Loop on the columns m - 1 to 0 in the lattice

       // Normalization
    	 if(trellis_of_lattices_b[i][0][0][half_corridor_width][0] != 0 && trellis_of_lattices_b[i][0][0][half_corridor_width][0] < 1)
       	 {
       	    for(int j = 0; j < vpc; j++)
       	        for(int k = 0; k < nb_branches; k++)
       	    	    for(int l = 0; l < corridor_width; l++)
       			trellis_of_lattices_b[i][j][k][l][0] *= 10e100;
       	 }

     } // Loop on time steps

  // **************************
  // A posteriori probabilities
  // **************************

    Matrix3 <double> sumprobmergeback(input_length_symb, vpc, corridor_width);
  
  // Sum all merging branches, backward recursion
    for(int i = 1; i < input_length_symb; i++)
    {
        temp_position = b * i;
      
        for(int j = 0; j < vpc; j++)
	    for(int k = 0; k < nb_branches; k++)
		for(int l = corridor_top[temp_position]; l <= corridor_bottom[temp_position]; l++)
		    sumprobmergeback[i][j][l] += trellis_of_lattices_b[i][j][k][l][0]; 
    }
  
  // Computes the posterior probabilities of the branches in the trellis
  
  // Vector containing the posterior probabilities of the branches in the trellis
    Matrix3 <double> app(input_length_symb, vpc, nb_branches);
  
  // Loop on all the time steps in the trellis except the last one
    for(int i = 0; i < (input_length_symb - 1); i++)
    {
        temp_position = b * (i + 1);
      
        for(int j = i * vpc; j < (i + 1) * vpc; j++)
	    for(int k = 0; k < nb_branches; k++)
		for(int l = corridor_top[temp_position]; l <= corridor_bottom[temp_position]; l++)
		    app[i][j % vpc][k] += trellis_of_lattices_f[i][j % vpc][k][l][b] * sumprobmergeback[i + 1][(state_diagram[j].transitions[k]) % vpc][l]; 
    }

  // Last time step in the trellis
  // We use the fact that the code is terminated
     for(int j = (input_length_symb - 1) * vpc; 
        j < input_length_symb * vpc; j++)
    
     for(int k = 0; k < nb_branches; k++)
        if(state_diagram[j].transitions[k] == 0)
	    app[input_length_symb - 1][j % vpc][k] = trellis_of_lattices_f[input_length_symb - 1][j % vpc][k][half_corridor_width][b];

  // Posterior probabilities on the input bits (not normalized)
    vector < vector <double> > sum_app(input_length_symb, vector<double> (nb_branches, 0));

    for(int i = 0; i < input_length_symb; i++)
	for(int j = 0; j < vpc; j++)
    	    for(int k = 0; k < nb_branches; k++)
		sum_app[i][k] += app[i][j][k];
  
  // Posterior probabilities on the input bits (normalized)
    double sum_temp;
  
    for(int i = 0; i < input_length_symb; i++)
    {
        sum_temp = 0;

        for(int k = 0; k < nb_branches; k++)
	    sum_temp += sum_app[i][k];
        for(int k = 0; k < nb_branches; k++)
	    decoded_app[i][k] = sum_app[i][k] / sum_temp;
    } 
}
