# include "definitions.h"
# include "basic_functions.h"
#include <gsl/gsl_matrix.h>
// ******************************************************************************
// Stores the transitions and the bits of a feedforward convolutional encoder with rate R = a : b in a state diagram
// It is assumed that the overall constraint length of the encoder is "a" times the memory order.
// Convention: register 1101 = 13 (from left to right) the bits are pushed from left to rignt
// nb_states: number of states in the state diagram
// memory_order: the memory order of the encoder
// conv_array: convolutions
//             example: for the R = 2 / 3 encoder with
//             [[1 + D   D   1 + D] [D  1  1]],
//             conv_array = {1,1,0,1,1,1,0,1,1,0,1,0}
// Output: modifies state_diagram

void create_trellis(const int a, const int b, const int memory_order, Node state_diagram[], const int nb_states){

    bool output_temp[b];
    bool temp;
    bool init_data[(memory_order + 1) * a * b];
    int n = 1;

//Generate the inner trellis code(This must not the destroy the Property of the Marcov chain)
    double rate_code = (double)a / (double)b;
    int max = (memory_order + 1) * a * b;
    int first_num = round((double)max * rate_code);
    for(int i = 0; i < first_num; i++){
	init_data[i] = 1;
    }

    for(int j = first_num; j < max; j++){
	init_data[j] = 0;
    }

    gsl_matrix_int * matrix_trellis = gsl_matrix_int_alloc(max + 1, max + 1);
    
    for(int n = 0; n < max + 1; n++){
	for(int m = 0; m < max + 1; m ++){
	    gsl_matrix_int_set(matrix_trellis, n, m, init_data[m]);
	}
    }

    gsl_matrix_int_transpose(matrix_trellis);

    for(int i = 0; i < max; i++){
	init_data[i] = gsl_matrix_int_get(matrix_trellis, 1, i);
    }
    
    
    int input_plus_state; // concatenation of the input its and the state of the registers 

    if(1 << (a * memory_order) != nb_states)
    {
        cout << endl 
	    << "ERROR: The memory order and the number of states in the state "
	    << "diagram do not match." << endl;
        exit(1);
    }

  // loop on all the possible states of the registers
    for(int cur_state = 0; cur_state < nb_states; cur_state++)
    {
      // Loop on all the possible transtions from the current state
        for(int cur_tran = 0; cur_tran < (1 << a); cur_tran++)
	{
	  // next state of the transition
	    state_diagram[cur_state].transitions[cur_tran] = (cur_state >> a) | (cur_tran << (a * (memory_order - 1)));
      	    input_plus_state = cur_state | cur_tran << (a * memory_order);
	 
	  // Initianlize temp outputs
	    for (int i = 0; i < b; i++) output_temp[i] = 0;
	
	    for(int i = 0; i < memory_order + 1; i++)
	        for(int j = 0; j < a; j++)
	        {
		    temp = (input_plus_state) & (1 << ((a * (memory_order + 1)) - 1 - i * a - j));
		     
		    for (int k = 0; k < b; k++)
		    {
			output_temp[k] = output_temp[k] ^ (temp & init_data[(j * b + k) * (memory_order + 1) + i]);
		    }
	        }

	    for(int i = 0; i < b; i ++){
		state_diagram[cur_state].bits[cur_tran][i] = output_temp[i];
	    }
 	} // Loop on all the possible transtions from the current state

    }
    return;
}

void conv_to_trellis(const int a, const int b, const int memory_order, 
		     const bool conv_array[], Node state_diagram[], 
		     const int nb_states) 
{
    bool output_temp[b];
    bool temp;
    int input_plus_state; // concatenation of the input its and the state 
                        // of the registers 

    if(1 << (a * memory_order) != nb_states)
    {
      cout << endl 
	   << "ERROR: The memory order and the number of states in the state "
	   << "diagram do not match." << endl;
      exit(1);
    }

  // loop on all the possible states of the registers
    for(int cur_state = 0; cur_state < nb_states; cur_state++)
    {

      // Loop on all the possible transtions from the current state
        for(int cur_tran = 0; cur_tran < (1 << a); cur_tran++)
	{
	
	  // next state of the transition
	  state_diagram[cur_state].transitions[cur_tran] = (cur_state >> a) | (cur_tran << (a * (memory_order - 1)));
      	  input_plus_state = cur_state | cur_tran << (a * memory_order);
	 
	  // Initianlize temp outputs
	    for (int i = 0; i < b; i++)
		output_temp[i] = 0;
	
	    for(int i = 0; i < memory_order + 1; i++)
	        for(int j = 0; j < a; j++)
	        {

		    temp = (input_plus_state)& (1 << 
		      ( 	// I ADDED PARENTHESES TO REMOVE C++ WARNING!!!
		      (a * (memory_order + 1)) - 1 - i * a - j
			)  // I ADDED PARENTHESES TO REMOVE C++ WARNING!!!
		     );
		     
		
		  for (int k = 0; k < b; k++)
		  {
			output_temp[k] = output_temp[k] ^ (temp & conv_array[(j * b + k) * (memory_order + 1) + i]);
		  }
	      }

	    for(int i = 0; i < b; i++)
	    {
	      // cout << output_temp[i];
	      state_diagram[cur_state].bits[cur_tran][i] = output_temp[i];
	    }

 	} // Loop on all the possible transtions from the current state

    } // loop on all the possible states of the registers
}

// ******************************************************************************
// Prints a boolean array of length n
void print_array(bool * input_array, int n)
{
  for(int i = 0; i < n; i++)
    {
      cout << input_array[i];
    }

  cout << endl;
}


// ******************************************************************************
// This procedure prints a trellis. The input arguments are the trellis and 
// the number of nodes.
void print_trellis(vector < vector <DecodingNode> > trellis, int n)
{
    if(trellis.size()==0)
    {
        cout << "Empty trellis. Cannot print.";
        return;
    }
  
    cout << endl << "***** Trellis *****" << endl;
  
    for(int i = 0; i < n; i++)
    {
      
        for(unsigned int j = 0; j < trellis.size(); j++)
	{
	    cout << setw(3) << trellis[j][i].metric;
	    cout << setw(3) << trellis[j][i].ancestor;
	    cout << setw(3) << trellis[j][i].nb_opt_branches;
	    cout << "|";
	}
	cout << endl;
    }
    cout << endl;
};


// ******************************************************************************
// Computes the Hamming distance of two binary strings of length n
int hamming(const bool * str1, const bool * str2, int length)
{
  int hamming_dist = 0;

  for(int i = 0; i < length; i++)
    {
      if(str1[i] != str2[i])
	hamming_dist++;
    }

  return hamming_dist;
}

// ******************************************************************************
// Generates a uniform random binary string of length n

// OLDER FUNCTION. IT USES THE DEFAULT GENERATOR AND REDUCES THE NUMBER OF 
// FUNCTION CALLS TO SPEED UP RUNNING TIME


// ******************************************************************************
// Generates a uniform random k-ary string of length str_length

void generate_rand_string(int * random_string, const int str_length, const int k, gsl_rng * r)
{
    for(int i = 0; i < str_length; i++)
    {
//	random_string[i] = gsl_rng_uniform_int(r, k);
	random_string[i] = i % 3;
    }
}


// ******************************************************************************
// Takes an input string and encodes using a bit distribution over a state 
// diagram.
// Assumes binary alphapbet.

// Note: it can be used for convolutional codes with feedforward shift 
// registers, but if the code is terminated, then one must include enough zeros 
// at the end of the input sequences in order to put the register back in the 
// all-zero state. 

// input: the input string
// length: the length of the input string
// trellis: the state diagram
// n: number of nodes

void encoder(const int b, const int * input_str, const int input_length, const Node * state_diagram, bool * output)
{
    int state = 0;  // Initial state of the encoder
    int tran ;

    for(int i = 0; i < input_length; i++) 
    // loop on all the input bits
    {
        tran = input_str[i];  // Transition taken
        for(int j = 0; j < b; j++) // outputs b bits for each a input bits
	    output[b * i + j] = state_diagram[state].bits[tran][j];
        state = state_diagram[state].transitions[tran];
    }

  // Sanity test to check if the encoder finishes in the all-zero state.
  // This will not necessarily occur if the code is not terminated.
    if(state != 0)
    {
        cout << "NO STATE 0" << endl;
        exit(1);
    }

};


// ******************************************************************************
// Computes the number of mismatched blocks between two binary strings
int mismatched_blocks(const bool * str1, const bool * str2, int length, 
		      int block_size)
{
  
  int mismatches = 0;
  int nb_blocks = length / block_size; // The number of complete blocks
  
  for(int i = 0; i < nb_blocks; i++)
    {
      for(int j = 0; j < block_size; j++)
	{
	  if(str1[i * block_size + j] != str2[i * block_size + j])
	    {
	      mismatches++;
	      break;
	    }
	}
    }
  
  // Last block
  for(int i = 0; i < length % block_size; i++)
    {
      if(str1[length - i - 1] != str2[length - i - 1])
	{
	  mismatches++;
	  break;
	}
      
    }
  
  return mismatches;

}


// ******************************************************************************
// Computes the number of successes before a failure occurs. 
// The prob. of success is given as an input argument
int geom(const double prob)
{
  unsigned long int temp = rand();

  if (temp > double(prob * RAND_MAX))
    return 0;
  else
      return geom(prob) + 1;
}


// ******************************************************************************
// Computes the average binary entropy of a vector of probabilities
double entropy_vector(double * prob_vector, int vector_length)
{
  double entropy = 0;
  double p;

  for(int i = 0; i < vector_length; i++)
    {
      p = prob_vector[i]; 
      
      if(p > 0 && p < 1)
	entropy = entropy - p * log(p) / log(2) - (1 - p) * log(1 - p) / log(2);
    }
  
  return entropy / vector_length;

}

// ******************************************************************************
// Pseudorandom interleaver based on a quadratic congruence 
// (see Lin-Costello, P.778)
//
// k is the parameter of the congruence
// 
// condition 0 : interleaving
// condition 1 : interleaving
//
// NOTE: Only works when the length of the vector is a power of two
 
void qc_interleaver_bool(const int length_vector, const bool * vector1, 
			 bool * vector2, const int condition) 
{
  int k = 1;  // Parameter of the quadratic congruence (see Lin-Costello, P.778)
  int prev = 0;
  int next;

  for(int m = 1; m < length_vector; m++)
    {
      next = ((m * (m + 1) * k) / 2) % length_vector; 
      
      // Interleaving
      if(condition == 0)
	vector2[prev] = vector1[next];
      // Deinterleaving
      else 
	vector2[next] = vector1[prev];
      
      prev = next; 
    }
  
  // Interleaving
  if(condition == 0)
    vector2[prev] = vector1[0];
  // Deinterleaving
  else 
    vector2[0] = vector1[prev];
}

void qc_interleaver_double(const int length_vector, const double * vector1, 
			   double * vector2, const int condition) 
{
  int k = 1;  // Parameter of the quadratic congruence (see Lin-Costello, P.778)
  int prev = 0;
  int next;

  for(int m = 1; m < length_vector; m++)
    {
      next = ((m * (m + 1) * k) / 2) % length_vector; 
 
      // Interleaving
      if(condition == 0)
	vector2[prev] = vector1[next];
      // Deinterleaving
      else 
	vector2[next] = vector1[prev];
      
      prev = next; 
    }
  
  // Interleaving
  if(condition == 0)
    vector2[prev] = vector1[0];
  // Deinterleaving
  else 
    vector2[0] = vector1[prev];
}

void qc_interleaver_vec(const int length_vector, 
			const double * vector1, 
			vector <double>  vector2, const int condition) 
{
  int k = 1;  // Parameter of the quadratic congruence (see Lin-Costello, P.778)
  int prev = 0;
  int next;

  for(int m = 1; m < length_vector; m++)
    {
      next = ((m * (m + 1) * k) / 2) % length_vector; 
 
      // Interleaving
      if(condition == 0)
	vector2[prev] = vector1[next];
      // Deinterleaving
      else 
	vector2[next] = vector1[prev];
      
      prev = next; 
    }
  
  // Interleaving
  if(condition == 0)
    vector2[prev] = vector1[0];
  // Deinterleaving
  else 
    vector2[0] = vector1[prev];
}


// *******************************************************************************
// Kronecker Delta

bool kd(const int b)
{
	return !b;
}


// *******************************************************************************
// Reverse the binary expression of an unsigned integer
// 
// v : the integer to reverse
// length: the number of digits in the decimal representation. This avoids
// reversing dummy zeroes. For instance, reversing 1011 with length=4 gives 1101.

// c : the reverse of v

// Uses a look-up table and reverses one byte at a time.

// Important note: assumes 4 bytes per integer

unsigned int bit_reverse(unsigned int v, int length = -1)
{
	static const unsigned char BitReverseTable256[256] = 
	{
    	    #define R2(n)     n,     n + 2*64,     n + 1*64,     n + 3*64
    	    #define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
    	    #define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
                  R6(0), R6(2), R6(1), R6(3)
	};
	unsigned int answer = ((BitReverseTable256[v & 0xff] << 24) | 
		(BitReverseTable256[(v >> 8) & 0xff] << 16) | 
		(BitReverseTable256[(v >> 16) & 0xff] << 8) |
		(BitReverseTable256[(v >> 24) & 0xff]));
  
	if(length > 0 || length < 32)
	    answer = answer >> (32 - length);
	return answer;
}


// *******************************************************************************
// This function generates rectangular graphs and returns them in an output file.
// Both vertices of each edge e=(va,vb) of the graph are included in the file:

// va1 vb1
// va2 vb2
// va3 vb3
//   ...
 
// Input arguments
// k : number of outgoing edges per vertex
// vpc : number of vertices per column
// s: number of columns in the graph (minimum shortcut distance)

// Output file names have the form r_k_vpc_s.txt

void gen_rect(const int k, const int vpc, const int s)
{
	if(k < 2 || vpc < 1 || s < 1)
	{
    	    cout << "Invalid input parameters. " << endl;
    	    exit(1);
	}
  
	int b = round(log(vpc) / log(k)) + 1; // Minimum bypass distance
  
	if(vpc != pow(k, b-1))
	{
    	    cout << "Invalid input parameter" << endl
		<< "The number of vertices per column must be a power of " << k 
		<< endl; 
    	    exit(1);
	}

  // Output files have the form r_k_vpc_s.txt
	stringstream sstemp;
	sstemp << "r_" << k << "_" << vpc << "_" << s << ".txt";
	string fname = sstemp.str();

  // Output file for the graph
	ofstream output_file(fname.c_str(), ios::out);
	output_file << k << " " << vpc << " " << s << endl << endl;
 
  // Graph transitions
	for(int col = 0; col < s; col++)
	    for(int v = 0; v < vpc; v++)
    		for(int t = 0; t < k; t++)
		    output_file << col * vpc + v << " " << (k * v + t) % vpc + ((col + 1) % s) * vpc  << endl;
 
	output_file.close();
}

// ******************************************************************************

void gen_markov_seq_init(const int length,  bool codeword[], const int order, 
			 const long double p[], const int init_state,
			 gsl_rng * r)
  
// This function generates a random binary sequence following a Markov chain 
// whose initial state is known.

// - length: the length of the sequence
// - codeword: the table containing the sequence
// - order: the order of the Markov chain
// - p: a vector containing the probability distributions for the 
//      transitions of the markov chains
// - init_state: the initial state of the Markov chain
{
    int state; // Temporary state
  
  // Markov chain of order 0 - A random string is generated 
    if(order == 0)
    {
        for(int i = 0; i < length; i++)
	    codeword[i] = gsl_rng_uniform_int(r, 2);
    }
  // Markov chain of higher order
    else
    {
      // Initial state from the input argument
        state = init_state;

      // String generation
    	for(int i = 0; i < length; i++)
	{
	  // State starts with a 0
	    if((state << 1) < (1 << order))
	    {
	        if(gsl_rng_uniform(r) <= p[state])
		{
		    state = (state << 1);
		    codeword[i] = 0;
		}
	        else
		{
		    state = (state << 1) + 1;
		    codeword[i] = 1;
		}
	    }
	  // State starts with a 1
	    else
	    {
	    	if(gsl_rng_uniform(r) <= p[state ^ ((1 << order) - 1)])
		{
		    state = ((state << 1) + 1) & ((1 << order) - 1);
		    codeword[i] = 1;
		}
	    	else
		{
		    state = ((state << 1)) & ((1 << order) - 1);
		    codeword[i] = 0;
		}
	    }
	} // Loop on the bits of the sequence
    } // Markov chain of higher order
}

// ******************************************************************************

// This function generates a code over a rectangular graph

// Input arguments:
// option : main option
// k : number of branches per vertex (input alphabet size)
// vpc : number of vertices per column in the graph (related to the memory order)
// input_length_symb : the number of columns in the graph (# of input symbols)
// b : number of output bits per input symbol (code of rate a:b)
// conv_array: The array containing the convolutions (for conv. codes only)
// markov_tran : the transition of the Markov process modeling the input source

Node * rect_code(const int option, const int k, const int vpc, const int input_length_symb, const int b, 
		const bool conv_array[], const long double markov_tran[], gsl_rng * r)
{
	gsl_rng_default_seed = time(NULL);  		// Random Seed

	if(1 << int(round(log(k) / log(2))) != k)
	{
    	    cout << "k is not a power of 2. " 
		<< "Only binary alphabets are curently supported" << endl;
    	    exit(1);
	}
 
  // Encoder of rate a:b (IN BITS)
	int a = int(round(log(k) / log(2)));

	if(a > MAX_a || b > MAX_b) 
	{
    	    cout << endl 
		<< "ERROR: Adjust MAX_a and/or MAX_b" << endl;
    	    exit(1);
	}

	int nb_states = vpc * input_length_symb;  // Nb of states in the encoder

  // Endoder (structure and bits)
	Node * state_diagram = new Node[nb_states];

  // Memory order required to put the encoder back in the all-zero state
  // (The number of NONBINARY symbols required to terminate the code)
	int memory_order = round(log(vpc) / log(k)); // Memory order of the code

  // Watermark codes
	if(option == 1)
	{
      // File containing the trellis encoder
      // Its name has the form r_k_vpc_$(input_length_symb).txt
    	    stringstream sstemp;
    	    sstemp << "r_" << k << "_" << vpc << "_" << input_length_symb << ".txt";
    	    string fname = sstemp.str();
    	    ifstream tempf(fname.c_str(), ios::in);

      // If the graph does not exist, it is created.
    	    if (!tempf) gen_rect(k, vpc, input_length_symb);
	    else tempf.close();
      
      // File containing the rectangular graph
    	    ifstream graph_file(fname.c_str(), ios::in);
    	    if(!graph_file)
	    {
		cout << "The file containing the state diagram was not found." << endl;
		exit(1);
	    }
      
    	    stringstream sstemp2;
      
      // File containing the bits of the watermark code 
      // Its name has the form watermark_a_b.txt
    	    sstemp2 << "watermark_" << a << "_" << b << ".txt";
    	    fname = sstemp2.str();
      
    	    ifstream watermark_file(fname.c_str(), ios::in);
    	    if(!watermark_file)
	    {
		cout << "The file containing the watermark bits was not found." << endl;
		exit(1);
	    }
      
      // Read the graph info at the beginning of the file. 
      // It allows the buffer to go directly to the graph itself
    	    int temp;
    	    graph_file >> temp;
    	    graph_file >> temp;
    	    graph_file >> temp;
           
      // Get structure of the encoder
    	    for(int i = 0; i < nb_states; i++)
	    {
		for(int j = 0; j < k; j++)
		{
		    graph_file >> temp;
		    graph_file >> state_diagram[i].transitions[j];
		}
	    }
    	    graph_file.close();
      
      // Bits of the watermark code are coming from the file
    	    for(int i = 0; i < k; i++)
	    {
		for(int j = 0; j < b; j++){
		    watermark_file >> state_diagram[0].bits[i][j];
		}      
	    }
    	    watermark_file.close();
      
      // Copy bits of the 1st column to the other columns
    	    for(int c = 1; c < input_length_symb; c++)
		for(int i = 0; i < k; i++)
		    for(int j = 0; j < b; j++)
			state_diagram[c].bits[i][j] = state_diagram[0].bits[i][j] ;
	}

  // Convolutional codes
	else if(option == 2)
	{
      // Generate the first column (code bits and transitions) of the trellis from the convolutional encoder
    	    conv_to_trellis(a, b, memory_order, conv_array, state_diagram, vpc);
    
      // Change the transitions of the first column so that they point to the second column
    	    for(int v = 0; v < vpc; v++)  
   		for(int i = 0; i < k; i++)
   		    state_diagram[v].transitions[i] += vpc;

       // Copy the first column (bits and transitions) into the other columns
    	    for(int c = 1 ; c < input_length_symb; c++)
   		for(int v = 0; v < vpc; v++)
   		    for(int i = 0; i < k; i++)
   		    {
   	    		state_diagram[c * vpc + v].transitions[i] = (state_diagram[v].transitions[i] + c * vpc) % nb_states;
   	    		for(int j = 0; j < b; j++)
   			    state_diagram[c * vpc + v].bits[i][j] = state_diagram[v].bits[i][j];
   		    }
	}

    //Inner trellis codes
	else if(option == 3){
       // Generate the first column (code bits and transitions) of the trellis
    	    create_trellis(a, b, memory_order, state_diagram, vpc);
    
      // Change the transitions of the first column so that they point to the second column
    	    for(int v = 0; v < vpc; v++)  
   		for(int i = 0; i < k; i++)
   		    state_diagram[v].transitions[i] += vpc;

       // Copy the first column (bits and transitions) into the other columns
    	    for(int c = 1 ; c < input_length_symb; c++)
   		for(int v = 0; v < vpc; v++)
   		    for(int i = 0; i < k; i++)
   		    {
   	    		state_diagram[c * vpc + v].transitions[i] = (state_diagram[v].transitions[i] + c * vpc) % nb_states;
   	    		for(int j = 0; j < b; j++){
   			    state_diagram[c * vpc + v].bits[i][j] = state_diagram[v].bits[i][j];
			}
   		    }
	}

    //Outter ldpc codes
	else if(option == 4){
      // File containing the trellis encoder
      // Its name has the form r_k_vpc_$(input_length_symb).txt
    	    stringstream sstemp;
    	    sstemp << "r_" << k << "_" << vpc << "_" << input_length_symb << ".txt";
    	    string fname = sstemp.str();
    	    ifstream tempf(fname.c_str(), ios::in);

      // If the graph does not exist, it is created.
    	    if (!tempf) gen_rect(k, vpc, input_length_symb);
	    else tempf.close();
      
      // File containing the rectangular graph
    	    ifstream graph_file(fname.c_str(), ios::in);
    	    if(!graph_file)
	    {
		cout << "The file containing the state diagram was not found." << endl;
		exit(1);
	    }

	//Generate the outter ldpc code file(FileName : outter_ldpc_code.txt)      
	    stringstream gen_outter_ldpc;
	    gen_outter_ldpc << "./gen_outter_ldpc " << nb_states * k << " " << b << " " << gsl_rng_default_seed;
	    cout << gen_outter_ldpc.str().c_str() << endl;
	    
	    // Execute the make-ldpc for the generate the outter ldpc code.
	    system((const char*)gen_outter_ldpc.str().c_str());
	    
	    // This is the outter ldpc code file name. 

    	    ifstream outter_ldpc_file("./outter_ldpc_code.txt", ios::in);
    	    if(!outter_ldpc_file)
	    {
		cout << "The file containing the outter ldpc  bits was not found." << endl;
		exit(1);
	    }
      
      // Read the graph info at the beginning of the file.
      // It allows the buffer to go directly to the graph itself
    	    int temp;
    	    graph_file >> temp;
    	    graph_file >> temp;
    	    graph_file >> temp;
           
      // Get structure of the encoder
    	    for(int i = 0; i < nb_states; i++)
	    {
		for(int j = 0; j < k; j++)
		{
		    graph_file >> temp;
		    graph_file >> state_diagram[i].transitions[j];
		}
	    }
    	    graph_file.close();
      
    	    for(int c = 0; c < nb_states; c++)
		for(int i = 0; i < k; i++)
		    for(int j = 0; j < b; j++)
			outter_ldpc_file >> state_diagram[c].bits[i][j];
	    outter_ldpc_file.close();
	}

      // Rectangular graph with random code bits
	else if (option == 5)
	{

      // Output files have the form r_k_vpc_$(input_length_symb).txt
    	    stringstream sstemp;
    	    sstemp << "r_" << k << "_" << vpc << "_" << input_length_symb << ".txt";
    	    string fname = sstemp.str();

      // File containing the trellis encoder
    	    ifstream tempf(fname.c_str(), ios::in);  
  
      // If the graph does not exist, it is created.
    	    if (!tempf) gen_rect(k, vpc, input_length_symb);
  	    else tempf.close();
      
      // File containing the rectangular graph
    	    ifstream graph_file(fname.c_str(), ios::in);
    	    if(!graph_file)
  	    {
  		cout << "The file containing the state diagram was not found." << endl;
  		exit(1);
  	    }
      
      // Read the graph info at the beginning of the file. 
      // It allows the buffer to go directly to the graph itself
    	    int temp;
    	    graph_file >> temp;
    	    graph_file >> temp;
    	    graph_file >> temp;
      
    	    nb_states = vpc * input_length_symb;
    	    memory_order = round(log(vpc) / log(k));
    	    state_diagram = new Node[nb_states];
      
      // Get structure of the encoder
    	    for(int i = 0; i < nb_states; i++)
  		for(int j = 0; j < k; j++)
  		{
  		    graph_file >> temp;
  		    graph_file >> state_diagram[i].transitions[j];
  		}
      
    	    graph_file.close();

      // ****************************
      // Generate bits of the encoder
      // ****************************
      // Bits of the encoder are initially generated randomly from a Markov chain. 
      // All the additional machinery is used so that the outgoing edges of each vertex have different output bits.

      // Initial state for each node in a column
    	    int init_state[vpc];

      // Maximum decimal value ( + 1) of the output bits of a branch in the rect. graph
    	    int max_range = 1 << b; 
      
      // Used to avoid identical bits for multiple branches from the same node
    	    int mult_branch[max_range];

      // Loop on the columns of the rectangular graph
    	    for(int i = 0; i < input_length_symb; i++)
  	    {
  	  // Initial states
	  
  	  // First column: initial states chosen randomly
  		if(i == 0)
  		{
  	    	    for(int j = 0; j < vpc; j++)
  			init_state[j] = gsl_rng_uniform_int(r, 2);
  		}
	  
  	  // Other columns: initial states chosen from the last bits of the
  	  // previous column
  		else
  		{
  	    	    for(int j = 0; j < vpc; j++)
  			init_state[j] = 0;

  	      // Count the number of branches to a state whose last code bit is 
  	      // a zero
  	    	    for(int j = (i - 1) * vpc; j < i * vpc; j++)
  			for(int v = 0; v < k; v++)
  			    if(state_diagram[j].bits[v][b - 1] == 0)
  				init_state[(state_diagram[j].transitions[v]) % vpc]++;

  	    	    for(int j = 0; j < vpc; j++)
  		    {
  			if((double(init_state[j]) / double(k)) <= gsl_rng_uniform(r)) init_state[j] = 1;
  			else init_state[j] = 0;
  		    }
  		}

  	  // Loop on the vertices of the column
  		for(int v = i * vpc; v < (i + 1) * vpc; v++)
  		{
  	    	    for(int l = 0; l < max_range; l++)
  			mult_branch[l] = 0;
	  
  	      // Loop on the branches of the node
  	    	    for(int j = 0; j < k; j++)
  		    {
  			bool done = false;
	      
  			while(!done)
  			{
  		      // Generate bits of the branch
  		    	    gen_markov_seq_init(b, state_diagram[v].bits[j], 1, markov_tran, init_state[v % vpc], r);
  		    	    temp = 0;
		      
  		    	    for(int l = 0; l < b; l++)
  				temp = 2 * temp + state_diagram[v].bits[j][l];
		      
  		    	    if(mult_branch[temp] == 0)
  			    {
  			  // Branch is different! 
  				mult_branch[temp] = 1;
  				done = true;
  			    }
  			}
  		    }
  		}
  	    }
	}
  return state_diagram;
}
