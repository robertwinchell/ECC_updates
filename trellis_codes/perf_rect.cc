// ******************************************************************************embedded_str

#include "definitions.h"
#include "basic_functions.h"
#include "sync_functions.h"
#include "decoders.h"

// ******************************************************************************

// Evaluates the performance of rectangular codes over the NSC with 
// our MAP decoding algorith.

// IMPORTANT NOTES: 

// Uses nonbinary input alphabets of the form 2^k and binary output alphabets

// We use the following shortcut: to minimize decoding complexity
// and memory requirements, the width of the decoding corridor is set at 2 times 
// the maximum drift of the channel for each specific codeword. 
// The decoder, of course, has no way to know this maximum drift in practice.

// The codes are terminated in State 0

// In all cases, the encoders / decoders use a rectangular graph.

int main(int argc, char * argv[])
{

  // **********************
  // Pseudorandom Generator
  // **********************

  // We use the Mersenne Twister Generator from the GNU Scientific Library
  // WHEN USED IN PARALLEL, IT IS FUNDAMENTAL TO USE ONE SEED PER NODE!
  gsl_rng_default_seed = time(NULL);  // Random Seed
  gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // Generator

  // *********
  // Variables
  // *********
  
  // Number of states in the state diagram
  int nb_states;

  // Memory order required to put the encoder back in the all-zero state
  // (The number of NONBINARY symbols required to terminate the code)
  int memory_order;

  int k;     // Number of outgoing edges per vertex
  int vpc;   // Number of vertices per column

  // The state diagram
  Node * state_diagram; 

  int input_length_symb;      
  // Number of columns in the graph (minimum shortcut distance)
  // Also corresponds to the number of nonbinary input symbols, including the
  // termination symbols
      
  // Encoder of rate a:b (IN BITS)
  int a, b;
  
  double p_d;   // Parameters of the NSC
  double p_i;
  double p_t;
  double p_s;
      
  bool emb;   // Embedded pseudorandom sequence (if 1)

  // ****************
  // Input parameters
  // ****************

  if(argc != 7)
    {
      cout << "Performance of rectangular graphs over the NSC" << endl
	   << "The input arguments are:" << endl
	   << "1 : The input file containing the rectangular graph and the "
	   << "code bits " << endl
	   << "2 - 5 : Parameters of the NSC" << endl 
	   << "6 : Embedded pseudoransom sequences [0 for no or 1 yes] " 
	   << endl << endl;

      exit(1);
    }

  char * input_file_name = argv[1];
  p_d = atof(argv[2]);
  p_i = atof(argv[3]);
  p_t = atof(argv[4]);
  p_s = atof(argv[5]);
  emb = atoi(argv[6]);

  // *****************************
  // Graph structure and code bits
  // *****************************
  
  // File containing the rectangular graph
  ifstream graph_file(input_file_name, ios::in);
  if(!graph_file)
    {
      cout << "The file containing the encoder and code bits was not found." 
	   << endl;
      exit(1);
    }

  graph_file >> k;
  graph_file >> vpc;
  graph_file >> input_length_symb;
  graph_file >> b;
      
  a = int(round(log(k) / log(2)));
  
  nb_states = vpc * input_length_symb;
  memory_order = round(log(vpc) / log(k));
  state_diagram = new Node[nb_states];
  
  // Get structure of the encoder
  int temp;
  
  for(int i = 0; i < nb_states; i++)
    for(int j = 0; j < k; j++)
      {
	graph_file >> temp;
	graph_file >> state_diagram[i].transitions[j];
	
	for(int bit = 0; bit < b; bit++)
	  graph_file >> state_diagram[i].bits[j][bit];
      }
  
  graph_file.close();
      
  cout << "Performance test - rectangular graphs over the NSC" << endl
       << "Encoder file : " << argv[1] << endl
       << "k = " << k << endl << vpc << " vertices per columns" << endl
       << input_length_symb << " input symbols " << endl
       << b << " output bits per input symbol" << endl;
  if(emb)
    cout << "Embedded pseudorandom sequences" << endl;
  else
    cout << "No embedded pseudorandom sequences" << endl;
  cout << "p_d = " << p_d << "; p_i = " << p_i
       << "; p_t = " << p_t << "; p_s = " << p_s << endl << endl;
  
  // *************************
  // Variables for simulations
  // *************************

  // Maximum distance from the expected synchronization drift
  int drift;
  int * temp_drift = new int;
  
  // Input sequence (IN NONBINARY SYMBOLS)
  int input_seq[input_length_symb];
  
  // Length of the codewords (IN BITS)
  const int codeword_length = input_length_symb * b;
  bool initial_codeword[codeword_length];
  bool corrupted_codeword[codeword_length];

  int corrupted_length; // Length of the corrupted codewords

  // Embedded sequence
  int embedded_str[input_length_symb * b];
  
  // Number of symbol errors (in the input sequence) 
  unsigned long int total_error_bit = 0;

  // Prior propabilities for the input symbols
  vector < vector <double> > prior_prob(input_length_symb, 
					vector<double> (k, 1 / double(k)));
  // App for the received symbols
  vector < vector <double> > decoded_app(input_length_symb, 
					 vector<double> (k, 0)); 
  // App for the received bits
  vector < vector <double> > decoded_app_bits(input_length_symb * a, 
					      vector<double> (2, 0));

  // **********************
  // Performance simulation
  // **********************  

  for(int nb_sim = 1; nb_sim <= 50000; nb_sim++)
    {

      // ******************************
      // Embedded pseudorandom sequence
      // ******************************
      
      // Only required for watermark codes and specific convolutional encoders
      // We use a different embedded sequence for each codeword
        if(emb)
      	{
      	  // generate random string
	  generate_rand_string(embedded_str, input_length_symb * b, 2, r);

      	  for(int c = 0; c < input_length_symb; c++)
      	    for(int v = 0; v < vpc; v++)
      	      for(int i = 0; i < k; i++)
      		for(int j = 0; j < b; j++)
      		  state_diagram[c * vpc + v].bits[i][j] = state_diagram[c * vpc + v].bits[i][j] ^ embedded_str[c * b + j];
      	}
      
      // Generate a random input sequence
      generate_rand_string(input_seq, input_length_symb, k, r);

      // Termination SYMBOLS
      for(int i = 0; i < memory_order; i++)
	input_seq[input_length_symb - i - 1] = 0;
      
      // Encode the input sequence
      encoder(b, input_seq, input_length_symb, state_diagram, initial_codeword);
      
      // Copy the initial codeword
      for(int i = 0; i < codeword_length; i++)
	corrupted_codeword[i] = initial_codeword[i];
      
      // The codeword is corrupted by the NSC
      corrupted_length = binary_NSC_drift(initial_codeword, corrupted_codeword, 
					  codeword_length, p_d, p_i, p_t, p_s, 
					  temp_drift, r); 
      // Copy max drift
      drift = * temp_drift;
      
      // Decode received sequence
      MAP_NSC_rect(k, input_length_symb, 2 * drift, b, p_d, p_i, p_t, p_s, 
		   state_diagram, vpc, corrupted_codeword, corrupted_length, 
		   prior_prob, decoded_app);

      // **************
      // Bit Error Rate
      // **************

      for(int i = 0; i < input_length_symb * a; i++)
	{
	  decoded_app_bits[i][0] = 0;
	  decoded_app_bits[i][1] = 0;
	}

      // Change the symbol apps from symbols to bits.
      for(int i = 0; i < input_length_symb; i++)
	{
	  for(int kk = 0; kk < k; kk++)
	    {
	      for(int j = 0; j < a; j++)
		{
		  if((kk & (1 << (a - j - 1))) == 0)
		    decoded_app_bits[i * a + j][0] += decoded_app[i][kk];
		  else
		    decoded_app_bits[i * a + j][1] += decoded_app[i][kk];
		}
	    }
	}

      // Calculate total number of errors
      for(int i = 0; i < input_length_symb; i++)
	for(int j = 0; j < a; j++)
	  if((((input_seq[i] & (1 << (a - j - 1))) == 0) && (decoded_app_bits[i * a + j][0] < 0.5)) || (((input_seq[i] & (1 << (a - j - 1))) != 0) 
		 && (decoded_app_bits[i * a + j][1] <= 0.5)))
	    total_error_bit++;
      
      if(nb_sim % 10 == 0)
	cout << nb_sim << " codewords; " << "Bit error rate = " 
	     << setprecision(10) << (double(total_error_bit) / double(a * (input_length_symb - memory_order))) / double(nb_sim) << endl;
    }

    delete temp_drift;
    delete state_diagram;
}
