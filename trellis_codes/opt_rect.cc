// ******************************************************************************

#include "definitions.h"
#include "basic_functions.h"
#include "sync_functions.h"
#include "decoders.h"

// ******************************************************************************

// Optimizes rectangular graphs over the NSC with our MAP decoding algorith.

// IMPORTANT NOTES: 

// Uses nonbinary input aplphabets of the form 2^k and binary output alphabets

// We use the following shortcut: to minimize decoding complexity
// and memory requirements, the width of the decoding corridor is set at 2 times 
// the maximum drift of the channel for each specific codeword. 
// The decoder, of course, has no way to know this maximum drift in practice.

// The codes are terminated in State 0

int main(int argc, char * argv[])
{

  // **********************
  // Pseudorandom Generator
  // **********************

  // We use the Mersenne Twister Generator from the GNU Scientific Library
  // WHEN USED IN PARALLEL, IT IS FUNDAMENTAL TO USE ONE SEED PER NODE!
    gsl_rng_default_seed = time(NULL);  // Random Seed
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // Generator


    if(argc != 10 && argc != 11)
    {
        cout << "Optimization of rectangular graphs over the NSC" << endl
	    << "The input arguments are:" << endl
	    << "1 : Number of outgoing edges per vertex" << endl
	    << "2 : Number of vertices per column" << endl
	    << "3 : Number of nonbinary input symbols" << endl
	    << "4 : b (encoder of rate a:b (IN BITS))" << endl
	    << "5 - 8 : Parameters of the NSC" << endl 
	    << "9 : transition of the input Markov chain of order 1" << endl
	    << "10 - The input file with the code (optional)" << endl << endl;
        exit(1);
    }
 
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
      
  // Transition probability for a Markov chain of order 1
    long double markov_tran[1];
    int nb_sim_max = 5000;  // Number of simulations per encoder

  // ****************
  // Input parameters
  // ****************
  
    k = atoi(argv[1]); 
    vpc = atoi(argv[2]);
    input_length_symb = atoi(argv[3]);      
    b = atoi(argv[4]);
    p_d = atof(argv[5]);
    p_i = atof(argv[6]);
    p_t = atof(argv[7]);
    p_s = atof(argv[8]);
    markov_tran[0] = atof(argv[9]);

    a = int(round(log(k) / log(2)));
    cout << a << endl;

    nb_states = vpc * input_length_symb;
    memory_order = round(log(vpc) / log(k));
    state_diagram = new Node[nb_states];


    if(1 << int(round(log(k)/log(2))) != k)
    {
        cout << "k is not a power of 2. "
	    << "Only binary alphabets are curently supported" << endl;
        exit(1);
    }
      
    if(a >= b)
    {
        cout << "a >= b!!!" << endl;
        exit(1);
    }
  
    if(a > MAX_a || b > MAX_b) 
    {
        cout << endl << "ERROR: Adjust MAX_a and/or MAX_b" << endl;
        cout << MAX_b;
        exit(1);
    }
  
  // *******************************
  // Rectangular graph and code bits
  // *******************************

  // ****************
  // Random code bits
  // ****************
    if(argc == 10)
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
      // Bits of the encoder are initially generated randomly from a Markov 
      // chain. All the additional machinery is used so that the outgoing edges 
      // of each vertex have different output bits.

      // Initial state for each node in a column
        int init_state[vpc];

      // Maximum decimal value ( + 1) of the output bits of a 
      // branch in the rect. graph
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

  // Code from an input file
    else if (argc == 11)
    {
      // The input file with the code is the last input argument
        char * input_file_name = argv[10];
     
      // File containing the rectangular graph 
        ifstream graph_file(input_file_name, ios::in);
        if(!graph_file)
	{
	    cout << "The file containing the encoder and code bits was not found." 
	        << endl;
	    exit(1);
	}

      // Sanity check
        int temp;
        graph_file >> temp;
        if(temp != k)
	{
	    cout << "The values for k from the input argument and the input file "
	        << "do not match" << endl;
	    exit(1);
	}
        
	graph_file >> temp;
        if(temp != vpc)
	{
	    cout << "The values for vpc from the input argument and the input file "
	        << "do not match" << endl;
	    exit(1);
	}
        
	graph_file >> temp;
        if(temp != input_length_symb)
	{
	    cout << "The values for input_length_symb from the input argument "
	        << "and the input file do not match" << endl;
	    exit(1);
	}

        graph_file >> temp;
        if(temp != b)
	{
	    cout << "The values for b from the input argument "
	        << "and the input file do not match" << endl;
	    exit(1);
	}

      // Get structure of the encoder as well as code bits
        for(int i = 0; i < nb_states; i++)
	    for(int j = 0; j < k; j++)
	    {
		graph_file >> temp;
		graph_file >> state_diagram[i].transitions[j];
	    
		for(int bit = 0; bit < b; bit++)
	    	    graph_file >> state_diagram[i].bits[j][bit];
	    }
      
        graph_file.close();
    }

    cout << "Optimization of rectangular graphs for the NSC" << endl
        << "k = " << k << endl << vpc << " vertices per columns" << endl
        << input_length_symb << " input symbols " << endl
        << b << " output bits per input symbol" << endl
        << "p_d = " << p_d << "; p_i = " << p_i
        << "; p_t = " << p_t << "; p_s = " << p_s << endl 
        << "Markov chain with transition probability " << markov_tran[0] << endl;
    if(argc == 11)
	cout << "Encoder file : " << argv[10] << endl;
    cout << nb_sim_max << " simulated codewords per encoder" << endl << endl;

  // *********************************
  // Print state diagram and code bits
  // *********************************
 
    cout << "********************" << endl << "Initial Encoder" << endl;
    for(int i = 0; i < nb_states; i++)
	for(int j = 0; j < k; j++)
        {
  	    cout << i << " -> " << state_diagram[i].transitions[j] << " : ";
  	    for(int k = 0; k < b; k++)
  		cout << state_diagram[i].bits[j][k];
  	    cout << endl;
        }

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
  
  // Number of symbol errors (in the input sequence) 
    unsigned long int total_error_bit, total_error_bit_new;

  // Prior propabilities for the input symbols
    vector < vector <double> > prior_prob(input_length_symb, vector<double> (k, 1 / double(k)));

  // App for the received symbols
    vector < vector <double> > decoded_app(input_length_symb, vector<double> (k, 0));

  // App for the received bits
    vector < vector <double> > decoded_app_bits(input_length_symb * a, vector<double> (2, 0));

  // Table containing the number of bit errors for each column
    vector <unsigned long int> error_table(input_length_symb, 0);

    bool old_column[vpc * k * b]; // Store the bits of the worst column of the code

    unsigned long int max; // Find maximum error in table
    int max_pos;           // Worst column in the original graph
    int max_pos_new;       // Worst column in the updated graph

    bool reset = true;
    int temp;

  // **********************
  // Performance simulation
  // **********************  

    for(int opt = 1; opt <= 10000; opt++)
    {
      // No need to simulate a new code for the first time
        if(reset == true)
	{
	    total_error_bit = 0;
	    for(int i = 0; i < input_length_symb; i++)
		error_table[i] = 0;
	  
	    for(int nb_sim = 1; nb_sim <= nb_sim_max; nb_sim++)
	    {
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
	        corrupted_length = binary_NSC_drift(initial_codeword, corrupted_codeword, codeword_length, p_d, p_i, p_t, p_s, temp_drift, r); 

	      // Copy max drift
	        drift = * temp_drift;
	      
	      // Decode received sequence
	        MAP_NSC_rect(k, input_length_symb, 2 * drift, b, p_d, p_i, p_t, p_s, state_diagram, vpc, corrupted_codeword, corrupted_length, prior_prob, decoded_app);
	  	  
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
			    if((kk & (1 << (a - j - 1))) == 0) decoded_app_bits[i * a + j][0] += decoded_app[i][kk];
			    else decoded_app_bits[i * a + j][1] += decoded_app[i][kk];
			}
		    }
		}
	      
	      // Calculate total number of errors
	        for(int i = 0; i < input_length_symb; i++)
		    for(int j = 0; j < a; j++)
			if((((input_seq[i] & (1 << (a - j - 1))) == 0) && (decoded_app_bits[i * a + j][0] < 0.5)) || 
		    	    (((input_seq[i] & (1 << (a - j - 1))) != 0) && (decoded_app_bits[i * a + j][1] <= 0.5)))
			{
		    	    error_table[i]++;
		    	    total_error_bit++;
			}
	      
	    } // Loop on the simulations for one code
      
	  // Find maximum error in table
	    max = 0;
	    max_pos = 0;
	  
	    for(int i = 0; i < input_length_symb; i++)
	    {
	        if(error_table[i] > max)
		{
		    max = error_table[i];
		    max_pos = i;
		}
	    }
	  
	    cout << "**************************************************" 
	        << endl << "Original code:" << endl
	        << "Maximum # of errors in column " << max_pos << " (" << max 
	        << " errors)" << endl << total_error_bit << " total errors" 
	        << endl << "Bit error rate = " << setprecision(10) 
	        << (double(total_error_bit) / double(a * (input_length_symb - memory_order))) / double(nb_sim_max) << endl;
	}

      // ***********************************
      // Change the worst column of the code
      // ***********************************

      // Store worst column
        for(int i = 0; i < vpc; i++)
	    for(int j = 0; j < k; j++)
		for(int l = 0; l < b; l++)
		    old_column[i * vpc + j * k + l] = state_diagram[max_pos * vpc + i].bits[j][l];
      
      // Generate random bits to replace the bits from the worst column 
      
      // Initial state for each node in a column
        int init_state[vpc];
      
      // Maximum decimal value ( + 1) of the output bits of a 
      // branch in the rect. graph
        int max_range = 1 << b; 
      
      // Used to avoid identical bits for multiple branches from the same node
        int mult_branch[max_range];
      
      // Initial states
      
      // First column: initial states chosen randomly
        if(max_pos == 0)
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
	    for(int j = (max_pos - 1) * vpc; j < max_pos * vpc; j++)
		for(int v = 0; v < k; v++)
	    	    if(state_diagram[j].bits[v][b - 1] == 0) init_state[(state_diagram[j].transitions[v]) % vpc]++;
	  
	    for(int j = 0; j < vpc; j++)
	    {
	        if((double(init_state[j]) / double(k)) <= gsl_rng_uniform(r))init_state[j] = 1;
	        else init_state[j] = 0;
	    }
	}
      
      // Loop on the vertices of the column
        for(int v = max_pos * vpc; v < (max_pos + 1) * vpc; v++)
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
      
      // Run simulations one more time
      
        total_error_bit_new = 0;
        for(int i = 0; i < input_length_symb; i++)
	    error_table[i] = 0;
      
        for(int nb_sim = 1; nb_sim <= nb_sim_max; nb_sim++)
	{
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
	    corrupted_length = binary_NSC_drift(initial_codeword, corrupted_codeword, codeword_length, p_d, p_i, p_t, p_s, temp_drift, r); 

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
		        if((kk & (1 << (a - j - 1))) == 0) decoded_app_bits[i * a + j][0] += decoded_app[i][kk];
		        else decoded_app_bits[i * a + j][1] += decoded_app[i][kk];
		    }
		}
	    }
	  
	  // Calculate total number of errors
	    for(int i = 0; i < input_length_symb; i++)
		for(int j = 0; j < a; j++)
	    	    if((((input_seq[i] & (1 << (a - j - 1))) == 0) && (decoded_app_bits[i * a + j][0] < 0.5)) || 
			(((input_seq[i] & (1 << (a - j - 1))) != 0) 
			&& (decoded_app_bits[i * a + j][1] <= 0.5)))
		    {
			error_table[i]++;
			total_error_bit_new++;
		    }

	} // Loop on the simulations for one code
      
      // Find maximum error in table
        max = 0;
        max_pos_new = 0;
	  
        for(int i = 0; i < input_length_symb; i++)
	{
	    if(error_table[i] > max)
	    {
	        max = error_table[i];
	        max_pos_new = i;
	    }
	}

        cout << "-------------------------" << endl 
	    << "Code with new column " << max_pos << endl
	    << "Maximum # of errors in column " << max_pos_new << " (" << max 
	    << " errors)" << endl << total_error_bit_new << " total errors" 
	    << endl 
	    << "Bit error rate = " 
	    << setprecision(10) << (double(total_error_bit_new) / double(a * (input_length_symb - memory_order))) / double(nb_sim_max) << endl;
      
      // Go back to previous graph!!!
        if(total_error_bit < total_error_bit_new)
	{
	    cout << "GO BACK TO PREVIOUS CODE" << endl;
	    reset = true;

	  // Change code in the old place
	    for(int i = 0; i < vpc; i++)
		for(int j = 0; j < k; j++)
	    	    for(int l = 0; l < b; l++)
			state_diagram[max_pos * vpc + i].bits[j][l] = old_column[i * vpc + j * k + l]; 
	}
        else
	{
	    max_pos = max_pos_new;
	    total_error_bit = total_error_bit_new;
	    reset = false;
	}
      
      // Print state diagram and code bits
        if(opt % 10 == 0)
	{
	    cout << "********************" << endl << "Best code so far" << endl;
	    for(int i = 0; i < nb_states; i++)
		for(int j = 0; j < k; j++)
	        {
		    cout << i << " -> " << state_diagram[i].transitions[j] << " : ";
		    for(int k = 0; k < b; k++)
			cout << state_diagram[i].bits[j][k] << " ";
		    cout << endl;
	        }
	}

    }  // Optimization loop
    delete temp_drift;
    delete state_diagram;
}

