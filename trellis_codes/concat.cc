// ******************************************************************************

#include "definitions.h"
#include "basic_functions.h"
#include "sync_functions.h"
#include "decoders.h"

// ******************************************************************************

// Concatenated codes for the NSC.

// Several inner codes are possible:
// - watermark codes
// - convolutional codes
// - random codes over rectangular graphs
// - predefined codes over rectangular graphs

// We use the LDPC code routines from Radford Neal as outer LDPC codes.

// IMPORTANT NOTES: 

// Uses nonbinary input aplphabets of the form 2^k and binary output alphabets

// We use the following shortcut: to minimize decoding complexity
// and memory requirements, the width of the decoding corridor is set at 2 times 
// the maximum drift of the channel for each specific codeword. 
// The decoder, of course, has no way to know this maximum drift in practice.

// The codes are terminated in State 0

// An all cases, the encoders / decoders use a rectangular graph.

int main(int argc, char * argv[])
{

  // **********************
  // Pseudorandom Generator
  // **********************

  // We use the Mersenne Twister Generator from the GNU Scientific Library
  // WHEN USED IN PARALLEL, IT IS FUNDAMENTAL TO USE ONE SEED PER NODE!
    gsl_rng_default_seed = time(NULL);  // Random Seed
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // Generator

  // ************
  // Main options
  // ************

    if(argc == 1) 
    {
        cout << "Performance of rectangular codes over the NSC" << endl 
	    << "The options are : " << endl
	    << "1 : Watermark codes" << endl
	    << "2 : Convolutional codes" << endl
	    << "3 : Rectangular graph with random code bits" << endl
	    << "4 : Rectangular graph including code bits" << endl;
        exit(1);
    }

  // *********
  // Variables
  // *********
  
    const int option = atoi(argv[1]); // Main option

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

  // Transition probability for a Markov chain of order 1
    long double markov_tran[1];  

  // ***********************************************
  // Input parameters, graph structure and code bits
  // ***********************************************
  
    if(option < 1 || option > 4)
    {
        cout << "Invalid main option" << endl;
        exit(1);
    }

  // Watermark codes
    if(option == 1)
    {
        if(argc != 9)
	{
	    cout << "Performance of watermark codes over the NSC" << endl
	        << "The input arguments are:" << endl
	        << "1 : Option 1" << endl
	        << "2 : Number of outgoing edges per vertex" << endl
	        << "3 : Number of nonbinary input symbols" << endl
	        << "4 : b (encoder of rate a:b (IN BITS))" << endl
	        << "5 - 8 : Parameters of the NSC" << endl << endl;
	    exit(1);
	}

        k = atoi(argv[2]);  
        vpc = 1;            
        input_length_symb = atoi(argv[3]);      
        b = atoi(argv[4]); 
        p_d = atof(argv[5]);
        p_i = atof(argv[6]);
        p_t = atof(argv[7]);
        p_s = atof(argv[8]);

        if(1 << int(round(log(k) / log(2))) != k)
	{
	    cout << "k is not a power of 2. "
	        << "Only binary alphabets are curently supported" << endl;
	    exit(1);
	}
      
      // Encoder of rate a:b (IN BITS)
      // Codes of rate 1:1 are accepted.
        a = int(round(log(k) / log(2)));

        if(a > MAX_a || b > MAX_b) 
	{
	    cout << endl << "ERROR: Adjust MAX_a and/or MAX_b" << endl;
	    exit(1);
	}
      
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
	  cout << "The file containing the watermark bits was not found." 
	       << endl;
	  exit(1);
	}

        cout << "Performance test - watermark codes" << endl
	    << "a:b = " << a << ":" << b << endl
	    << input_length_symb << " input symbols " << endl
	    << "p_d = " << p_d << "; p_i = " << p_i
	    << "; p_t = " << p_t << "; p_s = " << p_s << endl << endl;
    
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

      // Bits of the watermark code are coming from the file
        for(int i = 0; i < k; i++)
	    for(int j = 0; j < b; j++)
		watermark_file >> state_diagram[0].bits[i][j];
      
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
    	    if(argc < 11)
	    {
		cout << "Performance of convolutional codes over the NSC" << endl
	    	    << "The input arguments are:" << endl
	    	    << "1 : Option 2" << endl
	    	    << "2 : Number of outgoing edges per vertex" << endl
	    	    << "3 : Number of vertices per column" << endl
	    	    << "4 : Number of nonbinary input symbols" << endl
	    	    << "5 : b (encoder of rate a:b (IN BITS))" << endl
	    	    << "6 - 9 : Parameters of the NSC" << endl 
	    	    << "10 : Embedded pseudo-random sequence (yes / no)" << endl
	    	    << "11 - more : parameters of the convolutional encoder" << endl
		    << endl;
	    exit(1);
	}
 
        k = atoi(argv[2]); 
        vpc = atoi(argv[3]); 
        input_length_symb = atoi(argv[4]);
        b = atoi(argv[5]); 
        p_d = atof(argv[6]);
        p_i = atof(argv[7]);
        p_t = atof(argv[8]);
        p_s = atof(argv[9]);

        emb = atoi(argv[10]);  // Embedded pseudorandom sequence (if 1)

        if(1 << int(round(log(k) / log(2))) != k)
	{
	    cout << "k is not a power of 2. "
	    	<< "Only binary alphabets are curently supported" << endl;
	    exit(1);
	}

        a = int(round(log(k) / log(2)));
  
        if(a >= b)
	{
	    cout << "a >= b!!!" << endl;
	    exit(1);
	}
  
        if(a > MAX_a || b > MAX_b) 
	{
	    cout << endl << "ERROR: Adjust MAX_a and/or MAX_b" << endl;
	    exit(1);
	}
  
        nb_states = vpc * input_length_symb;
        memory_order = round(log(vpc) / log(k));
    	state_diagram = new Node[nb_states];

        if(argc != (11 + (memory_order + 1) * a * b))
	{
	    cout << "Incorrect number of input arguments for the " << "convolutional encoder." << endl;
	    exit(1);
	}

        bool conv_array[(memory_order + 1) * a * b];

        for(int i = 0; i < (memory_order + 1) * a * b; i++)
	    conv_array[i] = atoi(argv[11 + i]);

      // Generate the first column (code bits and transitions) of the trellis from the convolutional encoder
        conv_to_trellis(a, b, memory_order, conv_array, state_diagram, vpc);
    
      // Change the transitions of the first column so that they point 
      // to the second column
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

        cout << "Performance test - convolutional codes over the NSC" << endl
	    << "k = " << k << endl << vpc << " vertices per column" << endl
	    << input_length_symb << " input symbols " << endl
	    << b << " output bits per input symbol" << endl
	    << "p_d = " << p_d << "; p_i = " << p_i
	    << "; p_t = " << p_t << "; p_s = " << p_s << endl;
        if(emb) cout << "Embedded pseudorandom sequences" << endl;
        else cout << "NO embedded pseudorandom sequences" << endl;
        cout << endl;

        cout << "Convolutional encoder : " << endl;
        for(int i = 0; i < a * b; i++)
	{
	    for(int j = 0; j < memory_order + 1; j++)
		cout << conv_array[i * (memory_order + 1) + j] << " ";
	    cout << endl;
	}
    }

  // Rectangular graph with random code bits
    else if (option == 3)
    {
        if(argc != 11)
	{
	    cout << "Performance of rectangular graphs over the NSC"
	        << endl << "Random code bits" << endl
	        << "The input arguments are:" << endl
	        << "1 : Option 3" << endl
	        << "2 : Number of outgoing edges per vertex" << endl
	        << "3 : Number of vertices per column" << endl
	        << "4 : Number of nonbinary input symbols" << endl
	        << "5 : b (encoder of rate a:b (IN BITS))" << endl
	        << "6 - 9 : Parameters of the NSC" << endl 
	        << "10 - The transition of the Markov chain modeling the code" 
	        << endl << endl;
	    exit(1);
	}

        k = atoi(argv[2]); 
        vpc = atoi(argv[3]);
        input_length_symb = atoi(argv[4]);      
        b = atoi(argv[5]);
        p_d = atof(argv[6]);
        p_i = atof(argv[7]);
        p_t = atof(argv[8]);
        p_s = atof(argv[9]);
        markov_tran[0] = atof(argv[10]);

	if(1 << int(round(log(k)/log(2))) != k)
	{
	    cout << "k is not a power of 2. "
	        << "Only binary alphabets are curently supported" << endl;
	    exit(1);
	}

      // Encoder of rate a:b (IN BITS)
      // Codes of rate 1:1 are accepted
        a = int(round(log(k) / log(2)));
        if(a > MAX_a || b > MAX_b) 
	{
	    cout << endl << "ERROR: Adjust MAX_a and/or MAX_b" << endl;
	    cout << MAX_b;
	    exit(1);
	}

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

        cout << "Performance test - rectangular graph with random code bits" << endl
	    << "k = " << k << endl << vpc << " vertices per columns" << endl
	    << input_length_symb << " input symbols " << endl
	    << b << " output bits per input symbol" << endl
	    << "p_d = " << p_d << "; p_i = " << p_i
	    << "; p_t = " << p_t << "; p_s = " << p_s << endl 
	    << "Markov chain with transition probability "<< markov_tran[0]
	    << endl << endl;

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
      //All the additional machinery is used so that the outgoing edges 
      // of each vertex have different output bits.

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

  // Rectangular graph and code bits from an output file
    else if (option == 4)
    {
        if(argc != 7)
	{
	    cout << "Performance of rectangular graphs over the NSC" << endl
		<< "The input arguments are:" << endl
		<< "1 : Option 4" << endl
	        << "2 : The output file containing the rectangular graph and the "
	        << "code bits " << endl
	        << "3 - 6 : Parameters of the NSC" << endl << endl;
	  exit(1);
	}

        p_d = atof(argv[3]);
        p_i = atof(argv[4]);
        p_t = atof(argv[5]);
        p_s = atof(argv[6]);

        char * input_file_name = argv[2];
     
      // File containing the rectangular graph
        ifstream graph_file(input_file_name, ios::in);
        if(!graph_file)
	{
	    cout << "The file containing the encoder and code bits was not found." << endl;
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
        cout << "Performance test - rectangular graph with code bits" << endl 
	    << "Encoder file : " << argv[2] << endl
	    << "k = " << k << endl << vpc << " vertices per columns" << endl
	    << input_length_symb << " input symbols " << endl
	    << b << " output bits per input symbol" << endl
	    << "p_d = " << p_d << "; p_i = " << p_i
	    << "; p_t = " << p_t << "; p_s = " << p_s << endl << endl;
    }

  // *************************
  // Variables for simulations
  // *************************

  // Maximum distance from the expected synchronization drift
    int drift;
    int * temp_drift = new int;

  // Outer LDPC code
    int input_length_outer = 100;   // The number of 


  //int input_seq_outer[]
    stringstream gen_outter_ldpc;
    gen_outter_ldpc << "./gen_outter_ldpc " << vpc << " " << b << " " << gsl_rng_default_seed;
    cout << gen_outter_ldpc.str().c_str() << endl;
    system((const char*)gen_outter_ldpc.str().c_str());

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
    vector < vector <double> > prior_prob(input_length_symb, vector<double> (k, 1 / double(k)));

  // App for the received symbols
    vector < vector <double> > decoded_app(input_length_symb, vector<double> (k, 0)); 

  // App for the received bits
    vector < vector <double> > decoded_app_bits(input_length_symb * a, vector<double> (2, 0));

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
        if(option == 1 || (option == 2 && emb))
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
		    total_error_bit++;
      
        if(nb_sim % 10 == 0)
	    cout << nb_sim << " codewords; " << "Bit error rate = " 
	        << setprecision(10) << (double(total_error_bit) / double(a * (input_length_symb - memory_order))) / double(nb_sim) << endl;
    }
    delete temp_drift;
    delete state_diagram;
}
