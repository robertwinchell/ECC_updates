// ******************************************************************************

#include "definitions.h"
#include "basic_functions.h"

// ******************************************************************************

// Generates a code over a rectangular graph.
// Several options are possible

int main(int argc, char * argv[])
{

  // **********************
  // Pseudorandom Generator
  // **********************

  // We use the Mersenne Twister Generator from the GNU Scientific Library
  // WHEN USED IN PARALLEL, IT IS FUNDAMENTAL TO USE ONE SEED PER NODE!
	gsl_rng_default_seed = time(NULL);  			// Random Seed
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); 		// Generator

  // ************
  // Main options
  // ************

	if(argc == 1) 
	{
    	    cout << "Performance of rectangular codes over the NSC" << endl 
		<< "The options are : " << endl
		<< "1 : Watermark codes" << endl
		<< "2 : Convolutional codes" << endl
		<< "3 : Inner Trellis codes" << endl
		<< "4 : Outter Ldpc codes" << endl
		<< "5 : Rectangular graph with random code bits" << endl;
	    exit(1);
	}

  // *********
  // Variables
  // *********
  
	const int option = atoi(argv[1]); 		// Main option
	char * output_file_name = NULL; 		// Output file
  
  // Number of states in the state diagram
	int nb_states;
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
    
  // Memory order of the code
	int memory_order;

  // Contains the convolutions
	bool * conv_array;

  // Transition probability for a Markov chain of order 1
	long double markov_tran[1];  

  // ***********************************************
  // Input parameters, graph structure and code bits
  // ***********************************************
  
	if(option < 1 || option > 5)
	{
    	    cout << "Invalid main option" << endl;
    	    exit(1);
	}

  // Watermark codes
	if(option == 1)
	{
    	    if(argc != 6)
	    {
		cout << "Generate a watermark code over a rectangular graph" << endl
	    	    << "The input arguments are:" << endl
	    	    << "1 : Option 1" << endl
	    	    << "2 : Number of outgoing edges per vertex" << endl
	    	    << "3 : Number of nonbinary input symbols" << endl
	    	    << "4 : b (encoder of rate a:b (IN BITS))" << endl 
	    	    << "5 : Output file name" << endl << endl;
		exit(1);
	    }
      
    	    k = atoi(argv[2]);  
	    vpc = 1;
//    	    vpc = pow(k, 2);            
    	    input_length_symb = atoi(argv[3]);      
    	    b = atoi(argv[4]); 
    	    output_file_name = argv[5];

	}
    // Convolutional codes
	else if(option == 2)
	{

    	    if(argc < 8)
  	    {
  		cout << "Performance of convolutional codes over the NSC" << endl
  	    	    << "The input arguments are:" << endl
  	    	    << "1 : Option 2" << endl
  	    	    << "2 : Number of outgoing edges per vertex" << endl
  	    	    << "3 : Number of vertices per column" << endl
  	    	    << "4 : Number of nonbinary input symbols" << endl
  	    	    << "5 : b (encoder of rate a:b (IN BITS))" << endl
	    	    << "6 : Output file name" << endl
	    	    << "7 - more : parameters of the convolutional encoder" << endl
	    	    << endl;
  		exit(1);
  	    }
 
    	    k = atoi(argv[2]); 
    	    vpc = atoi(argv[3]); 
    	    input_length_symb = atoi(argv[4]);
    	    b = atoi(argv[5]);
    	    output_file_name = argv[6];
      
    	    a = int(round(log(k) / log(2)));
 
    	    memory_order = round(log(vpc) / log(k)); 

    	    if(argc != (7 + (memory_order + 1) * a * b))
   	    {
		cout << "argc:" << (memory_order + 1) * a * b << endl;
   		cout << "Incorrect number of input arguments for the "
   	    	    << "convolutional encoder." << endl;
   		exit(1);
   	    }

    	    conv_array = new bool[(memory_order + 1) * a * b];

    	    for(int i = 0; i < (memory_order + 1) * a * b; i++)
		conv_array[i] = atoi(argv[7 + i]); 

    	    for(int i = 0; i < (memory_order + 1) * a * b; i++)
		cout << conv_array[i];
	}
    // Rectangle graph with inner trellis codes
	else if(option == 3){
	    if(argc != 7){
		cout << "Performance of inner trellis codes over the NSC" << endl
  	    	    << "The input arguments are:" << endl
  	    	    << "1 : Option 3" << endl
  	    	    << "2 : Number of outgoing edges per vertex" << endl
  	    	    << "3 : Number of vertices per column" << endl
  	    	    << "4 : Number of nonbinary input symbols" << endl
  	    	    << "5 : b (encoder of rate a:b (IN BITS))" << endl
	    	    << "6 : Output file name" << endl
	    	    << endl;
  		exit(1);
	    }

	    // Input the initial parmeters.	    
    	    k = atoi(argv[2]);  
	    vpc = atoi(argv[3]);    
	    input_length_symb = atoi(argv[4]);      
    	    b = atoi(argv[5]); 
    	    output_file_name = argv[6];
	}
    // Rectangle graph with outter ldpc codes
	else if(option  == 4){
	    if(argc != 7){
		cout << "Performance of outter ldpc codes over the NSC" << endl
  	    	    << "The input arguments are:" << endl
  	    	    << "1 : Option 4" << endl
  	    	    << "2 : Number of outgoing edges per vertex" << endl
  	    	    << "3 : Number of vertices per column" << endl
  	    	    << "4 : Number of nonbinary input symbols" << endl
  	    	    << "5 : b (encoder of rate a:b (IN BITS))" << endl
	    	    << "6 : Output file name" << endl
	    	    << endl;
  		exit(1);
	    }

	    //Input the initiaj parameters.
	    k = atoi(argv[2]);  
	    vpc = atoi(argv[3]);
    	    input_length_symb = atoi(argv[4]);      
    	    b = atoi(argv[5]); 
    	    output_file_name = argv[6];
	}
  // Rectangular graph with random code bits
	else if (option == 5)
	{
    	    if(argc != 8)
  	    {
  		cout << "Performance of rectangular graphs over the NSC"
  	    	    << endl << "Random code bits" << endl
  	    	    << "The input arguments are:" << endl
  	    	    << "1 : Option 3" << endl
  	    	    << "2 : Number of outgoing edges per vertex" << endl
  	    	    << "3 : Number of vertices per column" << endl
  	    	    << "4 : Number of nonbinary input symbols" << endl
  	    	    << "5 : b (encoder of rate a:b (IN BITS))" << endl
  	    	    << "6 : The transition of the Markov chain modeling the code" 
	    	    << endl
	    	    << "7 : Output file name" << endl << endl;
	 
  		exit(1);
  	    }
 
    	    k = atoi(argv[2]); 
    	    vpc = atoi(argv[3]);
    	    input_length_symb = atoi(argv[4]);
    	    b = atoi(argv[5]);
    	    markov_tran[0] = atof(argv[6]);
    	    output_file_name = argv[7];
	}

  // Generate code bits over a rectangular graph
	state_diagram = rect_code(option, k, vpc, input_length_symb, b, conv_array, markov_tran, r);

  // Output file for the graph
	ofstream output_file(output_file_name, ios::out);

	output_file << k << " " << vpc << " " << input_length_symb << " " << b << endl << endl;

	nb_states = input_length_symb * vpc;

	for(int i = 0; i < nb_states; i++)
	{
	    for(int j = 0; j < k; j++)
    	    {
  		output_file << i << " " << state_diagram[i].transitions[j] << " ";
  		for(int kk = 0; kk < b; kk++)
  		    output_file << state_diagram[i].bits[j][kk] << " ";
  		output_file << endl;
    	    }
	}

	output_file.close();  

	if(option == 2) delete [] conv_array;
	delete [] state_diagram;
}
