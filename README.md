# ECC_updates
additional features to new ECC.  Adds bit manipulation to a trellis based error correction code.


This is code for a trellis based ECC. Here is a breakdown of the code changes.

1.	In the basic_function.cc
Follow function was added.
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

2.	In the rect_code function of the basic_functions.cc
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
3. In the main function of the perf_rect_code.cc
Follow part was added.
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
