#include "definitions.h"
#include "basic_functions.h"

// *******************************************************************************

// This program generates rectangular graphs and returns them in an output file.
// Both vertices of each edge e=(va,vb) of the graph are included in the file:

// va1 vb1
// va2 vb2
// va3 vb3
//   ...
 
// Input arguments
// k : alphabet size
// vpc : number of vertices per column
// s: number of columns in the graph (minimum shortcut distance)

// Output file names have the form r_k_vpc_s.txt

// *******************************************************************************

int main(int argc, char * argv[])
{
    if(argc != 4)
    {
        cout << "The input arguments are " << endl
	    << "1 - Number of outgoing edges per vertex " << endl
	    << "2 - Number of vertices per column " << endl
	    << "3 - Minimum shortcut distance (number of columns in the graph)"
	    << endl;

        exit(1);
    }

    int k   = atoi(argv[1]);  // Number of outgoing edges per vertex
    int vpc = atoi(argv[2]);  // Number of vertices per column
    int s = atoi(argv[3]);    // Minimum shortcut distance (number of columns in the graph)
    gen_rect(k,vpc,s);
}
