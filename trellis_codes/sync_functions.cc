// ******************************************************************************

#include "definitions.h"
#include "basic_functions.h"

// ******************************************************************************
// This function computes how many bits have to be skipped to match a 
// first string into a second string. For example, mapping 010 in 00011000 
// requires to skip 3 bits. If no mapping is possible, it returns the negative 
// value of the number of bits in the first string which have not been matched. 
// For example, mapping 110 in 0001000 returns -2.
int gaps(const bool str1[], const bool str2[], int length_str1, int length_str2)
{

  int temp = 0;  // Moves in the substring from left to right
  int i = 0;     // Moves in the string from left to right 

  for (i = 0; i < length_str2; i++)
    {
      if(str1[temp] == str2[i])
	{
	  temp++;
	  
	  if(temp == length_str1)
	    return i - length_str1 + 1;
	}
    };

  return -(length_str1 - temp);

};


// ******************************************************************************
// This function deletes bits in a codeword of length initial_length. The 
// position of the bits that have to be deleted are given a table.  

// input: - the codeword, its length, the position of the bits to be deleted, 
//          and the number of deletions
// output: - the final length of the codeword
//         - the updated codeword at the beginning of the array, the remaining
//           bits are junk.
//
// It is assumed that the position of the deletions are stored in the table 
// in increasing order. 

// Although the bits are stored in positions [0 .. length-1] in the array, 
// they are counted starting from position 1.
void delete_bits_table(bool codeword[], int initial_length, int positions[], 
		       int nb_deletions)
{

  for(int i = 0; i < nb_deletions - 1; i++)
    {

      for(int j = 0; j < positions[i+1] - positions[i] - 1; j++)
	{
	  codeword[positions[i]-i+j-1] = codeword[positions[i]+j];
	}
    };

  // End of the codeword the bits after the last deletion must be shifted to 
  // the left 
  for(int j = positions[nb_deletions-1]-nb_deletions; 
      j < initial_length - nb_deletions; j++)
    {
      codeword[j] = codeword[j+nb_deletions];
    }
}


// ******************************************************************************
// This function corrupts a codeword following the binary memoryless noisy 
// synchronization channel.

// The probabilities are defined over each bit of the codeword 

// p_d: probability that the bit is deleted
// p_i: probability that a random bit is inserted before the next bit 
//      of the codeword 
// p_t: probability that the buffered input bit is duplicated
// p_s: probability that the bit is complemented

// The insertions and duplications occur before the input bits. Moreover, it is 
// assumed that no insertion or duplication can occur after the last input bit

// Note: the inserted and duplicated bits CANNOT be further corrupted 
// by substitution errors.

// output: - the length of the corrupted codeword
//         - the corrupted codeword is modified

int binary_NSC(const bool codeword[], bool corrupted_codeword[], 
	       const int codeword_length, const double p_del, const double p_ins, 
	       const double p_dup, const double p_sub, gsl_rng * r)
{
  
  int final_length = 0;
  
  // Verifies if the transition probabilities make sense
  if(p_del + p_ins + p_dup + p_sub > 1)
    {    
      cout << endl << "ERROR: p_d + p_i + p_t + p_s > 1" << endl;
      exit(1);
    }
  
  // Loop on all the bits of the codeword
  for(int input_bit = 0; input_bit < codeword_length; input_bit++) 
    {

      double temp = gsl_rng_uniform(r);
	cout << "error: " << temp << endl;

      // Correct transmission
      if(temp < (1 - p_del - p_ins - p_dup) * (1 - p_sub))
	{
	  corrupted_codeword[final_length] = codeword[input_bit];
	  final_length++;
	}
      // Substitution error
      else if (temp < (1 - p_del - p_ins - p_dup))
	{
	  corrupted_codeword[final_length] = !codeword[input_bit];
	  final_length++;
	}
      // Deletion error
      else if (temp < (1 - p_ins - p_dup))
	{
	  // Nothing to do
	}
      // Insertion error
      else if (temp < (1 - p_dup))
	{
	  corrupted_codeword[final_length] = gsl_rng_uniform_int(r,2);
	  final_length++;
	  input_bit--;
	}
      // Duplication error
      else
	{
	  corrupted_codeword[final_length] = codeword[input_bit];
	  final_length++;
	  input_bit--;
	}      
      
    } // Loop on all the bits of the codeword

  return final_length;

}


// ******************************************************************************
// This function corrupts a codeword following the binary memoryless noisy 
// synchronization channel.

// The probabilities are defined over each bit of the codeword 

// p_d: probability that the bit is deleted
// p_i: probability that a random bit is inserted before the next bit 
//      of the codeword 
// p_t: probability that the buffered input bit is duplicated
// p_s: probability that the bit is complemented

// The insertions and duplications occur before the input bits. Moreover, it is 
// assumed that no insertion or duplication can occur after the last input bit

// Note: the inserted and duplicated bits CANNOT be further corrupted 
// by substitution errors.

// The additional use of this function is to calculate the maximum distance from 
// the expected synchronization drift.

// output: - the length of the corrupted codeword
//         - the corrupted codeword is modified
//         - final_drift is dereferenced to contain the final maximum drift

int binary_NSC_drift(const bool codeword[], bool corrupted_codeword[], 
		     const int codeword_length, const double p_del, 
		     const double p_ins, const double p_dup, const double p_sub, 
		     int * final_drift, gsl_rng * r)
{
  
  int final_length = 0;
  int max_drift = 0;   // Maximum distance from the expected drift
  int temp_drift;

  // Verifies if the transition probabilities make sense
  if(p_del + p_ins + p_dup + p_sub > 1)
    {    
      cout << endl << "ERROR: p_d + p_i + p_t + p_s > 1" << endl;
      exit(1);
    }
  
  // The average number of output bits per input bit.
  const double N = (1 - p_del) / (1 - p_ins - p_dup); 

  // Loop on all the bits of the codeword
  for(int input_bit = 0; input_bit < codeword_length; input_bit++) 
    {

      double temp = gsl_rng_uniform(r);

      // Correct transmission
      if(temp < (1 - p_del - p_ins - p_dup) * (1 - p_sub))
	{
	  corrupted_codeword[final_length] = codeword[input_bit];
	  final_length++;
	}
      // Substitution error
      else if (temp < (1 - p_del - p_ins - p_dup))
	{
	  corrupted_codeword[final_length] = !codeword[input_bit];
	  final_length++;
	}
      // Deletion error
      else if (temp < (1 - p_ins - p_dup))
	{
	  // Nothing to do
	}
      // Insertion error
      else if (temp < (1 - p_dup))
	{
	  corrupted_codeword[final_length] = gsl_rng_uniform_int(r,2);
	  final_length++;
	  input_bit--;
	}
      // Duplication error
      else
	{
	  corrupted_codeword[final_length] = codeword[input_bit];
	  final_length++;
	  input_bit--;
	}
      
      // Calcluate drift and update maximum drift
      temp_drift = ceil(abs(N * (input_bit + 1) - final_length));
      if(temp_drift > max_drift)
	max_drift = temp_drift;

      // cout << "Drifts " << temp_drift << " " << max_drift << endl;
      
    } // Loop on all the bits of the codeword

  // Updates final value of the drift
  *final_drift = max_drift;

  return final_length;

}
