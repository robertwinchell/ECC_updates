#include <iomanip>
#include <iostream>
#include <fstream> //Needed to manipulate files
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib> 
#include <cstring>
#include <values.h>

#include <gsl/gsl_rng.h> // Required to use the good random generators from gsl

using namespace std; // New definition reference for all standard functions.
 
// ******************************************************************************
// Considers trellis codes with rate R = a / b

#define MAX_a 4   // The maximum number of bits per transition
#define MAX_b 20  // COULD BE IMPROVED WITH DYNAMIC ARRAYS

// ******************************************************************************
// Structure for the nodes in the state diagram

struct Node
{  
  int bits_per_tran;              // The number of bits in the transition
  int transitions[1 << MAX_a];    // Transitions to subsequent nodes
  bool bits[1 << MAX_a][MAX_b];   // Bits of the transitions  
};

// ******************************************************************************
// Structure for the nodes in the decoding trellis

struct DecodingNode
{
  int metric;             // The metric of the current node
  int ancestor;           // The ancestor which optimizes the metric of the
			  // current node (survivor)

  int nb_opt_branches;    // Contains the number of branches
			  // optimizing the metric of the survivor. This is only
			  // used to make sure that one survivor is chosen
			  // uniformly at random among them.

  DecodingNode();         // Constructor
};

// ******************************************************************************
// Matrix classes used to define multi-dimensional arrays
// 
// The array is created on the heap and accessed with rows of pointers 
// corresponding to its dimensions.
//
// The elements of the array initialized at zero

// IMPORTANT NOTE: The user must be VERY careful with the indices, indiced
// outside their boundaried when may well remain inside the main array.

template <typename T>
class Matrix5 {
 public:
  Matrix5 (size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5)
    {
      // data
      T* pdata = new T [dim1 * dim2 * dim3 * dim4 * dim5];
      
      // Initialization
      for(unsigned int i = 0; i < dim1 * dim2 * dim3 * dim4 * dim5; ++i)
	pdata[i] = 0;
      
      // Pointers to the 4th dimension
      ptr_dim4 = new T* [dim1 * dim2 * dim3 * dim4];
      T * temp = pdata;
      
      for (unsigned int i = 0; i < dim1 * dim2 * dim3 * dim4; ++i)
	{
	  ptr_dim4[i] = temp; 
	temp += dim5;      
	}
    
      // Pointers to the 3rd dimension
      ptr_dim3 = new T** [dim1 * dim2 * dim3];
      T ** temp1 = ptr_dim4;
      for (unsigned int i = 0; i < dim1 * dim2 * dim3; ++i)
	{
	  ptr_dim3[i] = temp1; 
	  temp1 += dim4;      
	}
      
      // Pointers to the 2nd dimension
      ptr_dim2 = new T*** [dim1 * dim2];
      T *** temp2 = ptr_dim3;
      for (unsigned int i = 0; i < dim1 * dim2; ++i)
	{
	  ptr_dim2[i] = temp2; 
	  temp2 += dim3;      
	}
      
      // Pointers to the 1st dimension
      ptr_dim1 = new T**** [dim1];
      T **** temp3 = ptr_dim2;
      for (unsigned int i = 0; i < dim1; ++i)
	{
	  ptr_dim1[i] = temp3; 
	  temp3 += dim2;      
	}
    }
  
  ~Matrix5()
    {
      delete[] *ptr_dim4; // free data
      delete[] ptr_dim4;  // free pointers dim 4
      delete[] ptr_dim3;  // free pointers dim 3
      delete[] ptr_dim2;  // free pointers dim 2
      delete[] ptr_dim1;  // free pointers dim 1
    }
  
  // Overloading []
  T**** operator[](size_t dim1)
    {
      return ptr_dim1[dim1];
    }
  
 private:
  T** ptr_dim4;
  T*** ptr_dim3;
  T**** ptr_dim2;
  T***** ptr_dim1;
  
};


template <typename T>
class Matrix3 {
 public:
  Matrix3 (size_t dim3, size_t dim4, size_t dim5)
    {
      // data
      T* pdata = new T [dim3 * dim4 * dim5];
      
      // Initialization
      for(unsigned int i = 0; i < dim3 * dim4 * dim5; ++i)
	pdata[i] = 0;
      
      // Pointers to the 4th dimension
      ptr_dim4 = new T* [dim3 * dim4];
      T * temp = pdata;
      
      for (unsigned int i = 0; i < dim3 * dim4; ++i)
	{
	  ptr_dim4[i] = temp; 
	temp += dim5;      
	}
    
      // Pointers to the 3rd dimension
      ptr_dim3 = new T** [dim3];
      T ** temp1 = ptr_dim4;
      for (unsigned int i = 0; i < dim3; ++i)
	{
	  ptr_dim3[i] = temp1; 
	  temp1 += dim4;      
	}
    };
  
  ~Matrix3()
    {
      delete[] *ptr_dim4; // free data
      delete[] ptr_dim4;  // free pointers dim 4
      delete[] ptr_dim3;  // free pointers dim 3
    }
  
  // Overloading []
  T** operator[](size_t dim3)
    {
      return ptr_dim3[dim3];
    }
  
 private:
  T** ptr_dim4;
  T*** ptr_dim3;
   
};
