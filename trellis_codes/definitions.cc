#include "definitions.h"

// ******************************************************************************
// Constructor of a node in the decoding trellis

DecodingNode::DecodingNode()
{
  metric = -1;      // They are set to -1, so it is obvious that no path has 
  ancestor = -1;    // yet reached the node.
  nb_opt_branches = 0; 
};

