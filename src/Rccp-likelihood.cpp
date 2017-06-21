src <- '
/* Arguments: 
tE: edge matrix of the ape object phylo
nIN: number of internal nodes (not necessary, but convinient to get it from phylo)
tS: integer vector of tip states (-1=missing value)
nS: state space size (e.g. 2 for a binary trait)
pM: array of probability matrices for each edge of the tree

Two important assumptions:
1. edges in the edge matrix and probability matrices are in the "pruningwise" order; 
see ?reorder.phylo for more details
2. tip state vector is ordered according to the tip numbering in the edge matrix

 return(partLike(my.tree$edge, my.tree$Nnode, my.data, 2, prob.array))
 */

using namespace Rcpp;


IntegerMatrix treeEdges(tE);      /* tE: edge matrix of the ape object phylo */

// get the number of edges
int numEdges = treeEdges.nrow(); /* tE: edge matrix of the ape object phylo */
int numIntNodes = as<int>(nIN); /* nIN: number of internal nodes  */
IntegerVector tipStates(tS);    /* tS: integer vector of tip states (-1=missing value) */
int numStates = as<int>(nS);    /* nS: state space size (e.g. 2 for a binary trait) */
NumericVector vecProbMat(pM);   /* pM: array of probability matrices for each edge of the tree */
arma::cube cubeProbMat(vecProbMat.begin(), numStates, numStates, numEdges, false);

//Rcpp::Rcout << "tE "<<treeEdges<< std::endl;
//Rcpp::Rcout << "numEdges "<<numEdges << std::endl;
//Rcpp::Rcout << "numIntNodes "<<numIntNodes << std::endl;
//Rcpp::Rcout << "tipStates "<< tipStates<< std::endl;
//Rcpp::Rcout << "numStates "<< numStates<< std::endl;
//Rcpp::Rcout << "vecProbMat "<< vecProbMat<< std::endl;
//Rcpp::Rcout << "cubeProbMat "<< cubeProbMat<< std::endl;


// get the number of tips in the tree
int numTips = tipStates.size();

// prepare a matrix for storing regular (backward) partial likelihoods
arma::mat partialLike = arma::zeros<arma::mat>(numTips + numIntNodes, numStates);

//Rcpp::Rcout << "partialLike1 "<< partialLike<< std::endl;


for (int i=0; i < numTips; i++){
  if (tipStates[i] == -1){// -1 denotes a missing value
    partialLike.row(i) = arma::ones<arma::rowvec>(numStates);
  }else{
    partialLike(i, tipStates[i]) = 1.0;
  }
}

//Rcpp::Rcout << "partialLike2 "<< partialLike<< std::endl;

// compute regular partial likelihoods for all internal nodes
for (int i=0; i < numEdges; i+=2){      
  // parent=treeEdges(i,0) or treeEdges(i+1,0); treeEdges indices should be shifted by one
  partialLike.row(treeEdges(i,0)-1) = (partialLike.row(treeEdges(i,1)-1)*cubeProbMat.slice(i).t())%(partialLike.row(treeEdges(i+1,1)-1)*cubeProbMat.slice(i+1).t());            
}

return wrap(partialLike);
'  


partLike = cxxfunction(signature(tE="integer", nIN="integer", tS="integer", nS="integer", pM="numeric"),body=src, plugin="RcppArmadillo")
  