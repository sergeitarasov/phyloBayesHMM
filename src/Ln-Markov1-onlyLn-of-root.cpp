#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double likelihood_root(
    IntegerMatrix treeEdges, //tE: edge matrix of the ape object phylo
    int numIntNodes, //IN: number of internal nodes (not necessary, but convinient to get it from phylo)
    IntegerVector tipStates, //tS: integer vector of tip states (-1=missing value)
    NumericVector edges, // vector of edge lengths
    arma::mat M, //Rate matrix
    arma::rowvec base// base distribution at root
) {

  int numEdges = treeEdges.nrow();
  //my block
  // arma::vec edges = Rcpp::as<arma::vec>(edges_in);
  //arma::mat M = Rcpp::as<arma::mat>(M_in);

  //int numStates, nS: state space size (e.g. 2 for a binary trait)
  //NumericVector pM, //pM: array of probability matrices for each edge of the tree

  int numStates=M.n_rows;
  int edge_len=edges.size();

  arma::cube cubeProbMat(numStates,numStates, edge_len);
  cubeProbMat.zeros();

  for (int i=0; i<edge_len; i++){
    cubeProbMat.slice(i)=arma::expmat(M*edges(i));
  }
  //my block


  //NumericVector vecProbMat=pM;
  //arma::cube cubeProbMat(vecProbMat.begin(), numStates, numStates, numEdges, false);




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

  int root_id=treeEdges(treeEdges.nrow()-1, 0);
  //arma::mat root_likelihood(1,2);
  //arma::rowvec base=arma::ones<arma::rowvec>(2);
 double root_likelihood=log(accu((partialLike.row(root_id-1)%base)));


  // Rcpp::Rcout << "root vals "<<root_likelihood<< std::endl;

  return root_likelihood;}
