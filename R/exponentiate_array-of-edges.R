src <- '
  using namespace Rcpp;

  NumericVector vecArray(myArray);
  IntegerVector arrayDims = vecArray.attr("dim");

  arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

  //change one element in the array/cube
  cubeArray(0,0,0) = 518;

  return(wrap(cubeArray));
'
sdcghn
src <- '
  using namespace Rcpp;

  NumericVector vecArray(myArray);
  IntegerVector arrayDims = vecArray.attr("dim");

  arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);



  return(wrap(cubeArray.slice(1)));
'



readCube = cxxfunction(signature(myArray="numeric"),body=src, plugin="RcppArmadillo")


testArray = array(rnorm(18), dim=c(3,3,2))
 print(testArray)

readCube(testArray)


cppFunction("arma::cube getCube(int dim, int slice) {
            arma::cube a(dim,dim,slice);
            a.zeros();
            return a; }",
            depends="RcppArmadillo")

getCube(2, 3)

prob.array = array(0, dim=c(2,2,length(my.tree$edge.length)))
for (i in 1:length(my.tree$edge.length)){
  prob.array[,,i] = two.state.trans.prob(forward.rate, backward.rate, my.tree$edge.length[i])
}

arma::mat my_matexp(arma::mat M) {
  arma::mat Ou=arma::expmat(M);
  return Ou;
}

cppFunction("arma::mat my_matexp(arma::mat M, NumericVector edges) {
  arma::mat Ou=arma::expmat(M*edges(0));
  return Ou;}",

            depends="RcppArmadillo")

my_matexp(M, my.tree$edge.length)

###
cppFunction('int m_dimen(arma::mat M) {
int Ou=M.n_rows;
return Ou;}',
  depends="RcppArmadillo")
m_dimen(M)
#########################
cppFunction("arma::cube my_cube(NumericVector edges, arma::mat M) {
            int matrix_row=M.n_rows;
            int edge_len=edges.size();
            arma::cube a(matrix_row,matrix_row,edge_len);
            a.zeros();


//arma::mat prob=arma::expmat(M*edges(0));
//a.slice(0)=arma::expmat(M*edges(0));

for (int i=0; i<edge_len; i++){
          a.slice(i)=arma::expmat(M*edges(i));
 }
            return a; }",
            depends="RcppArmadillo")

my_cube(my.tree$edge.length, M)
########################################### WORKING VERSION

src <- '
using namespace Rcpp;
using namespace arma;

arma::vec edges = Rcpp::as<arma::vec>(edges_in);
arma::mat M = Rcpp::as<arma::mat>(M_in);

int matrix_row=M.n_rows;
int edge_len=edges.size();
arma::cube aa(matrix_row, matrix_row, edge_len);
aa.zeros();

for (int i=0; i<edge_len; i++){
          aa.slice(i)=arma::expmat(M*edges(i));
 }
return(wrap(aa));'

exp_array_of_edges = cxxfunction(signature(edges_in="numeric", M_in="matrix"),body=src, plugin="RcppArmadillo")

#################
src <- '
using namespace Rcpp;
using namespace arma;

arma::vec edges = Rcpp::as<arma::vec>(edges_in);
arma::mat M = Rcpp::as<arma::mat>(M_in);

Rcpp::Rcout << "M_val"<<M << std::endl;

int matrix_row=M.n_rows;
int edge_len=edges.size();
arma::cube aa(matrix_row, matrix_row, edge_len);
aa.zeros();

for (int i=0; i<edge_len; i++){
          aa.slice(i)=arma::expmat(M*edges(i));
 }
return(wrap(aa));'

exp_array_of_edges = cxxfunction(signature(edges_in="numeric", M_in="matrix"),body=src, plugin="RcppArmadillo")

exp_array_of_edges(my.tree$edge.length, M)




