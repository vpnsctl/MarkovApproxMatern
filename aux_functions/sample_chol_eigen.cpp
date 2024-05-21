// [[Rcpp::depends(RcppEigen)]]
// Import necessary Eigen and Rcpp namespaces
#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd sampleCholSparse(Eigen::SparseMatrix<double> Q, Eigen::MatrixXd stdGauss) {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > chol_Q(Q);
    Eigen::SparseMatrix<double> L = chol_Q.matrixL();
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> Pinv = chol_Q.permutationPinv();
    Eigen::MatrixXd result = L.transpose().template triangularView<Eigen::Upper>().solve(stdGauss);
    result = Pinv * result;
    return result;
}