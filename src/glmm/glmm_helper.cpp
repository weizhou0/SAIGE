#include <armadillo>

struct SA_NULL {
    arma::mat XV;
    arma::mat XVX;
    arma::mat XXVX_inv;
    arma::mat XVX_inv;
    arma::vec S_a;
    arma::mat XVX_inv_XV;
    arma::vec V;
};

SA_NULL ScoreTest_NULL_Model(const arma::vec& mu, const arma::vec& mu2, const arma::vec& y, const arma::mat& X) {
    SA_NULL result;

    result.V = mu2;
    arma::vec res = y - mu;
    result.XV = X.each_col() % result.V;
    result.XVX = X.t() * result.XV;
    result.XVX_inv = arma::inv(result.XVX);
    result.XXVX_inv = X * result.XVX_inv;
    result.XVX_inv_XV = result.XXVX_inv.each_col() % result.V;
    result.S_a = arma::sum(X.each_col() % res, 0).t();

    return result;
}
