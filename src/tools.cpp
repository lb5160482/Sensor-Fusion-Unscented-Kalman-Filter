#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse.fill(0.0);
    for (int i = 0; i < estimations.size(); ++i) {
        VectorXd diff = estimations[i] - ground_truth[i];
        // element wise product
        VectorXd diffELementProd = diff.array() * diff.array();
        rmse += diffELementProd;
    }
    rmse = rmse / estimations.size();
    // element wise sqrt
    rmse = rmse.array().sqrt();

    return rmse;
}
