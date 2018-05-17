#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;

    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;

    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;

    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;
    
    ///* predicted state vector
    VectorXd x_pred;

    ///* state covariance matrix
    MatrixXd P_;
    
    ///* predicted state covariance matrix
    MatrixXd P_pred;

    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;

    ///* time when the state is true, in us
    long long time_us_;

    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;

    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;

    ///* Laser measurement noise standard deviation position1 in m
    double std_laspx_;

    ///* Laser measurement noise standard deviation position2 in m
    double std_laspy_;

    ///* Radar measurement noise standard deviation radius in m
    double std_radr_;

    ///* Radar measurement noise standard deviation angle in rad
    double std_radphi_;

    ///* Radar measurement noise standard deviation radius change in m/s
    double std_radrd_ ;

    ///* Weights of sigma points
    VectorXd weights_;

    ///* State dimension
    int n_x_;

    ///* Augmented state dimension
    int n_aug_;
    
    ///* radar measurement dimension
    int n_z_radar_;
    
    ///* lidar measurement dimension
    int n_z_lidar_;

    ///* Sigma point spreading parameter
    double lambda_;


    /**
   * Constructor
   */
    UKF();

    /**
   * Destructor
   */
    virtual ~UKF();

    /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
    void ProcessMeasurement(MeasurementPackage meas_package);

    /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
    void Prediction(double dt);

    /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
    void UpdateLidar(MeasurementPackage meas_package);

    /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
    void UpdateRadar(MeasurementPackage meas_package);

    /**
    * Generate sigma points using current state vector and state covariance
    * @return MatrixXd where each column represents a sigma point in the state space
    */
    MatrixXd GenerateSigmaAug();
    
    /**
     * Predict sigma points in state space using state transfer funcgtion
     * @param dt delta time from last timestamp to current timestamp
     * @return MatrixXd where each column represents a sigma point prediction in state space
     */
    void SigPrediction(MatrixXd& sigAug, double dt);
    
    /**
     * Using predidcted sigma points to predict state mean and covariance
     * @param xPred state mean to be predicted
     * @param PPred state covariance to be predicted
     */
    void MeanAndCovPrediction();
};

#endif /* UKF_H */
