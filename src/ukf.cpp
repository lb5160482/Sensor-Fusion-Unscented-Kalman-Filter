#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
    is_initialized_ = false;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state and vector and predicted state vector
    x_ = VectorXd(5);
    x_pred = VectorXd(5);

    // initial covariance matrix and predicted state covariance matris
    P_ = MatrixXd::Identity(5, 5);
    P_pred = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.5;

    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

    // state dimension
    n_x_ = 5;

    // augmented state dimension
    n_aug_ = 7;
    
    ///* radar measurement dimension
    n_z_radar_ = 3;
    
    ///* lidar measurement dimension
    n_z_lidar_ = 2;

    // sigma point spreading point
    lambda_ = 3 - n_x_;
    
    // weights vector for computing mean and covariance from predicted sigma points in state space
    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
        weights_(i) = 1 / (2 * (lambda_ + n_aug_));
    }
    
    //
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    // check sensor usage, ignore if specified not to use
    if ((!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) ||
            (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
        return;
    }

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        x_.fill(0.0);
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            float px = meas_package.raw_measurements_(0);
            float py = meas_package.raw_measurements_(1);
            x_(0) = px;
            x_(1) = py;
            cout << "Initialization finished with Laser data!" << endl;
        }
        else {
            float rho = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
            float px = cos(phi) * rho;
            float py = sin(phi) * rho;
            x_(0) = px;
            x_(1) = py;
            cout << "Initialization finished with Radar data!" << endl;
        }

        // Finish initialization
        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;

        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    
    Prediction(dt);
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    }
    else {
        UpdateLidar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
    // get aumented sigma points
    MatrixXd sigAug = GenerateSigmaAug();
    // convert augmented sigma points to state space as prediction
     SigPrediction(sigAug, dt);
    // using predicted sigma points in state space to approximate predicted state mean and state covariance
    MeanAndCovPrediction();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    // matrix for sigma points in measurement space
    MatrixXd ZSig = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);
    // mean predicted measurement
    VectorXd zPred = VectorXd(n_z_lidar_);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_lidar_,n_z_lidar_);
    
    // transform predicted sigma points into measurement space
    ZSig.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd xPred = Xsig_pred_.col(i);
        double px = xPred(0);
        double py = xPred(1);
        
        ZSig(0, i) = px;
        ZSig(1, i) = py;
    }
    
    // calculate mean predicted measurement
    zPred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        zPred += weights_(i) * ZSig.col(i);
    }
    
    // calculate measurement covariance matrix
    MatrixXd R = MatrixXd(2, 2);
    R << std_laspx_ * std_laspx_, 0,
         0, std_laspy_ * std_laspy_;
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diff = ZSig.col(i) - zPred;
        S += weights_(i) * diff * diff.transpose();
    }
    S += R;
    
    /************* Update states *************/
    // create and computecross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z_lidar_);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diffX = Xsig_pred_.col(i) - x_pred;
        while (diffX(3) < -M_PI) {
            diffX(3) += 2 * M_PI;
        }
        while (diffX(3) > M_PI) {
            diffX(3) -= 2 * M_PI;
        }
        VectorXd diffZ = ZSig.col(i) - zPred;
        Tc += weights_(i) * diffX * diffZ.transpose();
    }
    
    // compute Kalman gain
    MatrixXd K = Tc * S.inverse();
    
    // compute diffZ(between measurement and predicted, note different defination from the above one)
    VectorXd diffZ = meas_package.raw_measurements_ - zPred;
    
    this->x_ = this->x_pred + K * diffZ;
    this->P_ = this->P_pred - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    // matrix for sigma points in measurement space
    MatrixXd ZSig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
    // mean predicted measurement
    VectorXd zPred = VectorXd(n_z_radar_);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
    
    // transform predicted sigma points into measurement space
    ZSig.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd xPred = Xsig_pred_.col(i);
        double px = xPred(0);
        double py = xPred(1);
        double v = xPred(2);
        double yaw = xPred(3);
        
        double c1 = sqrt(px * px + py * py);
        ZSig(0, i) = c1;
        ZSig(1, i) = atan2(py, px);
        ZSig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / c1;
    }

    // calculate mean predicted measurement
    zPred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        zPred += weights_(i) * ZSig.col(i);
    }
    
    // calculate measurement covariance matrix
    MatrixXd R = MatrixXd(3, 3);
    R << std_radr_ * std_radr_, 0, 0,
         0, std_radphi_ * std_radphi_, 0,
         0, 0, std_radrd_ * std_radrd_ ;
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diff = ZSig.col(i) - zPred;
        while (diff(1) < -M_PI) {
            diff(1) += 2 * M_PI;
        }
        while (diff(1) > M_PI) {
            diff(1) -= 2 * M_PI;
        }
        S += weights_(i) * diff * diff.transpose();
    }
    S += R;

    /************* Update states *************/
    // create and computecross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diffX = Xsig_pred_.col(i) - x_pred;
        while (diffX(3) < -M_PI) {
            diffX(3) += 2 * M_PI;
        }
        while (diffX(3) > M_PI) {
            diffX(3) -= 2 * M_PI;
        }
        VectorXd diffZ = ZSig.col(i) - zPred;
        while (diffZ(1) < -M_PI) {
            diffZ(1) += 2 * M_PI;
        }
        while (diffZ(1) > M_PI) {
            diffZ(1) -= 2 * M_PI;
        }
        Tc += weights_(i) * diffX * diffZ.transpose();
    }
    
    // compute Kalman gain
    MatrixXd K = Tc * S.inverse();
    
    // compute diffZ(between measurement and predicted, note different defination from the above one)
    VectorXd diffZ = meas_package.raw_measurements_ - zPred;
    while (diffZ(1) < -M_PI) {
        diffZ(1) += 2 * M_PI;
    }
    while (diffZ(1) > M_PI) {
        diffZ(1) -= 2 * M_PI;
    }
    
    this->x_ = this->x_pred + K * diffZ;
    this->P_ = this->P_pred - K * S * K.transpose();
}

/**
* Generate sigma points using current state vector and state covariance
* @return MatrixXd where each column represents a sigma point in the state space
*/
MatrixXd UKF::GenerateSigmaAug() {
    // Create augmented mean vector
    VectorXd xAug = VectorXd(n_aug_);
    xAug.head(n_x_) = this->x_;
    xAug[5] = 0;
    xAug[6] = 0;

    // Create angmented state covariance
    MatrixXd PAug = MatrixXd::Zero(n_aug_, n_aug_);
    PAug.fill(0.0);
    PAug.topLeftCorner(n_x_, n_x_) = this->P_;
    PAug(5, 5) = std_a_ * std_a_;
    PAug(6, 6) = std_yawdd_ * std_yawdd_;

    // Create sigma points matrix
    MatrixXd xSigAUg = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    // Calculate the square root of PAug
    MatrixXd sqrtPAug = PAug.llt().matrixL();

    // Update sigma points matrix
    xSigAUg.col(0) = xAug;
    for (int i = 0; i < n_aug_; ++i) {
        xSigAUg.col(i + 1) = xAug + sqrt(lambda_ + n_aug_) * sqrtPAug.col(i);
        xSigAUg.col(i + 1 + n_aug_) = xAug - sqrt(lambda_ + n_aug_) * sqrtPAug.col(i);
    }

    return xSigAUg;
}

/**
 * Predict sigma points in state space using state transfer funcgtion
 * @param dt delta time from last timestamp to current timestamp
 * @return MatrixXd where each column represents a sigma point prediction in state space
 */
void UKF::SigPrediction(MatrixXd& sigAug, double dt) {
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd xSigAug = sigAug.col(i);
        double px = xSigAug(0);
        double py = xSigAug(1);
        double v = xSigAug(2);
        double yaw = xSigAug(3);
        double yawRate = xSigAug(4);
        double aV = xSigAug(5);
        double aYaw = xSigAug(6);
        // avoid dividing by zero
        if (fabs(yawRate) < 0.001) {
            this->Xsig_pred_(0, i) = px + v * cos(yaw) * dt + 0.5 * pow(dt, 2) * cos(yaw) * aV;
            this->Xsig_pred_(1, i) = py + v * sin(yaw) * dt + 0.5 * pow(dt, 2) * sin(yaw) * aV;
        }
        else {
            this->Xsig_pred_(0, i) = px + v / yawRate * (sin(yaw + yawRate * dt) - sin(yaw)) + 0.5 * pow(dt, 2) * cos(yaw) * aV;
            this->Xsig_pred_(1, i) = py + v / yawRate * (-cos(yaw + yawRate * dt) + cos(yaw)) + 0.5 * pow(dt, 2) * sin(yaw) * aV;
        }
        this->Xsig_pred_(2, i) = v + dt * aV;
        this->Xsig_pred_(3, i) = yaw + yawRate * dt + 0.5 * pow(dt, 2) * aYaw;
        this->Xsig_pred_(4, i) = yawRate + dt * aYaw;
    }

}

/**
 * Using predidcted sigma points to predict state mean and covariance
 * @param xPred state mean to be predicted
 * @param PPred state covariance to be predicted
 */
void UKF::MeanAndCovPrediction() {
    this->x_pred.fill(0.0);
    this->P_pred.fill(0.0);
    // predict x
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        this->x_pred += weights_(i) * Xsig_pred_.col(i);
    }
    // predict P
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diff = Xsig_pred_.col(i) - this->x_pred;
        while (diff(3) < -M_PI) {
            diff(3) += 2 * M_PI;
        }
        while (diff(3) > M_PI) {
            diff(3) -= 2 * M_PI;
        }
        this->P_pred += weights_(i) * diff * diff.transpose();
    }
}
