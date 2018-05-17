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
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;
    
    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;
    
    // initial state vector
    x_ = VectorXd(5);
    
    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.3;
    
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
    
    /**
     TODO:
     
     Complete the initialization. See ukf.h for other member properties.
     
     Hint: one or more values initialized above might be wildly off...
     */
    is_initialized_ = false;
    
    n_x_ = 5;
    
    n_aug_ = 7;
    
    Xsig_pred_ = MatrixXd(5, 2 * n_aug_ + 1);
    
    weights_ = VectorXd(2 * n_aug_ + 1);
    
    lambda_ = 3 - n_x_;
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
        
        x_ << 1, 1, 1, 1, 0.1;
        P_ << 0.15,    0, 0, 0, 0,
        0, 0.15, 0, 0, 0,
        0,    0, 1, 0, 0,
        0,    0, 0, 1, 0,
        0,    0, 0, 0, 1;
        
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
        // set weights
        weights_(0) = lambda_ / (lambda_ + n_aug_);
        for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
            weights_(i) = 0.5 / (lambda_ + n_aug_);
        }
        
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    
    Prediction(delta_t);
    
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
void UKF::Prediction(double delta_t) {
    /*****************************************************************************
     *  Generate augmented sigma points
     ****************************************************************************/
    // Create augmented mean vector
    VectorXd xAug = VectorXd(n_aug_);
    xAug.head(n_x_) = x_;
    xAug[5] = 0;
    xAug[6] = 0;
    
    // Create angmented state covariance
    MatrixXd PAug = MatrixXd(n_aug_, n_aug_);
    PAug.fill(0.0);
    PAug.topLeftCorner(5, 5) = P_;
    PAug(5, 5) = std_a_ * std_a_;
    PAug(6, 6) = std_yawdd_ * std_yawdd_;
    
    // Calculate the square root of PAug
    MatrixXd L = PAug.llt().matrixL();
    
    // Create augmented sigma points matrix
    MatrixXd XSigAug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    XSigAug.col(0) = xAug;
    for (int i = 0; i < n_aug_; ++i) {
        XSigAug.col(i + 1) = xAug + sqrt(lambda_ + n_aug_) * L.col(i);
        XSigAug.col(i + 1 + n_aug_) = xAug - sqrt(lambda_ + n_aug_) * L.col(i);
    }
    
    /*****************************************************************************
     *  Convert augmented sigma points to state space as prediction(update Xsig_pred_)
     ****************************************************************************/
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd xSigAug = XSigAug.col(i);
        double px = xSigAug(0);
        double py = xSigAug(1);
        double v = xSigAug(2);
        double yaw = xSigAug(3);
        double yawRate = xSigAug(4);
        double aV = xSigAug(5);
        double aYaw = xSigAug(6);
        // avoid dividing by zero
        if (fabs(yawRate) < 0.001) {
            Xsig_pred_(0, i) = px + v * cos(yaw) * delta_t + 0.5 * pow(delta_t, 2) * cos(yaw) * aV;
            Xsig_pred_(1, i) = py + v * sin(yaw) * delta_t + 0.5 * pow(delta_t, 2) * sin(yaw) * aV;
        }
        else {
            Xsig_pred_(0, i) = px + v / yawRate * (sin(yaw + yawRate * delta_t) - sin(yaw)) + 0.5 * pow(delta_t, 2) * cos(yaw) * aV;
            Xsig_pred_(1, i) = py + v / yawRate * (-cos(yaw + yawRate * delta_t) + cos(yaw)) + 0.5 * pow(delta_t, 2) * sin(yaw) * aV;
        }
        Xsig_pred_(2, i) = v + delta_t * aV;
        Xsig_pred_(3, i) = yaw + yawRate * delta_t + 0.5 * pow(delta_t, 2) * aYaw;
        Xsig_pred_(4, i) = yawRate + delta_t * aYaw;
    }
    
    /*****************************************************************************
     *  Using predidcted sigma points to approximate predicted state mean and covariance
     ****************************************************************************/
    // predict x
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        x_ += weights_(i) * Xsig_pred_.col(i);
    }
    // predict P
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diff = Xsig_pred_.col(i) - x_;
        // normalization
        while (diff(3) < -M_PI) {
            diff(3) += 2 * M_PI;
        }
        while (diff(3) > M_PI) {
            diff(3) -= 2 * M_PI;
        }
        P_ += weights_(i) * diff * diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    int n_z = 2;
    
    // compute matrix for sigma points in measurement space
    MatrixXd ZSig = MatrixXd(n_z, 2 * n_aug_ + 1);
    // transform predicted sigma points into measurement space
    ZSig.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd xsigPred = Xsig_pred_.col(i);
        double px = xsigPred(0);
        double py = xsigPred(1);
        
        ZSig(0, i) = px;
        ZSig(1, i) = py;
    }
    
    // mean predicted measurement
    VectorXd zPred = VectorXd(n_z);
    zPred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        zPred += weights_(i) * ZSig.col(i);
    }
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diff = ZSig.col(i) - zPred;
        // normalization
        while (diff(1) < -M_PI) {
            diff(1) += 2 * M_PI;
        }
        while (diff(1) > M_PI) {
            diff(1) -= 2 * M_PI;
        }
        S += weights_(i) * diff * diff.transpose();
    }
    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_laspx_ * std_laspx_,                         0,
                               0,   std_laspy_ * std_laspy_;
    S += R;
    
    /************* Update states *************/
    // create and computecross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // state difference
        VectorXd diffX = Xsig_pred_.col(i) - x_;
        while (diffX(3) < -M_PI) {
            diffX(3) += 2 * M_PI;
        }
        while (diffX(3) > M_PI) {
            diffX(3) -= 2 * M_PI;
        }
        
        // measurement difference
        VectorXd diffZ = ZSig.col(i) - zPred;
        
        Tc += weights_(i) * diffX * diffZ.transpose();
    }
    
    // compute Kalman gain
    MatrixXd K = Tc * S.inverse();
    
    // compute diffZ(between measurement and predicted, note different defination from the above one)
    VectorXd diffZ = meas_package.raw_measurements_ - zPred;
    
    // udpate state and state covariance
    x_ = x_ + K * diffZ;
    P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    int n_z = 3;
    
    // compute matrix for sigma points in measurement space
    MatrixXd ZSig = MatrixXd(n_z, 2 * n_aug_ + 1);
    // transform predicted sigma points into measurement space
    ZSig.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd xsigPred = Xsig_pred_.col(i);
        double px = xsigPred(0);
        double py = xsigPred(1);
        double v = xsigPred(2);
        double yaw = xsigPred(3);
        
        double c1 = sqrt(px * px + py * py);
        ZSig(0, i) = c1;
        ZSig(1, i) = atan2(py, px);
        ZSig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / c1;
    }
    
    // mean predicted measurement
    VectorXd zPred = VectorXd(n_z);
    zPred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        zPred += weights_(i) * ZSig.col(i);
    }
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        VectorXd diff = ZSig.col(i) - zPred;
        // normalization
        while (diff(1) < -M_PI) {
            diff(1) += 2 * M_PI;
        }
        while (diff(1) > M_PI) {
            diff(1) -= 2 * M_PI;
        }
        S += weights_(i) * diff * diff.transpose();
    }
    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(3, 3);
    R << std_radr_ * std_radr_,                         0,                       0,
                             0, std_radphi_ * std_radphi_,                       0,
                             0,                         0, std_radrd_ * std_radrd_;
    S += R;
    
    /************* Update states *************/
    // create and computecross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // state difference
        VectorXd diffX = Xsig_pred_.col(i) - x_;
        while (diffX(3) < -M_PI) {
            diffX(3) += 2 * M_PI;
        }
        while (diffX(3) > M_PI) {
            diffX(3) -= 2 * M_PI;
        }
        
        // measurement difference
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
    
    // udpate state and state covariance
    x_ = x_ + K * diffZ;
    P_ = P_ - K * S * K.transpose();
}
