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
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

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

    //set state dimension
    n_x_ = 5;
    n_aug_ = n_x_ + 2;
    lambda_ = 3 - n_x_;

    //set weights
    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_.fill(0.5 / (lambda_ + n_aug_));
    weights_(0) = lambda_ / (lambda_ + n_aug_);

    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        x_.fill(0);

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            // Convert radar from polar to cartesian coordinates and initialize state.
            float ro = meas_package.raw_measurements_(0);
            float theta = meas_package.raw_measurements_(1);
            float ro_dot = meas_package.raw_measurements_(2);
            x_(0) = ro * sin(theta);
            x_(1) = ro * cos(theta);
            x_(2) = ro_dot * sin(theta);
            x_(3) = ro_dot * cos(theta);
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            // LASER only measures location leave velocity default.
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }

        time_us_ = meas_package.timestamp_;

        //state covariance matrix P
        P_.fill(0);
        P_.diagonal().fill(1);

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    // - Time is measured in seconds.
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    Prediction(dt);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    // Generate sigma points

    //create augmented mean vector
    VectorXd x_aug = VectorXd(7);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }

    //predict sigma points
    //avoid division by zero
    //write predicted sigma points into right column
    //create matrix with predicted sigma points as columns

    for (int i = 0; i<Xsig_aug.cols(); i++)
    {
        float px = Xsig_aug(0, i);
        float py = Xsig_aug(1, i);
        float v = Xsig_aug(2, i);
        float psi = Xsig_aug(3, i);
        float psi_dot = Xsig_aug(4, i);
        float mu_a = Xsig_aug(5, i);
        float mu_psi_dot_dot = Xsig_aug(6, i);

        //   std::cout << abs(psi_dot) << std::endl;
        if (fabsf(psi_dot) < 0.000001) {
            float a = delta_t * cos(psi);
            float b = delta_t * sin(psi);

            Xsig_pred_.col(i) << px + v * a + 0.5 * delta_t * a * mu_a,
                                py + v * b + 0.5 * delta_t * b * mu_a,
                                v + delta_t * mu_a,
                                psi + 0.5 * delta_t * delta_t * mu_psi_dot_dot,
                                psi_dot + delta_t * mu_psi_dot_dot;
        } else {
            Xsig_pred_.col(i) << px + v * (sin(psi + delta_t * psi_dot) - sin(psi)) / psi_dot + 0.5 * delta_t * delta_t * cos(psi) * mu_a,
                                py + v * (-cos(psi + delta_t * psi_dot) + cos(psi)) / psi_dot + 0.5 * delta_t * delta_t * sin(psi) * mu_a,
                                v + delta_t * mu_a,
                                psi + psi_dot * delta_t + 0.5 * delta_t * delta_t * mu_psi_dot_dot,
                                psi_dot + delta_t * mu_psi_dot_dot;
        }
    }

    //predict state mean
    x_ = Xsig_pred_ * weights_;

    //predict state covariance matrix
    P_.fill(0);
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        VectorXd a = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (a(3) > M_PI) a(3) -= 2. * M_PI;
        while (a(3) < -M_PI) a(3) += 2. * M_PI;
        P_ += weights_(i) * a * a.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    const int n_z = 3; //  radar measurement dimentions

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    //transform sigma points into measurement space
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        float px  = Xsig_pred_(0, i);
        float py  = Xsig_pred_(1, i);
        float v   = Xsig_pred_(2, i);
        float psi = Xsig_pred_(3, i);

        Zsig(0, i) = sqrt(px * px + py * py); // rho
        Zsig(1, i) = atan2(py, px); // phi
        Zsig(2, i) = v * (px * cos(psi) + py * sin(psi)) / Zsig(0, i); // rho_dot
    }

    //calculate mean predicted measurement
    z_pred = Zsig * weights_;

    //calculate innovation covariance matrix S
    S.fill(0);
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        MatrixXd a = Zsig.col(i) - z_pred;
        //angle normalization
        while (a(1) > M_PI) a(1) -= 2. * M_PI;
        while (a(1) < -M_PI) a(1) += 2. * M_PI;
        S += weights_(i) * a * a.transpose();
    }

    S(0, 0) += std_radr_ * std_radr_;
    S(1, 1) += std_radphi_ * std_radphi_;
    S(2, 2) += std_radrd_ * std_radrd_;

    //calculate cross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    Tc.fill(0);
    for (int i = 0; i < Zsig.cols(); i++) {
        Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
    }

    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //update state mean and covariance matrix
    x_ = x_ + K * (meas_package.raw_measurements_ - z_pred);
    P_ = P_ - K * S * K.transpose();
}
