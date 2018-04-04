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
    x_ = VectorXd::Zero(5);

    // initial covariance matrix
    P_ = MatrixXd::Identity(5, 5);

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

    is_initialized_ = false;

    //set state dimension
    n_x_ = 5;
    n_aug_ = n_x_ + 2;
    lambda_ = 3 - n_x_;

    //set weights
    weights_ = VectorXd::Constant(2 * n_aug_ + 1, 0.5 / (lambda_ + n_aug_));
    weights_(0) = lambda_ / (lambda_ + n_aug_);

    // Initialise matrix to hold sigma points in state space
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    // create a file for nis logging
    nis_file.open("nis.csv");
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    // Update timestamp. Time is measured in seconds.
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            // Convert radar from polar to cartesian coordinates and initialize state.
            float ro = meas_package.raw_measurements_(0);
            float theta = meas_package.raw_measurements_(1);
            float ro_dot = meas_package.raw_measurements_(2);
            x_(0) = ro * sin(theta);
            x_(1) = ro * cos(theta);
            x_(2) = ro_dot * sin(theta);
            x_(3) = ro_dot * cos(theta);
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            // LASER only measures location leave velocity default.
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    cout << "Prediction\n";
    Prediction(dt);

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    cout << "Update\n";
    Update(meas_package);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    //create augmented mean stat e
    VectorXd x_aug = VectorXd::Zero(n_aug_);
    x_aug.head(n_x_) = x_;

    //create augmented covariance matrix
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_aug_-2, n_aug_-2) = std_a_*std_a_;
    P_aug(n_aug_-1, n_aug_-1) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }

    // Predict sigma points
    for (int i = 0; i<Xsig_aug.cols(); i++)
    {
        float px             = Xsig_aug(0, i);
        float py             = Xsig_aug(1, i);
        float v              = Xsig_aug(2, i);
        float psi            = Xsig_aug(3, i);
        float psi_dot        = Xsig_aug(4, i);
        float mu_a           = Xsig_aug(5, i);
        float mu_psi_dot_dot = Xsig_aug(6, i);

        if (fabsf(psi_dot) < 0.000001) { //avoid division by zero
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
        a(3) = fmod(a(3), 2*M_PI); // angle normalization
        P_ += weights_(i) * a * a.transpose();
    }
}

/**
* Updates the state and the state covariance matrix given a measurement
* @param meas_package The measurement at k+1
*/
void UKF::Update(MeasurementPackage meas_package)
{
    // measurement dimentions
    int n_z = 0;

    // process noise matrix
    MatrixXd R;

    // Initialise measurement dimentions and process noise matrix depending on sensor type
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        n_z = 3;
        R = MatrixXd::Zero(n_z, n_z);
        R.diagonal() << std_radr_ * std_radr_,
                        std_radphi_ * std_radphi_,
                        std_radrd_ * std_radrd_;
    } else {
        n_z = 2;
        R = MatrixXd::Zero(n_z, n_z);
        R.diagonal() << std_laspx_ * std_laspx_,
                        std_laspy_ * std_laspy_;
    }

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space depending on sensor type
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            float px  = Xsig_pred_(0, i);
            float py  = Xsig_pred_(1, i);
            float v   = Xsig_pred_(2, i);
            float psi = Xsig_pred_(3, i);

            Zsig(0, i) = sqrt(px * px + py * py); // rho
            Zsig(1, i) = atan2(py, px); // phi
            Zsig(2, i) = v * (px * cos(psi) + py * sin(psi)) / Zsig(0, i); // rho_dot
        } else {
            Zsig(0, i) = Xsig_pred_(0, i);
            Zsig(1, i) = Xsig_pred_(1, i);
        }
    }

    //calculate mean predicted measurement
    VectorXd z_pred = Zsig * weights_;

    //calculate innovation covariance matrix S
    MatrixXd S = R;
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        MatrixXd a = Zsig.col(i) - z_pred; // residual
        //angle normalization
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            a(1) = fmod(a(1), 2*M_PI);
        }
        S += weights_(i) * a * a.transpose();
    }

    //calculate cross correlation matrix
    MatrixXd Tc = MatrixXd::Zero(n_x_, Zsig.rows());
    for (int i = 0; i < Zsig.cols(); i++) {
        Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
    }

    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    //angle normalization
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        z_diff(1) = fmod(z_diff(1), 2*M_PI);
    }

    //update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();

    //angle normalization
    x_(3) = fmod(x_(3), 2*M_PI);

    // NIS analysis
    static unsigned int radar_count = 0;
    static unsigned int radar_over_threshold_count = 0;
    static unsigned int lidar_count = 0;
    static unsigned int lidar_over_threshold_count = 0;
    static const float radar_threshold = 7.815;
    static const float lidar_threshold = 5.991;

    float NIS = (z_diff.transpose() * S.inverse() * z_diff)(0,0);

    // log the NIS value
    nis_file << meas_package.sensor_type_ << ',' << NIS << endl;

    // NIS threshold
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        ++radar_count;
        if (NIS > radar_threshold)
            ++radar_over_threshold_count;
    } else {
        ++lidar_count;
        if (NIS > lidar_threshold)
            ++lidar_over_threshold_count;
    }
    printf("nis over threshold rate: radar %.03f lidar %.03f\n",
            radar_over_threshold_count / float(radar_count),
            lidar_over_threshold_count / float(lidar_count));
}
