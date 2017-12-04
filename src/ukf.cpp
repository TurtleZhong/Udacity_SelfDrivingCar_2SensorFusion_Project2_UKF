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
    std_a_ = 1.2;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.6;

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

    weights_ = VectorXd(2 * n_aug_ + 1);

    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    time_us_ = 0;

    lambda_ = 3 - n_x_;

    radar_NIS = 0.;

    lidar_NIS = 0.;

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

    if(!is_initialized_)
    {

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            double rho = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double rhodot = meas_package.raw_measurements_(2);

            x_ << rho * cos(phi), rho * sin(phi), 3, rhodot * cos(phi), rhodot * sin(phi);

            P_ << std_radr_*std_radr_, 0, 0, 0, 0,
                    0, std_radr_*std_radr_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, std_radphi_ * std_radphi_, 0,
                    0, 0, 0, 0, std_radphi_*std_radphi_;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            double px = meas_package.raw_measurements_(0);
            double py = meas_package.raw_measurements_(1);
            x_ << px, py, 3, 0.4, 0.0;

            P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
                    0, std_laspy_*std_laspy_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1;
        }

        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;

    }
    else
    {
        double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
        time_us_ = meas_package.timestamp_;

        Prediction(delta_t);

        /*measurements*/

        if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            UpdateRadar(meas_package);
        }
        else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            UpdateLidar(meas_package);
        }

    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    /*generated the sig_points*/
    //create sigma point matrix
    MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

    //calculate square root of P
    MatrixXd A = P_.llt().matrixL();
    /*******************************************************************************
   * generated the sig_points
   ******************************************************************************/
    double lambda_sqrt = sqrt(lambda_ + n_x_);
    Xsig.block<5,1>(0,0) = x_;
    MatrixXd state = MatrixXd::Zero(5,5);
    for(int i = 0; i <5; i++)
    {
        state.block<5,1>(0,i) = x_;
    }
    MatrixXd first = state + lambda_sqrt * A;
    MatrixXd second = state - lambda_sqrt * A;
    Xsig.block<5,5>(0,1) = first;
    Xsig.block<5,5>(0,6) = second;


    //create augmented mean vector
    VectorXd x_aug = VectorXd(7);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    /*******************************************************************************
    * Augmentation
    ******************************************************************************/

    x_aug.segment<5>(0) = x_;
    x_aug.segment<2>(5) = Eigen::Vector2d(0,0);
    P_aug = MatrixXd::Zero(7,7);
    P_aug.topLeftCorner(5,5) = P_;
    MatrixXd Q = MatrixXd(2,2);
    Q << std_a_*std_a_, 0,
            0, std_yawdd_*std_yawdd_;
    P_aug.block<2,2>(5,5) = Q;

    MatrixXd A_aug = P_aug.llt().matrixL();

    Xsig_aug.block<7,1>(0,0) = x_aug;
    MatrixXd state_aug = MatrixXd(7,7);
    for(int i = 0; i < 7; i++)
    {
        state_aug.block<7,1>(0,i) = x_aug;
    }
    MatrixXd aug_first = state_aug + sqrt(3) * A_aug;
    MatrixXd aug_second = state_aug - sqrt(3) * A_aug;
    Xsig_aug.block<7,7>(0,1) = aug_first;
    Xsig_aug.block<7,7>(0,n_aug_+1) = aug_second;


    /*sigma point prediction*/
    //create matrix with predicted sigma points as columns
    //    MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    /*double delta_t = 0.1;*/ //time diff in sec
    /*******************************************************************************
   * Prediction Points
   ******************************************************************************/

    for(int col = 0; col< 2 * n_aug_ + 1; col++)
    {
        VectorXd cols = Xsig_aug.block<7,1>(0,col);
        VectorXd cols_5 = cols.segment<5>(0);
        double niu_a = cols(5);
        double niu_phi = cols(6);
        if(cols_5(4)==0)
        {
            VectorXd part1(5);
            part1 << cols_5(2) * cos(cols_5(3)) *delta_t,
                    cols_5(2) * sin(cols_5(3)) *delta_t,
                    0,
                    0,
                    0;
            VectorXd part2(5);
            part2 << 0.5 * delta_t * delta_t * cos(cols_5(3)) * niu_a,
                    0.5 * delta_t * delta_t * sin(cols_5(3)) * niu_a,
                    delta_t * niu_a,
                    0.5 * delta_t * delta_t * niu_phi,
                    delta_t * niu_phi;

            Xsig_pred_.block<5,1>(0,col) = Xsig_aug.block<5,1>(0,col) + part1 + part2;
        }
        else
        {
            double phi_dot = cols_5(4);
            VectorXd part3(5);
            part3 << cols_5(2)*(sin(cols_5(3) + phi_dot * delta_t) - sin(cols_5(3))) / phi_dot,
                    cols_5(2)*(-cos(cols_5(3) + phi_dot * delta_t) + cos(cols_5(3))) / phi_dot,
                    0,
                    phi_dot * delta_t,
                    0;
            VectorXd part4(5);
            part4 << 0.5 * delta_t * delta_t * cos(cols_5(3)) * niu_a,
                    0.5 * delta_t * delta_t * sin(cols_5(3)) * niu_a,
                    delta_t * niu_a,
                    0.5 * delta_t * delta_t * niu_phi,
                    delta_t * niu_phi;

            //            std::cout << "part2 = \n" << part3 << "part3\n" << part4 << std::endl;
            Xsig_pred_.block<5,1>(0,col) = Xsig_aug.block<5,1>(0,col) + part3 + part4;
        }
    }

    //      //create vector for weights
    //      VectorXd weights = VectorXd(2*n_aug+1);

    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x_, n_x_);


    /*******************************************************************************
     * predicted mean and covariance
     ******************************************************************************/

    //set weights
    //predict state mean
    //predict state covariance matrix

    for(int i = 0; i < 2*n_aug_ +1; i++)
    {
        if(i==0)
        {
            weights_(i) = lambda_ / (lambda_ +  n_aug_);
        }
        else
        {
            weights_(i) = 0.5 / (lambda_ + n_aug_);
        }
    }

    for(int i = 0; i < 2*n_aug_ +1; i++)
    {
        x += weights_(i) * Xsig_pred_.block<5,1>(0,i);
    }

    for(int i = 0; i < 2*n_aug_ +1; i++)
    {
        P += weights_(i) * (Xsig_pred_.block<5,1>(0,i) - x) * (Xsig_pred_.block<5,1>(0,i) - x).transpose();
    }


    this->x_ = x;
    this->P_ = P;
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

    /*******************************************************************************
     * Predict measurement
     ******************************************************************************/
    //create matrix for sigma points in measurement space
    int n_z = 2;
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    //transform sigma points into measurement space
    //calculate mean predicted measurement
    //calculate measurement covariance matrix S

    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        VectorXd xk = Xsig_pred_.block<5,1>(0,i);
        double px = xk(0);
        double py = xk(1);
        Zsig.block<2,1>(0,i) = Eigen::Vector2d(px,py);

    }

    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        z_pred += weights_(i)*Zsig.block<2,1>(0,i);
    }

    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        S += weights_(i)*(Zsig.block<2,1>(0,i) - z_pred)*(Zsig.block<2,1>(0,i) - z_pred).transpose();
    }

    MatrixXd R = MatrixXd::Zero(2,2);
    R(0,0) = std_laspx_ *std_laspx_;
    R(1,1) = std_laspy_ * std_laspy_;
    S+=R;



    /*Update*/
    VectorXd z = Eigen::Vector2d(meas_package.raw_measurements_(0),
                                 meas_package.raw_measurements_(1));
    MatrixXd Tc = MatrixXd(n_x_,n_z);
    //calculate cross correlation matrix
    //calculate Kalman gain K;
    //update state mean and covariance matrix
    for(int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        Tc += weights_(i) * (Xsig_pred_.block<5,1>(0,i) - x_) * (Zsig.block<2,1>(0,i)-z_pred).transpose();
    }

    MatrixXd K = Tc * S.inverse();
    x_ = x_ + K*(z - z_pred);
    P_-= K*S*K.transpose();

    VectorXd diff_z = z - z_pred;
    while (diff_z(1)> M_PI) diff_z(1)-=2.*M_PI;
    while (diff_z(1)<-M_PI) diff_z(1)+=2.*M_PI;

    lidar_NIS = diff_z.transpose() * S.inverse() * diff_z;


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

    /*******************************************************************************
     * Predict measurement
     ******************************************************************************/
    //create matrix for sigma points in measurement space
    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    //transform sigma points into measurement space
    //calculate mean predicted measurement
    //calculate measurement covariance matrix S

    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        VectorXd xk = Xsig_pred_.block<5,1>(0,i);
        double px = xk(0);
        double py = xk(1);
        double v  = xk(2);
        double xphi = xk(3);
        double xphi_dot = xk(4);

        double rho = sqrt(px*px + py*py);
        double zphi = atan2(py,px);
        double zphi_dot = (px*cos(xphi)*v + py*sin(xphi)*v)/rho;
        Zsig.block<3,1>(0,i) = Eigen::Vector3d(rho,zphi,zphi_dot);

    }

    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        z_pred += weights_(i)*Zsig.block<3,1>(0,i);
    }

    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        S += weights_(i)*(Zsig.block<3,1>(0,i) - z_pred)*(Zsig.block<3,1>(0,i) - z_pred).transpose();
    }

    MatrixXd R = MatrixXd::Zero(3,3);
    R(0,0) = std_radr_ *std_radr_;
    R(1,1) = std_radphi_ * std_radphi_;
    R(2,2) = std_radrd_ * std_radrd_;
    S+=R;



    /*Update*/
    VectorXd z = Eigen::Vector3d(meas_package.raw_measurements_(0),
                                 meas_package.raw_measurements_(1),
                                 meas_package.raw_measurements_(2));
    MatrixXd Tc = MatrixXd(n_x_,n_z);
    //calculate cross correlation matrix
    //calculate Kalman gain K;
    //update state mean and covariance matrix
    for(int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        Tc += weights_(i) * (Xsig_pred_.block<5,1>(0,i) - x_) * (Zsig.block<3,1>(0,i)-z_pred).transpose();
    }

    MatrixXd K = Tc * S.inverse();
    x_ = x_ + K*(z - z_pred);
    P_-= K*S*K.transpose();

    VectorXd diff_z = z - z_pred;
    while (diff_z(1)> M_PI) diff_z(1)-=2.*M_PI;
    while (diff_z(1)<-M_PI) diff_z(1)+=2.*M_PI;

    radar_NIS = diff_z.transpose() * S.inverse() * diff_z;


}
