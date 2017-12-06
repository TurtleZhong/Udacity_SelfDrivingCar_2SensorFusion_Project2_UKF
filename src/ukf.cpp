#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // state vector dimension
    n_x_ = 5;

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1.8;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.3;

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

    /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */



    ///* Augmented state dimension
    n_aug_ = n_x_ + 2;

    ///* Sigma points dimension
    n_sig_ = 2 * n_aug_ + 1;

    lambda_ = 3 - n_aug_;

    // Initialize weights.
    weights_ = VectorXd(n_sig_);
    for(int i = 0; i < n_sig_; i++)
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

    Xsig_pred_ = MatrixXd(n_x_, n_sig_);

    /*measurement covariance*/

    R_lidar_ = MatrixXd(2, 2);
    R_lidar_ << std_laspx_*std_laspx_,0,
            0,std_laspy_*std_laspy_;

    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;

    NIS_laser_ = 0.;
    NIS_radar_ = 0.;



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

    // skip predict/update if sensor type is ignored
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
        (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_))
    {

      /*****************************************************************************
      *  Initialization
      ****************************************************************************/
      if (!is_initialized_) {

        x_ << 1, 1, 1, 1, 0.1;

        P_ = MatrixXd::Identity(5,5);

        time_us_ = meas_package.timestamp_;

        if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
        {

          x_(0) = meas_package.raw_measurements_(0);
          x_(1) = meas_package.raw_measurements_(1);

        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
        {

          double rho = meas_package.raw_measurements_(0);
          double phi = meas_package.raw_measurements_(1);
          x_(0) = rho * cos(phi);
          x_(1) = rho * sin(phi);
        }

        is_initialized_ = true;

        return;
      }

      //compute the time elapsed between the current and previous measurements
      float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
      time_us_ = meas_package.timestamp_;

      Prediction(delta_t);


      if (meas_package.sensor_type_ == MeasurementPackage::LASER)
      {
        UpdateLidar(meas_package);
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
      {
        UpdateRadar(meas_package);
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

    // step1. generate sigma points.
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd( n_aug_, n_sig_ );

    //calculate square root of P
    MatrixXd A = P_aug.llt().matrixL();
    Xsig_aug.col(0) = x_aug;
    double sqrt_lambda = sqrt(lambda_ + n_aug_);
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_aug.col( i + 1 ) = x_aug + sqrt_lambda * A.col(i);
        Xsig_aug.col( i + 1 + n_aug_ ) = x_aug - sqrt_lambda * A.col(i);
    }

    // step2: Predict Sigma Points.

    for (int i = 0; i< n_sig_; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawdot = Xsig_aug(4,i);
        double niu_a = Xsig_aug(5,i);
        double niu_yawdd = Xsig_aug(6,i);

        double px_p, py_p;

        if (fabs(yawdot) > 0.001)
        {
            px_p = p_x + v/yawdot * ( sin (yaw + yawdot*delta_t) - sin(yaw));
            py_p = p_y + v/yawdot * ( cos(yaw) - cos(yaw+yawdot*delta_t) );
        }
        else
        {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawdot*delta_t;
        double yawdot_p = yawdot;

        //add noise
        px_p = px_p + 0.5*niu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*niu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + niu_a*delta_t;

        yaw_p = yaw_p + 0.5*niu_yawdd*delta_t*delta_t;
        yawdot_p = yawdot_p + niu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawdot_p;
    }


    // step3 Predict Mean and Covariance

    //create vector for predicted state
    VectorXd x = VectorXd::Zero(n_x_);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd::Zero(n_x_, n_x_);

    for(int i = 0; i < n_sig_; i++)
    {
        x += weights_(i) * Xsig_pred_.block<5,1>(0,i);
    }

    for(int i = 0; i < n_sig_; i++)
    {
        //        P += weights_(i) * (Xsig_pred_.block<5,1>(0,i) - x) * (Xsig_pred_.block<5,1>(0,i) - x).transpose();
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        //angle normalization
        NormalizeAngle(x_diff,3);

        P = P + weights_(i) * x_diff * x_diff.transpose();
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
    // 1. Predit measurement
    int n_z = 2;

    MatrixXd Zsig = MatrixXd(n_z, n_sig_);
    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        VectorXd xk = Xsig_pred_.block<5,1>(0,i);
        double px = xk(0);
        double py = xk(1);
        Zsig.block<2,1>(0,i) = Eigen::Vector2d(px,py);

    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        z_pred += weights_(i)*Zsig.block<2,1>(0,i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sig_; i++)
    {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    S += R_lidar_;
    // step2 Update state
    VectorXd z = meas_package.raw_measurements_;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    Tc.fill(0.0);
    for(int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        Tc += weights_(i) * (Xsig_pred_.block<5,1>(0,i) - x_) * (Zsig.block<2,1>(0,i)-z_pred).transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //update state mean and covariance matrix
    x_ +=  K * z_diff;
    P_ -= K*S*K.transpose();

    //NIS Lidar Update
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
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
    // Radar measument dimension
    int n_z = 3;
    // 1. Predict measurement
    MatrixXd Zsig = MatrixXd(n_z, n_sig_);
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

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for(int i =0; i < 2 * n_aug_ + 1; i++)
    {
        z_pred += weights_(i)*Zsig.block<3,1>(0,i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sig_; i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // normalization
        NormalizeAngle(z_diff, 1);

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    S += R_radar_;

    // 2. Update state
    // Incoming radar measurement
    VectorXd z = meas_package.raw_measurements_;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    Tc.fill(0.0);
    for (int i = 0; i < n_sig_; i++)
    {

        VectorXd z_diff = Zsig.col(i) - z_pred;
        //normalization
        NormalizeAngle(z_diff, 1);
        // state diff
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        NormalizeAngle(x_diff, 3);
        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //normalization
    NormalizeAngle(z_diff, 1);
    //update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K*S*K.transpose();

    //NIS Update
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 *  Normalized the component `index` of the vector `vector` to be inside [-M_PI, M_PI] interval.
 */
void UKF::NormalizeAngle(VectorXd &vector, int index)
{
    while (vector(index)> M_PI) vector(index)-=2.*M_PI;
    while (vector(index)<-M_PI) vector(index)+=2.*M_PI;
}
