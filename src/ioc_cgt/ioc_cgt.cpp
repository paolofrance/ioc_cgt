#include <lqr_cart_impedance/lqr_cart_impedance.h>
#include <state_space_filters/filtered_values.h>
#include <eigen_matrix_utils/overloads.h>
#include <geometry_msgs/PoseStamped.h>
#include <pluginlib/class_list_macros.h>
#include <Eigen/Dense>
#include <eigen_conversions/eigen_msg.h>
#include <std_msgs/Float32.h>
#include <sensor_msgs/JointState.h>
#include <rosdyn_core/primitives.h>
#include <name_sorting/name_sorting.h>
#include <tf2_eigen/tf2_eigen.h>
#include <lqr_cart_impedance/utils.h>

PLUGINLIB_EXPORT_CLASS(cnr::control::LQRCartImpedance  , controller_interface::ControllerBase)








namespace cnr
{
namespace control
{


/**
 * @brief LQRCartImpedance::LQRCartImpedance
 */
LQRCartImpedance::LQRCartImpedance()
{
}

/**
 * @brief LQRCartImpedance::doInit
 * @return
 */
bool LQRCartImpedance::doInit()
{
  //INIT PUB/SUB
  std::string external_wrench_topic ;
  GET_AND_RETURN( this->getControllerNh(), "external_wrench_topic"  , external_wrench_topic );
  this->template add_subscriber<geometry_msgs::WrenchStamped>(
        external_wrench_topic,5,boost::bind(&LQRCartImpedance::callback,this,_1), false);
  
  GET_AND_DEFAULT( this->getControllerNh(), "use_cartesian_reference", use_cartesian_reference_, false);
  if(use_cartesian_reference_)
  {
    std::string pose_target;
    GET_AND_RETURN( this->getControllerNh(), "pose_target"  , pose_target);
    this->template add_subscriber<geometry_msgs::PoseStamped>(pose_target,5,boost::bind(&LQRCartImpedance::setTargetPoseCallback,this,_1), false);
  }
  else
  {
    std::string joint_target;
    GET_AND_RETURN( this->getControllerNh(), "joint_target_topic"  , joint_target);
    this->template add_subscriber<sensor_msgs::JointState>(joint_target,5,boost::bind(&LQRCartImpedance::setTargetJointsCallback,this,_1), false);
  }
  
  this->setPriority(this->QD_PRIORITY);
  {
      ect::FilteredVectorXd::Value dead_band;
      ect::FilteredVectorXd::Value saturation;
      ect::FilteredVectorXd::Value init_value;

      dead_band = 0.0 * this->chain().getDQMax();
      saturation = this->chain().getDQMax();
      init_value = dead_band;
      if(!vel_fitler_sp_.activateFilter ( dead_band, saturation, (10.0 / 2.0 / M_PI), this->m_sampling_period, init_value ))
      {
        CNR_RETURN_FALSE(this->logger());
      }
  }
  dq_sp_ = vel_fitler_sp_.getUpdatedValue();

  wrench_deadband_.setZero();
  w_b_            .setZero();


  urdf::Model urdf_model;
  if ( !urdf_model.initParam ( "/robot_description" ) ) {
      ROS_ERROR ( "Urdf robot_description '%s' does not exist", (  this->getControllerNamespace()+"/robot_description" ).c_str() );
      return false;
  }
  Eigen::Vector3d gravity;
  gravity << 0, 0, -9.806;

  std::string robot_base_frame;
  GET_AND_RETURN( this->getControllerNh(), "robot_base_frame"  , robot_base_frame);
  std::string robot_tip_frame;
  GET_AND_RETURN( this->getControllerNh(), "robot_tip_frame"   , robot_tip_frame);
  std::string force_sensor_frame;
  GET_AND_RETURN( this->getControllerNh(), "force_sensor_frame", force_sensor_frame);


  chain_bs_ = rosdyn::createChain ( urdf_model,robot_base_frame, force_sensor_frame, gravity );
  chain_bt_ = rosdyn::createChain ( urdf_model,robot_base_frame, robot_tip_frame   , gravity );
  
  // =========================== INIT VAR END

  std::vector<double> wrench_deadband(6,0);
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "wrench_deadband", wrench_deadband, 6, "<=" );
  wrench_deadband_   = Eigen::Vector6d( wrench_deadband.data() );

  GET_AND_DEFAULT(this->getControllerNh(),"use_filtered_wrench",use_filtered_wrench_,false);
  {
      double omega;
      GET_AND_DEFAULT(this->getControllerNh(),"omega_wrench",omega,10.0);
      ect::FilteredVectorXd::Value dead_band;
      ect::FilteredVectorXd::Value saturation;
      ect::FilteredVectorXd::Value init_value;

      dead_band  = wrench_deadband_;
      saturation = 1000.0 * dead_band;
      init_value = 0.0 * dead_band;
      if(!wrench_fitler_.activateFilter ( dead_band, saturation, (omega / (2 * M_PI)), this->m_sampling_period, init_value ))
      {
        CNR_RETURN_FALSE(this->logger());
      }
  }
  w_b_filt_ = wrench_fitler_.getUpdatedValue();

  std::vector<double> M_r(6,0), D_r(6,0), K_r(6,0);
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "M_r", M_r, 6 , "<=" );
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "K_r", K_r, 6 , "<"  );
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "D_r", D_r, 6 , "<"  );

  bool is_damping_ratio;
  GET_AND_RETURN( this->getControllerNh(), "damping_is_ratio", is_damping_ratio);


  M_ = Eigen::Vector6d( M_r.data() ).asDiagonal();
  K_ = Eigen::Vector6d( K_r.data() ).asDiagonal();

  if (is_damping_ratio)
  {
    Eigen::Vector6d d_tmp;
      for (int i=0; i<D_r.size();i++)
        d_tmp(i) = D_r.data()[i] * 2 * std::sqrt( M_r.data()[i] * K_r.data()[i] );
      D_ = d_tmp.asDiagonal();
  }
  else
      D_ = Eigen::Vector6d( D_r.data() ).asDiagonal();

  CNR_INFO(this->logger(),"mass: \n"<<M_);
  CNR_INFO(this->logger(),"spring: \n"<<K_);
  CNR_INFO(this->logger(),"damping: \n"<<D_);

  std::vector<double> Q(12,0), R(6,0);
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "Q", Q, 12 , "<" );
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "R", R, 6 , "<" );

  // system params initi
  
  A_.resize(12,12);
  B_.resize(12,6);
  Q_.resize(12,12);
  R_.resize(6,6);

  A_.setZero();
  B_.setZero();
  Q_.setZero();
  R_.setZero();
  
  Eigen::VectorXd Q_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Q.data(), Q.size());
  Eigen::VectorXd R_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(R.data(), R.size());
  
  A_.block(0,0,6,6) = Eigen::MatrixXd::Zero(6,6);
  A_.block(0,6,6,6) = Eigen::MatrixXd::Identity(6,6);
  A_.block(6,0,6,6) = -1*K_*M_.inverse();
  A_.block(6,6,6,6) = -1*D_*M_.inverse();
  
  B_.block(0,0,6,6) = Eigen::MatrixXd::Zero(6,6);
  B_.block(6,0,6,6) = M_.inverse();
  
  Q_ = Q_vec.asDiagonal();
  R_ = R_vec.asDiagonal();
  
  CNR_FATAL(this->logger(),"A\n"<<A_);
  CNR_FATAL(this->logger(),"B\n"<<B_);
  CNR_FATAL(this->logger(),"Q\n"<<Q_);
  CNR_FATAL(this->logger(),"R\n"<<R_);

  w_b_init_ = false;

  first_cycle_ = true;  

  q_sp_ .setZero();
  dq_sp_.setZero();
  q_    .setZero();
  dq_   .setZero();
  ddq_  .setZero();
  
  P_.setZero();
  LQR_gain_.setZero();
  
  filtered_wrench_base_pub_ = this->template add_publisher<geometry_msgs::WrenchStamped>("filtered_wrench_base",5);
  wrench_base_pub_          = this->template add_publisher<geometry_msgs::WrenchStamped>("wrench_base",5);
  wrench_tool_pub_          = this->template add_publisher<geometry_msgs::WrenchStamped>("wrench_tool",5);
  robot_wrench_pub_         = this->template add_publisher<geometry_msgs::WrenchStamped>("/robot_wrench",5);
  current_pose_pub_         = this->template add_publisher<geometry_msgs::PoseStamped>  ("/current_pose",5);
  
  CNR_RETURN_TRUE(this->logger());
}

/**
 * @brief LQRCartImpedance::doStarting
 * @param time
 */
bool LQRCartImpedance::doStarting(const ros::Time& time)
{
  CNR_TRACE_START(this->logger(),"Starting Controller");

  q_sp_  = this->getPosition();
  dq_sp_ = this->getVelocity();
  q_  = q_sp_;
  dq_ = dq_sp_;
  this->setCommandPosition(q_);
  this->setCommandVelocity(dq_);
  T_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  pose_sp_.pose = tf2::toMsg (T_base_targetpose_);
  
  CNR_WARN(this->logger(),pose_sp_);

  dq_sp_ = 0 * this->getVelocity();
  count_update_ = 0;

  solveRiccati(A_,B_,Q_,R_,P_);
  LQR_gain_ = R_.inverse()*B_.transpose()*P_;
  CNR_WARN(this->logger(),"gain matrix: \n"<<LQR_gain_);

  CNR_RETURN_TRUE(this->logger());
}

/**
 * @brief LQRCartImpedance::stopping
 * @param time
 */
bool LQRCartImpedance::doStopping(const ros::Time& time)
{
  CNR_TRACE_START(this->logger(),"Stopping Controller");
  CNR_RETURN_TRUE(this->logger());
}

/**
 * @brief LQRCartImpedance::doUpdate
 * @param time
 * @param period
 * @return
 */
bool LQRCartImpedance::doUpdate(const ros::Time& time, const ros::Duration& period)
{
  auto start = std::chrono::steady_clock::now();

  CNR_TRACE_START_THROTTLE_DEFAULT(this->logger());
  std::stringstream report;
  std::lock_guard<std::mutex> lock(m_mtx);

  count_update_++;

  if (first_cycle_)
  {
    Eigen::Vector3d x = this->chainState().toolPose().translation();
    q_sp_ = this->getPosition();
    dq_sp_ = this->getVelocity();
    T_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
    pose_sp_.pose = tf2::toMsg (T_base_targetpose_);
    first_cycle_ = false;
  }
  
  if(use_cartesian_reference_)
    tf2::fromMsg (pose_sp_.pose, T_base_targetpose_);
  else
    T_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  
  Eigen::Vector6d cart_target_vel_of_t_in_b  = chain_bt_->getJacobian(q_sp_)*dq_sp_;
  Eigen::Affine3d T_b_t = chain_bt_->getTransformation(q_);
  Eigen::Matrix6Xd J_of_t_in_b  = chain_bt_->getJacobian(q_);
  Eigen::Vector6d cart_vel_of_t_in_b  = J_of_t_in_b*dq_;
  Eigen::Vector6d cart_acc_nl_of_t_in_b  = chain_bt_->getDTwistNonLinearPartTool(q_,dq_); // DJ*Dq
  Eigen::Vector6d cart_acc_of_t_in_b; cart_acc_of_t_in_b.setZero();
//   Eigen::Vector6d cart_vel_of_t_in_b; cart_vel_of_t_in_b.setZero();
  Eigen::Matrix<double,6,1> cartesian_error_actual_target_in_b;
  rosdyn::getFrameDistance(T_b_t, T_base_targetpose_ , cartesian_error_actual_target_in_b);
  
  rosdyn::VectorXd state(cart_vel_of_t_in_b.size() + cartesian_error_actual_target_in_b.size());
  state << cartesian_error_actual_target_in_b, cart_vel_of_t_in_b;
  
  Eigen::Vector6d wrench;
  if (use_filtered_wrench_)
    wrench = w_b_filt_;
  else
    wrench = w_b_;
  
  rosdyn::VectorXd ur = -LQR_gain_ * state;
  
  Eigen::VectorXd dX = A_*state + B_*wrench + B_*ur;//
  
  cart_vel_of_t_in_b = dX.head(6);
  cart_acc_of_t_in_b = dX.tail(6); 

  Eigen::Matrix6Xd J_b = chain_bt_->getJacobian(q_);

  Eigen::FullPivLU<Eigen::MatrixXd> pinv_J ( J_b );
  pinv_J.setThreshold ( 1e-2 );
  if(pinv_J.rank()<6)
  {
    CNR_FATAL_THROTTLE(this->logger(),1.0,"count update: "<<count_update_<<", rank: "<<pinv_J.rank());
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(J_b, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if (svd.singularValues()(svd.cols()-1)==0)
    ROS_WARN_THROTTLE(1,"SINGULARITY POINT");
  else if (svd.singularValues()(0)/svd.singularValues()(svd.cols()-1) > 1e2)
    ROS_WARN_THROTTLE(1,"SINGULARITY POINT");

  ddq_ = svd.solve(cart_acc_of_t_in_b-cart_acc_nl_of_t_in_b);
  q_  += dq_  * period.toSec() + ddq_*std::pow(period.toSec(),2.0)*0.5;
  dq_ += ddq_ * period.toSec();
  
  Eigen::VectorXd qd = svd.solve(cart_vel_of_t_in_b);

  this->setCommandPosition( q_ );
  this->setCommandVelocity( dq_);
  
  Eigen::Affine3d T_bt = chain_bt_->getTransformation(q_);
  geometry_msgs::Pose cp = tf2::toMsg (T_bt);
  
  geometry_msgs::PoseStamped ps;
  ps.header.stamp = ros::Time::now();
  ps.pose = cp;
  
  this->publish(current_pose_pub_,ps);

  geometry_msgs::WrenchStamped robot_w;

  robot_w.header.frame_id = "ur5_base_link";
  robot_w.header.stamp = ros::Time::now();
  robot_w.wrench.force.x  = ur( 0 );
  robot_w.wrench.force.y  = ur( 1 );
  robot_w.wrench.force.z  = ur( 2 );
  robot_w.wrench.torque.x = ur( 3 );
  robot_w.wrench.torque.y = ur( 4 );
  robot_w.wrench.torque.z = ur( 5 );


  this->publish(robot_wrench_pub_,robot_w);
  
  
  auto mid = std::chrono::steady_clock::now();
  CNR_FATAL_COND(this->logger(),std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count()>=8000
                 ,"too much time to command: "<<std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count());

  CNR_RETURN_TRUE_THROTTLE_DEFAULT(this->logger());

  }

  /**
   * @brief LQRCartImpedance::callback
   * @param msg
   */
  void LQRCartImpedance::callback(const geometry_msgs::WrenchStampedConstPtr& msg )
  {
    if(!w_b_init_)
    {
      w_b_0_ ( 0 ) = msg->wrench.force.x;
      w_b_0_ ( 1 ) = msg->wrench.force.y;
      w_b_0_ ( 2 ) = msg->wrench.force.z;
      w_b_0_ ( 3 ) = msg->wrench.torque.x;
      w_b_0_ ( 4 ) = msg->wrench.torque.y;
      w_b_0_ ( 5 ) = msg->wrench.torque.z;

      w_b_init_ = true;
    }

    Eigen::Vector6d wrench_s;
    wrench_s( 0 ) = msg->wrench.force.x  - w_b_0_ ( 0 );
    wrench_s( 1 ) = msg->wrench.force.y  - w_b_0_ ( 1 );
    wrench_s( 2 ) = msg->wrench.force.z  - w_b_0_ ( 2 );
    wrench_s( 3 ) = msg->wrench.torque.x - w_b_0_ ( 3 );
    wrench_s( 4 ) = msg->wrench.torque.y - w_b_0_ ( 4 );
    wrench_s( 5 ) = msg->wrench.torque.z - w_b_0_ ( 5 );

    Eigen::Affine3d T_bs = chain_bs_->getTransformation ( this->getPosition() );
    Eigen::Affine3d T_bt = chain_bt_->getTransformation ( this->getPosition() );
    Eigen::Affine3d T_ts = T_bt.inverse() * T_bs;
    Eigen::Vector6d w_t = rosdyn::spatialDualTranformation ( wrench_s , T_ts );
    Eigen::Vector6d wrench;
    wrench = rosdyn::spatialRotation ( w_t, T_bt.linear() );

    for ( unsigned int idx=0; idx<6; idx++ )
    {
      if ( ( wrench ( idx ) >wrench_deadband_ ( idx ) ) )
      {
          w_b_ ( idx ) = wrench ( idx )-wrench_deadband_ ( idx );
      }
      else if ( ( wrench ( idx ) <-wrench_deadband_ ( idx ) ) )
      {
          w_b_ ( idx ) = wrench ( idx )+wrench_deadband_ ( idx );
      }
      else
      {
          w_b_ ( idx ) =0;
      }
    }

    geometry_msgs::WrenchStamped tool_w;

    tool_w.header.frame_id = "robotiq_ft_frame_id";
    tool_w.header.stamp = ros::Time::now();
    tool_w.wrench.force.x  = wrench_s( 0 );
    tool_w.wrench.force.y  = wrench_s( 1 );
    tool_w.wrench.force.z  = wrench_s( 2 );
    tool_w.wrench.torque.x = wrench_s( 3 );
    tool_w.wrench.torque.y = wrench_s( 4 );
    tool_w.wrench.torque.z = wrench_s( 5 );

    geometry_msgs::WrenchStamped base_w;

    base_w.header.frame_id = "ur5_base_link";
    base_w.header.stamp = ros::Time::now();
    base_w.wrench.force.x  = w_b_( 0 );
    base_w.wrench.force.y  = w_b_( 1 );
    base_w.wrench.force.z  = w_b_( 2 );
    base_w.wrench.torque.x = w_b_( 3 );
    base_w.wrench.torque.y = w_b_( 4 );
    base_w.wrench.torque.z = w_b_( 5 );

    geometry_msgs::WrenchStamped filter_base_w;

    wrench_fitler_.update(wrench);
    w_b_filt_ = wrench_fitler_.getUpdatedValue();

    filter_base_w.header.frame_id = "ur5_base_link";
    filter_base_w.header.stamp = ros::Time::now();
    filter_base_w.wrench.force.x  = w_b_filt_( 0 );
    filter_base_w.wrench.force.y  = w_b_filt_( 1 );
    filter_base_w.wrench.force.z  = w_b_filt_( 2 );
    filter_base_w.wrench.torque.x = w_b_filt_( 3 );
    filter_base_w.wrench.torque.y = w_b_filt_( 4 );
    filter_base_w.wrench.torque.z = w_b_filt_( 5 );


    this->publish(wrench_base_pub_,base_w);
    this->publish(filtered_wrench_base_pub_ ,filter_base_w);
    this->publish(wrench_tool_pub_,tool_w);

    if (wrench_s(5)>1.5)
    {
        control_msgs::GripperCommandActionGoal goal;
        goal.goal.command.position = -0.01;
        goal.goal.command.max_effort = 30.0;
//        this->publish(activate_gripper_pub,goal);
    }
    else if (wrench_s(5)<-1.5)
    {
        control_msgs::GripperCommandActionGoal goal;
        goal.goal.command.position = 0.08;
        goal.goal.command.max_effort = 30.0;
//        this->publish(activate_gripper_pub,goal);
    }

  }

  bool LQRCartImpedance::solveRiccati(const Eigen::MatrixXd &A,
                                 const Eigen::MatrixXd &B,
                                 const Eigen::MatrixXd &Q,
                                 const Eigen::MatrixXd &R, Eigen::MatrixXd &P)
  {

    const uint dim_x = A.rows();
    const uint dim_u = B.cols();

    Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
    Ham << A, -B * R.inverse() * B.transpose(), -Q, -A.transpose();

    Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);

    Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);
    int j = 0;
    for (int i = 0; i < 2 * dim_x; ++i) {
      if (Eigs.eigenvalues()[i].real() < 0.) {
        eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
        ++j;
      }
    }

    Eigen::MatrixXcd Vs_1, Vs_2;
    Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
    Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
    P = (Vs_2 * Vs_1.inverse()).real();

    return true;
  }

  void LQRCartImpedance::setTargetPoseCallback(const geometry_msgs::PoseStampedConstPtr& msg)
  {
    try
    {
      pose_sp_ = *msg;
//       CNR_FATAL(this->logger(),"new cartesian setpoint recived:"<<pose_sp_);
    }
    catch(...)
    {
      ROS_ERROR("Something wrong in target callback");
    }
  }

void LQRCartImpedance::setTargetJointsCallback(const sensor_msgs::JointStateConstPtr& msg)
{
  try
  {
    sensor_msgs::JointState tmp_msg=*msg;
    if (!name_sorting::permutationName(this->jointNames(),tmp_msg.name,tmp_msg.position,tmp_msg.velocity,tmp_msg.effort))
    {
      CNR_ERROR(this->logger(),"joints not found");
      return;
    }
    for (unsigned int iAx=0;iAx<q_sp_.rows();iAx++)
    {
//       CNR_INFO(this->logger(),"new joint setpoint recived");
      q_sp_(iAx)=tmp_msg.position.at(iAx);
      dq_sp_(iAx)=tmp_msg.velocity.at(iAx);
    }

  }
  catch(...)
  {
    CNR_ERROR(this->logger(),"Something wrong in target callback");
  }
}

  }
  }
