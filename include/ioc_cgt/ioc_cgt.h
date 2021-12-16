#pragma once

#include <cmath>
#include <Eigen/Core>
#include <ros/time.h>
#include <geometry_msgs/WrenchStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <sensor_msgs/JointState.h>

#include <state_space_filters/filtered_values.h>
#include <cnr_controller_interface/cnr_joint_command_controller_interface.h>
#include <cnr_hardware_interface/posveleff_command_interface.h>
#include <cnr_hardware_interface/veleff_command_interface.h>
#include <control_msgs/GripperCommandAction.h>
#include <actionlib/client/simple_action_client.h>

namespace ect = eigen_control_toolbox;

namespace cnr
{
namespace control
{


/**
 * @brief The LQRCartImpedance class
 */
class LQRCartImpedance: public cnr::control::JointCommandController<
        hardware_interface::JointHandle, hardware_interface::VelocityJointInterface>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  LQRCartImpedance();
  bool doInit();
  bool doUpdate  (const ros::Time& time, const ros::Duration& period);
  bool doStarting(const ros::Time& time);
  bool doStopping(const ros::Time& time);

protected:

  std::mutex m_mtx;


  ect::FilteredVectorXd vel_fitler_sp_;
  ect::FilteredVectorXd wrench_fitler_;

  rosdyn::VectorXd dq_sp_;
  rosdyn::VectorXd q_sp_;
  rosdyn::VectorXd ddq_;
  rosdyn::VectorXd dq_;
  rosdyn::VectorXd q_;

  geometry_msgs::PoseStamped pose_sp_;
  Eigen::Affine3d T_base_targetpose_;
  
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;
  
  Eigen::MatrixXd Q_;
  Eigen::MatrixXd R_;

  Eigen::MatrixXd P_;
  Eigen::MatrixXd LQR_gain_;

  bool use_cartesian_reference_;

  bool first_cycle_;
  int count_update_;

  bool w_b_init_;
  bool use_filtered_wrench_;
  
  Eigen::Vector6d w_b_filt_;
  Eigen::Vector6d w_b_;
  Eigen::Vector6d w_b_0_;
  Eigen::Vector6d wrench_deadband_;

  rosdyn::ChainPtr chain_bs_;
  rosdyn::ChainPtr chain_bt_;

  size_t filtered_wrench_base_pub_;
  size_t wrench_base_pub_;
  size_t wrench_tool_pub_;
  size_t robot_wrench_pub_;
  size_t current_pose_pub_;

  double exponential_;
  double min_val_;

  Eigen::Matrix66d M_;
  Eigen::Matrix66d D_;
  Eigen::Matrix66d K_;

  void callback                (const geometry_msgs::WrenchStampedConstPtr&  msg );
  void setTargetPoseCallback   (const geometry_msgs::PoseStampedConstPtr&    msg );
  void setTargetJointsCallback (const sensor_msgs::JointStateConstPtr&       msg );
  bool solveRiccati(const Eigen::MatrixXd &A,
                    const Eigen::MatrixXd &B,
                    const Eigen::MatrixXd &Q,
                    const Eigen::MatrixXd &R,
                          Eigen::MatrixXd &P) ;

  
  
  
};


}
}
