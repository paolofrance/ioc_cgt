# state state cartesian impedance
  lqr_impedance:
    type: cnr/control/LQRCartImpedance
    controlled_joints: all
    kin_update_period: 0.008
    
    robot_base_frame  : base_link
    robot_tip_frame   : tip
    force_sensor_frame: robotiq_ft_frame_id
    
    use_cartesian_reference: false                 # are you going to use a reference Pose or JointState?  default: false
    joint_target_topic:     "/joint_pos_target"   # incoming joints setpoint topic 
    pose_target:            "/target_cart_pose"   # incoming pose setpoint topic

    external_wrench_topic: "/robotiq_ft_wrench"   # topic for reading the wrench
    wrench_deadband: [5,5,5,.1,.1,.1]             # deadbnd in which the force is ignored
    use_filtered_wrench: true                     # do you want to filter the wrench?
    omega_wrench: 500                             # omega filtering

    M_r: [10,10,10,10,10,10]          # diagonal inertia matrix values
    K_r: [100,100,100,10,10,10]       # diagonal stiffness matrix values
    D_r: [50,50,50,10,10,10]          # diagonal damping matrix values
    damping_is_ratio: false           # is the damping defined above absolute or ratio?

    Q: [10,10,10,10,10,10,1,1,1,1,1,1]    # weight state matrix
    R: [1,1,1,1,1,1]                      # weight control matrix
