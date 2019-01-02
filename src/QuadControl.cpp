#include "Common.h"
#include "QuadControl.h"

#include <iostream>

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    
  // ****************************************************************************
  //  cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
  //  cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
  //  cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
  //  cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right
  // *****************************************************************************
 
  float l = L/sqrt(2);
  float total_thrust = collThrustCmd;
  float tau_x_over_l = momentCmd.x / l;
  float tau_y_over_l = momentCmd.y / l;
  float tau_z_over_K = -momentCmd.z / kappa;
    
  cmd.desiredThrustsN[0] = (total_thrust + tau_x_over_l + tau_y_over_l + tau_z_over_K) / 4.f; // front left, cw
  cmd.desiredThrustsN[1] = (total_thrust - tau_x_over_l + tau_y_over_l - tau_z_over_K) / 4.f; // front right, ccw
  cmd.desiredThrustsN[2] = (total_thrust + tau_x_over_l - tau_y_over_l - tau_z_over_K) / 4.f; // rear left, ccw
  cmd.desiredThrustsN[3] = (total_thrust - tau_x_over_l - tau_y_over_l + tau_z_over_K) / 4.f; // rear right, cw

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;
  
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
  V3F body_rates_error = pqrCmd - pqr;
    
  float p_dot = kpPQR.x * body_rates_error[0]; // u_bar_p.
  float q_dot = kpPQR.y * body_rates_error[1]; // u_bar_q.
  float r_dot = kpPQR.z * body_rates_error[2]; // u_bar_r.
  
  momentCmd.x = Ixx * p_dot;
  momentCmd.y = Iyy * q_dot;
  momentCmd.z = Izz * r_dot;
    
    
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
    float target_net_acc = -collThrustCmd/mass;
    float b_x_c_target = CONSTRAIN(accelCmd.x/target_net_acc, -maxTiltAngle, maxTiltAngle);
    float b_y_c_target = CONSTRAIN(accelCmd.y/target_net_acc, -maxTiltAngle, maxTiltAngle);

    float b_x_actual = R(0,2);
    float b_y_actual = R(1,2);

    float b_dot_x_c = kpBank * (b_x_c_target - b_x_actual);
    float b_dot_y_c = kpBank * (b_y_c_target - b_y_actual);

    float p_c = 1/R(2,2) * (R(1,0)*b_dot_x_c - R(0,0)*b_dot_y_c);
    float q_c = 1/R(2,2) * (R(1,1)*b_dot_x_c - R(0,1)*b_dot_y_c);

    pqrCmd.x = p_c;
    pqrCmd.y = q_c;
    pqrCmd.z = 0;
    
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  float z_error = posZCmd - posZ;
  float z_dot_error = velZCmd - velZ;
  integratedAltitudeError += z_error*dt;

  float P_term = kpPosZ * z_error;
  float D_term = kpVelZ * z_dot_error;
  float I_term = KiPosZ * integratedAltitudeError;

  // This should be the controller commanded global_acceleration to the Quad to make it move with accelZCmd.
  float global_acc_cmd = P_term + D_term + I_term + accelZCmd;

  float b_z_actual = R(2,2);
  float g = 9.81f;
  // This should be the controller body frame acceleration to the Quad to make it move with accelZCmd in the world frame.
  float body_acc_cmd = (global_acc_cmd - g)/b_z_actual;

  // We choose to contrain the acc of the quadcopter to a maximum and miminum acceleration.
  float constrained_body_acc_cmd = CONSTRAIN(body_acc_cmd, -maxAscentRate / dt, maxAscentRate / dt);

  // This is the controller commanded thrust to the Quad to achieve accelZCmd in the world frame.
  thrust = -mass * constrained_body_acc_cmd;
    
  /////////////////////////////// END STUDENT CODE ////////////////////////////
  
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
  // This block ensures that the quadcopter maintains a speed in the x,y axis is less than or equal to the set maximum speed for the quadcopter.
  if( velCmd.mag() > maxSpeedXY ){
      V3F unit_velCmd_vector = velCmd/velCmd.mag();
      velCmd = maxSpeedXY * unit_velCmd_vector;
  }

  V3F XY_error = posCmd - pos;
  V3F velXY_error = velCmd - vel;

  XY_error.z = 0.f; // XY_error.z = 0
  velXY_error.z = 0; // velXY.z = 0

  V3F P_term_xy = kpPosXY * XY_error;
  V3F D_term_xy = kpVelXY * velXY_error;

  accelCmd += (P_term_xy + D_term_xy);
  
  // This block ensures that the quadcopter maintains an acceleration in the x,y axis is less than or equal to the set maximum acceleration for the quadcopter.
  if( accelCmd.mag() > maxAccelXY ){
      V3F unit_accelCmd_vector = accelCmd/accelCmd.mag();
      accelCmd = maxAccelXY * unit_accelCmd_vector;
   }
  
  // z values for acc, velocity & positon were all set to zero because the Altitude Controller already has the responsibility of
  // controlling the quadcopter's altitude.
    
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b]. 
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
  if (yawCmd > 0){
     yaw = fmodf(yaw, 2*F_PI);
  }
  else {
      yawCmd = -fmodf(-1.f*yawCmd, 2*F_PI);
  }
  
  float yaw_err = yawCmd - yaw;

  if (yaw_err > F_PI) {
     yaw_err -= 2*F_PI;
  }
  else if (yaw_err < -F_PI) {
        yaw_err += 2*F_PI;
  }

  yawRateCmd = kpYaw * yaw_err;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
