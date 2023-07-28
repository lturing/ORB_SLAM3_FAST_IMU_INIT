/**************************************************************************
* VIG-Init
*
* Copyright SenseTime. All Rights Reserved.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
**************************************************************************/

#ifndef IMUINITIALIZER_H
#define IMUINITIALIZER_H

// Recipie of ORB-SLAM
// init_recipie(
//     solve_gyro_bias,
//     solve_scale_given_gravity,
//     solve_gravity_scale,
//     skip,
//     refine_scale_ba_via_gravity);

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "ImuTypes.h"
#include "Frame.h"
#include "MapPoint.h"
#include "ORBVocabulary.h"
#include "KeyFrameDatabase.h"
#include "ORBextractor.h"
#include "ORBmatcher.h"
#include "Settings.h"
#include "GeometricCamera.h"

namespace ORB_SLAM3
{

class Frame;
class MapPoint;
class Settings;

class ImuInitializer {
  public:
    ImuInitializer(std::vector<Frame>& frames, Settings *settings);
    ~ImuInitializer();

    //bool init_sfm();
    bool sfm_check(Frame& F1, Frame& F2);
    bool init_imu();
    void preintegrate();

    void solve_gyro_bias();
    void solve_gravity_scale_velocity();
    void refine_scale_velocity_via_gravity();

    void solve_gravity_scale();
    void solve_scale_ba_given_gravity(); 
    void refine_scale_ba_via_gravity(); 

    void reset_states();

    Eigen::Matrix<float, 3, 2> s2_tangential_basis(Eigen::Vector3f &x);

    //void apply_init();
    

    // IMU bias
    //IMU::Bias imuBias;
    Eigen::Vector3f ba, bg, gravity;
    float scale;
    std::vector<Eigen::Vector3f> velocities;
    std::vector<Frame> frames;
    Settings *settings;
    GeometricCamera* camera;
    std::vector<cv::Point3f> iniP3D;
    //Eigen::Vector3f gravity;
    //Eigen::Vector3f mVw;
    float PVIO_GRAVITY_NOMINAL;

};

}

#endif
