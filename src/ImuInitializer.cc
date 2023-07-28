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

#include "ImuInitializer.h"
#include "ImuTypes.h"
#include <bits/stdc++.h>
#include <iostream>

using namespace std;

namespace ORB_SLAM3
{
using namespace Eigen;

ImuInitializer::ImuInitializer(std::vector<Frame>& frames, Settings *settings): frames(frames),
settings(settings) {
    camera = settings->camera1();
    PVIO_GRAVITY_NOMINAL = 9.80665;
}

ImuInitializer::~ImuInitializer() {};

void ImuInitializer::reset_states() {
    ba.setZero();
    bg.setZero();
    gravity.setZero();
    scale = 1.0;
    velocities.resize(frames.size(), Eigen::Vector3f::Zero());

}

/*
bool ImuInitializer::sfm_check_v1() {
    // [1] try initializing using raw_map
    for (int i = 0; i < frames.size(); i++)
    {
        frames[i].ComputeBoW();
    }
    
    std::vector<cv::Point2f> vMatched;
    std::vector<int> iniMatches;

    vMatched.resize(frames[0].mvKeysUn.size());
    for(size_t i=0; i<frames[0].mvKeysUn.size(); i++)
        vMatched[i]=frames[0].mvKeysUn[i].pt;

    fill(iniMatches.begin(),iniMatches.end(),-1);

    // Find correspondences
    ORBmatcher matcher(0.9,true);
    int num_frames = frames.size();

    int nmatches = matcher.SearchForInitialization(frames[0],frames[num_frames-1],vMatched,iniMatches,100);


    double total_parallax = 0;
    int track_num = 0;
    for(size_t i = 0; i < iniMatches.size(); i++)
    {
        if (iniMatches[i] >= 0)
        {
            int j = iniMatches[i];
            total_parallax += (frames[0].mvKeysUn[i] - frames[num_frames -1].mvKeysUn[j]).norm();
            track_num += 1;
        }
    }

    //if (track_num < (int)config->initializer_min_matches()) return false;
    if (track_num < initializer_min_matches) 
        return false;

    total_parallax /= track_num;

    if (total_parallax < initializer_min_parallax)
        return false;

    Sophus::SE3f Tcw;
    vector<bool> vbTriangulated; // Triangulated Correspondences (mvIniMatches)
    set<MapPoint*> vpoints;

    //vpoints.find(vp) != vpoints.end()

    if(camera->ReconstructWithTwoViews(frames[0].mvKeysUn,frames[num_frames-1].mvKeysUn,iniMatches,Tcw,iniP3D,vbTriangulated))
    {
        for(size_t i=0, iend=iniMatches.size(); i<iend;i++)
        {
            if(iniMatches[i]>=0 && !vbTriangulated[i])
            {
                iniMatches[i]=-1;
                nmatches--;
            }
        }

        if (nmatches < initializer_min_triangulation)
            return false;
        
        // Set Frame Poses
        //frames[0].SetPose(Sophus::SE3f());
        //frames[num_frames-1].SetPose(Tcw);

        return true;
    }

    return false;

}
*/


bool ImuInitializer::sfm_check(Frame& F1, Frame& F2) {
    F1.ComputeBoW();
    F2.ComputeBoW();

    int initializer_min_matches = 10;
    int initializer_min_parallax = 6;
    int initializer_min_triangulation = 6;

    // Find correspondences
    //ORBmatcher matcher(0.9,true);
    ORBmatcher matcher(0.9,true);
    vector<int> vpMatches12;
    int nmatches = matcher.SearchByBoW(F1, F2, vpMatches12);


    double total_parallax = 0;
    int track_num = 0;
    for(size_t i = 0; i < vpMatches12.size(); i++)
    {
        if (vpMatches12[i] >= 0)
        {
            int j = vpMatches12[i];
            total_parallax += sqrt(pow(F1.mvKeysUn[i].pt.x - F2.mvKeysUn[j].pt.x, 2) + pow(F1.mvKeysUn[i].pt.y - F2.mvKeysUn[j].pt.y, 2));

            track_num += 1;
        }
    }

    //if (track_num < (int)config->initializer_min_matches()) return false;
    if (track_num < initializer_min_matches) 
    {
        cout << "rejected by initializer_min_matches, track_num=" << track_num << endl;
        return false;
    }

    total_parallax /= track_num;

    if (total_parallax < initializer_min_parallax)
    {
        cout << "rejected by initializer_min_parallax, total_parallax=" << total_parallax << endl;
        return false;
    }

    return true;

    Sophus::SE3f Tcw;
    vector<bool> vbTriangulated; // Triangulated Correspondences (mvIniMatches)
    std::vector<cv::Point3f> iniP3D;

    //vpoints.find(vp) != vpoints.end()


    
    if(camera->ReconstructWithTwoViews_v2(F1.mvKeysUn, F2.mvKeysUn,vpMatches12,Tcw,iniP3D,vbTriangulated))
    {
        for(size_t i=0, iend=vpMatches12.size(); i<iend;i++)
        {
            if(vpMatches12[i]>=0 && !vbTriangulated[i])
            {
                vpMatches12[i]=-1;
                nmatches--;
            }
        }

        if (nmatches < initializer_min_triangulation)
        {
            cout << "rejected by initializer_min_triangulation, nmatches=" << nmatches << endl;
            return false;
        }
            
        
        // Set Frame Poses
        //frames[0].SetPose(Sophus::SE3f());
        //frames[num_frames-1].SetPose(Tcw);

        return true;
    }

    
    cout << "rejected by initializer_min_triangulation(failed), nmatches=" << nmatches << endl;
    return false;
    
}




bool ImuInitializer::init_imu() {
    
    //check imu observibility
    
    Eigen::Vector3f sum_g;
    for (int i = 1; i < frames.size(); i++)
    {
        float dt = frames[i].mpImuPreintegratedFrame->delta.t;
        Eigen::Vector3f tmp_g = frames[i].mpImuPreintegratedFrame->delta.v / dt; 
        sum_g += tmp_g;
    }

    Eigen::Vector3f aver_g;
    aver_g = sum_g * 1.0 / ((int)frames.size() - 1);
    double var = 0;
    for (int i = 1; i < frames.size(); i++)
    {
        float dt = frames[i].mpImuPreintegratedFrame->delta.t;
        Eigen::Vector3f tmp_g = frames[i].mpImuPreintegratedFrame->delta.v / dt; 
        var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
    }

    var = sqrt(var / ((int)frames.size() - 1));
    //ROS_WARN("IMU variation %f!", var);
    if(var < 0.4)
    {
        return false;
    }
    

    float scale_min_thre = 0.002;
    float scale_max_thre = 1.0;

    float bg_min_thre = -0.2;
    float bg_max_thre = 0.2;

    float bg_min, bg_max;
    float ba_min, ba_max;
    

    reset_states();
    solve_gyro_bias();

    bg_max = bg.maxCoeff();
    bg_min = bg.minCoeff();
    std::cout << std::fixed << std::setprecision(5) << "bg=" << bg.transpose() << std::endl;

    if (bg_min < bg_min_thre || bg_max > bg_max_thre)
        return false;
    
    return true;

    solve_gravity_scale_velocity();

    if (scale < scale_min_thre || scale > scale_max_thre) 
    {
        //std::cout << "scale=" << scale << ", failed in solve_gravity_scale_velocity" << std::endl;
        return false;
    }


    refine_scale_velocity_via_gravity();
    if (scale < scale_min_thre || scale > scale_max_thre)
    {
        //std::cout << "scale=" << scale << ", failed in refine_scale_velocity_via_gravity" << std::endl;
        return false;
    } 

    // static const vector<3> gravity_nominal{0, 0, -PVIO_GRAVITY_NOMINAL};

    //Eigen::Vector3f gravity_nominal{0, 0, -PVIO_GRAVITY_NOMINAL};
    //Eigen::Quaternionf q = quaternion::FromTwoVectors(gravity, gravity_nominal);
    
    //Sophus::SE3f Twg(mRwg.cast<float>().transpose(), Eigen::Vector3f::Zero());

    std::cout << std::fixed << std::setprecision(5) << "scale=" << scale << ", bg=" << bg.transpose() << ", ba=" << ba.transpose() << ", gravity=" << gravity.transpose() << std::endl;
    
    preintegrate();
    return true;
    //return apply_init();
}







/*
bool ImuInitializer::init_imu() {
    
    float scale_min_thre = 0.002;
    float scale_max_thre = 1.0;

    float bg_min_thre = -0.1;
    float bg_max_thre = 0.1;

    float bg_min, bg_max;
    float ba_min, ba_max;
    

    reset_states();
    solve_gyro_bias();

    bg_max = bg.maxCoeff();
    bg_min = bg.minCoeff();
    if (bg_min < bg_min_thre || bg_max > bg_max_thre)
        return false;
    
    //solve_gravity_scale_velocity();

    solve_gravity_scale();
    if (scale < scale_min_thre || scale > scale_max_thre) 
    {
        //std::cout << "scale=" << scale << ", failed in solve_gravity_scale_velocity" << std::endl;
        return false;
    }

    solve_scale_ba_given_gravity();
    ba_max = ba.maxCoeff();
    ba_min = ba.minCoeff();
    if (scale < scale_min_thre || scale > scale_max_thre || ba_max > 1.0 || ba_min < -1.0) 
    {
        //std::cout << "scale=" << scale << ", failed in solve_gravity_scale_velocity" << std::endl;
        return false;
    }

    //refine_scale_velocity_via_gravity();
    refine_scale_ba_via_gravity();
    ba_max = ba.maxCoeff();
    ba_min = ba.minCoeff();
    if (scale < scale_min_thre || scale > scale_max_thre || ba_max > 1.0 || ba_min < -1.0)
    {
        //std::cout << "scale=" << scale << ", failed in refine_scale_velocity_via_gravity" << std::endl;
        return false;
    } 

    // static const vector<3> gravity_nominal{0, 0, -PVIO_GRAVITY_NOMINAL};

    //Eigen::Vector3f gravity_nominal{0, 0, -PVIO_GRAVITY_NOMINAL};
    //Eigen::Quaternionf q = quaternion::FromTwoVectors(gravity, gravity_nominal);
    
    //Sophus::SE3f Twg(mRwg.cast<float>().transpose(), Eigen::Vector3f::Zero());

    std::cout << std::fixed << std::setprecision(5) << "scale=" << scale << ", bg=" << bg.transpose() << ", ba=" << ba.transpose() << ", gravity=" << gravity.transpose() << std::endl;
    
    preintegrate();
    return true;
    //return apply_init();
}

*/

// IMU::Bias b(0,0,0,0,0,0); frames[i].mpImuPreintegratedFrame->Initialize(b);
void ImuInitializer::preintegrate() {
    IMU::Bias b(ba[0],ba[1],ba[2],bg[0],bg[1],bg[2]);
    for (int i = 1; i < frames.size(); i++)
    {
        frames[i].mpImuPreintegratedFrame->Reintegrate_v2(b);
    }

}


/*
void ImuInitializer::preintegrate_old() {
    IMU::Bias b(ba[0],ba[1],ba[2],bg[0],bg[1],bg[2]);
    for (int i = 1; i < frames.size(); i++)
    {
        int n = frames[i].mvImus.size()-1;
        for (int j = 0; j < n; j++)
        {
            float tstep;
            Eigen::Vector3f acc, angVel;
            if ((j == 0) && (j < (n-1)))
            {
                float tab = frames[i].mvImus[j+1].t-frames[i].mvImus[j].t;
                float tini = frames[i].mvImus[j].t-frames[i-1].mTimeStamp;

                acc = (frames[i].mvImus[j].a+frames[i].mvImus[j+1].a-
                    (frames[i].mvImus[j+1].a-frames[i].mvImus[j].a)*(tini/tab))*0.5f;
                angVel = (frames[i].mvImus[j].w+frames[i].mvImus[j+1].w-
                    (frames[i].mvImus[j+1].w-frames[i].mvImus[j].w)*(tini/tab))*0.5f;
                tstep = frames[i].mvImus[j+1].t-frames[i-1].mTimeStamp;
            }
            else if (j < (n-1))
            {
                acc = (frames[i].mvImus[j].a+frames[i].mvImus[j+1].a)*0.5f;
                angVel = (frames[i].mvImus[j].w+frames[i].mvImus[j+1].w)*0.5f;
                tstep = frames[i].mvImus[j+1].t-frames[i].mvImus[j].t;
            }
            else if ( (j > 0) && (j == (n-1)))
            {
                float tab = frames[i].mvImus[j+1].t-frames[i].mvImus[j].t;
                float tend = frames[i].mvImus[j+1].t-frames[i].mTimeStamp;
                acc = (frames[i].mvImus[j].a+frames[i].mvImus[j+1].a-
                        (frames[i].mvImus[j+1].a-frames[i].mvImus[j].a)*(tend/tab))*0.5f;
                angVel = (frames[i].mvImus[j].w+frames[i].mvImus[j+1].w-
                        (frames[i].mvImus[j+1].w-frames[i].mvImus[j].w)*(tend/tab))*0.5f;
                tstep = frames[i].mTimeStamp-frames[i].mvImus[j].t;
            }
            else if ( (j == 0) && (j == (n-1)))
            {
                acc = frames[i].mvImus[j].a;
                angVel = frames[i].mvImus[j].w;
                tstep = frames[i].mTimeStamp-frames[i-1].mTimeStamp;
            }

            if (tstep > 0)
            {
                frames[i].mpImuPreintegratedFrame->IntegrateNewMeasurement_v2(acc,angVel,tstep, true, false);
            }

        }
    }

}
*/


void ImuInitializer::solve_gyro_bias()
{
    preintegrate();
    Eigen::Matrix3f A;
    Eigen::Vector3f b;

    A.setZero();
    b.setZero();
    //std::cout << "frames[0].mnId=" << frames[0].mnId << std::endl;

    for (size_t j = 1; j < frames.size(); j++){
        const size_t i = j - 1;

        Sophus::SE3<float> pose_i = frames[i].GetImuPose();
        Sophus::SE3<float> pose_j = frames[j].GetImuPose();

        Eigen::Quaternionf& dq = frames[j].mpImuPreintegratedFrame->delta.q;
        Eigen::Matrix3f& dq_dbg = frames[j].mpImuPreintegratedFrame->jacobian.dq_dbg;
        if(0) //j == 1 || true)
        {
            std::cout << "pose_i=" << pose_i.matrix() << std::endl;
            std::cout << "pose_j=" << pose_j.matrix() << std::endl;
            std::cout << "dq(w, x, y, z)=" << dq.w() << ", " << dq.x() << ", " << dq.y() << ", " << dq.z() << std::endl;
            std::cout << "dq_dbg=" << dq_dbg << std::endl;
        }
        
        A += dq_dbg.transpose() * dq_dbg;
        b += dq_dbg.transpose() * ORB_SLAM3::IMU::logmap((pose_i.unit_quaternion() * dq).conjugate() * pose_j.unit_quaternion());
    }
    
    //std::cout << std::fixed << std::setprecision(5) << "A=" << A << std::endl;
    //std::cout << "b=" << b.transpose() << std::endl;

    Eigen::JacobiSVD<Eigen::Matrix3f> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    bg = svd.solve(b);
    //std::cout << "bg=" << bg.transpose() << std::endl;
    //std::cout << "A=" << A << std::endl;
    //std::cout << "b=" << b << std::endl;

    //std::cout << "bg=" << bg.transpose() << std::endl;
}


void ImuInitializer::solve_gravity_scale_velocity()
{
    preintegrate();
    int N = frames.size();
    Eigen::MatrixXf A;
    Eigen::VectorXf b;
    A.resize((N - 1) * 6, 3 + 1 + 3 * N);
    b.resize((N - 1) * 6);
    A.setZero();
    b.setZero();

    for (size_t j = 1; j < N; ++j) {
        const size_t i = j - 1;

        ORB_SLAM3::IMU::Delta &delta = frames[j].mpImuPreintegratedFrame->delta;
        Sophus::SE3<float> pose_i = frames[i].GetPose().inverse(); // Twc
        Sophus::SE3<float> pose_j = frames[j].GetPose().inverse();

        Sophus::SE3<float> pose_i_w_i = frames[i].GetImuPose(); // Twb
        Sophus::SE3<float> pose_j_w_i = frames[j].GetImuPose();

        Sophus::SE3<float> Tbc = frames[i].mImuCalib.mTcb.inverse();


        // translation()
        // rotationMatrix()

        A.block<3, 3>(i * 6, 0) += -0.5 * delta.t * delta.t * Eigen::Matrix3f::Identity();
        A.block<3, 1>(i * 6, 3) += pose_j.translation() - pose_i.translation();
        A.block<3, 3>(i * 6, 4 + i * 3) += -delta.t * Eigen::Matrix3f::Identity();

        b.segment<3>(i * 6) += pose_i_w_i.unit_quaternion() * delta.p + (pose_j_w_i.unit_quaternion() * Tbc.translation() - pose_i_w_i.unit_quaternion() * Tbc.translation());
        A.block<3, 3>(i * 6 + 3, 0) += -delta.t * Eigen::Matrix3f::Identity();
        A.block<3, 3>(i * 6 + 3, 4 + i * 3) += -Eigen::Matrix3f::Identity();
        A.block<3, 3>(i * 6 + 3, 4 + j * 3) += Eigen::Matrix3f::Identity();
        b.segment<3>(i * 6 + 3) += pose_i_w_i.unit_quaternion() * delta.v;
        
    }

    Eigen::VectorXf x = A.fullPivHouseholderQr().solve(b);
    
    gravity = x.segment<3>(0).normalized() * PVIO_GRAVITY_NOMINAL;
    scale = x(3);
    for (size_t i = 0; i < frames.size(); ++i) {
        velocities[i] = x.segment<3>(4 + i * 3);
    }
}

void ImuInitializer::refine_scale_velocity_via_gravity()
{
    static const double damp = 0.1;
    preintegrate();
    int N = frames.size();
    Eigen::MatrixXf A;
    Eigen::VectorXf b;
    Eigen::VectorXf x;

    A.resize((N - 1) * 6, 2 + 1 + 3 * N);
    b.resize((N - 1) * 6);
    x.resize(2 + 1 + 3 * N);

    for (size_t iter = 0; iter < 5; ++iter) {
        A.setZero();
        b.setZero();
        x.setZero();
        Eigen::Matrix<float, 3, 2> Tg = s2_tangential_basis(gravity);

        for (size_t j = 1; j < frames.size(); j++){
            const size_t i = j - 1;

            const ORB_SLAM3::IMU::Delta &delta = frames[j].mpImuPreintegratedFrame->delta;
            Sophus::SE3<float> pose_i = frames[i].GetPose().inverse(); // Twc
            Sophus::SE3<float> pose_j = frames[j].GetPose().inverse(); // Twc

            Sophus::SE3<float> pose_i_w_i = frames[i].GetImuPose();
            Sophus::SE3<float> pose_j_w_i = frames[j].GetImuPose();

            Sophus::SE3<float> Tbc = frames[i].mImuCalib.mTcb.inverse();;

            A.block<3, 2>(i * 6, 0) += -0.5 * delta.t * delta.t * Tg;
            A.block<3, 1>(i * 6, 2) += pose_j.translation() - pose_i.translation();
            
            A.block<3, 3>(i * 6, 3 + i * 3) += -delta.t * Eigen::Matrix3f::Identity();
            b.segment<3>(i * 6) += 0.5 * delta.t * delta.t * gravity + pose_i_w_i.unit_quaternion() *  delta.p + (pose_j_w_i.unit_quaternion() * Tbc.translation() - pose_i_w_i.unit_quaternion() * Tbc.translation());


            A.block<3, 2>(i * 6 + 3, 0) += -delta.t * Tg;
            A.block<3, 3>(i * 6 + 3, 3 + i * 3) += -Eigen::Matrix3f::Identity();
            
            A.block<3, 3>(i * 6 + 3, 3 + j * 3) += Eigen::Matrix3f::Identity();
            b.segment<3>(i * 6 + 3) += delta.t * gravity + pose_i_w_i.unit_quaternion() * delta.v;
        }

        x = A.fullPivHouseholderQr().solve(b);
        Eigen::Vector2f dg = x.segment<2>(0);
        gravity = (gravity + damp * Tg * dg).normalized() * PVIO_GRAVITY_NOMINAL;
    }

    scale = x(2);
    for (size_t i = 0; i < frames.size(); ++i) {
        velocities[i] = x.segment<3>(3 + i * 3);
    }
}




void ImuInitializer::solve_gravity_scale() {
    preintegrate();
    Eigen::Matrix4f A;
    Eigen::Vector4f b;
    A.setZero();
    b.setZero();

    for (size_t j = 1; j + 1 < frames.size(); ++j) {
        const size_t i = j - 1;
        const size_t k = j + 1;

        ORB_SLAM3::IMU::Delta &delta_ij = frames[j].mpImuPreintegratedFrame->delta;
        ORB_SLAM3::IMU::Delta &delta_jk = frames[k].mpImuPreintegratedFrame->delta;
        ORB_SLAM3::IMU::Jacobian &jacobian_ij = frames[j].mpImuPreintegratedFrame->jacobian;
        ORB_SLAM3::IMU::Jacobian &jacobian_jk = frames[k].mpImuPreintegratedFrame->jacobian;

        Sophus::SE3<float> pose_i = frames[i].GetImuPose();
        Sophus::SE3<float> pose_j = frames[j].GetImuPose();
        Sophus::SE3<float> pose_k = frames[k].GetImuPose();

        Eigen::Matrix<float, 3, 4> C;
        C.setZero();
        C.block<3, 3>(0, 0) = -0.5 * delta_ij.t * delta_jk.t * (delta_ij.t + delta_jk.t) * Eigen::Matrix3f::Identity();
        C.block<3, 1>(0, 3) = delta_ij.t * (pose_k.translation() - pose_j.translation()) - delta_jk.t * (pose_j.translation() - pose_i.translation());
        Eigen::Vector3f d = delta_ij.t * (pose_j.unit_quaternion() * delta_jk.p) + delta_ij.t * delta_jk.t * (pose_i.unit_quaternion() * delta_ij.v) - delta_jk.t * (pose_i.unit_quaternion() * delta_ij.p);
        A += C.transpose() * C;
        b += C.transpose() * d;
    }

    JacobiSVD<Eigen::Matrix4f> svd(A, ComputeFullU | ComputeFullV);
    Eigen::Vector4f x = svd.solve(b);
    gravity = x.segment<3>(0).normalized() * PVIO_GRAVITY_NOMINAL;
    scale = x(3);
}



void ImuInitializer::solve_scale_ba_given_gravity() {
    preintegrate();
    Eigen::Matrix4f A; 
    Eigen::Vector4f b; 
    A.setZero();
    b.setZero();

    for (size_t j = 1; j + 1 < frames.size(); ++j) {
        const size_t i = j - 1;
        const size_t k = j + 1;

        ORB_SLAM3::IMU::Delta &delta_ij = frames[j].mpImuPreintegratedFrame->delta;
        ORB_SLAM3::IMU::Delta &delta_jk = frames[k].mpImuPreintegratedFrame->delta;
        ORB_SLAM3::IMU::Jacobian &jacobian_ij = frames[j].mpImuPreintegratedFrame->jacobian;
        ORB_SLAM3::IMU::Jacobian &jacobian_jk = frames[k].mpImuPreintegratedFrame->jacobian;

        Sophus::SE3<float> pose_i = frames[i].GetImuPose();
        Sophus::SE3<float> pose_j = frames[j].GetImuPose();
        Sophus::SE3<float> pose_k = frames[k].GetImuPose();

        Eigen::Matrix<float, 3, 4> C;
        C.setZero();
        C.block<3, 1>(0, 0) = delta_ij.t * (pose_k.translation() - pose_j.translation()) - delta_jk.t * (pose_j.translation() - pose_i.translation());
        C.block<3, 3>(0, 1) = -(pose_j.unit_quaternion() * jacobian_jk.dp_dba * delta_ij.t + pose_i.unit_quaternion() * jacobian_ij.dv_dba * delta_ij.t * delta_jk.t - pose_i.unit_quaternion() * jacobian_ij.dp_dba * delta_jk.t);
        Eigen::Vector3f d = 0.5 * delta_ij.t * delta_jk.t * (delta_ij.t + delta_jk.t) * gravity + delta_ij.t * (pose_j.unit_quaternion() * delta_jk.p) + delta_ij.t * delta_jk.t * (pose_i.unit_quaternion() * delta_ij.v) - delta_jk.t * (pose_i.unit_quaternion() * delta_ij.p);
        A += C.transpose() * C;
        b += C.transpose() * d;
    }

    JacobiSVD<Eigen::Matrix4f> svd(A, ComputeFullU | ComputeFullV);
    Eigen::Vector4f x = svd.solve(b);
    scale = x(0);
    ba = x.segment<3>(1);
}




void ImuInitializer::refine_scale_ba_via_gravity() {
    static const double damp = 0.1;
    Eigen::Matrix<float, 6, 6> A;
    Eigen::Matrix<float, 6, 1> b;
    for (size_t iter = 0; iter < 2; ++iter) {
        preintegrate();
        A.setZero();
        b.setZero();
        Eigen::Matrix<float, 3, 2> Tg = s2_tangential_basis(gravity);

        for (size_t j = 1; j + 1 < frames.size(); ++j) {
            const size_t i = j - 1;
            const size_t k = j + 1;

            ORB_SLAM3::IMU::Delta &delta_ij = frames[j].mpImuPreintegratedFrame->delta;
            ORB_SLAM3::IMU::Delta &delta_jk = frames[k].mpImuPreintegratedFrame->delta;
            ORB_SLAM3::IMU::Jacobian &jacobian_ij = frames[j].mpImuPreintegratedFrame->jacobian;
            ORB_SLAM3::IMU::Jacobian &jacobian_jk = frames[k].mpImuPreintegratedFrame->jacobian;

            Sophus::SE3<float> pose_i = frames[i].GetImuPose();
            Sophus::SE3<float> pose_j = frames[j].GetImuPose();
            Sophus::SE3<float> pose_k = frames[k].GetImuPose();

            Eigen::Matrix<float, 3, 6> C;
            C.setZero();
            C.block<3, 1>(0, 0) = delta_ij.t * (pose_k.translation() - pose_j.translation()) - delta_jk.t * (pose_j.translation() - pose_i.translation());
            C.block<3, 3>(0, 1) = -(pose_j.unit_quaternion() * jacobian_jk.dp_dba * delta_ij.t + pose_i.unit_quaternion() * jacobian_ij.dv_dba * delta_ij.t * delta_jk.t - pose_i.unit_quaternion() * jacobian_ij.dp_dba * delta_jk.t);
            C.block<3, 2>(0, 4) = -0.5 * delta_ij.t * delta_jk.t * (delta_ij.t + delta_jk.t) * Tg;
            Eigen::Vector3f d = 0.5 * delta_ij.t * delta_jk.t * (delta_ij.t + delta_jk.t) * gravity + delta_ij.t * (pose_j.unit_quaternion() * delta_jk.p) + delta_ij.t * delta_jk.t * (pose_i.unit_quaternion() * delta_ij.v) - delta_jk.t * (pose_i.unit_quaternion() * delta_ij.p);
            A += C.transpose() * C;
            b += C.transpose() * d;

        }

        JacobiSVD<Eigen::Matrix<float, 6, 6>> svd(A, ComputeFullU | ComputeFullV);
        Eigen::Matrix<float, 6, 1> x = svd.solve(b);
        scale = x(0);
        ba += damp * x.segment<3>(1);
        gravity = (gravity + damp * Tg * x.segment<2>(4)).normalized() * PVIO_GRAVITY_NOMINAL;
    }
}



Eigen::Matrix<float, 3, 2> ImuInitializer::s2_tangential_basis(Eigen::Vector3f &x) {
    int d = 0;
    for (int i = 1; i < 3; ++i) {
        if (abs(x[i]) > abs(x[d])) d = i;
    }
    Eigen::Vector3f b1 = x.cross(Eigen::Vector3f::Unit((d + 1) % 3)).normalized();
    Eigen::Vector3f b2 = x.cross(b1).normalized();
    return (Eigen::Matrix<float, 3, 2>() << b1, b2).finished();
}




}
