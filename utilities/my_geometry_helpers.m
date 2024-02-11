source "./utilities/geometry_helpers_3d.m"



# This function, given a pose of two robots X_1 and X_2
# with a camera attached on both (in pose T_cam w.r.t. robot 
# pose), returns the pose of the camera attached on robot 2
# w.r.t. the camera attached on robot 1.
# Inputs:
#   X_1: pose of the robot 1 w.r.t. world
#   X_2: pose of the robot 2 w.r.t. world
#   T_cam: pose of the camera c1 w.r.t. robot r1, that is 
#          assumed to be also the pose of the camera c2
#          w.r.t. robot r2
# Outputs:
#   cam2_wrt_cam1: pose of the camera mounted on robot 2
#                  w.r.t. camera mounted on robot 1
function cam2_wrt_cam1 = pose2CameraRelativeTransform(X_1, X_2, T_cam)
    cam2_wrt_cam1 = inv(T_cam)*(inv(X_1)*X_2)*T_cam;
endfunction



# This function, given the pose of a robot 1 with a camera
# mounted on it (at pose T_cam w.r.t. its pose) and a 
# pose of a camera 2 mounted on a robot 2 (at pose T_cam
# w.r.t. pose of robot 2) expressed w.r.t. the camera 
# mounted on the robot 1, returns the pose of the robot 2
# w.r.t. world.
# Inputs:
#   X_1: pose of the robot 1 w.r.t. world
#   cam2_wrt_cam1: pose of the camera 2 w.r.t. camera 1
#   T_cam: pose of the camera c1 w.r.t. robot r1, that is 
#          assumed to be also the pose of the camera c2
#          w.r.t. robot r2
# Outputs:
#   X_2: pose of the robot 2 w.r.t. world
function X_2 = cameraRelativeTransform2Pose(X_1, cam2_wrt_cam1, T_cam)
    X_2 = X_1 * T_cam * cam2_wrt_cam1 * inv(T_cam);
endfunction



# Function that, given a pose of the robot w.r.t. world
# return the pose of the world, w.r.t. camera attached
# on it.
# Inputs:
#   X_rob: pose of the robot w.r.t. world
#   T_cam: pose of the camera w.r.t. robot
# Outputs:
#   wrld_wrt_cam: pose of the world w.r.t. camera
function wrld_wrt_cam = rw2wcam(X_rob, T_cam)
    wrld_wrt_cam = inv(X_rob*T_cam);
endfunction



# Function that, given a pose of the world w.r.t.
# camera, return the pose of the robot on which the 
# camera is attached w.r.t. world
# Inputs:
#   wrld_wrt_cam: pose of the world w.r.t. camera
#   T_cam: pose of the camera w.r.t. robot
# Outputs:
#   X_rob: pose of the robot w.r.t. world
function X_rob = wcam2rw(wrld_wrt_cam, T_cam)
    X_rob = inv(wrld_wrt_cam)*inv(T_cam);
endfunction



# Function that, given a set of points expressed in the camera
# frame of a robot, return the points expressed in world frame
# Inputs:
#   X_l = points expressed w.r.t. the camera
#   X_r = pose of the robot w.r.t. world
#   T_cam = pose of the camera w.r.t. the robot
# Outputs:
#   X_l_new = points expressed w.r.t. the world
function X_l_new = pointsCam2World(X_l, X_r, T_cam)
    
    # output vector initialization
    X_l_new = zeros(size(X_l));

    # rotation matrix and translation vector of X_r
    R_r = X_r(1:3, 1:3);
    t_r = X_r(1:3, 4);

    # rotation matrix and translation vector of T_cam
    R_cam = T_cam(1:3, 1:3);
    t_cam = T_cam(1:3, 4);

    # transformations to apply at each point
    R = R_r*R_cam;
    t = t_r+t_cam;

    # apply transformation at each point and fill the output vector
    for point_idx = 1:size(X_l, 2)
        X_l_new(:,point_idx) = R*X_l(:,point_idx) + t;
    endfor

endfunction



# Function that, given a pose of a robot expressed as a 4x4
# transformation matrix, returns the pose of the robot as
# a 3D vector (x, y, theta) assuming that the robot is
# on a plane with z = 0
# Inputs:
#   X_pose: 4x4 transformation matrix that represent the pose
#           of the robot w.r.t. the world
# Outputs:
#   plane_pose = 3D vector (x, y, theta)
function plane_pose = transform2planePose(X_pose)
    plane_pose = zeros(3,1);
    plane_pose(1:2) = X_pose(1:2,4);

    # because the rotational part of X_pose is assumed to be around z-axis:
    plane_pose(3) = atan2(X_pose(2,1), X_pose(1,1));    
endfunction



# Function that, given a set of poses of a robot expressed as a 4x4
# transformation matrix, returns the poses of the robot as
# a 3D vector (x, y, theta) assuming that the robot is
# on a plane with z = 0
# Inputs:
#   X_r: 4x4 transformation matrices that represent the poses
#           of the robot w.r.t. the world
# Outputs:
#   plane_poses = 3D vectors (x, y, theta)
function plane_poses = transforms2planePoses(X_r)

    plane_poses = [];
    for idx = 1:size(X_r, 3)
        plane_poses(:,end+1) = transform2planePose(X_r(:,:,idx));
    endfor
   
endfunction



# Function that, given a pose of the robot as a 3D vector
# (x, y, theta) assuming that the robot is on a plane with
# z = 0, returns the same robot pose expressed as a 4x4
# transformation matrix
# Inputs:
#   plane_pose = 3D vector (x, y, theta)
# Outputs:
#   X_pose: 4x4 transformation matrix that represent the pose
#           of the robot w.r.t. the world
function X_pose = planePose2Transform(plane_pose)
    
    # initialize output
    X_pose = eye(4);
    
    # get translation vector
    t = [plane_pose(1); plane_pose(2); 0];

    # get rotation matrix
    R = Rz(plane_pose(3));

    # ensemble the output
    X_pose(1:3, 4) = t;
    X_pose(1:3, 1:3) = R;

endfunction



# Function that, given a set of 4x4 transformation matrices
# representing robot poses w.r.t. world, returns a set of 4x4 
# transformation matrices representing camera poses w.r.t. world.
# Inputs:
#   X_r: set of robot poses w.r.t. world
#   T_cam: pose of the camera w.r.t. robot
# Output:
#   X_c: set of camera poses w.r.t. world
function X_c = robotPoses2CamPoses(X_r, T_cam)

    # get the number of poses
    num_poses = size(X_r, 3);

    # initialize the output matrix
    X_c = zeros(4, 4, num_poses);

    # for each pose
    for pose_idx = 1:num_poses

        # compute the camera pose attached on the robot in the current pose
        current_camera_pose = X_r(:, :, pose_idx) * T_cam;

        # fill the output matrix
        X_c(:,:,pose_idx) = current_camera_pose;

    endfor

endfunction



# Inverse of the function before
function X_r = camPoses2RobotPoses(X_c, T_cam)

    # get the number of poses
    num_poses = size(X_c, 3);

    # initialize the output matrix
    X_r = zeros(4, 4, num_poses);

    # for each pose
    for pose_idx = 1:num_poses

        # compute the camera pose attached on the robot in the current pose
        current_robot_pose = X_c(:, :, pose_idx) * inv(T_cam);

        # fill the output matrix
        X_r(:,:,pose_idx) = current_robot_pose;

    endfor

endfunction
