source "./utilities/loading_helpers.m"
source "./utilities/visualization_helpers.m"
source "./utilities/measurement_and_state_helpers.m"
source "./utilities/my_geometry_helpers.m"
source "./my_functions/my_main_functions.m"
source "./my_functions/my_epipolar_functions.m"
source "./my_functions/my_projective_icp_functions.m"
source "./my_functions/my_total_ls_functions.m"



##### INITIALIZE PARAMETERS ####################################################
data_folder = "data/";
camera_data_filename = strcat(data_folder, "camera.dat");
world_data_filename = strcat(data_folder, "world.dat");
trajectory_data_filename = strcat(data_folder, "trajectory.dat");

scale_factor = 1/18;        # scale factor to apply to 8-point algorithm
scale_factor_triangulation = 1;
num_iterations_posit = 100;
num_iterations_final_refinement = 50;
kernel_threshold_posit = 10;
threshold_to_ignore = 2000;     # error threshold that determine if an outlier is too outlier to be considered
kernel_threshold_final_refinement = 10;
damping = 1;

associations = zeros(1,2000)-1; # landmark id to state position
associations_inv = [];          # state position to landmark id
pending_landmarks = [];         # [landmark_id, pose_id, coord_x, coord_y]  # redundant!


##### LOAD DATA ################################################################
[K, T_cam, z_near, z_far, cam_width, cam_height] = loadCameraParams(camera_data_filename);
[landmark_positions, _] = loadWorld(world_data_filename);
measurements = loadAllMeasurements(data_folder);
trajectory = loadTrajectory(trajectory_data_filename);



##### STATE INITIALIZATION ####################################################
"START STATE INITIALIZATION"

# take the first two measurements
first_measurement = measurements{1};
second_measurement = measurements{2};

# Compute the 8-point algorithm to get the first two robot poses
X_r = initilizePoses(first_measurement, second_measurement, K, T_cam, scale_factor);

# Initialize the X_l state vector to empty
X_l = [];

# Compute the first landmark positions, that are the landmarks seen from both
# the first and the second measurement
[X_l, associations, associations_inv, pending_landmarks] = augmentConstellation(X_r, first_measurement, second_measurement, getPoseFromId(X_r, 1), getPoseFromId(X_r, 2), K, T_cam, X_l, associations, associations_inv, pending_landmarks, 2, scale_factor_triangulation);

# Get all the landmarks seen from the first measure that are not seen from the second
# and add them to the pending vector
pending_landmarks = getPendingFromFirstMeasurement(pending_landmarks, first_measurement, associations_inv);



#### ODOMETRY #################################################################
"START EXPLORING TRAJECTORY"

# For each measurement, starting from the 3th
for meas_idx = 3:30#size(measurements, 2)

    strcat("PREDICT POSE: ", int2str(meas_idx), "/", int2str(size(measurements, 2)))
    
    # Compute the current pose given the measurements
    [X_guess, chi_stats, num_inliers] = doProjectiveICP(measurements{meas_idx}, X_l, associations, ...
                                                                num_iterations_posit, kernel_threshold_posit, damping, ...
                                                                K, T_cam, z_near, z_far, cam_width, cam_height, threshold_to_ignore);

    # Add the new pose to the state
    X_r(:,:,meas_idx) = X_guess; 
    
    # Triangulate new points given the new measurement and the new pose
    [X_l, associations, associations_inv, pending_landmarks] = augmentConstellation(X_r, measurements{meas_idx-1}, measurements{meas_idx}, getPoseFromId(X_r, meas_idx-1), getPoseFromId(X_r, meas_idx), K, T_cam, X_l, associations, associations_inv, pending_landmarks, meas_idx, scale_factor_triangulation);

endfor



#### EXTRACT TRAJECTORIES ####
ground_truth_trajectory = getTrajectory(trajectory);
predicted_trajectory = transforms2planePoses(X_r);


#### VISUALIZE ERRORS ####
trajectory_error = computeTrajectoryPredictionError(ground_truth_trajectory, predicted_trajectory);
landmarks_error = computeLandmarksPredictionError(landmark_positions, X_l, associations_inv);
trajectory_error
landmarks_error


#### VISUALIZE TRAJECTORIES ####
figure(1)
plot3(ground_truth_trajectory(1,:), ground_truth_trajectory(2,:), zeros(1,size(ground_truth_trajectory, 2)));
hold on;
plot3(predicted_trajectory(1,:), predicted_trajectory(2,:), zeros(1,size(predicted_trajectory, 2)));
legend("ground truth trajectory", "predicted trajectory");
xlabel('x'); ylabel('y'); zlabel('z');
title("Trajectory");
grid;

#### VISUALIZE LANDMARKS ####
figure(2);
plot3(landmark_positions(1,:), landmark_positions(2,:), landmark_positions(3,:), 'b*', "linewidth", 2);
hold on;
plot3(X_l(1,:), X_l(2,:), X_l(3,:), 'ro', "linewidth", 2);
legend("ground truth landmarks", "predicted landmarks");
xlabel('x'); ylabel('y'); zlabel('z');
title("Landmarks");
grid;
waitfor(gcf);

