source "./utilities/loading_helpers.m"
source "./utilities/visualization_helpers.m"
source "./utilities/measurement_and_state_helpers.m"
source "./utilities/my_geometry_helpers.m"
source "./my_functions/my_main_functions.m"
source "./my_functions/my_epipolar_functions.m"
source "./my_functions/my_projective_icp_functions.m"
source "./my_functions/my_total_ls_functions.m"



function [X_l, predicted_trajectory, landmarks_error, trajectory_error] = performVisualOdometry(data_folder, camera_data_filename, world_data_filename, trajectory_data_filename, ...
                                            scale_factor, scale_factor_triangulation, ...
                                            num_iterations_posit, num_iterations_final_refinement, kernel_threshold_posit, ...
                                            threshold_to_ignore, kernel_threshold_final_refinement, damping);


    # given a landmark id i this vector in position i+1 will have the position
    # in the state where we can find the prediction of the landmark i
    # we don't know how many landmark do we have so we set it's dimension
    # "big enough" and we initialize it with -1 in each position, where
    # -1 means that we have no prediction in the state for the landmark i
    associations = zeros(1,2000)-1;  

    # given a landmark in position j in the state this vector in position j will have
    # the id i+1 of this landmark
    # this vector will grow when we meet a new landmark 
    associations_inv = [];  

    # this vector will contain a list of:
    #   [landmark_id, pose_id, coord_x, coord_y]'
    # to memorize that in the measurement maked in pose_id pose we have
    # found the landmark_id for the first time, in this way, when we'll encounter
    # again such landmark in a measurement we have to triangulate it
    # between the pose pose_id and the new pose
    # and add the result in the state, eliminating the corresponding
    # element in this vector
    # coord_x and coord_y represent the two measured value for that landmark
    # from the pose pose_id
    pending_landmarks = [];



    ##### LOAD DATA ################################################################
    [K, T_cam, z_near, z_far, cam_width, cam_height] = loadCameraParams(camera_data_filename);
    [landmark_positions, landmark_appearances] = loadWorld(world_data_filename);
    measurements = loadAllMeasurements(data_folder);
    trajectory = loadTrajectory(trajectory_data_filename);



    ##### STATE INITIALIZATION ####################################################
    # Our state is composed by two components, one component contains all the
    # camera poses (trajectory) and the other component contains all the poses of the 
    # observed landmarks; all expressed in the world frame.
    # In order to start to fill the second component we need
    # first to triangulate some measurement from the first two positions.
    # So, assuming that the camera in the initial timestep is in the origin
    # we can compute, taking the first two measurements, through the 8-point
    # algorithm, the position of the second camera w.r.t. the first one
    # and so the position of the second pose w.r.t. the world frame.
    # Then, using the computed poses we can triangulate the measurement to 
    # obtain an initial constellation for our landmarks, that is the landmark
    # predicted poses seen from both the camera.

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



    ##### EXPLORE TRAJECTORY ######################################################
    # Now it's time to explore all other poses in the trajectory. Here basically
    # we take each remaining measurement and, from each of them, we extract first
    # a new pose with Projective ICP and add it to the X_r part of the state and then
    # from this new pose and the previous one, we extract all the landmark seen from 
    # both that we have not added to the state before (because we observe them now
    # for the first time with two cameras) through triangulation.

    "START EXPLORING TRAJECTORY"

    # For each measurement, starting from the 3th
    for meas_idx = 3:size(measurements, 2)

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



    ##### FINAL REFINEMENT #########################################################
    # Here we basically take all the camera poses and all the measurements and do
    # a final refinement considering them as a factor graph.

    "START FINAL REFINEMENT"

    # Get the camera poses on the predicted trajectory
    X_c = robotPoses2CamPoses(X_r, T_cam);

    # Compute the projection association vector and Zp
    [Zp, projection_associations] = formatProjectiveMeasurement(measurements, associations);

    # Compute the pose association vector and Zr
    [Zr, pose_associations] = formatPoseMeasurements(X_c);

    # Execute the final refinement
    [X_c, X_l, chi_stats_p, num_inliers_p, chi_stats_r, num_inliers_r, H, b] = doFinalOptimization(X_c, X_l,
            Zp, projection_associations,
            Zr, pose_associations,
            num_iterations_final_refinement,
            damping,
            kernel_threshold_final_refinement,
            K, cam_width, cam_height, threshold_to_ignore);

    # Compute predicted trajectory from camera trajectory
    predicted_trajectory = camPoses2RobotPoses(X_c, T_cam);

    # Compute predicted trajectory as 3D vectors
    predicted_trajectory = transforms2planePoses(predicted_trajectory);

    # Compute errors
    landmarks_error = computeLandmarksPredictionError(landmark_positions, X_l, associations_inv);
    ground_truth_trajectory = getTrajectory(trajectory);
    trajectory_error = computeTrajectoryPredictionError(ground_truth_trajectory, predicted_trajectory);

endfunction
