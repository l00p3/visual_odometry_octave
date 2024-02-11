source "visual_odometry.m"
source "./utilities/loading_helpers.m"
source "./utilities/visualization_helpers.m"
source "./utilities/measurement_and_state_helpers.m"
source "./my_functions/my_main_functions.m"



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



##### PERFORM VISUAL ODOMETRY ##################################################
[predicted_landmark_poses, predicted_trajectory, landmarks_error, trajectory_error] = performVisualOdometry(data_folder, camera_data_filename, world_data_filename, trajectory_data_filename, ...
                                            scale_factor, scale_factor_triangulation, ...
                                            num_iterations_posit, num_iterations_final_refinement, kernel_threshold_posit, ...
                                            threshold_to_ignore, kernel_threshold_final_refinement, damping);



##### GET GROUND TRUTH #########################################################
trajectory = loadTrajectory(trajectory_data_filename);
ground_truth_trajectory = getTrajectory(trajectory);
[landmark_positions, _] = loadWorld(world_data_filename);



##### PLOT #####################################################################

"RESULTS:"

# VISUALIZE ERRORS
landmarks_error
trajectory_error

# PLOT LANDMARKS AND POSES
figure(1);
title("Landmarks");
visualize_points(landmark_positions, "real landmarks", predicted_landmark_poses, "predicted landmarks", "Landmarks");
figure(2);
title("Poses");
visualize_poses(ground_truth_trajectory, "real poses", predicted_trajectory, "predicted poses", "Poses");

waitfor(gcf);
