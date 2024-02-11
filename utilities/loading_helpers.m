# In this file there are all the utilities to load the 
# informations contained in the files in the data folder


# Load the camera parameters from the file filename in the given format
# Inputs:
#   filename: file containing all the camera parameters:
# Outputs:
#   K: camera matrix
#   T_cam: pose of the camera w.r.t. robot
#   z_near: min distance of the space where the camera can perceive stuff
#   z_far: max distance of the space where the camera can perceive stuff
#   width: width of the camera projective plane
#   height: height of the camera projective plane
function [K, T_cam, z_near, z_far, cam_width, cam_height] = loadCameraParams(filename)
    # load raw data
    data_cam = dlmread(filename);

    # camera matrix
    K = data_cam(2:4, 1:3);

    # pose of the camera w.r.t. the world
    T_cam = data_cam(6:9, 1:4);

    # camera frustum parameters
    z_near = data_cam(10, 2);
    z_far = data_cam(11, 2);
    cam_width = data_cam(12, 2);
    cam_height = data_cam(13, 2);
endfunction



# Load the landmarks gorundtruth
# Inputs:
#   filename: file containing all the world info in the format:
#               LANDMARK_ID POSITION APPEARANCE
#             for each row
# Outputs:
#   landmark_positions: vector containing for each column a landmark
#                       position, in particular in position i there is
#                       the landmark with id i-1
#   landmark_appearances: vector containing for each column a landmark
#                         appearance, in particular in position i there
#                         is the landmark with id i-1
function [landmark_positions, landmark_appearances] = loadWorld(filename)
    # load raw data
    data_wrld = dlmread(filename);
    num_landmark = size(data_wrld, 1);

    # initialize empty vectors
    landmark_positions = zeros(3, num_landmark);
    landmark_appearances = zeros(10, num_landmark);

    # fill vectors
    for(data_idx = 1:size(data_wrld, 1))
        id_landmark = data_wrld(data_idx, 1);
        landmark_positions(:,id_landmark+1) = data_wrld(data_idx, 2:4)';
        landmark_appearances(:,id_landmark+1) = data_wrld(data_idx, 5:14)';
    endfor

endfunction



# Load a single measurement, that is a list of image points seen from the actual
# position of the camera in the trajectory
# Inputs:
#   filename: file containing a single measurement. Every measurement contains a 
#             sequence number, ground truth (of the robot) and odometry pose 
#             and measurement information:
#               point POINT_ID_CURRENT_MESUREMENT ACTUAL_POINT_ID IMAGE_POINT APPEARANCE
# Outputs:
#   measurement : a struct with these fields:
#                       - id_pose: id of the current pose
#                       - ground_truth_pose: the real pose of the camera that makes 
#                                               these observations
#                       - odometry_pose: the read pose of the camera that makers 
#                                               these observations
#                       - measurements: a vector containing, in each column, a 2D 
#                                                image point measured
#                       - appearances: a vector containing, in each column, a 10-dim 
#                                                appearances vector of the measured point
#                       - landmarks_ids: a vector of ids, each id in the position j 
#                                                corresponds to the id of the landmark 
#                                                measured in the position j of the vector 
#                                                measurements (same for appearances)
function measurement = loadSingleMeasurement(filename)

    # load the id of the pose
    data_meas = dlmread(filename, '', [0, 0, 0, 1]);
    id_pose = data_meas(1, 2);

    # load the ground truth pose and the odometry pose
    # here we use a different approach because of an error
    # in reading the second line of the file
    data_meas = fileread(filename);
    data_meas = strsplit(data_meas, "\n");
    ground_truth_pose = str2num(data_meas{2}(9:end))';
    odometry_pose = str2num(data_meas{3}(11:end))';

    # load remaining raw data (measurements)
    data_meas = dlmread(filename, '', 3, 0);
    num_measured_landmark = size(data_meas, 1); 

    # initialize empty vectors
    measurements = zeros(2, num_measured_landmark);
    appearances = zeros(10, num_measured_landmark);
    landmarks_ids = zeros(1,num_measured_landmark);

    # fill vectors
    for(data_idx = 1:size(data_meas, 1))
        id_meas = data_meas(data_idx, 2);    # id of the actual point in the measurement
        id_measured_landmark = data_meas(data_idx, 3);
        measurements(:,id_meas+1) = data_meas(data_idx, 4:5)';
        appearances(:,id_meas+1) = data_meas(data_idx, 6:15)';
        landmarks_ids(id_meas+1) = id_measured_landmark;
    endfor

    # create the data structure
    measurement = struct("id_pose", id_pose,
                        "ground_truth_pose", ground_truth_pose,
                        "odometry_pose", odometry_pose,
                        "measurements", measurements,
                        "appearances", appearances,
                        "landmarks_ids", landmarks_ids);

endfunction



# Load all the measurements from the given folder
# Inputs:
#   data_folder: folder containing a list of files in the format "meas-#####.dat"
#                where each file contains a set of measured landmarks
# Outputs:
#   measurements : a list of measurements, each measurement is a struct with these
#                  fields:
#                       - id_pose: id of the current pose
#                       - ground_truth_pose: the real pose of the camera that makes 
#                                               these observations
#                       - odometry_pose: the read pose of the camera that makers 
#                                               these observations
#                       - measurements: a vector containing, in each column, a 2D 
#                                                image point measured
#                       - appearances: a vector containing, in each column, a 10-dim 
#                                                appearances vector of the measured point
#                       - landmarks_ids: a vector of ids, each id in the position j 
#                                                corresponds to the id of the landmark 
#                                                measured in the position j of the vector 
#                                                measurements (same for appearances)
function measurements = loadAllMeasurements(data_folder)

    # get the list of files in the folder
    files_list = ls(data_folder);

    # initialize the output list
    measurements = {};

    # load the measurements
    for file_idx = 1:size(files_list, 1)

        if(strcmp(files_list(file_idx, 1:4), "meas"))

            # load the corresponding measurement
            actual_measurement = loadSingleMeasurement(strcat(data_folder, "/", files_list(file_idx, :)));

            # store measurement
            measurements(end+1) = actual_measurement;

        endif

    endfor

endfunction



# Load the trajectory
# Inputs:
#   filename: a file containing info about the trajectory with the format:
#                   pose: POSE_ID ODOMETRY_POSE GROUND_TRUTH_POSE
#             the ODOMETRY_POSE is obtained by adding Gaussian Noise (0; 0.001) 
#             to the actual robot commands
# Outputs:
#   trajectory: a list of structs, each of them containing the following fields:
#                   - pose_id: id of the pose
#                   - odometry_pose: read pose
#                   - ground_truth_pose: actual pose 
function trajectory = loadTrajectory(filename)

    # Load raw data
    data_trj = fileread(filename);
    data_trj = strsplit(data_trj, "\n");

    # initialie the output list
    trajectory = {};

    for element_idx = 1:size(data_trj, 2)
        # take the current element from the raw data
        current_row = strsplit(data_trj{element_idx}, " ");
        
        # retrieve useful informations
        pose_id = str2num(current_row{2});
        odometry_pose = [str2num(current_row{3}); str2num(current_row{4}); str2num(current_row{5})];
        ground_truth_pose = [str2num(current_row{6}); str2num(current_row{7}); str2num(current_row{8})];

        # add the new element to the trajectory
        trajectory(end+1) = struct("pose_id", pose_id,
                                   "odometry_pose", odometry_pose,
                                   "ground_truth_pose", ground_truth_pose);
    endfor

endfunction