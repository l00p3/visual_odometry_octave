# This function, given a pose id, return the corresponding pose
# in the state vector (the id goeas from 1 to number of poses)
# Inputs:
#   X_r: poses state
#   id: pose id
# Outputs:
#   retrieved_pose: the pose with id retrieved from the state
function retrieved_pose = getPoseFromId(X_r, id)
    retrieved_pose = X_r(:,:,id);
endfunction




# Function that, given an id of tha landmark in the state
# (to not confuse with the real id of the landmark, this is
# just the position in the state), return the corresponding 
# landmark from the state.
# Inputs:
#   X_l: landmarks state
#   id: position of the landmark in the state
# Outputs:
#   retrieved_landmark: the pose of the landmark in position id
function retrieved_landmark = getLandmarkFromStateId(X_l, id)
    retrieved_landmark = X_l(:, id);
endfunction



# This function takes two measurements and return two set of points
# that contain the common points seen by the first measurement and
# from the second, respectively and two set of points containing all
# the points that are not in common.
# Inputs:
#   first_measurement
#   second_measurement
# Outputs:
#   common_points1: set of points containing all the projective points
#                   seen from the first measurement
#   common_points2: same set of points contained in common_points1, also
#                   ordered in the same way, but observed from the second
#                   measurements
#   common_ids: list of actual landmark ids, in the same order of common_points vectors
#   uncommon_points: set of points containing all the projective points seen from
#                    the second measurement but not seen in the first
#   uncommon_ids: list of landmark ids of the uncommon_points
function [common_points1, common_points2, common_ids, uncommon_points, uncommon_ids] = getIntersectingMeasurements(measure1, measure2)

    # list of ids of landmarks from the two measurements
    ids1 = measure1.landmarks_ids;
    ids2 = measure2.landmarks_ids; 

    # initialize them to empty
    common_points1 = [];
    common_points2 = [];
    common_ids = [];
    uncommon_points = [];
    uncommon_ids = [];
    
    # for each point id in the first measurement
    for ids1_element_idx = 1:size(ids1,2)

        # check if the current point is seen also in the second measurement
        current_id = ids1(ids1_element_idx);
        ids2_element_idx = find(ids2 == ids1(ids1_element_idx));

        # if it is visible in both measurements add them
        if(!isempty(ids2_element_idx))
            # else add it as common point
            common_points1(:,end+1) = measure1.measurements(:,ids1_element_idx);
            common_points2(:,end+1) = measure2.measurements(:,ids2_element_idx);
            common_ids(end+1) = current_id;

            # keep track that this point is in common
            ids2(ids2_element_idx) = -1;
        endif

    endfor

    # for each point id in the second measurement
    for ids2_element_idx = 1:size(ids2,2)

        # if this point is not in common
        if(ids2(ids2_element_idx) != -1)

            # memorize it as uncommon
            uncommon_points(:,end+1) = measure2.measurements(:,ids2_element_idx);
            uncommon_ids(end+1) = ids2(ids2_element_idx);

        endif

    endfor

endfunction



# This function, given a set of common points (seen by two camera)
# and the associations vector return only the common points that are 
# not in the associations vector
# Inputs:
#   common_points1: set of common points seen from the camera 1
#   common_points2: set of common points seen from the camera 2
#   common_ids: list of the ids of the common points
#   associations_inv: a map between landmark id in the state vector
#                     to actual landmark id
# Outputs:
#   new_common_points1: set of only new common points seen from cam 1
#   new_common_points2: set of only new common points seen from cam 2 
#   new_common_ids: set of only new common points ids
function [new_common_points1, new_common_points2, new_common_ids] = getNewCommons(common_points1, common_points2, common_ids, associations_inv)

    # initialize output vector
    new_common_points1 = [];
    new_common_points2 = [];
    new_common_ids = [];

    # for each landmark id seen
    for idx = 1:size(common_ids, 2)

        # get current landmark id
        current_landmark_id = common_ids(idx);

        # check if the current landmark id is inside the associations vector
        res = find(associations_inv == current_landmark_id+1);

        # if the current seen landmark is not yet in the associations
        # vector is new so add it to the output vectors
        if(isempty(res))

            # add new measurement to the output vectors
            new_common_points1(:,end+1) = common_points1(:,idx);
            new_common_points2(:,end+1) = common_points2(:,idx);
            new_common_ids(:,end+1) = current_landmark_id;

        endif

    endfor

endfunction



# This function computes, from all the measurements and the 
# set of known landmarks X_l, two vector that represent the 
# positions of the measured landmarks in the world (where we
# think they are) and the measurements, ordered in the same way
# Inputs:
#   measurements
#   X_l: set of landmark poses (where we think they are located)
#   associations: a map between actual landmark id to landmark id 
#                 in the state vector X_l
function [X_l_new, z] = getPosesAndMeasurements(measurement, X_l, associations)

    # initialize output vectors
    X_l_new = [];
    z = [];

    # get the measured points list
    measured_points = measurement.measurements;

    # get the list of the ids of the measured landmarks
    measurement_ids = measurement.landmarks_ids;

    # for each id in measurement_ids
    for idx = 1:size(measurement_ids, 2)

        # get the actual id
        actual_landmark_id = measurement_ids(idx);

        # check if we have already seen this landmark
        state_landmark_id = associations(actual_landmark_id+1);

        # if it is not -1 we have seen again this landmark before
        # so we have a pose in the state vector X_l
        if(state_landmark_id != -1)

            # add it to the outputs vectors
            X_l_new(:,end+1) = X_l(:, state_landmark_id);
            z(:,end+1) = measured_points(:,idx);

        endif

    endfor

endfunction



# Function that compute the projection measurements in a format
# that is compatible with the final refinement function
# Inputs:
#   measurements: all the measurements 
#   associations: a map between actual landmark id to landmark id 
#                 in the state vector X_l
# Outputs:
#   Zp: the projective measurements (2xnum_measurements)
#   projection_associations: 2xnum_measurements. 
#                 projection_associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
function [Zp, projection_associations] = formatProjectiveMeasurement(measurements, associations)

    # initialize to empty the output vectors
    projection_associations = [];
    Zp = [];

    # compute the number of measurements in the measurement vector
    num_measurements = size(measurements, 2);

    for meas_idx = 1:num_measurements

        # compute the pose that makes this measurement
        pose_idx = meas_idx;

        # get the current measurement
        current_measurement = measurements{meas_idx};

        # compute the number of measured landmark in the current measurement
        num_measured_landmark = size(current_measurement.measurements, 2);

        # for each measured landmark
        for land_meas_idx = 1:num_measured_landmark

            # get the position in the X_l vector of the measured landmark
            landmark_idx = associations(current_measurement.landmarks_ids(land_meas_idx)+1);

            # check if we have a registered prediction for this landmark
            # if not just ignore it
            if(landmark_idx == -1)
                continue;
            endif

            # get the measured position of this landmark and add it to Zp
            Zp(:,end+1) = current_measurement.measurements(:,land_meas_idx);

            # add the new element to the associations vector
            projection_associations(:,end+1) = [pose_idx; landmark_idx];

        endfor

    endfor

endfunction



# Function that compute the pose measurements in a format
# that is compatible with the final refinement function 
# Inputs:
#   X_r: list of camera poses
# Outputs:
#   Zr: the poses measurements, where from each pose the next one is observed
#   pose_associations: 2xnum_measurements. 
#                 pose_associations(:,k)=[i_idx, j_idx]' means the kth measurement
#                 refers to an observation made from pose i_idx, that
#                 observed the pose j_idx
function [Zr, pose_associations] = formatPoseMeasurements(X_r)

    # compute the number of poses
    num_poses = size(X_r, 3);

    # compute the number of measurements
    num_pose_measurements = num_poses-1;

    # initialize the outputs
    Zr = zeros(4, 4, num_pose_measurements);
    pose_associations = zeros(2, num_pose_measurements);

    # generate measurements
    measurement_num=1;
    for (pose_num=1:num_poses-1)

        # observer pose
        Xi = X_r(:, :, pose_num);

        # observed pose
        Xj = X_r(:, :, pose_num+1);

        # generate a new element for the association vector
        pose_associations(:,measurement_num) = [pose_num, pose_num+1]';

        # how the pose Xj is seen from the Xi
        # the pose of Xj w.r.t. Xi
        Zr(:, :, measurement_num) = inv(Xi)*Xj;

        measurement_num++;
    endfor

endfunction



# Function that delete all the pending landmarks using a list
# as referement
# Inputs:
#   pending_landmarks: this vector will contain a list of:
#                           [landmark_id, pose_id, coord_x, coord_y]'
#                      to memorize that in the measurement maked in pose_id pose we have
#                      found the landmark_id for the first time, in this way, when we'll encounter
#                      again such landmark in a measurement we have to triangulate it
#                      between the pose pose_id and the new pose
#                      and add the result in the state, eliminating the corresponding
#                      element in this vector
#                      coord_x and coord_y represent the two measured value for that landmark
#                      from the pose pose_id
#   no_more_pending_ids: list of landmark ids that are not pending anymore
# Outputs:
#   clean_pending_landmarks: the cleaned version of the vector pending_landmarks
function clean_pending_landmarks = cleanPendingLandmarks(pending_landmarks, no_more_pending_ids)

    # initialize output vector
    clean_pending_landmarks = [];

    # for each pending_landmark element
    for idx = 1:size(pending_landmarks, 2)

        # get the id of the current pending landmark
        currend_landmark_id = pending_landmarks(1,idx);

        # if this id is not in the no_more_pending_ids list add 
        # the corresponding pending element to the output vector
        res = find(no_more_pending_ids == currend_landmark_id);
        if(isempty(res))
            clean_pending_landmarks(:,end+1) = pending_landmarks(:,idx);
        endif

    endfor

endfunction



# Function that, given a list of landmark ids makes two things:
#   - check if one of this landmark ids correspond to a pending
#     landmark and in that case add it to the vectors common_pendingi
#     and remove it from the pending landmarks
#   - if not consider it as a new pending landmark
# Inputs:
#   pending_landmarks: this vector will contain a list of:
#                           [landmark_id, pose_id, coord_x, coord_y]'
#                      to memorize that in the measurement maked in pose_id pose we have
#                      found the landmark_id for the first time, in this way, when we'll encounter
#                      again such landmark in a measurement we have to triangulate it
#                      between the pose pose_id and the new pose
#                      and add the result in the state, eliminating the corresponding
#                      element in this vector
#                      coord_x and coord_y represent the two measured value for that landmark
#                      from the pose pose_id
#   landmark_measures: a list of landmark 2D points measuring the landmarks corresponding to the 
#                      ids in the landmark_ids vector
#   landmark_ids: a list of landmark ids
#   pose_id: pose from which the measurements contained in landmark_measures are taken
#   associations_inv: vector that, given a landmark id in the state vector
#                     return the corresponding actual landmark id
# Output:
#   new_pending_landmarks: the updated version of pending_landmarks
#   common_pending1: a vector containing a list of 2D measured points that was pending
#   common_pending2: a vector containing a list of 2D measured points, same points in 
#                    common_points1 but estracted from landmark_measures vector
#   common_pending_ids: a vector with the elements:
#                           [landmark_id, pose_id]'
#                       where in position i we have the landmark id of the point measured
#                       in common_pendingi and the pose_id of the pose from which the point
#                       in position i in common_pending1 is measured
function [new_pending_landmarks, common_pending1, common_pending2, common_pending_ids] = updatePendingLandmarks(pending_landmarks, ...
                                                                                                    landmark_measures, landmark_ids,
                                                                                                    pose_id, associations_inv)

    # initialize output vectors
    new_pending_landmarks = [];
    common_pending1 = [];
    common_pending2 = [];
    common_pending_ids = [];

    # for each pending_landmark element
    for idx = 1:size(pending_landmarks, 2)

        # get the id of the current pending landmark
        currend_landmark_id = pending_landmarks(1,idx);  

        # if this id is not in the landmark_ids list add 
        # the corresponding pending element to the vector
        # new_pending_landmarks
        res = find(landmark_ids == currend_landmark_id);
        if(isempty(res))
            new_pending_landmarks(:,end+1) = pending_landmarks(:,idx);
        
        # else it is no more pending and so add it to the common_pendingi 
        # vectors
        else
            common_pending1(:,end+1) = pending_landmarks(3:4,idx);
            common_pending2(:,end+1) = landmark_measures(:,res);
            common_pending_ids(:,end+1) = pending_landmarks(1:2,idx);

            # remember the already used ids
            landmark_ids(res) = -1;
        endif

    endfor

    # now consider all the not used landmark ids in the landmark_ids vector
    # as new pending landmarks, only if they are not in the associations_inv 
    # vector
    for idx = 1:size(landmark_ids,2)
        l_id = landmark_ids(idx);
        res = find(associations_inv == l_id+1);
        if(l_id != -1 && isempty(res))
            new_pending_landmarks(:,end+1) = [l_id; pose_id; landmark_measures(:,idx)];
        endif
    endfor

endfunction



# Function that deal with the first measurement in order to add
# all the landmarks that are not seen also by the second measurement
# in the pending landmarks vector
# Inputs:
#   pending_landmarks: this vector will contain a list of:
#                           [landmark_id, pose_id, coord_x, coord_y]'
#                      to memorize that in the measurement maked in pose_id pose we have
#                      found the landmark_id for the first time, in this way, when we'll encounter
#                      again such landmark in a measurement we have to triangulate it
#                      between the pose pose_id and the new pose
#                      and add the result in the state, eliminating the corresponding
#                      element in this vector
#                      coord_x and coord_y represent the two measured value for that landmark
#                      from the pose pose_id
#   first_measurement: the first measurement
#   associations_inv: vector that, given a landmark id in the state vector
#                     return the corresponding actual landmark id
# Outputs:
#   pending_landmarks: the updated version of pending_landmarks
function pending_landmarks = getPendingFromFirstMeasurement(pending_landmarks, first_measurement, associations_inv)

    # get all the ids seen from the first measurement
    first_measurement_landmark_ids = first_measurement.landmarks_ids;

    # get all the pending ids
    pending_landmarks_ids = [];
    for idx = 1:size(pending_landmarks, 2)
        pending_landmarks_ids(end+1) = pending_landmarks(1, idx);
    endfor

    # for each id
    for idx = 1:size(first_measurement_landmark_ids, 2)

        # if it is not present in the pending_landmark vector already
        # and we have not yet triangulated it
        # then add it to as a pending_landmark
        current_landmark_id = first_measurement_landmark_ids(idx);
        res1 = find(pending_landmarks_ids == current_landmark_id);
        res2 = find(associations_inv == current_landmark_id+1);
        if(isempty(res1) && isempty(res2))
            pending_landmarks(:,end+1) = [current_landmark_id; 1; first_measurement.measurements(:,idx)];
        endif

    endfor

endfunction



# Function that simply extract the ground truth trajectory
# from the loaded struct
function ground_truth_trajectory = getTrajectory(trajectory)
    ground_truth_trajectory = [];
    for idx = 1:size(trajectory, 2)
        ground_truth_trajectory(:,end+1) = trajectory(idx){1}.ground_truth_pose;
    endfor
endfunction