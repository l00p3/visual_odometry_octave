# This function, given the first two measurement return an initial guess
# for the first two poses, using the eight-point algorithm.
# Inputs:
#   measurement1: measurement obtained by the camera of the first pose
#   measurement2: measurement obtained by the camera of the second pose
#   K: camera matrix
#   T_cam: pose of the camera w.r.t. the robot pose
#   scale_factor: scale factor for the computed pose, needed because the pose
#                   predicted with the 8-point algorithm is up to scale
function X_r = initilizePoses(measurement1, measurement2, K, T_cam, scale_factor)

    # first pose is assumed to be at the origin
    first_pose = eye(4);

    # take the common seen landmarks, ordered in the same way
    [common_points1, common_points2, _] = getIntersectingMeasurements(measurement1, measurement2);

    # compute the 8-point algorithm
    pose_second_camera_wrt_first = eightPointAlgorithm(common_points1, common_points2, K, true);

    # compute the second pose in the trajectory w.r.t. world
    second_pose = cameraRelativeTransform2Pose(first_pose, pose_second_camera_wrt_first, T_cam);

    # apply the scale factor to the translational part
    second_pose(1:3,4) = second_pose(1:3,4)*scale_factor;

    # initialize the poses state
    X_r(:,:,1) = first_pose;
    X_r(:,:,2) = second_pose;

endfunction



# This function, given the associations vector, two measurements with the
# corresponding robot poses from which the measurements are taken and
# the actual landmark constellation, return the augmented landmark constellation
# (the part of the state that contains the landmark poses) togheter with the
# augmented associations_vector.
# Inputs:
#   measurement1
#   measurement2
#   X_1: pose of the robot that performs measurement1
#   X_2: pose of the robot that performs measurement2
#   pose2
#   K: camera matrix
#   T_cam: pose of the camera w.r.t. the robot pose
#   X_l: part of the state that contains all the believed landmark poses
#   associations: vector that, given a landmark id returns the corresponding
#                 id in the state vector
#   associations_inv: vector that, given a landmark id in the state vector
#                     return the corresponding actual landmark id
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
#   X_2_id: id of the pose X_2
#   scale_factor: scale factor to apply to the triangulated points
# Outputs:
#   X_l: part of the state that contains all the believed landmark poses
#   associations: vector that, given a landmark id returns the corresponding
#                 id in the state vector
#   associations_inv: vector that, given a landmark id in the state vector
#                     return the corresponding actual landmark id
#   pending_landmarks: updated version of the input pending_landmarks
function [X_l, associations, associations_inv, pending_landmarks] = augmentConstellation(X_r, measurement1, measurement2, X_1, X_2, K, T_cam, X_l, associations, associations_inv, pending_landmarks, X_2_id, scale_factor)

    # take the common seen landmarks, ordered in the same way
    [common_points1, common_points2, common_ids, uncommon_points, uncommon_ids] = getIntersectingMeasurements(measurement1, measurement2);

    # check if there is some new landmark seen in common
    [common_points1, common_points2, common_ids] = getNewCommons(common_points1, common_points2, common_ids, associations_inv);

    # if there is some new landmark seen, add it to the state with triangulation
    if(!isempty(common_ids))
        
        # compute the pose of the second camera w.r.t. the first
        pose_second_camera_wrt_first = pose2CameraRelativeTransform(X_1, X_2, T_cam);

        # triangulate common points
        # X_l_new is expressed w.r.t. the first camera
        [n_success, X_l_new, errors] = triangulatePoints(inv(pose_second_camera_wrt_first), common_points1, common_points2, K);

        # express points w.r.t the world frame
        X_l_new = pointsCam2World(X_l_new, X_1, T_cam);

        # apply the scale
        X_l_new *= scale_factor;

        # augment landmark state
        X_l = [X_l, X_l_new];

        # update the associations vectors
        for idx = 1:size(common_ids,2)

            strcat("ADDED LANDMARK WITH id: ", int2str(common_ids(idx)))
            
            # compute the new landmark id in the state
            # the current dimension of associations_inv is equal to 
            # the number of landmarks added in the state, so the next
            # will have id equal to its dimension + 1
            new_state_id = size(associations_inv, 2)+1;

            # just update
            associations(common_ids(idx)+1) = new_state_id;
            associations_inv(end+1) = common_ids(idx)+1;

        endfor
    
    endif

    # clean the mutual pending landmarks deleting all the added landmarks
    pending_landmarks = cleanPendingLandmarks(pending_landmarks, common_ids);

    # check if there is some pending landmark that can be solved now
    # and delete them in case or add all the other pending landmark
    # finded
    [pending_landmarks, common_pending1, common_pending2, common_pending_ids] = updatePendingLandmarks(pending_landmarks, ...
                                                                                                    uncommon_points, uncommon_ids, 
                                                                                                    X_2_id, associations_inv);

    # if there is some pending that can be solved, add it to the state with triangulation
    if(!isempty(common_pending1))

        # for each pending landmark
        for pending_l_idx = 1:size(common_pending_ids, 2)

            # get the first pose that saw the landmark
            X_old = X_r(:,:,common_pending_ids(2,pending_l_idx));

            # compute the pose of the second camera w.r.t. the first
            pose_second_camera_wrt_first = pose2CameraRelativeTransform(X_old, X_2, T_cam);

            # triangulate 
            # X_l_new is expressed w.r.t. the first camera
            [n_success, X_l_new, errors] = triangulatePoints(inv(pose_second_camera_wrt_first), [common_pending1(:,pending_l_idx)], [common_pending2(:,pending_l_idx)], K);
            
            # express points w.r.t the world frame
            X_l_new = pointsCam2World(X_l_new, X_old, T_cam);

            # apply the scale
            X_l_new *= scale_factor;

            # augment landmark state
            X_l = [X_l, X_l_new];

            # update the associations vectors
            new_state_id = size(associations_inv, 2)+1;
            associations(common_pending_ids(1,pending_l_idx)+1) = new_state_id;
            associations_inv(end+1) = common_pending_ids(1,pending_l_idx)+1;

        endfor

    endif

endfunction



# This function, given the predicted landmarks and the true landmarks
# pose in the world, compute the error of the predictions (obviosuly only
# on the seen landmarks) as the average distance between the real pose and
# the predicted pose
# Inputs:
#   real_landmark_poses: ground truth of the landmark poses
#   predicted_landmark_poses: predicted landmark poses
#   associations_inv: vector that maps the position of the landmark in 
#                     the state to its landmark id
# Outputs:
#   error
function error = computeLandmarksPredictionError(real_landmark_poses, predicted_landmarks, associations_inv)

    num_predicted_landmarks = size(associations_inv, 2);
    error = 0;
    for idx = 1:num_predicted_landmarks
        predicted_pose = predicted_landmarks(:,idx);
        real_pose = real_landmark_poses(:,associations_inv(idx));
        error += norm(predicted_pose-real_pose);
    endfor
    error /= num_predicted_landmarks;

endfunction



# This function, given the predicted poses and the true trajectory
# poses in the world, compute the error of the predictions as 
# the average distance between the real poses and
# the predicted poses
# Inputs:
#   real_trajectory: ground truth of the trajectory poses
#   predicted_trajectory: predicted trajectory poses
# Outputs:
#   error
function error = computeTrajectoryPredictionError(real_trajectory, predicted_trajectory)

    num_poses = size(predicted_trajectory, 2);
    error = 0;
    for idx = 1:num_poses
        predicted_pose = predicted_trajectory(:,idx);
        real_pose = real_trajectory(:,idx);
        error += norm(predicted_pose-real_pose);
    endfor
    error /= num_poses;

endfunction
