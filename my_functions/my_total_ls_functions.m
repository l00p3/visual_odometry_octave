source "./utilities/geometry_helpers_3d.m"



# retrieves the index in the perturbation vector, that corresponds to
# a certain pose
# input:
#   pose_index:     the index of the pose for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the sub-vector corrsponding to 
#          pose_index, in the array of perturbations  (-1 if error)
function v_idx = poseMatrixIndex(pose_index, num_poses, num_landmarks)
  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*6;
endfunction;



# retrieves the index in the perturbation vector, that corresponds to
# a certain landmark
# input:
#   landmark_index:     the index of the landmark for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the perturnation corrsponding to the
#           landmark_index, in the array of perturbations
function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*6 + (landmark_index-1) * 3;
endfunction;



# Error and jacobian implementation for the projection measurements
# input:
#   Xr: the robot pose in world frame (4x4 homogeneous matrix)
#   Xl: the landmark pose (3x1 vector, 3d pose in world frame)
#   z: projection of the landmark on the image plane
#   K: camera matrix
#   image_rows: number of rows in the projection plane
#   image_cols: number of columns in the projection plane
# output:
#   is_valid: true if projection ok
#   e: 2x1 is the difference between prediction and measurement
#   Jr: 2x6 derivative w.r.t a the error and a perturbation on the
#       pose
#   Jl: 2x3 derivative w.r.t a the error and a perturbation on the
#       landmark
function [is_valid, e, Jr, Jl] = projectionErrorAndJacobian(Xr, Xl, z, K, image_rows, image_cols)

  # initializations
  is_valid = false; # initialize the projection as invalid
  e = [0; 0];
  Jr = zeros(2,6);  # Jacobian w.r.t. the pose involved in the measurement
  Jl = zeros(2,3);  # Jacobian w.r.t. the landmark involved in the measurement
  
  # inverse transform, world w.r.t. to Xr
  iR = Xr(1:3,1:3)';
  it = -iR*Xr(1:3,4);

  # point prediction, in world coordinates
  pw = iR*Xl + it;

  # check if the prediction is in front of camera
  if (pw(3)<0)
     return;
  endif

  # Jacobian of landmark prediction for the ICP part
  # (same as multi-point registration)
  # w.r.t. the part of the state that refers to the
  # robot pose Xr and then w.r.t. the part of the state
  # that refers to the landmark pose
  Jwr = zeros(3,6);
  Jwr(1:3,1:3) = -iR;
  Jwr(1:3,4:6) = iR*skew(Xl);
  Jwl = iR;
  
  # point prediction, in camera frame
  p_cam = K*pw;
  
  # point prediction, on projection plane
  iz = 1./p_cam(3);
  z_hat = p_cam(1:2)*iz;
  
  # check if the point prediction on projection plane
  # is inside the camera frustum
  if (z_hat(1)<0 || 
      z_hat(1)>image_cols ||
      z_hat(2)<0 || 
      z_hat(2)>image_rows)
    return;
  endif;

  # compute the derivative of the projection function
  iz2 = iz*iz;
  Jp = [iz, 0, -p_cam(1)*iz2;
        0, iz, -p_cam(2)*iz2];
  
  # compute the error
  e = z_hat-z;
  
  # compute the final jacobians
  Jr = Jp*K*Jwr;
  Jl = Jp*K*Jwl;
  
  # set to valid the projection
  is_valid = true;

endfunction;



# Linearizes the robot-landmark measurements
#   Xr: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
#   Xl: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
#   Zl:  the measurements (2xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   kernel_threshod: robust kernel threshold
#   K: camera matrix
#   image_rows: number of rows in the projection plane
#   image_cols: number of columns in the projection plane
#   threshold_to_ignore: error threshold that determine if an outlier is too outlier to be considered
# output:
#   H: the resulting H matrix
#   b: the resulting b vector
#   chi_tot: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers
function [H, b, chi_tot, num_inliers] = buildLinearSystemProjections(Xr, Xl, Zl, associations, kernel_threshold, K, image_rows, image_cols, threshold_to_ignore)

    # initilization
    num_poses = size(Xr, 3);
    num_landmarks = size(Xl,2);
    system_size = 6*num_poses + 3*num_landmarks; 
    H = zeros(system_size, system_size);
    b = zeros(system_size, 1);
    chi_tot = 0;
    num_inliers = 0;
    
    # for each measurement
    for (measurement_num = 1:size(Zl,2))
        
        # get the indices
        pose_index = associations(1, measurement_num);
        landmark_index = associations(2, measurement_num);
        
        # get the actual emasurement, pose and landmark
        current_z = Zl(:,measurement_num);
        current_Xr = Xr(:,:,pose_index);
        current_Xl = Xl(:,landmark_index);

        # compute the error and the Jacobian
        [is_valid, e, Jr, Jl] = projectionErrorAndJacobian(current_Xr, current_Xl, current_z, K, image_rows, image_cols);
        
        # if the measurement is not valid ignore it
        if (!is_valid)
            continue;
        endif;
        
        # compute the chi error
        chi = e'*e;

        # robust kernel
        if(chi > threshold_to_ignore)
            continue
        endif
        if (chi>kernel_threshold)
            e *= sqrt(kernel_threshold/chi);
            chi = kernel_threshold;
        else
            num_inliers++;
        endif;

        # update error evolution
        chi_tot += chi;

        # retrieve the positions in the H matrix and b vector
        pose_matrix_index = poseMatrixIndex(pose_index, num_poses, num_landmarks);
        landmark_matrix_index = landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

        # update H and b
        H(pose_matrix_index:pose_matrix_index+6-1,
            pose_matrix_index:pose_matrix_index+6-1) += Jr'*Jr;

        H(pose_matrix_index:pose_matrix_index+6-1,
            landmark_matrix_index:landmark_matrix_index+3-1) += Jr'*Jl;

        H(landmark_matrix_index:landmark_matrix_index+3-1,
            landmark_matrix_index:landmark_matrix_index+3-1) += Jl'*Jl;

        H(landmark_matrix_index:landmark_matrix_index+3-1,
            pose_matrix_index:pose_matrix_index+6-1) += Jl'*Jr;

        b(pose_matrix_index:pose_matrix_index+6-1) += Jr'*e;
        b(landmark_matrix_index:landmark_matrix_index+3-1) += Jl'*e;

    endfor

endfunction



# Error and jacobian of a measured pose, all poses are in world frame
# Inputs:
#   Xi: the observing robot pose (4x4 homogeneous matrix)
#   Xj: the observed robot pose (4x4 homogeneous matrix)
#   Z: the relative transform measured between Xr1 and Xr2
# Outputs:
#   e: 12x1 is the difference between prediction, and measurement, vectorized
#   Ji: 12x6 derivative w.r.t a the error and a perturbation of the
#       first pose
#   Jj: 12x6 derivative w.r.t a the error and a perturbation of the
#       second pose
function [e, Ji, Jj] = poseErrorAndJacobian(Xi, Xj, Z)
    global Rx0;
    global Ry0;
    global Rz0;
    
    Ri=Xi(1:3,1:3);
    Rj=Xj(1:3,1:3);
    ti=Xi(1:3,4);
    tj=Xj(1:3,4);
    tij=tj-ti;
    Ri_transpose=Ri';
    Ji=zeros(12,6);
    Jj=zeros(12,6);
  
    dR_dax=Ri_transpose*Rx0*Rj;
    dR_day=Ri_transpose*Ry0*Rj;
    dR_daz=Ri_transpose*Rz0*Rj;
  
    Jj(1:9,4)=reshape(dR_dax, 9, 1);
    Jj(1:9,5)=reshape(dR_day, 9, 1);
    Jj(1:9,6)=reshape(dR_daz, 9, 1);
    Jj(10:12,1:3)=Ri_transpose;
  
    Jj(10:12,4:6)=-Ri_transpose*skew(tj);
    Ji=-Jj;

    Z_hat=eye(4);
    Z_hat(1:3,1:3)=Ri_transpose*Rj;
    Z_hat(1:3,4)=Ri_transpose*tij;
    e=flattenIsometryByColumns(Z_hat-Z);
endfunction;



# Linearizes the robot-robot measurements
# Inputs:
#   Xr: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
#   Xl: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
#   Zr: the robot_robot measuremenrs (4x4xnum_measurements: array of homogeneous matrices)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=[i_idx, j_idx]' means the kth measurement
#                 refers to an observation made from pose i_idx, that
#                 observed the pose j_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   kernel_threshod: robust kernel threshold
# Outputs:
#   H: the H matrix, filled
#   b: the b vector, filled
#   chi_tot: the total chi2 of the current round
#   num_inliers: number of measurements whose error is below kernel_threshold

function [H, b, chi_tot, num_inliers]=buildLinearSystemPoses(Xr, Xl, Zr, associations, kernel_threshold)

    # initializations
    num_poses = size(Xr, 3);
    num_landmarks = size(Xl,2);
    system_size = 6*num_poses+3*num_landmarks; 
    H = zeros(system_size, system_size);
    b = zeros(system_size,1);
    chi_tot=0;
    num_inliers=0;
    
    # for each measurement
    for (measurement_num = 1:size(Zr,3))

        # create the Omega matrix because we need 
        # to pimp the rotation part a little
        Omega = eye(12);
        Omega(1:9,1:9) *= 1e3;

        # get the indices of the current measurement
        pose_i_index = associations(1,measurement_num);
        pose_j_index = associations(2,measurement_num);

        current_z = Zr(:,:,measurement_num);
        Xi = Xr(:,:,pose_i_index);
        Xj = Xr(:,:,pose_j_index);

        # compute the error and jacobian
        [e,Ji,Jj] = poseErrorAndJacobian(Xi, Xj, current_z);

        # compute the chi error
        chi = e'*Omega*e;

        # robust kernel
        if (chi > kernel_threshold)
            Omega *= sqrt(kernel_threshold/chi);
            chi = kernel_threshold;
        else
            num_inliers ++;
        endif;

        # update stats
        chi_tot += chi;

        # retrieve the positions in the H matrix and b vector
        pose_i_matrix_index = poseMatrixIndex(pose_i_index, num_poses, num_landmarks);
        pose_j_matrix_index = poseMatrixIndex(pose_j_index, num_poses, num_landmarks);

        # update H and b
        H(pose_i_matrix_index:pose_i_matrix_index+6-1,
            pose_i_matrix_index:pose_i_matrix_index+6-1) += Ji'*Omega*Ji;

        H(pose_i_matrix_index:pose_i_matrix_index+6-1,
            pose_j_matrix_index:pose_j_matrix_index+6-1) += Ji'*Omega*Jj;

        H(pose_j_matrix_index:pose_j_matrix_index+6-1,
            pose_i_matrix_index:pose_i_matrix_index+6-1) += Jj'*Omega*Ji;

        H(pose_j_matrix_index:pose_j_matrix_index+6-1,
            pose_j_matrix_index:pose_j_matrix_index+6-1) += Jj'*Omega*Jj;

        b(pose_i_matrix_index:pose_i_matrix_index+6-1) += Ji'*Omega*e;
        b(pose_j_matrix_index:pose_j_matrix_index+6-1) += Jj'*Omega*e;

    endfor

endfunction




# Implementation of the boxplus; pplies a perturbation to a set 
# of landmarks and robot poses
# Inputs:
#   Xr: the robot poses (4x4xnum_poses: array of homogeneous matrices)
#   Xl: the landmark pose (3xnum_landmarks matrix of landmarks)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# Outputs:
#   Xr: the robot poses obtained by applying the perturbation
#   Xl: the landmarks obtained by applying the perturbation
function [Xr, Xl] = boxPlus(Xr, Xl, dx)
  
  # compute the number of poses and landmarks in the state
  num_poses = size(Xr, 3);
  num_landmarks = size(Xl, 2);
  
  # for each pose
  for(pose_index = 1:num_poses)
    # take the positions in the state of the current pose
    pose_matrix_index = poseMatrixIndex(pose_index, num_poses, num_landmarks);

    # take the corresponding perturbation
    dxr = dx(pose_matrix_index:pose_matrix_index+6-1);

    # apply the perturbation
    Xr(:,:,pose_index) = v2t(dxr)*Xr(:,:,pose_index);
  endfor;

  # for each landmark
  for(landmark_index=1:num_landmarks)
    # take the positions in the state of the current landmark
    landmark_matrix_index = landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

    # take the corresponding perturbation
    dxl = dx(landmark_matrix_index:landmark_matrix_index+3-1,:);
    
    # apply the perturbation
    Xl(:,landmark_index) += dxl;
  endfor;

endfunction;



# Implementation of the  final optimization loop with factor
# graph
# Inputs:
#   Xr: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
#   Xl: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
#   Zp: the projective measurements (2xnum_measurements)
#   projection_associations: 2xnum_measurements. 
#                 projection_associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   Zr: the poses measurements, where from each pose the next one is observed
#   num_iterations: the number of iterations of least squares
#   pose_associations: 2xnum_measurements. 
#                 pose_associations(:,k)=[i_idx, j_idx]' means the kth measurement
#                 refers to an observation made from pose i_idx, that
#                 observed the pose j_idx
#   damping: damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold
#   K: camera matrix
#   image_rows: number of rows in the projection plane
#   image_cols: number of columns in the projection plan
#   threshold_to_ignore: error threshold that determine if an outlier is too outlier to be considered
# Outputs:
#   Xr: the robot poses after optimization
#   Xl: the landmarks after optimization
#   chi_stats_{p,r}: array 1:num_iterations, containing evolution of chi2 for projections and poses
#   num_inliers{p,r}: array 1:num_iterations, containing evolution of inliers projections and poses
function [Xr, Xl, chi_stats_p, num_inliers_p,chi_stats_r, num_inliers_r, H, b] = doFinalOptimization(Xr, Xl,
	     Zp, projection_associations,
	     Zr, pose_associations,
	     num_iterations,
	     damping,
	     kernel_threshold,
         K, image_rows, image_cols, threshold_to_ignore)

    # initialization
    num_poses = size(Xr, 3);
    num_landmarks = size(Xl, 2);
    chi_stats_p = zeros(1, num_iterations);
    num_inliers_p = zeros(1, num_iterations);
    chi_stats_r = zeros(1, num_iterations);
    num_inliers_r = zeros(1, num_iterations);

    # size of the linear system
    system_size = 6*num_poses + 3*num_landmarks; 
  
    # start optimization process
    for (iteration=1:num_iterations)

        strcat("FINAL REFINEMENT ITERATION: ", int2str(iteration), "/", int2str(num_iterations))
        
        # initialize H and b
        H = zeros(system_size, system_size);
        b = zeros(system_size, 1);
   
        # check if there are landmarks in the state
        if (num_landmarks)    

            # build the H matrix and the b vector from the projections 
            [H_new, b_new, chi_, num_inliers_] = buildLinearSystemProjections(Xr, Xl, Zp, projection_associations, kernel_threshold, K, image_rows, image_cols, threshold_to_ignore);
            H += H_new;
            b += b_new;
            
            # update stats
            chi_stats_p(iteration) += chi_;
            num_inliers_p(iteration) = num_inliers_;

        endif;

    # build the H matrix and the b vector from the pose graph 
    [H_new, b_new, chi_, num_inliers_] = buildLinearSystemPoses(Xr, Xl, Zr, pose_associations, kernel_threshold);
    H += H_new;
    b += b_new;

    # update the stats
    chi_stats_r(iteration) += chi_;
    num_inliers_r(iteration) = num_inliers_;

    # damping of the H matrix
    H += eye(system_size)*damping;

    # initialize perturbation vector
    dx = zeros(system_size, 1);

    # the linear system is underdetermined so we need
    # to eliminate redundant variables. For this reasor we
    # solve the linear system blocking the first pose,
    # this corresponds to "remove" from H and b the blocks
    # of the 1st pose, while solving the system
    dx(6+1:end) = -(H(6+1:end,6+1:end)\b(6+1:end,1));

    # apply the perturbation
    [Xr, Xl] = boxPlus(Xr, Xl, dx);

  endfor
endfunction