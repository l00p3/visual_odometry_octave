 source "./utilities/geometry_helpers_3d.m"



# Compute the error and Jacobian for projective ICP
# Inputs: 
#   X_guess: guessed pose of the world w.r.t. the camera
#   X_l: point in the space where we think is located the
#        measured landmark
#   z: measured point on the projective plane
#   K: camera matrix
#   z_near: min distance of the space where the camera can perceive stuff
#   z_far: max distance of the space where the camera can perceive stuff
#   width: width of the camera projective plane
#   height: height of the camera projective plane
# Outputs:
#   is_valid: true if projection is ok
#   e: error
#   J: jacobian matrix of the error w.r.t. state X_guess 
function [is_valid, e, J] = errorAndJacobianPosit(X_guess, p, z, K, z_near, z_far, cam_width, cam_height)

    # initialization
    is_valid = false;   # initialize the projection as invalid
    e = [0; 0];
    J = zeros(2,6);

    # compute the position of the point w.r.t. camera frame
    p_cam = X_guess(1:3,1:3)*p + X_guess(1:3,4);

    # check if the prediction is in front of camera
    # or not too far
    if (p_cam(3) < z_near || p_cam(3) > z_far)
        return;
    endif

    # compute the prediction (projection)
    p_cam_k = K*p_cam;
    iz = 1./p_cam_k(3);
    z_hat = p_cam_k(1:2)*iz;

    # check if the point prediction on projection plane
    # is inside the camera frustum
    if (z_hat(1)<0 || 
            z_hat(1)>cam_width ||
            z_hat(2)<0 || 
            z_hat(2)>cam_height)
        return;
    endif;

    # compute the error
    e = z_hat - z;

    # Jacobian of landmark prediction for the ICP part
    # (same as multi-point registration)
    Jr = zeros(3,6);
    Jr(1:3,1:3) = eye(3);
    Jr(1:3,4:6) = skew(-p_cam);
  
    # other component of Jacobian chain rule
    iz2 = iz*iz;
    Jp = [iz, 0,  -p_cam_k(1)*iz2;
          0,  iz, -p_cam_k(2)*iz2];

    # compute the final Jacobian
    J = Jp*K*Jr;

    # set to valid the projection
    is_valid=true;

endfunction




# Function that perform, given a measurement taken from a camera,
# the projective ICP to get the pose of the world w.r.t. the camera
# from which such measurement are taken. There will be taken in consideration
# only such measurements for which the pose of the landmark is already
# triangulated (so guessed)
# Inputs:
#   measurement
#   X_l: set of already known landmarks
#   associations: in position i has the position j in the state X_l of the 
#                 landmark with id i
#   num_iterations
#   kernel_threshold: threshold for the outliers
#   damping: damping factor
#   K: camera matrix
#   T_cam: pose of the camera w.r.t. robot
#   z_near: min distance of the space where the camera can perceive stuff
#   z_far: max distance of the space where the camera can perceive stuff
#   width: width of the camera projective plane
#   height: height of the camera projective plane
#   threshold_to_ignore: error threshold that determine if an outlier is too outlier to be considered
# Outputs:
#   X_guess: guessed pose of the camera that better explain the measurements
#            (actually the pose of the world w.r.t. the camera)
#   chi_stats: error evolution
#   num_inliers: num inliers for each iteration
function [X_guess, chi_stats, num_inliers] = doProjectiveICP(measurement, X_l, associations, ...
                                                             num_iterations, kernel_threshold, damping, ...
                                                             K, T_cam, z_near, z_far, cam_width, cam_height, ...
                                                             threshold_to_ignore)

    # compute the initial guess from the odometry measurements
    # X_guess will be the pose of the world w.r.t. robot
    X_guess = planePose2Transform(measurement.odometry_pose);   # pose of the robot from odometry
    X_guess = rw2wcam(X_guess, T_cam);  # pose of the world w.r.t. camera

    # get only the useful measurements
    [X_l_new, z] = getPosesAndMeasurements(measurement, X_l, associations);

    # compute how many measurement correnspond to already seen landmarks
    n_landmarks = size(X_l_new, 2);

    # if no measurement correspond to an already seen landmark just return
    if(n_landmarks == 0)
        return;
    endif

    # initialize some LS stats
    chi_stats = zeros(1, num_iterations);
    num_inliers = zeros(1, num_iterations);
  
    for (iteration = 1:num_iterations)
        
        # reset all to zero
        H = zeros(6,6);
        b = zeros(6,1);
        chi_stats(iteration) = 0;

        # for each landmark
        for (i = 1:n_landmarks)

            # compute error and jacobian
            [is_valid, e, J] = errorAndJacobianPosit(X_guess, X_l_new(:,i), z(:,i), K, z_near, z_far, cam_width, cam_height);
            
            # check if the projection is ok
            if(!is_valid)
                continue
            endif
            
            # compute chi
            chi = e'*e;
            
            # deal with outliers
            if(chi > threshold_to_ignore)
                continue
            endif
            if (chi > kernel_threshold)
	            e *= sqrt(kernel_threshold/chi);
	            chi = kernel_threshold;
            else
	            num_inliers(iteration)++;
            endif;
        
            # update stats
            chi_stats(iteration)+=chi;

            # update H matrix and b vector
            H += J'*J;
            b += J'*e;
    
        endfor
    
        # damping of the H matrix
        H += eye(6)*damping;

        # solve linear system to get the perturbation
        dx = -H\b;

        # apply the perturbation
        X_guess = v2t(dx)*X_guess;
  
    endfor

    # pose of the robot w.r.t. world
    X_guess = wcam2rw(X_guess, T_cam); 

endfunction
