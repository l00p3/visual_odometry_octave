# triangulates a point, passing through two lines
# one passing through the origin, and having
# direction vector d1 (for which p1 is on the origin so 0)
# one passing through a point p2, and having
# direction d2
function [success, p, e]=triangulatePoint(p2, d1, d2)
  
  p=zeros(3,1);
  success=false;
  e=-1;     # error: distance where the lines "meet"
                      
  D=[-d1, d2];                # assemble system matrix to find ascissa 
  s=-(D'*D)\(D'*p2);          # s: ascissa of closest point on p1 and p2
  
  # if the point is behind the camera, return with success false
  if (s(1)<0 || s(2)<0)
    return;
  endif;

  # otherwise we succeed
  success=true;
  
  p1_triangulated=d1*s(1);    # point on 1st line
  p2_triangulated=d2*s(2)+p2; # point on 2nd line

  # difference between the points (basically how far are the triangulated points p1 and p2)
  e=norm(p1_triangulated-p2_triangulated); 

  # midpoint, the triangulated point             
  p=0.5*(p1_triangulated+p2_triangulated); 
               
endfunction;



# Triangulates a batch of points in image coords,
# Inputs
#   X: the pose of the 1st camera w.r.t the 2nd camera
#   P1_img: list of 2D points seen on the projective plane
#           of the 1st camera
#   P2_img: same points in P1_img but seen on the projective
#           plane of the 2nd camera
#   K: camera matrix
# Outputs:
#   n_success: 
#   P: bunch of triangulated points, a triangulated point is a point 
#       in the space where, given two points in two projective planes 
#       (seen by two cameras), and considering such points as 
#       directions, the two directions intersect
#   errors: corresponding errors, that are the distances between the 
#             lines where they "intersect"
function [n_success, P, errors] = triangulatePoints(X, P1_img, P2_img, K)
  
  #initialize vars
  n_points=size(P1_img,2);
  P=zeros(3,n_points);        # here we accumulate the triangulated points
  errors=zeros(1,n_points);   # here we accumulate the errors
  n_success=0;
  
  # inverse transform
  iX=inv(X);

  # inversa camera matrix
  iK = inv(K);
  
  # inverse rotation * inverse camera matrix
  iRiK=iX(1:3,1:3)*iK;

  # express the points in camera coordinates (adding 1 after x and y coords)
  # and rotate the direction vector of P2 in world frame
  D1_cam=iK*[P1_img; ones(1,n_points)];
  D2_cam=iRiK*[P2_img; ones(1,n_points)];

  # position of the origin of camera 2
  p2=iX(1:3,4);

  # triangulate each couple of point 
  for (i=1:n_points)
    p1_cam=D1_cam(:,i);
    p2_cam=D2_cam(:,i);

    # triangulation
    [success, p, e]=triangulatePoint(p2, p1_cam, p2_cam);

    if (success==true)
      ++n_success;
      P(:,n_success)=p;
      errors(n_success)=e;
    endif;

  endfor;

endfunction



# Computes a preconditioning matrix, that scales the points around the
# origin
# A(1:2,1:2): inverse sigma of the points
# A(1:2,3)  : -inverse sigma*mean
function A=computePreconditioningMatrix(P_img)
  n_points=size(P_img,2);
  P_img=[P_img; ones(1,n_points)];
  s=sum(P_img,2);
  s2=P_img*P_img';
  mu=s/n_points;
  sigma=s2/n_points-mu*mu';
  A=eye(3);
  A(1:2,1:2)=inv(chol(sigma(1:2,1:2)));
  A(1:2,3)=-A(1:2,1:2)*mu(1:2);
endfunction



# Take a transformation matrix that represent a pose of a 
# camera w.r.t. to the origin of the world frame
# and return the corresponding foundamental matrix E
function F = transform2fundamental(K, X)
  iK = inv(K);
  F = iK'*transform2essential(X)*iK;
endfunction



# Take a transformation matrix that represent a pose of a 
# camera w.r.t. to the origin of the world frame
# and return the corresponding essential matrix E
function  E = transform2essential(X)
  E = X(1:3,1:3)'*skew(X(1:3,4));
endfunction;



# This function here compute the transformation of a camera
# w.r.t. another camera just considering a set of common points 
# seen on the projective image plane
# Inputs:
#   P1_img: a set of points seen from the camera 1
#   P2_img: the same set of points seen from the camera 2
#   preconditioning: if use preconditioning or not
# Outputs:
#   X_pred: the transformation matrix that represent the pose of 
#           the camera 2 w.r.t. the camera 1
function [X_pred] = eightPointAlgorithm(P1_img, P2_img, K, preconditioning)
    # --- ALGORITHM ---------------------------------------------
    # 1. Ensemble the H matrix
    # 2. Extract the minimal eigenvector from H
    # 3. Build the F matrix from it (important here we estimate F)
    # 4. Extract the essentiam matrix E from F
    # 5. Compute the SVD decomposition of the E matrix
    # 6. Extract the R matrix and t vector form the decomposition (2 solutions)
    # 7. Ensemble them in a transformation matrix
    # 8. Choose the best from the 2 solutions

    # some initialization
    n_points = size(P1_img, 2);
    W=[0, -1,  0;
       1,  0,  0;
       0,  0,  1];

    # Preconditioning (scale points around the origin)
    if(preconditioning)
        A1 = computePreconditioningMatrix(P1_img);
        A2 = computePreconditioningMatrix(P2_img);

        P1_img = A1(1:2,1:2)*P1_img+repmat(A1(1:2,3),1,n_points);
        P2_img = A2(1:2,1:2)*P2_img+repmat(A2(1:2,3),1,n_points);
    endif

    # Ensemble the H matrix
    H=zeros(9,9);
    for (i=1:n_points)
        p1_img=[P1_img(:,i); 1];
        p2_img=[P2_img(:,i); 1];
        A = reshape(p1_img*p2_img',1,9);
        H += A'*A;
    endfor;

    # Extract the eigenvector from H and build the F matrix
    [V,lambda] = eig(H);
    F_pred = reshape(V(:,1),3,3);

    # Preconditioning
    if(preconditioning)
        F_pred = A1'*F_pred*A2;
    endif

    # Extract the essential matrix
    E = K'*F_pred*K;

    # Compute the SVD decomposition
    [U,S,V]=svd(E);

    # Extract the R matrix (2 solutions)
    R1 = V*W*U';
    if (det(R1)<0) #right handed condition
        [U,S,V]=svd(-E);
        R1=V*W*U';
    endif;
    R2 = V*W'*U';

    # Extract the t vector (2 solutions)
    t_cross=R1*E;
    t1 = [t_cross(3,2)-t_cross(2,3);
          t_cross(1,3)-t_cross(3,1);
          t_cross(2,1)-t_cross(1,2)];
    t_cross=R2*E;
    t2 = [t_cross(3,2)-t_cross(2,3);
          t_cross(1,3)-t_cross(3,1);
          t_cross(2,1)-t_cross(1,2)];

    # Ensemble all in two transformation matrices
    X_pred1 = eye(4);
    X_pred1(1:3, 1:3) = R1;
    X_pred1(1:3, 4) = t1;
    X_pred2 = eye(4);
    X_pred2(1:3, 1:3) = R2;
    X_pred2(1:3, 4) = t2;

    # Choose the best prediction (the transformation that, 
    # triangulating all the points has the most points)
    n_in_front=0;
    X_test=X_pred1;
    [n_test, P]=triangulatePoints(X_test, P1_img, P2_img, K);
    if (n_test>n_in_front)
        X=X_test;
        n_in_front=n_test;
    endif;
    X_test(1:3,4)=-X_test(1:3,4);
    [n_test, P]=triangulatePoints(X_test, P1_img, P2_img, K);
    if (n_test>n_in_front)
        X=X_test;
        n_in_front=n_test;
    endif;

    X_test=X_pred2;
    [n_test, P]=triangulatePoints(X_test, P1_img, P2_img, K);
    if (n_test>n_in_front)
        X=X_test;
        n_in_front=n_test;
    endif;
    X_test(1:3,4)=-X_test(1:3,4);
    [n_test, P]=triangulatePoints(X_test, P1_img, P2_img, K);
    if (n_test>n_in_front)
        X=X_test;
        n_in_front=n_test;
    endif;

    X_pred = X;
endfunction