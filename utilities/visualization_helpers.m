# Here we have utility functions for visualization



# visualize two set of 3D points in the space
function visualize_points(points1, legend1, points2, legend2, title_txt)
    plot3(points1(1,:), points1(2,:), points1(3,:), 'b*', "linewidth", 2);
    hold on;
    plot3(points2(1,:), points2(2,:), points2(3,:), 'ro', "linewidth", 2);
    legend(legend1, legend2);
    title(title_txt);
    grid;
endfunction



# visualize two set of 2D poses
function visualize_poses(points1, legend1, points2, legend2, title_txt)
    plot(points1(1,:), points1(2,:), 'b*', "linewidth", 2);
    hold on;
    plot(points2(1,:), points2(2,:), 'ro', "linewidth", 2);
    legend(legend1, legend2);
    title(title_txt);
    grid;
endfunction



# visualize the evolution of the error
function visualize_chi_evolution(chi_stats)
    plot(chi_stats, 'r-', "linewidth", 2);
    legend("chi");
    title("chi evolution");
    grid; 
    final_chi = chi_stats(end);
    xlabel(strcat("iterations\nfinal chi value:\t", mat2str(final_chi)));
endfunction



# visualize the evolution of the inliers
function visualize_inliers_evolution(num_inliers)
    plot(num_inliers, 'b-', "linewidth", 2);
    legend("#inliers");
    title("inliers evolution");
    grid; 
    final_inliers_num = num_inliers(end);
    xlabel(strcat("iterations\nfinal #inliers:\t", mat2str(final_inliers_num)));
endfunction



# visualize the H matrix in a readible way
function visualize_H(H)
    H_ =  H./H;                      # NaN and 1 element
    H_(isnan(H_))=0;                 # Nan to Zero
    H_ = abs(ones(size(H_)) - H_);   # switch zero and one
    H_ = flipud(H_);                 # switch rows
    colormap(gray(64));
    hold on;
    image([0.5, size(H_,2)-0.5], [0.5, size(H_,1)-0.5], H_*64);
    hold off;
endfunction