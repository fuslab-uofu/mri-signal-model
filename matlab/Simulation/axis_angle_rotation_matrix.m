function R = axis_angle_rotation_matrix(ax, ang)
%% R = axis_angle_rotation_matrix(ax, ang)
% Generates a rotation matrix by the given angle about the specified axis
%
% ~ Input ~
% * ax: The axis about which to rotate. Expects column vector in R3 or 
% matrix of column vectors
% * ang: The angle (in radians) by which to rotate about each specified
% axis. Expects one entry in ang for every column in ax.
% ~ Output ~
% * R: 3x3 rotation matrix by ang radians about axis ax. If input contains
% multiple rotation axes/column vectors, the corresponding matrix for each
% rotation is stored in each page (dim3), i.e., R(:, :, 2) is the rotation
% matrix for the rotation by ang(2) radians about axis ax(:, 2)
%
% ~ Example ~
% >> in = [1, 1
%          0, 0
%          0, 0];
% >> ax = [0, 0
%          1, 0
%          0, 1]; 
% >> ang = [-pi/2, -pi/2];
% >> axis_angle_rotation_matrix(ax, ang)
% ans(:,:,1) =
%     0    0   -1
%     0    1    0
%     1    0    0
% ans(:,:,2) =
%     0    1    0
%    -1    0    0
%     0    0    1
%
%% 2023-05-04 Samuel Adams-Tew

ang(isnan(ang)) = 0;
u = ax./sqrt(sum(ax.*ax)); % Convert to unit vector(s)
u(isnan(u)) = 0;

c = cos(ang);
s = sin(ang);
C = 1 - c;

xyC = u(1,:,:,:,:,:).*u(2,:,:,:,:,:).*C;
xzC = u(1,:,:,:,:,:).*u(3,:,:,:,:,:).*C;
yzC = u(2,:,:,:,:,:).*u(3,:,:,:,:,:).*C;
xs = u(1,:,:,:,:,:).*s;
ys = u(2,:,:,:,:,:).*s;
zs = u(3,:,:,:,:,:).*s;

R = [u(1,:,:,:,:,:).^2.*C + c, xyC - zs, xzC + ys
     xyC + zs, u(2,:,:,:,:,:).^2.*C + c, yzC - xs
     xzC - ys, yzC + xs, u(3,:,:,:,:,:).^2.*C + c];

end