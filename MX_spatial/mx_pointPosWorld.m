function [Pos, R] = mx_pointPosWorld(bodyN, pos)

global XBase

X_0_b = XBase{bodyN};
% X_0_b = [R_0_b, 0; R_0_b * mx_skew(P_b_0)', R_0_b]

% R_b_0: bodyN连杆的姿态
R_0_b = X_0_b(1:3, 1:3);
R_b_0 = R_0_b';

% P_b_0: bodyN坐标系原点在World坐标系下的位置
temp1 = X_0_b(4:6, 1:3);
P_b_0 = mx_skewInv((R_b_0 * temp1)');

Pos = P_b_0 + R_b_0 * pos;

if nargout > 1
    R = R_b_0;
end

end