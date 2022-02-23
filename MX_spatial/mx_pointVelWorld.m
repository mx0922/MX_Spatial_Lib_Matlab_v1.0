function VEL = mx_pointVelWorld(bodyN, pos)

global XBase VBody

Pos = mx_pointPosWorld(bodyN, pos);

% 注意：这里要将bodyN在bodyN坐标系中的速度转换到世界坐标系中表示
Vb = XBase{bodyN} \ VBody{bodyN};

Omega = Vb(1:3);
% Vel = v + cross(w) * Pos —— 要乘以该点在世界坐标系中的位置
Vel = Vb(4:6) + mx_skew(Omega) * Pos;

VEL = [Omega; Vel];

end