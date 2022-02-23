function dJ = mx_pointJacDotWorld(bodyN, pos)

global tree VBody XBase S

Pos = mx_pointPosWorld(bodyN, pos);
POS = Xtrans(Pos);

Vel_6D = mx_pointVelWorld(bodyN, pos);
Vel = Vel_6D(4:6);
VEL = [zeros(3, 6); -mx_skew(Vel), zeros(3)]; % 注意！别搞错了，这是对POS求导
% VEL = Xtrans(Vel); % 这是错的！！！

NB = tree.NB;
e = zeros(1, NB);
body = bodyN;
while body ~= 0
    e(body) = 1;
    body = tree.parent(body);
end

dJ = zeros(6, NB);
for i = 1:NB
    if e(i)
        X_i_0 = XBase{i}^(-1);
        dX = X_i_0 * mx_vCross(VBody{i});      
        dJ(:, i) = VEL * X_i_0 * S{i} + POS * dX * S{i};
    end    
end

end