function J = mx_pointJacWorld(bodyN, pos)

global tree S XBase

Pos = mx_pointPosWorld(bodyN, pos);
POS = Xtrans(Pos);

NB = tree.NB;

e = zeros(1, NB);
body = bodyN;
while body ~= 0
    e(body) = 1;
    body = tree.parent(body);
end

J = zeros(6, NB);

for i = 1:NB
    if e(i)
        J(:, i) = POS * XBase{i}^(-1) * S{i};
    end    
end

end