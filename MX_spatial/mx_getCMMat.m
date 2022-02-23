function [cG, hG, AG, mass, vG, X_G_0, IG] = mx_getCMMat(q, qd)

global tree Xup S

I_0_C = zeros(6, 6);
Ic = tree.I;

NB = tree.NB;
for i = NB:-1:1
    j = tree.parent(i);
    if j == 0
        I_0_C = I_0_C + (Xup{i})' * Ic{i} * Xup{i};
    else
        Ic{j} = Ic{j} + (Xup{i})' * Ic{i} * Xup{i};  
    end
end

mass = I_0_C(6, 6);

cG = mx_skewInv(I_0_C(1:3, 4:6)/mass);

X_G_0 = [eye(3), zeros(3); mx_skew(cG), eye(3)];

hG = zeros(6, 1);

AG = zeros(6, NB);
X_G = cell(1, NB);

for i = 1:NB
    j = tree.parent(i);
    if j == 0
        X_G{i} = Xup{i} * X_G_0;
    else
        X_G{i} = Xup{i} * X_G{j};
    end
    
    AG(:, i) = (X_G{i})' * Ic{i} * S{i};
    
    hG = hG + AG(:, i) * qd(i);    
end

IG = X_G_0' * I_0_C * X_G_0;

vG = IG \ hG;

end