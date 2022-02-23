function [hG, cG, IG, vG] = mx_getHGandCG_IGvG(model, q, qd)

I_0_C = zeros(6, 6);

NB = model.NB;

Ic = cell(1, NB);
Xup = cell(1, NB);
S = cell(1, NB);

for i = 1:NB
    Ic{i} = model.I{i};  
    [XJ, S{i}] = jcalc( model.pitch(i), q(i) );
    Xup{i} = XJ * model.Xtree{i};
end

for i = NB:-1:1
    j = model.parent(i);
    if j == 0
        I_0_C = I_0_C + (Xup{i})' * Ic{i} * Xup{i};
    else
        Ic{j} = Ic{j} + (Xup{i})' * Ic{i} * Xup{i};  
    end
end

Mtotal = I_0_C(6, 6);

cG = mx_skewInv(I_0_C(1:3, 4:6)/Mtotal);

X_G_0 = [eye(3), zeros(3); mx_skew(cG), eye(3)];

hG = zeros(6, 1);

X_G = cell(1, NB);
AG = zeros(6, NB);
for i = 1:NB
    j = model.parent(i);
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