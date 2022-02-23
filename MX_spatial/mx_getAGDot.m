function [dAG, AG] = mx_getAGDot(q, qd)

global tree XBase

[~, ~, AG, ~, vG, X_G_0] = mx_getCMMat(q, qd);

NB = tree.NB;
dAG = zeros(6, NB);

for i = 1:NB
    % body jacobian
    Jb = mx_getJacBody(i);
    
     % body jacobian dot
    dJb = mx_getJacDotBody(i);
    
    % projector
    X_0_i = XBase{i};
    X_G_i = X_0_i * X_G_0;
    
    proj = X_G_i' * tree.I{i};
    
    vG = [zeros(3, 1); vG(4:6)];
    % projector dot
    projDot = mx_getProjectorDot(X_G_i, vG, i);
    
    dAG = dAG + proj * dJb + projDot * Jb;    
end

end