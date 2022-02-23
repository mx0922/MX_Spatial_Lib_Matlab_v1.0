function projDot = mx_getProjectorDot(X_G_i, vG, bodyN)

global tree VBody

I = tree.I{bodyN};

vi = VBody{bodyN};

% ∂‘X_G_i«Ûµº
dX = X_G_i * mx_vCross(vG) - mx_vCross(vi) * X_G_i;

projDot = dX' * I;

end