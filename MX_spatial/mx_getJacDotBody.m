function dJb = mx_getJacDotBody(bodyN)

global tree S XBase VBody

NB = tree.NB;

e = zeros(1, NB);
body = bodyN;
while body ~= 0
    e(body) = 1;
    body = tree.parent(body);
end

dJb = zeros(6, NB);

for i = 1:NB
    if e(i)
        if i == bodyN
            dJb(:, i) = zeros(6, 1);
        else
            Xtemp = XBase{bodyN} * (XBase{i})^(-1);
            dX = Xtemp * mx_vCross(VBody{i}) - mx_vCross(VBody{bodyN}) * Xtemp;
            dJb(:, i) = dX * S{i};
        end
    end    
end

end