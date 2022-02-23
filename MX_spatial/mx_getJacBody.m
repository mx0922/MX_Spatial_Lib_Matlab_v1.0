function Jb = mx_getJacBody(bodyN)

global tree S XBase

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
        J(:, i) = XBase{i}^(-1) * S{i};
    end    
end

Jb = XBase{bodyN} * J;

end