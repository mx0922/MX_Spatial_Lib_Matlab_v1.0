function  fbmodel = mx_floatBase( model , typeEuler)

if any( model.parent(2:model.NB) == 0 )
  error('only one connection to a fixed base allowed');
end

if ~isequal( model.Xtree{1}, Xtrans([0,0,0]) )
  warning('Xtree{1} not identity');
end

fbmodel.NB = model.NB + 5;

fbmodel.pitch(1:6) = [inf,inf,inf,0,0,0];
fbmodel.pitch(7:fbmodel.NB) = model.pitch(2:model.NB);

fbmodel.parent(1:6) = [0 1 2 3 4 5];
fbmodel.parent(7:fbmodel.NB) = model.parent(2:model.NB) + 5;

% 写成惯用的ZYX欧拉角
fbmodel.Xtree{1} = Xroty(pi/2);
fbmodel.Xtree{2} = Xrotx(-pi/2) * Xroty(-pi/2);
fbmodel.Xtree{3} = Xrotx(pi/2);

% fbmodel.Xtree{4} = Xroty(pi/2);
% fbmodel.Xtree{5} = Xrotx(-pi/2) * Xroty(-pi/2);
% fbmodel.Xtree{6} = Xrotx(pi/2);

if nargin < 2
    fbmodel.Xtree{4} = Xroty(0);
    fbmodel.Xtree{5} = Xrotx(-pi/2) * Xroty(0);
    fbmodel.Xtree{6} = Xroty(pi/2) * Xrotx(pi/2);
else
    if strcmp(typeEuler, 'ZYX')
        fbmodel.Xtree{4} = Xroty(0);
        fbmodel.Xtree{5} = Xrotx(-pi/2) * Xroty(0);
        fbmodel.Xtree{6} = Xroty(pi/2) * Xrotx(pi/2);
    else
        fbmodel.Xtree{4} = Xroty(pi/2);
        fbmodel.Xtree{5} = Xrotx(-pi/2) * Xroty(-pi/2);
        fbmodel.Xtree{6} = Xrotx(pi/2);
    end
end

for i = 7:fbmodel.NB
  fbmodel.Xtree{i} = model.Xtree{i-5};
end

for i = 1:fbmodel.NB
  if i < 6
    fbmodel.I{i} = mcI( 0, [0,0,0], zeros(3) );
  else
    fbmodel.I{i} = model.I{i-5};
  end
end

end