function VEL = mx_pointVelWorld(bodyN, pos)

global XBase VBody

Pos = mx_pointPosWorld(bodyN, pos);

% ע�⣺����Ҫ��bodyN��bodyN����ϵ�е��ٶ�ת������������ϵ�б�ʾ
Vb = XBase{bodyN} \ VBody{bodyN};

Omega = Vb(1:3);
% Vel = v + cross(w) * Pos ���� Ҫ���Ըõ�����������ϵ�е�λ��
Vel = Vb(4:6) + mx_skew(Omega) * Pos;

VEL = [Omega; Vel];

end