% derive of two links
clear; close all; clc
%% 初始化
i = sym([1; 0; 0]);
j = sym([0; 1; 0]);
k = sym([0; 0; 1]);

syms l1 c1 l2 c2 'real'
syms q1 q2 qd1 qd2 'real'

q = [q1; q2];
qd = [qd1; qd2];

derivative = @(in)( jacobian(in, q) * qd );
Rot_z = @(q)( [cos(q) * i + sin(q) * j, cos(q) * j - sin(q) * i, k] );

%% 运动学
% 点的位置
R1 = Rot_z(q1);
P_c1 = R1 * (c1 * j);
P_l1 = R1 * (l1 * j);

R2 = Rot_z(q2);
R12 = R1 * R2;
P_c2 = P_l1 + R12 * (c2 * j);
P_ee = P_l1 + R12 * (l2 * j);

% 末端jacobian
Jac_w = [k, R1 * k]; % 需要思考一下
Jac_p = jacobian(P_ee, q);
Jac_ee = [Jac_w; Jac_p];

% 末端速度
Vp_ee = simplify(derivative(P_ee));
% Vp_ee = Jac_p * qd;
Vw_ee = Jac_w * qd;
V_ee = [Vw_ee; Vp_ee];

% 末端jacobian的导数
dJac_ee = sym(zeros(6, 2));
for ii = 1:2
    dJac_ee(:, ii) = derivative(Jac_ee(:, ii));
end

%% 动力学
syms m1 I1x I1y I1z m2 I2x I2y I2z g 'real'
syms qdd1 qdd2 'real'

I1 = diag([I1x I1y I1z]);
I2 = diag([I2x I2y I2z]);

qdd = [qdd1; qdd2];

Vp_c1 = derivative(P_c1);
Vw_c1 = qd1 * k;

Vp_c2 = derivative(P_c2);
Vw_c2 = Vw_ee;

I1w = R1 * I1 * R1';
I2w = R12 * I2 * R12';

% 质心的位置
P_COM = (m1 * P_c1 + m2 * P_c2) / (m1 + m2);

% 连杆相对质心的动量
LM = m1 * Vp_c1 + m2 * Vp_c2;
AM = cross(P_c1 - P_COM, m1 * Vp_c1) + I1w * Vw_c1 + cross(P_c2 - P_COM, m2 * Vp_c2) + I2w * Vw_c2;
ALM = simplify([AM; LM]);

% 动能
KE = 0.5 * m1 * dot(Vp_c1, Vp_c1) + 0.5 * Vw_c1' * I1w * Vw_c1 + ...
    0.5 * m2 * dot(Vp_c2, Vp_c2) + 0.5 * Vw_c2' * I2w * Vw_c2;
% 势能
G = g * j;
PE = m1 * dot(P_c1, G) + m2 * dot(P_c2, G);

% lagrangian dynamics
DK_Dqd = jacobian(KE, qd);
Term1 = jacobian(DK_Dqd, q) * qd + jacobian(DK_Dqd, qd) * qdd;
Term2 = jacobian(KE, q) - jacobian(PE, q);

eqns = Term1 - Term2';
[MM, FF] = equationsToMatrix(eqns, qdd);

% 质心动量矩阵 AG
AG = equationsToMatrix(ALM, qd);

%% functions

matlabFunction(...
    P_ee, ...
    'file', 'autoGen_eePos.m', ...
    'vars', {'q1', 'q2', 'l1', 'l2'});

matlabFunction(...
    V_ee, ...
    'file', 'autoGen_eeVel.m', ...
    'vars', {'q1', 'q2', 'qd1', 'qd2', 'l1', 'l2'});

matlabFunction(...
    Jac_ee, ...
    'file', 'autoGen_eeJac.m', ...
    'vars', {'q1', 'q2', 'l1', 'l2'});

matlabFunction(...
    dJac_ee, ...
    'file', 'autoGen_eeJacDot.m', ...
    'vars', {'q1', 'q2', 'qd1', 'qd2', 'l1', 'l2'});

matlabFunction(...
    MM, FF, ...
    'file', 'autoGen_dynMat.m', ...
    'vars', {'q1', 'q2', 'qd1', 'qd2', 'l1', 'c1', 'c2', 'm1', 'm2', 'I1z', 'I2z', 'g'});

matlabFunction(...
    P_COM, ALM, AG, ...
    'file', 'autoGen_cmmMat.m', ...
    'vars', {'q1', 'q2', 'qd1', 'qd2', 'l1', 'c1', 'c2', 'm1', 'm2', 'I1z', 'I2z'});