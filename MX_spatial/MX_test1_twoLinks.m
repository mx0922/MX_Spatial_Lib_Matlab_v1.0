%% Two links planar example
% Written by Meng Xiang in BIT 2021/12/1.
clear; close all; clc
addpath('../Roy_spatial')

%% Set up models
global tree
NRJ = 2;
param = getModelParams();
tree = getTree(NRJ, param);

% joint config
q = [1.23; 0.756];
qd = [0.854; 2.581];
qdd = [2.14; 3.256];

%% 3D的方法求解运动学
% 求末端点的位置、速度、Jacobian及其导数
addpath('./twoLinksDerive_3D')

P_EE = autoGen_eePos(q(1),q(2),param.l1,param.l2);

V_EE = autoGen_eeVel(q(1),q(2),qd(1),qd(2),param.l1,param.l2);

J_EE = autoGen_eeJac(q(1),q(2),param.l1,param.l2);

dJ_EE = autoGen_eeJacDot(q(1),q(2),qd(1),qd(2),param.l1,param.l2);

%% 求解动力学矩阵
param.g = 9.81;
[MM,FF] = autoGen_dynMat(q(1),q(2),qd(1),qd(2),param.l1,param.c1,param.c2,param.m1,param.m2,param.I1z,param.I2z,param.g);

TAU = MM * qdd - FF;

[P_COM,ALM,AG] = autoGen_cmmMat(q(1),q(2),qd(1),qd(2),param.l1,param.c1,param.c2,param.m1,param.m2,param.I1z,param.I2z);

%% spatial的方法求解上述量
% 先求一些不变的量，用作全局变量
global S XBase VBody Xup

Xup = cell(1, NRJ);
S = cell(1, NRJ);
XBase = cell(1, NRJ);
VBody = cell(1, NRJ);
for i = 1:NRJ
    [XJ, S{i}] = jcalc(tree.pitch(i), q(i));
    vJ = S{i} * qd(i);
    Xup{i} = XJ * tree.Xtree{i};
    ip = tree.parent(i);
    if ip == 0
        XBase{i} = Xup{i};
        VBody{i} = vJ;
    else
        XBase{i} = Xup{i} * XBase{ip};
        VBody{i} = Xup{i} * VBody{ip} + vJ;
    end        
end

%% 如何求某一点在惯性系下的位置、速度、Jacobian及其导数
% 以末端为例: 在第二个body上，在body2的坐标系中的位置为pos
pos = [0; param.l2; 0];
bodyN = 2;

% 末端位置
P_ee = mx_pointPosWorld(bodyN, pos);
P_ee - P_EE

% 末端速度
V_ee = mx_pointVelWorld(bodyN, pos);
% V_ee - V_EE

% 末端Jacobian
J_ee = mx_pointJacWorld(bodyN, pos);
% J_ee - J_EE

% 末端Jacobian的导数 —— 求加速度会用到
dJ_ee = mx_pointJacDotWorld(bodyN, pos);
% dJ_ee - dJ_EE

%% spatial求解动力学矩阵
a_grav = [0; -9.81; 0];
f_ext = [];

[H,C] = mx_HandC( tree, q, qd, f_ext, a_grav);
% H - MM
% C + FF

tau = ID( tree, q, qd, qdd, f_ext, a_grav);
% tau - TAU

qdd_ab = FDab( tree, q, qd, tau, f_ext, a_grav );
% qdd_ab - qdd

% 动量相关
[P_COM_s, ALM_s, AG_s, MTOT] = mx_getCMMat(q, qd);
% P_COM_s - P_COM
% ALM_s - ALM
% AG_s - AG
% MTOT - param.m1 - param.m2

function model = getTree(nb, p)

model.NB = nb;
model.pitch = zeros(1, nb);

model.parent = [0, 1];

model.Xtree{1} =  Xtrans([0 0 0]);
model.Xtree{2} =  Xtrans([0 p.l1 0]);

CoM{1} = [0 p.c1 0];
CoM{2} = [0 p.c2 0];

mass = [p.m1, p.m2];

Icm{1} = diag([p.I1x, p.I1y, p.I1z]);
Icm{2} = diag([p.I2x, p.I2y, p.I2z]);

for i = 1:nb
    model.I{i} = mcI(mass(i), CoM{i}, Icm{i});    
end

end

function p = getModelParams()

p.l1 = 0.46;
p.c1 = 0.25;
p.m1 = 5.34;
p.I1x = 0.1025;
p.I1y = 0.024;
p.I1z = 0.1123;

p.l2 = 0.52;
p.c2 = 0.31;
p.m2 = 6.12;
p.I2x = 0.1209;
p.I2y = 0.019;
p.I2z = 0.0983;

end
