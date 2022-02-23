%% Float base example
% Written by Meng Xiang in BIT 2021/12/1.
clear; close all; clc
addpath('../Roy_spatial')

%% Set up models
global tree

NRJ = 13;
param = getModelParams();
model = getFloatBaseTree(NRJ, param);
% tree = mx_floatBase(model, 'XYZ');
tree = mx_floatBase_New(model, 'ZYX'); % options: 'XYZ'(default), 'ZYX'

q = [0.5000, 0.1000, 0.8000, 0.1745, 0.2618, 0.6109, 0.0873, -0.1745, -0.2618, 0.8727, -0.4363, 0.2618, 0.0524, 0.2618, -0.3491, 1.0472, -0.6109, 0.1745]';
qd = [0.5000, 0.1000, 0.8000, 0.1745, 0.2618, 0.6109, 0.0873, -0.1745, -0.2618, 0.8727, -0.4363, 0.2618, 0.0524, 0.2618, -0.3491, 1.0472, -0.6109, 0.1745]';

% 先求一些不变的量，用作全局变量
global S XBase VBody Xup

NB = tree.NB;
Xup = cell(1, NB);
S = cell(1, NB);
XBase = cell(1, NB);
VBody = cell(1, NB);
for i = 1:NB
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

%% 运动学相关
pos = [0.08; 0; 0];

% 注意：这里还可以输出姿态R，但这个姿态还需要乘以一个3D旋转矩阵
[P_ee, R_ee] = mx_pointPosWorld(NB, pos);
R_ee = R_ee * mx_Roty_3D(-pi/2);

V_ee = mx_pointVelWorld(NB, pos);

J_ee = mx_pointJacWorld(NB, pos);

dJ_ee = mx_pointJacDotWorld(NB, pos);

%% 动力学相关
a_grav = [0; 0; -9.81];
f_ext = [];

% M、H矩阵
[H,C] = mx_HandC( tree, q, qd, f_ext, a_grav);

% ID & FD
tau_mx = ones(12, 1);
tau = [0;0;0;0;0;0;tau_mx];
qdd = FDab( tree, q, qd, tau );

% X：X_0_fb; v, a: 浮动基的速度、加速度在世界坐标系下的表示
[X,v,a] = fbKin( q, qd, qdd );

[a1,tau1] = IDf( tree, X, v, q(7:end), qd(7:end), qdd(7:end) );
[a2,qdd2] = FDf( tree, X, v, q(7:end), qd(7:end), tau(7:end) );
% [ a1-a, a2-a ]
% [ tau1-tau(7:end), qdd2-qdd(7:end) ]

%% 动量相关，A、dA矩阵
[P_COM_s, ALM_s, AG_s, MTOT, vG, X_G_0, IG] = mx_getCMMat(q, qd);

% 求dA矩阵
dAG = mx_getAGDot(q, qd);


function model = getFloatBaseTree(nb, p)

model.NB = nb;
model.pitch = zeros(1, nb);

model.parent = [0,  1, 2, 3, 4, 5, 6,  1, 8, 9, 10, 11, 12]; % 先右后左

model.Xtree{1} =  Xtrans([0 0 0]);
% 右腿
model.Xtree{2} =  Xtrans([p.j2x p.j2y 0]);
model.Xtree{3} =  Xroty(pi/2);
model.Xtree{4} =  Xrotx(-pi/2) * Xroty(-pi/2);
model.Xtree{5} =  Xtrans([0 -p.j3z 0]);
model.Xtree{6} =  Xtrans([0 -p.j4z 0]);
model.Xtree{7} =  Xroty(pi/2) * Xrotx(pi/2);
% 左腿
model.Xtree{8} =   Xtrans([p.j5x p.j5y 0]);
model.Xtree{9} =   Xroty(pi/2);
model.Xtree{10} =  Xrotx(-pi/2) * Xroty(-pi/2);
model.Xtree{11} =  Xtrans([0 -p.j6z 0]);
model.Xtree{12} =  Xtrans([0 -p.j7z 0]);
model.Xtree{13} =  Xroty(pi/2) * Xrotx(pi/2);

% 质心在其坐标系下的表示
CoM{1} = [0 0 p.c1z];

CoM{2} = [0 0 0];
CoM{3} = [0 0 0];
CoM{4} = [p.c2x -p.c2z p.c2y];
CoM{5} = [p.c3x -p.c3z p.c3y];
CoM{6} = [0 0 0];
CoM{7} = [-p.c4z p.c4y p.c4x];

CoM{8} =  [0 0 0];
CoM{9} =  [0 0 0];
CoM{10} = [p.c5x -p.c5z p.c5y];
CoM{11} = [p.c6x -p.c6z p.c6y];
CoM{12} = [0 0 0];
CoM{13} = [-p.c7z p.c7y p.c7x];

model.CoM = CoM;

mass = [p.m1,  0, 0, p.m2, p.m3, 0, p.m4,  0, 0, p.m5, p.m6, 0, p.m7];
model.mass = mass;

Icm{1} = diag([p.I1x, p.I1y, p.I1z]);

Icm{2} = diag([0, 0, 0]);
Icm{3} = diag([0, 0, 0]);
Icm{4} = diag([p.I2x, p.I2z, p.I2y]);
Icm{5} = diag([p.I3x, p.I3z, p.I3y]);
Icm{6} = diag([0, 0, 0]);
Icm{7} = diag([p.I4z, p.I4y, p.I4x]);

Icm{8} =  diag([0, 0, 0]);
Icm{9} =  diag([0, 0, 0]);
Icm{10} = diag([p.I5x, p.I5z, p.I5y]);
Icm{11} = diag([p.I6x, p.I6z, p.I6y]);
Icm{12} = diag([0, 0, 0]);
Icm{13} = diag([p.I7z, p.I7y, p.I7x]);

model.Icm = Icm;

for i = 1:nb
    model.I{i} = mcI(mass(i), CoM{i}, Icm{i});    
end

end

function p = getModelParams()
p.c1x = 0;      p.c1y = 0;      p.c1z = 0.3;

p.c2x = 0.014;  p.c2y = -0.015; p.c2z = -0.106;
p.c3x = 0;      p.c3y = 0;      p.c3z = -0.092;
p.c4x = 0.015;  p.c4y = -0.010; p.c4z = -0.075;

p.c5x = 0.014;  p.c5y =  0.015; p.c5z = -0.106;
p.c6x = 0;      p.c6y = 0;      p.c6z = -0.092;
p.c7x = 0.015;  p.c7y =  0.010; p.c7z = -0.075;

p.j2x = 0.042;  p.j2y = -0.08;  p.j2z = 0;
p.j3x = 0;      p.j3y = 0;      p.j3z = -0.33;
p.j4x = 0;      p.j4y = 0;      p.j4z = -0.321;

p.j5x = 0.042;  p.j5y = 0.08;   p.j5z = 0;
p.j6x = 0;      p.j6y = 0;      p.j6z = -0.33;
p.j7x = 0;      p.j7y = 0;      p.j7z = -0.321;

p.m1 = 22.5;    p.g = 9.8;
p.m2 = 5.8;     p.m3 = 1.8;     p.m4 = 1.15;
p.m5 = 5.8;     p.m6 = 1.8;     p.m7 = 1.15;

p.I1x = 0.61;       p.I1y = 0.38;       p.I1z = 0.3;
p.I2x = 0.067;      p.I2y = 0.063;      p.I2z = 0.014;    
p.I3x = 0.03;       p.I3y = 0.03;       p.I3z = 0.0002;    
p.I4x = 0.0017;     p.I4y = 0.0047;     p.I4z = 0.0052;
p.I5x = 0.067;      p.I5y = 0.063;      p.I5z = 0.014;    
p.I6x = 0.03;       p.I6y = 0.03;       p.I6z = 0.0002;    
p.I7x = 0.0017;     p.I7y = 0.0047;     p.I7z = 0.0052;
end