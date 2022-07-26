function [x_k , P_k] = EKF(u, z, dt)
%%brief info：联合仿真可以仅参考这个function
% u:IMUb系下部分参数输入: x方向加速度， y方向加速度， 绕z轴角速度（3维）
% z:GPS在全局坐标系下测量值：x位置， y位置（2维）
% bias:模型原理上目前仅考虑游走噪声
% xk：过程状态：全局坐标系下的x位置， y位置， x方向速度， y方向速度，绕z轴角度（5维）
% 可调矩阵：两个持久变量xk、P的初始化值和过程噪声以及测量噪声矩阵数值

%% 初始化X,P
persistent xk
    if isempty(xk)
        xk = zeros(5, 1);
    end
persistent P
    if isempty(P)
        P = diag([0.1, 0.1, 0.2, 0.2, deg2rad(5.0)]);
    end
    
    
%% 噪声矩阵设计
%过程噪声
Q = diag([0.25  % variance of ax
    0.25 % variance of ay
    deg2rad(5.0)  % variance of wz
    ])^(2);

%测量噪声
R = diag([
    0.2 % GPS X量测噪声
    0.2 % GPS Y量测噪声
    ])^(2);


%% 运算系数矩阵
ax= u(1);
ay = u(2);
yaw = xk(5);

A =  [1, 0, dt, 0, 0;
    0, 1, 0, dt, 0;
    0, 0, 1, 0, 0;
    0, 0, 0, 1, 0;
    0, 0, 0, 0, 1];

B = [0.5*dt^2*cos(yaw), -0.5*dt^2*sin(yaw), 0;
    0.5*dt^2*sin(yaw), 0.5*dt^2*cos(yaw), 0;
    dt*cos(yaw), -dt*sin(yaw), 0;
    dt*sin(yaw), dt*cos(yaw), 0;
    0, 0, dt];

H = [1, 0, 0, 0, 0;  0, 1, 0, 0, 0];

JF_k_1 = [1, 0, dt, 0, 0.5*dt^2*(-ax*sin(yaw)-ay*cos(yaw));
    0, 1, 0, dt, 0.5*dt^2*(ax*cos(yaw)-ay*sin(yaw));
    0, 0, 1, 0, dt*(-ax*sin(yaw)-ay*cos(yaw));
    0, 0, 0, 1, dt*(ax*cos(yaw)-ay*sin(yaw));
    0, 0, 0, 0, 1;];

JG_k_1 = [0.5*dt^2*(-cos(yaw)), 0.5*dt^2*(sin(yaw)), 0;
    0.5*dt^2*(-sin(yaw)), 0.5*dt^2*(-cos(yaw)), 0;
    dt*(-cos(yaw)), dt*(sin(yaw)), 0;
    dt*(-sin(yaw)), dt*(-cos(yaw)), 0;
    0, 0, -dt;];

JH_k = [1, 0, 0, 0, 0 ; 0, 1, 0, 0, 0];

JL_k = [1, 0; 0, 1];


%% EKF

%predict
xk = A*xk + B*u;
P = JF_k_1*P*JF_k_1' + JG_k_1*Q*JG_k_1';
% correct
S = JH_k*P*JH_k' + JL_k*R*JL_k';
K = P*JH_k'* S^(-1);
y = z - H*xk;
xk = xk + K*y;
%update
P = (eye(length(xk)) - K*JH_k) * P;


%% output
x_k = xk;
P_k = P;
