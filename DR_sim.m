function [x_Rk,x_DR, x_GPS,u_IMU,x_t] = DR(ax,ay,wz,dt)
%%brief info：IMU无误差情况下的航迹推算、IMU误差情况下的航迹推算、理想IMU加噪、理想GPS加噪
% x、x_d：全局坐标系下的状态记录矩阵
% z:GPS在全局坐标系下测量值：x位置， y位置（2维）
% bias:模型原理上目前仅考虑游走噪声
% xk：过程状态：全局坐标系下的x位置， y位置， x方向速度， y方向速度，绕z轴角度（5维）
% 可调矩阵：两个持久变量xk、P的初始化值和过程噪声以及测量噪声矩阵数值
persistent x
    if isempty(x)
        x = zeros(5, 1);
    end
persistent x_d
    if isempty(x_d)
        x_d = zeros(5, 1);
    end
u = [ax,ay,wz]';
%% 运算系数矩阵
A =  [1, 0, dt, 0, 0;
    0, 1, 0, dt, 0;
    0, 0, 1, 0, 0;
    0, 0, 0, 1, 0;
    0, 0, 0, 0, 1];
yaw_RK = x(5);
B_RK = [0.5*dt^2*cos(yaw_RK), -0.5*dt^2*sin(yaw_RK), 0;
    0.5*dt^2*sin(yaw_RK), 0.5*dt^2*cos(yaw_RK), 0;
    dt*cos(yaw_RK), -dt*sin(yaw_RK), 0;
    dt*sin(yaw_RK), dt*cos(yaw_RK), 0;
    0, 0, dt];
yaw_DR = x_d(5);
B_DR = [0.5*dt^2*cos(yaw_DR), -0.5*dt^2*sin(yaw_DR), 0;
    0.5*dt^2*sin(yaw_DR), 0.5*dt^2*cos(yaw_DR), 0;
    dt*cos(yaw_DR), -dt*sin(yaw_DR), 0;
    dt*sin(yaw_DR), dt*cos(yaw_DR), 0;
    0, 0, dt];
H = [1, 0, 0, 0, 0;  0, 1, 0, 0, 0];

%% 理想航迹推算
x= A*x + B_RK*u;
x_Rk = x;

%% 加噪
GPS_NOISE = diag([0.3, 0.3]) ^(2);
INPUT_NOISE  = diag([0.1, 0.1, deg2rad(10.0)]) ^(2);
x_GPS = H*x + GPS_NOISE * randn(2,1);
u_IMU = u + INPUT_NOISE * randn(3,1);

%%误差航迹推算 
x_d = A*x_d + B_DR * u_IMU;
x_DR = x_d;
x_t = abs(x_DR-x_Rk);
