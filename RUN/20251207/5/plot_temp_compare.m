clc; clear all
% addpath('RES')
% load('res_N=20_rho=2.00_beta=1.0_epsilon=0.10_omega=4.36_dtScale = 1.mat')
% plot(T,X(:,1))
hold on
% load('res_N=20_rho=2.00_beta=1.0_epsilon=0.10_omega=4.36_dtScale = 0.5.mat')
% plot(T,X(:,1))
% load('res_N=20_rho=2.00_beta=1.0_epsilon=0.10_omega=4.36_dtScale = 0.1.mat')
% plot(T,X(:,1))
load('res_N=20_rho=2.00_beta=1.0_epsilon=0.10_omega=4.36_dtScale = 0.05.mat')
plot(T,X(:,1))
load('res_N=20_rho=2.00_beta=1.0_epsilon=0.10_omega=4.36_dtScale = 0.02.mat')
plot(T,X(:,1))
