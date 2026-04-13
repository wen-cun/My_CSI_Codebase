%% 这个脚本是gen_lightsource函数的demo
clear;
close all;
clc;
%% 定义输入变量和输出变量
type=4;%输入变量：光源类型
vk0=1./(0.3:0.001:0.8)';%输入变量：查询点波数

% vsk=nan*ones(size(vk0));%输出变量：查询点强度
%% 返回查询强度
vsk=gen_lightsource(vk0,type);
%% 绘图观看效果
figure(1);
plot(vk0,vsk,'LineWidth',1.5,'Color',GetColor(1,1))
defaultAxes(2);
xlabel('k/$\mu$ m','Interpreter','latex');