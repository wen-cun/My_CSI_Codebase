%% 这个脚本是gen_lightsource函数的demo
clear;
close all;
clc;
%% 定义输入变量和输出变量
type=1;%输入变量：光源类型
vk0=1./(0.3:0.001:0.8)';%输入变量：查询点波数

vsk=nan*ones(size(vk0));%输出变量：查询点强度
%%
switch type
    case 1
        vsk=expe_usual(vk0);
    otherwise
        error('光源类型不支持！！！');
end
%% expe_usual 函数
function vSk=expe_usual(vk0)

end