function [M_te,M_tm] = CalcTransferMatrix(vk0,theta,index,sample_stru)
% 这个函数用来计算膜层结构的传递矩阵

% 输入 vk0 为波数，double型1×1变量
% 输入 theta 为角度，double型1×1变量
% 输入 index 为每层结构的单波长复折射率，指定为layernum×2型cell元胞数组，每个子胞内都为double型1×1变量
% 输入 sample_stru 为该膜层结构，指定为layernum×2型cell元胞数字，第一列子胞为每层膜的介质，第二列子胞为每层膜的厚度

% 输出 M_te 为TE偏振下，该膜层结构的总传递矩阵，complex double型2×2型变量
% 输出 M_tm 为TM偏振下，该膜层结构的总传递矩阵，complex double型2×2型变量

%     例：
%        vk0=1/0.5;
%        theta=asin(0.15);
%        index{1,1}=1.5；
%        index{1,2}=0.02；
%        index{2,1}=1.75;
%        index{2,2}=0.1;
%        sample_stru={'SiO2',5;'Si',inf};
%        [M_te,M_tm]=CalcTransferMatrix(vk0,theta,index,samp_stru)

layernum=size(sample_stru,1); %介质的总个数
M_te=eye(2);
M_tm=eye(2);
if layernum==1
    %直接返回单位矩阵
else
    for kk=1:layernum-1 %内层循环计算每层介质传递矩阵
                %计算第kk层薄膜传递矩阵
                N = index{kk,1}+1i*index{kk,2};%介质复折射率
                q = Calcq(N,theta);
                p1_te = q; %介质的TE场导纳，这里要用到折射率的消光系数
                delta = 2*pi*vk0*sample_stru{kk,2}*q; %介质内部的TE场delta
                M1_te = CalcSingleTransferMatrix(p1_te,delta); %该介质的传递矩阵 
                M_te=M_te*M1_te;%更新TE场传递矩阵
                
                p1_tm = N^2/q; %介质的TM场导纳，这里要用到折射率的消光系数
                M1_tm = CalcSingleTransferMatrix(p1_tm,delta); %该介质的传递矩阵 
                M_tm=M_tm*M1_tm;%更新TM场传递矩阵
   end
end
end