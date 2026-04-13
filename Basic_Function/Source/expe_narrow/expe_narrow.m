function vsk=expe_narrow(vk0)
%生成指定查询波数点的窄带光源光谱强度
%输入vk0 为查询点的波数
%输出vsk 为查询点的强度
%  例： 
%    vk0=1./(0.3:0.001:0.8)';
%    vsk=expe_usual(vk0);

spectrm_data=readmatrix("C:\Users\lenovo\Desktop\MATLAB\My_CSI_Codebase\Basic_Function\Source\expe_narrow\LightSource.txt");%读取光谱数据

spectrm_data(:,2)=spectrm_data(:,2)-mean(spectrm_data(end-9:end,2),'omitnan');%强度列减去均值
spectrm_data(:,2)=spectrm_data(:,2)./max(spectrm_data(:,2));%归一化

spectrm_data(:,1)=1./(spectrm_data(:,1)/1000);%波长列先转为um制再取波数

vsk=interp1(spectrm_data(:,1),spectrm_data(:,2),vk0,'spline');%三次样条插值
end
