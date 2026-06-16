function vsk=expe_usual(vk0)
%生成指定查询波数点的实验常用光源光谱强度

%输入vk0 为查询点的波数,  变量类型： double型vector

%输出vsk 为查询点的强度，与输入的vk0同维度，  变量类型 double型vector

%  例： 
%    vk0=1./(0.3:0.001:0.8)';
%    vsk=expe_usual(vk0);

spectrm_data=readmatrix("C:\Users\lenovo\Desktop\MATLAB\My_CSI_Codebase\Basic_Function\Source\expe_usual\LightSource.txt");%读取光谱数据

spectrm_data(:,2)=spectrm_data(:,2)-mean(spectrm_data(end-9:end,2),'omitnan');%强度列减去均值
spectrm_data(:,2)=spectrm_data(:,2)./max(spectrm_data(:,2));%归一化
lam0=sum(spectrm_data(:,1).*spectrm_data(:,2))./sum(spectrm_data(:,2));

% Color=GetColor(2,1);
% figure();
% plot(spectrm_data(:,1),spectrm_data(:,2),'LineWidth',1.5,'Color',Color(1,:));
% hold on;

tau = 75;%超高斯窗标准差
order = 3;%超高斯窗阶数
window = exp(-((spectrm_data(:,1)-lam0).^2/(2*tau^2)).^order);

% plot(spectrm_data(:,1),window+1,'--','LineWidth',1.5,'Color',Color(2,:));

spectrm_data(:,2)=spectrm_data(:,2).*window; %超高斯窗滤波

spectrm_data(:,1)=1./(spectrm_data(:,1)/1000);%波长列先转为um制再取波数

vsk=interp1(spectrm_data(:,1),spectrm_data(:,2),vk0,'spline',0);%三次样条插值
end

