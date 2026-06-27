function signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol)
% 这个函数用来计算点式SDI信号

% 输入NA为系统的NA
% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量
% 输入system_pol为系统的偏振状态

%输出signal为CSI信号， double 型N×1维向量，等波长空间光谱分布


%     例：
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vsk=vsk./(max(abs(vsk)));
%     NA = 0.3; %系统NA
%     theta_max=asin(NA); %最大NA对应的空气中光线角度theta
%     theta_peri=0.01; %角度theta的采样周期
%     theta_array = 0:theta_peri:theta_max; %theta采样数组
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     system_pol='ideal';
%     signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol))

switch system_pol
    case 'ideal'  %理想情况
        A1 = 1; %样品振幅反射率
        A2 = 1; %参考镜振幅反射率
        if NA<eps
            %0NA
            signal = SDIPointIdeal0NA(vk0,vsk,A1,A2,sample_dis);
        else
            signal = SDIPointIdealNA(vk0,vsk,theta_array,A1,A2,sample_dis);
        end
    case 'unpolar' %自然光/非偏振模式
        if NA<eps
            %0NA
            signal=SDIPointUnpolar0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,sample_dis);
        else
            signal=SDIPointUnpolarNA(vk0,vsk,theta_array,r_Se,r_Sm,r_Me,r_Mm,sample_dis);
        end
    otherwise
        error('未包含该偏振状态,仅支持: ideal,unpolar');
end
end
