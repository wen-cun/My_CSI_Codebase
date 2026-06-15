function signal=CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol)
% 这个函数用来计算点式CSI信号

% 输入NA为系统的NA
% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入z_scan为扫描序列, double型1×K维向量，K为图像帧数
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量
% 输入system_pol为系统的偏振状态

%输出signal为CSI信号， double 1×K维向量

%     例：
%     z_peri=0.0719;%PZT单次移动间隔，um制
%     N_half = round(10 / z_peri); % 计算单侧点数
%     z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vk0=1./1.1:0.01:1./0.3;
%     NA=0.15；
%     theta_array=0:0.01:asin(NA);
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     system_pol='unpolar';
%     signal=CSIPointSingalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol)

switch system_pol
     case 'ideal'  %理想情况
        A1 = 1; %样品振幅反射率
        A2 = 1; %参考镜振幅反射率
        if NA<eps
            %0NA
            signal=CSIPointIdeal0NA(vk0,vsk,A1,A2,z_scan,sample_dis);
        else
            signal=CSIPointIdealNA(vk0,vsk,theta_array,A1,A2,z_scan,sample_dis);
        end
    case 'unpolar' %自然光/非偏振模式
        if NA<eps
            %0NA
            signal=CSIPointUnpolar0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis);
        else
            signal=CSIPointUnpolarNA(vk0,vsk,theta_array,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis);
        end
    otherwise
        error('未包含该偏振状态,仅支持: ideal,unpolar');
end
end

