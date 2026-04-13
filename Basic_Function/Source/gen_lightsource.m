function vsk = gen_lightsource(vk0,type)
%这个函数用来生成指定类型，指定波数的光源光谱强度

%输入vk0为查询点的波数，  变量类型：double型vector
%输入type为光源的类型，  变量类型:double

%输出vsk为查询点的强度，与输入的vk0同维度，  变量类型:doubkle型vector

%  例
%  vk0=1./(0.3:0.001:0.8)';
%  vsk=gen_lightsource(vk0,1);
%  返回指定第1中光源，在vk0位置的光谱强度

switch type
    case 1
        vsk=expe_usual(vk0);
    case 2
        vsk=expe_many_peaks(vk0);
    case 3
        vsk=expe_narrow(vk0);
    case 4
        %模拟高斯光源
        mu=1/0.55;
        sigma=0.05;
        SNR=inf;
        vsk=simu_gauss(vk0,mu,sigma,SNR);
    otherwise
        error('光源类型不支持！！！');
end

end

