function defaultColor(a)
%DEFAULECOLOR 显示默认的颜色
% a指定第几个
%20240807制定第一个
CM=[0.0941    0.1098    0.2627
    0.0431    0.3725    0.7451
    0.4588    0.6667    0.7451
    0.9451    0.9255    0.9235
    0.8157    0.5451    0.4510
    0.6510    0.1353    0.1431
    0.2353    0.0353    0.0706];
CMX=linspace(0,1,size(CM,1));
CMXX=linspace(0,1,256)';
CM=[interp1(CMX,CM(:,1),CMXX,'pchip'), ...
    interp1(CMX,CM(:,2),CMXX,'pchip'), ...
    interp1(CMX,CM(:,3),CMXX,'pchip')];
colormap(CM)
end

