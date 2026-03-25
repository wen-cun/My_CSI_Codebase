function defaultAxes(a)
%DEFAULTAXES该函数用来改变绘图的坐标区的参数
%如果a为2代表二维绘图，为3代表3为绘图
if abs(a-2)<1e-16
    ax=gca;
    hold on;
    box on;
    ax.XGrid='on';
    ax.YGrid='on';
    ax.XMinorTick='on';
    ax.YMinorTick='on';
    ax.LineWidth=0.8;
    ax.GridLineStyle='-.';
    ax.FontName='Cambria';
    ax.FontSize=12;
%     set(gcf,'position',[400,200,1.2800e+03,710]);%调整图大小
    hold off;
elseif abs(a-3)<1e-16
    ax=gca;
    hold on;
    box on;
    ax.XGrid='on';
    ax.YGrid='on';
    ax.ZGrid='on';
    ax.XMinorTick='on';
    ax.YMinorTick='on';
    ax.ZMinorTick='on';
    ax.LineWidth=.8;
    ax.GridLineStyle='-.';
    ax.FontName='Cambria';
    ax.FontSize=12;
%     set(gcf,'position',[400,200,1.2800e+03,710]);%调整图大小
    hold off;
else
    error('绘图指定的纬度不对！！！');
end
end

