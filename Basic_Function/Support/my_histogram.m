function my_histogram(data)
%这个函数用来绘制直方图
%只需输入要显示的数据即可
% 💡 核心操作 2：开启概率密度归一化 (PDF)，方便后续叠加密度曲线
h_hist = histogram(data, 'Normalization', 'pdf', 'NumBins', 150);

% 💡 核心操作 3：去掉黑边，开启透明度，换用高级配色
h_hist.FaceColor = [0.2, 0.6, 0.8]; % 高级灰蓝色 (类似浅海蓝)
h_hist.EdgeColor = 'none';          % 【绝对关键】去掉默认的丑陋黑边！
h_hist.FaceAlpha = 0.75;
defaultAxes(2);
end

