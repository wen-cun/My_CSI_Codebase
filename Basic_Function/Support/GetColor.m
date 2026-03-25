function Color= GetColor(ColorNum,ColorType)
%这个脚本用来选取指定个数的颜色
%输入ColorNum是颜色的个数，ColorType为指定的类型
%输出Color为颜色向量
%更新于2025年12月30日
switch ColorNum
    case 1
        switch ColorType
            case 1
                Color=[065 130 164]/255;%天蓝
            case 2
                Color=[166 064 054]/255;%朱红
            case 3
                Color=[250 215 151]/255;%树叶黄
            case 4
                Color=[106 97 141]/255;%暗紫
            otherwise
                error('没有指定类型的颜色');
        end
    case 2
        switch ColorType
            case 1
                Color=[065 130 164;166 064 054]/255;%雾霾蓝与朱红
            case 2
                Color=[214 145 173;250 215 151]/255;%胭脂红和树叶黄
            otherwise
                error('没有指定类型的颜色');
        end
    case 3
        switch ColorType
            case 1
                Color=[166 064 054;70 130 80;065 130 164]/255;
            case 2
                Color=[106 97 141;214 125 173;250 216 151]/255;
            otherwise
                error('没有指定类型的颜色');
        end
    case 4
        switch ColorType
            case 1
                Color=[166 064 054;70 130 80;065 130 164;053 078 107]/255;
            case 2
                Color=['#2C363F';'#2F6665';'#E5BD47';'#DCDCDD'];
            case 3
                Color=[106 97 141;214 125 173;250 216 151;217 217 217]/255;
            otherwise
                error('没有指定类型的颜色');
        end
    case 5
        switch ColorType
            case 1
                Color=[0.8745 0.8902    0.8078;0.7098    0.7843    0.3804;0.5412    0.6549    0.5373;...
                      0.7961    0.5216    0.4510;0.5922    0.3804    0.3255];
            case 2
               Color=[106 97 141;168 136 175;214 145 173;248 162 153;250 215 151]; 
            otherwise
                error('没有指定类型的颜色');
        end
    case 6
        switch ColorType
            case 1
                Color=[0.8627    0.7608    0.4784;0.6902    0.7255    0.7451;0.3882    0.3765    0.3725;...
                       0.5961    0.3686    0.3608;0.6824    0.7490    0.6588;0.9451    0.6078    0.2039];
            case 2
                Color=[106 97 141;168 136 175;214 145 173;248 162 153;250 215 151;216 216 217]; 
            otherwise
                error('没有指定类型的颜色');
        end        
    otherwise
        error('没有包含该数目的的颜色');
end
end

