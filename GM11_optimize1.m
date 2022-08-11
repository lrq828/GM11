function gm11 = GM11_optimize1(X,td)

% GM11_dome函数用来建立GM(1,1)模型，并做预测分析
% 输入原始数据X，非负,td表示预期未来td期数据,select表示是否使用基于背景值优化改进的算法模型
% 输出参数gm11表示一个结构体，包括模型本身的参数和模型精度检验值，以及预测值

    %% 数据的预处理
    format short
    n = length(X);     %数据的个数，数据量
    Ago = cumsum(X);    %一次累加
    Z = (Ago(2:end) - Ago(1:end-1))./...
        (log(Ago(2:end)) - log(Ago(1:end-1)));
        
    %% 构造B和Yn，并用最小二乘法求解a和u
    B = [-Z;ones(1,n-1)]';
    Yn = X(2:n);
    
    So1 = (B'*B)\(B'*Yn');   %最小二乘解
    a = So1(1);    %发展系数
    u = So1(2);    %灰度作用量
    
    %% 建立GM(1,1)白化形式的一阶一元微分方程
    preData1 = [X(1),(1 - exp(a))*(Ago(1) - u/a)*exp(-a*(1:n+td-1))];

    %% 可视化
    t = 1:n;
    plot(t,X,'ko--','MarkerFaceColor','k');     %原数据
    hold on
    grid on
    plot(t,preData1(1:n),'b*-','LineWidth',1.5);    %预测
    plot(n:n+td,preData1(n:n+td),'r*-','LineWidth',1.5);
    title('基于背景值优化的GM(1,1)人口预测')
    legend('实际2011-2020年总人口' ,' 预测2011-2020年总人口',...
        '未来10年预测人口','LOcation','best')
    legend('boxoff')
    year=2010:2:2030;
    set(gca,'XTickLabel',{year})
    xlabel('年份')
    ylabel('人口数(万)')
    
    %% 模型精度检验
    Error = X - preData1(1:n);    %误差
    Abs_Error = abs(Error);  %绝对误差
    RelE = Abs_Error./X;    %相对误差
    RelMean = mean(RelE);   %相对误差均值
    S1 = std(X,1);  %原数据标准差
    S2 = std(Abs_Error(2:end),1);   %误差标准差
    C = S2/S1;  %后验方差比
    P = sum(abs(Abs_Error - mean(Abs_Error)) < 0.6745*S1)/n;    %小误差概率
    R_k = (min(Abs_Error(2:end)) + 0.5*max(Abs_Error(2:end)))./...
        (Abs_Error(2:end) + 0.5*max(Abs_Error(2:end)));
    R = sum(R_k)/(n-1);     %关联度检验
    lamda = X(1:n-1)./X(2:n);   %级比检验
    descend_lamda = sort(lamda,'descend'); 
    if descend_lamda(end)<exp(2/(n+1)) && descend_lamda(end)>exp(-2/(n+1))
        disp('通过级比检验')
    else
        disp('没有通过级比检验，不能进行灰度预测')
    end
    
    %% 组合输出参数各种结果
    gm11.model = '基于背景值做改进的GM(1,1)模型';
    gm11.Coeff_a = a;
    gm11.Coeff_u = u;
    gm11.Predict_Value = preData1;
    gm11.Error = Error;
    gm11.Abs_Error = Abs_Error;
    gm11.Rel_Error = RelE;  
    gm11.Rel_Error_Mean = RelMean;
    gm11.C = C;
    gm11.P = P;
    gm11.R = R;
    
end