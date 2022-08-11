function gmV = GMVerhulst_demo(X,pre_number,select)
% GMVerhulst_demo函数用来建立GMGMVerhulst模型，并做预测分析
% 输入原始数据X，非负,pre_number表示预期未来td期数据,select表示是否使用基于背景值优化改进的算法模型
% 输出参数gmV表示一个结构体，包括模型本身的参数和模型精度检验值，以及预测值

    %% 数据的预处理
    format short
    n = length(X);     %数据的个数，数据量
%     Ago = cumsum(X);    %一次累加
    if select == 1
        Z = (X(1:n-1) + X(2:n))/2;  %邻均值序列
    else
         Z = (X(2:end) - X(1:end-1))./...
            (log(X(2:end)) - log(X(1:end-1)));
    end
    %% 构造B和Y，并用最小二乘法求解a和b
    B = [-Z;Z.^2]';
    Y = (diff(X))';
    
    Sol = (B'*B)\(B'*Y);
    a = Sol(1);
    b = Sol(2);
    
    %% 构造时间相应绪论，累减还原
    F = [X(1),a*X(1)./(b*X(1)+(a-b*X(1))*exp(a*(1:n+pre_number-1)))];
%     preData1 = [F(1),F(2:end) - F(1:end-1)];
    preData1 = F;
    
     %% 可视化
    subplot(1,2,2)
    t = 1:n;
    plot(t,X,'ko--','MarkerFaceColor','k');     %原数据
    hold on
    grid on
    plot(t,preData1(1:n),'b*-','LineWidth',1.5);    %预测
    plot(n:n+pre_number,preData1(n:n+pre_number),'r*-','LineWidth',1.5);
    title('GMVerhulst——Original VS Current And Futrue Predict')
    legend('Original Data' ,' Predict Current Value',...
        'Predict Future Value','LOcation','best')
    legend('boxoff')
    xlabel('Time')
    ylabel('Value')
    
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
    if select == 1
        gmV.model = ('传统GMVerhulst模型');
    else
        gmV.model = ('背景值优化的GMVerhulst模型');
    end
    gmV.Coeff_a = a;
    gmV.Coeff_b = b;
    gmV.Predict_Value = preData1;
    gmV.Error = Error;
    gmV.Abs_Error = Abs_Error;
    gmV.Rel_Error = RelE;
    gmV.Rel_Error_Mean = RelMean;
    gmV.C = C;
    gmV.P = P;
    gmV.R = R;
    gmV.lamda = lamda;
    
end