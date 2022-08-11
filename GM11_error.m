function [gm11,gme11] = GM11_error(X,td,select)
%GM(1,1)残差修正模型
    n = length(X);      %原数据长度
    gm11 = GM11_demo(X,td,select);      %对原数据建模
    F = gm11.Predict_Value;     %获取预测值
    Error = X - F(1:n);
    s = sign(Error);    %根据误差判断符号
    k0 = 1;
    for i =length(s):-1:2
        if s(i) == s(i-1)
            k0 = k0 + 1;
        else
            break
        end
    end
    if k0 >= 4      %获取可建模的残差尾部
        k0 = n- k0 + 1;
        Error0 = abs(Error(k0:end));
        Error_demo = GM11_demo(Error0,0,1);     %对残差尾部Error0进行建模 ,注意此处会在GM11_demo绘制出关于残差尾部Error0的GM(1,1)的曲线
        ea = Error_demo.Coeff_a;
        eu = Error_demo.Coeff_u;
        ErrorF = zeros(1,n+td);
        for k = k0:n+td
            ErrorF(k) = (Error0(1)-eu/ea)/exp(ea*(k-k0));      %对应的时间响应函数的离散形式
        end
        if s(end) == -1
            preData2 = F(k0:end) + ea*ErrorF(k0:end);
        elseif s(end) == 1
            preData2 = F(k0:end) - ea*ErrorF(k0:end);
        end
        preData2 = [F(1:k0-1),preData2];
        
        %% 可视化
        subplot(1,2,2)
        t2 = 1:n;
        plot(t2,X,':','MarkerFaceColor','y');     %原数据
        hold on
        grid on
        plot(t2,preData2(1:n),'g*-','LineWidth',1);    %修正预测
        plot(n:n+td,preData2(n:n+td),'m*-','LineWidth',1);
        title('GM(1,1) error_correction——Original VS Current And Futrue Predict')
        legend('Original Data' ,' Predict Current Value2',...
            'Predict Future Value2','LOcation','best')
        legend('boxoff')
        xlabel('Time')
        ylabel('Value')
        
        %% 修正模型精度检验
        Error2 = X - preData2(1:n);    %误差
        Abs_Error2 = abs(Error2);  %绝对误差
        RelE2 = Abs_Error2./X;    %相对误差
        RelMean2 = mean(RelE2);   %相对误差均值
        S1 = std(X,1);  %原数据标准差
        S2 = std(Abs_Error2(2:end),1);   %误差标准差
        C2 = S2/S1;  %后验方差比
        P2 = sum(abs(Abs_Error2 - mean(Abs_Error2)) < 0.6745*S1)/n;    %小误差概率
        R_k2 = (min(Abs_Error2(2:end)) + 0.5*max(Abs_Error2(2:end)))./...
            (Abs_Error2(2:end) + 0.5*max(Abs_Error2(2:end)));
        R2 = sum(R_k2)/(n-1);     %关联度检验
        lamda_error = X(1:n-1)./X(2:n);   %级比检验
        descend_lamda_error = sort(lamda_error,'descend');
        if descend_lamda_error(end)<exp(2/(n+1)) && descend_lamda_error(end)>exp(-2/(n+1))
            disp('误差序列通过级比检验')
        else
            disp('误差序列没有通过级比检验，不能进行灰度预测')
        end
        %% 组合输出参数各种结果
        gme11.Type = '误差修正GM(1,1)模型';
        gme11.Coeff_a = gm11.Coeff_a;
        gme11.Coeff_u = gm11.Coeff_u;
        gme11.Coeff_ea = ea;
        gme11.Coeff_eu = eu;
        gme11.Predict_Value2 = preData2;
        gme11.Error2 = Error2;
        gme11.Abs_Error2 = Abs_Error2;
        gme11.Rel_Error2 = RelE2;
        gme11.Rel_Error_Mean2 = RelMean2;
        gme11.C2 = C2;
        gme11.P2 = P2;
        gme11.R2 = R2;
        gme11.lamda_error = lamda_error;
    else
        disp('不能进行残差尾部修正')
        gme11 = [];
    end
end


