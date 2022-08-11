function [gm11,gme11] = GM11_error(X,td,select)
%GM(1,1)�в�����ģ��
    n = length(X);      %ԭ���ݳ���
    gm11 = GM11_demo(X,td,select);      %��ԭ���ݽ�ģ
    F = gm11.Predict_Value;     %��ȡԤ��ֵ
    Error = X - F(1:n);
    s = sign(Error);    %��������жϷ���
    k0 = 1;
    for i =length(s):-1:2
        if s(i) == s(i-1)
            k0 = k0 + 1;
        else
            break
        end
    end
    if k0 >= 4      %��ȡ�ɽ�ģ�Ĳв�β��
        k0 = n- k0 + 1;
        Error0 = abs(Error(k0:end));
        Error_demo = GM11_demo(Error0,0,1);     %�Բв�β��Error0���н�ģ ,ע��˴�����GM11_demo���Ƴ����ڲв�β��Error0��GM(1,1)������
        ea = Error_demo.Coeff_a;
        eu = Error_demo.Coeff_u;
        ErrorF = zeros(1,n+td);
        for k = k0:n+td
            ErrorF(k) = (Error0(1)-eu/ea)/exp(ea*(k-k0));      %��Ӧ��ʱ����Ӧ��������ɢ��ʽ
        end
        if s(end) == -1
            preData2 = F(k0:end) + ea*ErrorF(k0:end);
        elseif s(end) == 1
            preData2 = F(k0:end) - ea*ErrorF(k0:end);
        end
        preData2 = [F(1:k0-1),preData2];
        
        %% ���ӻ�
        subplot(1,2,2)
        t2 = 1:n;
        plot(t2,X,':','MarkerFaceColor','y');     %ԭ����
        hold on
        grid on
        plot(t2,preData2(1:n),'g*-','LineWidth',1);    %����Ԥ��
        plot(n:n+td,preData2(n:n+td),'m*-','LineWidth',1);
        title('GM(1,1) error_correction����Original VS Current And Futrue Predict')
        legend('Original Data' ,' Predict Current Value2',...
            'Predict Future Value2','LOcation','best')
        legend('boxoff')
        xlabel('Time')
        ylabel('Value')
        
        %% ����ģ�;��ȼ���
        Error2 = X - preData2(1:n);    %���
        Abs_Error2 = abs(Error2);  %�������
        RelE2 = Abs_Error2./X;    %������
        RelMean2 = mean(RelE2);   %�������ֵ
        S1 = std(X,1);  %ԭ���ݱ�׼��
        S2 = std(Abs_Error2(2:end),1);   %����׼��
        C2 = S2/S1;  %���鷽���
        P2 = sum(abs(Abs_Error2 - mean(Abs_Error2)) < 0.6745*S1)/n;    %С������
        R_k2 = (min(Abs_Error2(2:end)) + 0.5*max(Abs_Error2(2:end)))./...
            (Abs_Error2(2:end) + 0.5*max(Abs_Error2(2:end)));
        R2 = sum(R_k2)/(n-1);     %�����ȼ���
        lamda_error = X(1:n-1)./X(2:n);   %���ȼ���
        descend_lamda_error = sort(lamda_error,'descend');
        if descend_lamda_error(end)<exp(2/(n+1)) && descend_lamda_error(end)>exp(-2/(n+1))
            disp('�������ͨ�����ȼ���')
        else
            disp('�������û��ͨ�����ȼ��飬���ܽ��лҶ�Ԥ��')
        end
        %% �������������ֽ��
        gme11.Type = '�������GM(1,1)ģ��';
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
        disp('���ܽ��вв�β������')
        gme11 = [];
    end
end


