function gmV = GMVerhulst_demo(X,pre_number,select)
% GMVerhulst_demo������������GMGMVerhulstģ�ͣ�����Ԥ�����
% ����ԭʼ����X���Ǹ�,pre_number��ʾԤ��δ��td������,select��ʾ�Ƿ�ʹ�û��ڱ���ֵ�Ż��Ľ����㷨ģ��
% �������gmV��ʾһ���ṹ�壬����ģ�ͱ���Ĳ�����ģ�;��ȼ���ֵ���Լ�Ԥ��ֵ

    %% ���ݵ�Ԥ����
    format short
    n = length(X);     %���ݵĸ�����������
%     Ago = cumsum(X);    %һ���ۼ�
    if select == 1
        Z = (X(1:n-1) + X(2:n))/2;  %�ھ�ֵ����
    else
         Z = (X(2:end) - X(1:end-1))./...
            (log(X(2:end)) - log(X(1:end-1)));
    end
    %% ����B��Y��������С���˷����a��b
    B = [-Z;Z.^2]';
    Y = (diff(X))';
    
    Sol = (B'*B)\(B'*Y);
    a = Sol(1);
    b = Sol(2);
    
    %% ����ʱ����Ӧ���ۣ��ۼ���ԭ
    F = [X(1),a*X(1)./(b*X(1)+(a-b*X(1))*exp(a*(1:n+pre_number-1)))];
%     preData1 = [F(1),F(2:end) - F(1:end-1)];
    preData1 = F;
    
     %% ���ӻ�
    subplot(1,2,2)
    t = 1:n;
    plot(t,X,'ko--','MarkerFaceColor','k');     %ԭ����
    hold on
    grid on
    plot(t,preData1(1:n),'b*-','LineWidth',1.5);    %Ԥ��
    plot(n:n+pre_number,preData1(n:n+pre_number),'r*-','LineWidth',1.5);
    title('GMVerhulst����Original VS Current And Futrue Predict')
    legend('Original Data' ,' Predict Current Value',...
        'Predict Future Value','LOcation','best')
    legend('boxoff')
    xlabel('Time')
    ylabel('Value')
    
    %% ģ�;��ȼ���
    Error = X - preData1(1:n);    %���
    Abs_Error = abs(Error);  %�������
    RelE = Abs_Error./X;    %������
    RelMean = mean(RelE);   %�������ֵ
    S1 = std(X,1);  %ԭ���ݱ�׼��
    S2 = std(Abs_Error(2:end),1);   %����׼��
    C = S2/S1;  %���鷽���
    P = sum(abs(Abs_Error - mean(Abs_Error)) < 0.6745*S1)/n;    %С������
    R_k = (min(Abs_Error(2:end)) + 0.5*max(Abs_Error(2:end)))./...
        (Abs_Error(2:end) + 0.5*max(Abs_Error(2:end)));
    R = sum(R_k)/(n-1);     %�����ȼ���
    lamda = X(1:n-1)./X(2:n);   %���ȼ���
    descend_lamda = sort(lamda,'descend'); 
    if descend_lamda(end)<exp(2/(n+1)) && descend_lamda(end)>exp(-2/(n+1))
        disp('ͨ�����ȼ���')
    else
        disp('û��ͨ�����ȼ��飬���ܽ��лҶ�Ԥ��')
    end
    
    %% �������������ֽ��
    if select == 1
        gmV.model = ('��ͳGMVerhulstģ��');
    else
        gmV.model = ('����ֵ�Ż���GMVerhulstģ��');
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