function gm11 = GM11_optimize1(X,td)

% GM11_dome������������GM(1,1)ģ�ͣ�����Ԥ�����
% ����ԭʼ����X���Ǹ�,td��ʾԤ��δ��td������,select��ʾ�Ƿ�ʹ�û��ڱ���ֵ�Ż��Ľ����㷨ģ��
% �������gm11��ʾһ���ṹ�壬����ģ�ͱ���Ĳ�����ģ�;��ȼ���ֵ���Լ�Ԥ��ֵ

    %% ���ݵ�Ԥ����
    format short
    n = length(X);     %���ݵĸ�����������
    Ago = cumsum(X);    %һ���ۼ�
    Z = (Ago(2:end) - Ago(1:end-1))./...
        (log(Ago(2:end)) - log(Ago(1:end-1)));
        
    %% ����B��Yn��������С���˷����a��u
    B = [-Z;ones(1,n-1)]';
    Yn = X(2:n);
    
    So1 = (B'*B)\(B'*Yn');   %��С���˽�
    a = So1(1);    %��չϵ��
    u = So1(2);    %�Ҷ�������
    
    %% ����GM(1,1)�׻���ʽ��һ��һԪ΢�ַ���
    preData1 = [X(1),(1 - exp(a))*(Ago(1) - u/a)*exp(-a*(1:n+td-1))];

    %% ���ӻ�
    t = 1:n;
    plot(t,X,'ko--','MarkerFaceColor','k');     %ԭ����
    hold on
    grid on
    plot(t,preData1(1:n),'b*-','LineWidth',1.5);    %Ԥ��
    plot(n:n+td,preData1(n:n+td),'r*-','LineWidth',1.5);
    title('���ڱ���ֵ�Ż���GM(1,1)�˿�Ԥ��')
    legend('ʵ��2011-2020�����˿�' ,' Ԥ��2011-2020�����˿�',...
        'δ��10��Ԥ���˿�','LOcation','best')
    legend('boxoff')
    year=2010:2:2030;
    set(gca,'XTickLabel',{year})
    xlabel('���')
    ylabel('�˿���(��)')
    
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
    gm11.model = '���ڱ���ֵ���Ľ���GM(1,1)ģ��';
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