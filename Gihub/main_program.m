%% main����
% ������ʼ���������ݲɼ������ع������������
clc
clear all

%% =====================������ʼ��========================
global num_nodes;
global h ;                    %����ʱ��
h = 0.001;%0.01
num_meas =12;%45                 %�۲����(<34)64
d_dimension = 3;              %ÿ�����ӵ�ά��
Data_Collection =0;           %���ݲɼ�����
iterationTime =20000;    % 80000(5kaishi) %200000;%��������500000
ref_point=5000;         % 30000   %�ο����ѡ��(80000)
history=[];                   %��¼�켣�����е�
TPR1=[];
TPR2=[];
FPR1=[];
FPR2=[];
ACC1=[];
ACC2=[];
AUroc1=[];
AUroc2=[];
Nmse1 = [];
Nmse2 = [];
inner_mat=[0,1,0;0,0,1;1,0,0];%��������
%����һ����
Lamda = 0.1;                  %L1������ϵ��
%����������
p1=0.2;                       %0�ĸ���
%% =====================���ݲɼ�=====================
%������ȡ1��
for t=1:1
    %      p2 = t*0.1;
    %      num_meas=num_meas+5;
    str1='SW-DATA';str2=num2str(t);
    filename = [str1,str2];
    K = load(filename);%��ȡ����ṹ
    Adj = K.a;
    % K = importdata('SW-DATA1.mat');%��ȡ����ṹ
    % Adj = full(K);
    % Adj = importdata('swnet20-2.mat');%��ȡ��
    num_nodes = size(Adj,1);
    %������ȡ2��
    % K = load('A.txt');%��ȡ����ṹ
    %  num_nodes = max(max(K(:,1)),max(K(:,1)));
    K_Bin=(d_dimension+1)*num_nodes+1;%��Чֵ����ʼ
    n = d_dimension*(num_nodes+1)+1;%����һ��������ά��
    % Adj = zeros(num_nodes,num_nodes);%�ڽӾ���
    % for i=1:size(K)
    %     a=K(i,1);
    %     b=K(i,2);
    %     Adj(a,b)=K(i,3);
    % end
    outer_mat = Adj - diag(sum(Adj));
    if Data_Collection == 1
        Point_Init = random('Normal',0,1,3,num_nodes);
        Network_parameters.Point_Init = Point_Init;
        Network_parameters.iterationTime = iterationTime;
        Network_parameters.N = num_nodes;
        Network_parameters.AA = inner_mat;
        Network_parameters.C = outer_mat;
        history = Net_Generating(Network_parameters);
        str1='SW-DATA12';str2=num2str(t);
        filename = [str1,str2];
        save(filename,'history');
        %      save('SW-DATA11.mat','history');
    else
        str1='SW-DATA12';str2=num2str(t);
        filename = [str1,str2];
        his = load(filename);
        history = his.history;
        %        history = importdata('SW-DATA11.mat');
%         scatter3(history(1,1,1:20000),history(2,1,1:20000),history(3,1,1:20000));%���켣ͼ
    end
    
    %% =====================ѡ������뾶Ҫ������ݣ��ع��õ�C����=====================
    for P_S = 1:num_nodes %�Ե�P_S���ڵ�Ϊ���ģ�����ڵ�P_S�������ڵ�����ӹ�ϵ
        j=1;
        Rx=history(1,P_S,ref_point);
        Ry=history(2,P_S,ref_point);
        Rz=history(3,P_S,ref_point);%�ο���
        for i=1:iterationTime
            V_x=history(1,P_S,i)-Rx;
            V_y=history(2,P_S,i)-Ry;
            V_z=history(3,P_S,i)-Rz;
            %         if (power(V_x,2)+power(V_y,2)+power(V_z,2)) ~= 0
            if (i~=ref_point && i>1);
                if ((power(V_x,2)+power(V_y,2)+power(V_z,2)) <= 0.5^3)%0.8%0.6
                    Success(1:3,:,j)=history(1:3,:,i);
                    order_S(j) = i;
                    j=j+1;
                end
            end
        end
        %% ����Ԥ����ʹ֮����y=\Phi*x����ʽ--�˴���ʼ���з������Ż�
        Observations=Success(:,:,1:num_meas);
        %����һ��ȡfi
        for i=1:num_meas
            lo(1:3)=Observations(1:3,P_S,i);
            fi=[1,lo(1:3)-[Rx,Ry,Rz]];
            W=fi;
            for p=1:num_nodes
                R=W;
                W=[R,Observations(:,p,i)'];
            end
            FI1(i,:)=W;
        end
        %��������ȡfi
        for i=1:num_meas
            lo(1:3)=Observations(1:3,P_S,i);
            fi=[1,lo(1:3)-[Rx,Ry,Rz]];
            W=fi;
            for p=1:num_nodes
                R=W;
                Comb = 0;
                for k=1:d_dimension
                    Comb = Comb+inner_mat(k,:)*Observations(:,p,i);
                end
                W=[R,Comb];
            end
            FI2(i,:,P_S)=W;
        end
        
        %----------------------------------Y�ı�ʾ
        for i=1:num_meas
            j=order_S(i);
            Y(:,i,P_S)=(history(:,P_S,j)-history(:,P_S,j-1))/h;
        end
        %    Y
        %% ����һ��ѹ����֪���
        delta_Y1 = [];
        for k=1:d_dimension
            cvx_begin
            variable x(n)
            minimize (Lamda*norm(x,1) + square_pos(norm(Y(k,:,P_S)'- FI1*x,2)))
            cvx_end
            X(:,k)=x;
            Y_1 = FI1*x;
            delta_Y1 = [delta_Y1, Y(k,:,P_S)'-Y_1];
        end
        
        % ��ȡ���x�й���Cij����
        Xi_hat = X(d_dimension+2:end,:);
        Xi_hat_sum =sum(Xi_hat,2);%��Ci_hat�������
        Ci_hat1 = zeros(num_nodes,1);
        
        for i = 1:num_nodes
            Temp = 0;
            for j = 1:d_dimension
                index = d_dimension*(i-1)+j;
                Temp = Temp + Xi_hat_sum(index);
            end
            Ci_hat1(i) = Temp/sum(sum(inner_mat));
        end
        
        C_hat1(P_S,:) = Ci_hat1;
    end
    %% �����ع�
    FI_GR = FI2(:,:,1);%�����ع����Ϻ��FI
    for i = 2:num_nodes
        B = FI_GR;
        FI_GR = blkdiag(B,FI2(:,:,i));
    end
    Y_GR = [];%�����ع����Ϻ��Y��
    for i = 1:num_nodes
        G = Y_GR;
        H = sum(Y(:,:,i));
        Y_GR = [G,H];
    end
    
    % �������Ϻ��FI���н�һ���Ĵ���
    FI_GR_f = [];%FI_GR(:,1:d_dimension+1);
    FI_GR_b = [];
    for i = 1:num_nodes
        V = (i-1)*(num_nodes+1+d_dimension);
        S = FI_GR_f;
        FI_GR_f = [S,FI_GR(:,V+1:V+d_dimension+1)];
    end
    for i = 1:num_nodes-1
        for k=i+1:num_nodes
            A = FI_GR_b;
            FI_GR_b = [A,FI_GR(:,k+(i-1)*num_nodes+i*(d_dimension+1))+FI_GR(:,i+(k-1)*num_nodes+k*(d_dimension+1))-2*FI_GR(:,i*(num_nodes+d_dimension+2)-num_nodes)];
        end
    end
    FI_GR = [FI_GR_f,FI_GR_b];
    %% ���������
    x = BSSl0(FI_GR,Y_GR', p1,K_Bin, 0.1, 0.5, 2, 1000);
    delta_Y2 = [];
    Y_2 = FI_GR * x;
    delta_Y2 = [delta_Y2, Y_GR'-Y_2];
    %�������ó���C����
    j = (d_dimension+1)*num_nodes;
    for k = 1:num_nodes-1
        row = num_nodes-k;
        C_end(k,k+1:num_nodes) = x(j+1:j+row);
        j=j+row;
    end
    C_end(num_nodes,num_nodes)=0;
    C_hat2=C_end+C_end';
    
    
    %% =====================ָ����㼰��ͼ=====================
    
    for i=1:num_nodes
        outer_mat(i,i)=0;
    end
    for i=1:num_nodes
        C_hat1(i,i)=0;
    end
    %������ת��Ϊ��������ʽ
    x_orig = outer_mat(:);%ʵ�ʵ�C����
    C_1 = C_hat1(:);%����һ�Ľ��
    C_2 = C_hat2(:);%�������Ľ��
    
    x_ALL=[x_orig,C_1,C_2];
    %  prec_rec(C_1,x_orig,'holdFigure',1);
    prec_rec(C_2,x_orig,'holdFigure',2);
    %  [ AUROC1, AUPR1,prec1, tpr1, fpr1] = prec_rec(C_1,x_orig);
    %  [ AUROC2, AUPR2,prec2, tpr2, fpr2] = prec_rec(C_2,x_orig);
    [AUROC1, AUPR1, prec1, tpr1, fpr1, Accuarcy1,Tpr1, NMSE1] = EvaluationFcn( C_hat1,num_nodes,outer_mat);
    [AUROC2, AUPR2, prec2, tpr2, fpr2, Accuarcy2,Tpr2,NMSE2] = EvaluationFcn( C_hat2,num_nodes,outer_mat);
    TPR1=[TPR1,Tpr1];
    TPR2=[TPR2,Tpr2];
    %  FPR1=[FPR1,fpr1];
    %  FPR2=[FPR2,fpr2];
    ACC1=[ACC1,Accuarcy1];
    ACC2=[ACC2,Accuarcy2];
    AUroc1 = [AUroc1, AUROC1];
    AUroc2 = [AUroc2, AUROC2];
    Nmse1=[Nmse1,NMSE1];
    Nmse2=[Nmse2,NMSE2];
end

 
         