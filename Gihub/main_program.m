%% main函数
% 参数初始化——数据采集——重构——结果评价
clc
clear all

%% =====================参数初始化========================
global num_nodes;
global h ;                    %积分时间
h = 0.001;%0.01
num_meas =12;%45                 %观测点数(<34)64
d_dimension = 3;              %每个振子的维数
Data_Collection =0;           %数据采集开关
iterationTime =20000;    % 80000(5kaishi) %200000;%迭代次数500000
ref_point=5000;         % 30000   %参考点的选择(80000)
history=[];                   %记录轨迹上所有点
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
inner_mat=[0,1,0;0,0,1;1,0,0];%内联矩阵
%方法一参数
Lamda = 0.1;                  %L1范数的系数
%方法二参数
p1=0.2;                       %0的概率
%% =====================数据采集=====================
%数据提取1：
for t=1:1
    %      p2 = t*0.1;
    %      num_meas=num_meas+5;
    str1='SW-DATA';str2=num2str(t);
    filename = [str1,str2];
    K = load(filename);%读取网络结构
    Adj = K.a;
    % K = importdata('SW-DATA1.mat');%读取网络结构
    % Adj = full(K);
    % Adj = importdata('swnet20-2.mat');%读取网
    num_nodes = size(Adj,1);
    %数据提取2：
    % K = load('A.txt');%读取网络结构
    %  num_nodes = max(max(K(:,1)),max(K(:,1)));
    K_Bin=(d_dimension+1)*num_nodes+1;%有效值的起始
    n = d_dimension*(num_nodes+1)+1;%方法一中向量的维数
    % Adj = zeros(num_nodes,num_nodes);%邻接矩阵
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
%         scatter3(history(1,1,1:20000),history(2,1,1:20000),history(3,1,1:20000));%画轨迹图
    end
    
    %% =====================选择满足半径要求的数据，重构得到C矩阵=====================
    for P_S = 1:num_nodes %以第P_S个节点为中心，构造节点P_S与其他节点的连接关系
        j=1;
        Rx=history(1,P_S,ref_point);
        Ry=history(2,P_S,ref_point);
        Rz=history(3,P_S,ref_point);%参考点
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
        %% 数据预处理，使之满足y=\Phi*x的形式--此处开始进行方法的优化
        Observations=Success(:,:,1:num_meas);
        %方法一获取fi
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
        %方法二获取fi
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
        
        %----------------------------------Y的表示
        for i=1:num_meas
            j=order_S(i);
            Y(:,i,P_S)=(history(:,P_S,j)-history(:,P_S,j-1))/h;
        end
        %    Y
        %% 方法一：压缩感知求解
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
        
        % 提取求得x中关于Cij的项
        Xi_hat = X(d_dimension+2:end,:);
        Xi_hat_sum =sum(Xi_hat,2);%将Ci_hat按行求和
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
    %% 整体重构
    FI_GR = FI2(:,:,1);%整体重构整合后的FI
    for i = 2:num_nodes
        B = FI_GR;
        FI_GR = blkdiag(B,FI2(:,:,i));
    end
    Y_GR = [];%整体重构整合后的Y；
    for i = 1:num_nodes
        G = Y_GR;
        H = sum(Y(:,:,i));
        Y_GR = [G,H];
    end
    
    % 对于整合后的FI进行进一步的处理
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
    %% 方法二求解
    x = BSSl0(FI_GR,Y_GR', p1,K_Bin, 0.1, 0.5, 2, 1000);
    delta_Y2 = [];
    Y_2 = FI_GR * x;
    delta_Y2 = [delta_Y2, Y_GR'-Y_2];
    %方法二得出的C矩阵
    j = (d_dimension+1)*num_nodes;
    for k = 1:num_nodes-1
        row = num_nodes-k;
        C_end(k,k+1:num_nodes) = x(j+1:j+row);
        j=j+row;
    end
    C_end(num_nodes,num_nodes)=0;
    C_hat2=C_end+C_end';
    
    
    %% =====================指标计算及画图=====================
    
    for i=1:num_nodes
        outer_mat(i,i)=0;
    end
    for i=1:num_nodes
        C_hat1(i,i)=0;
    end
    %将矩阵转化为向量的形式
    x_orig = outer_mat(:);%实际的C矩阵
    C_1 = C_hat1(:);%方法一的结果
    C_2 = C_hat2(:);%方法二的结果
    
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

 
         