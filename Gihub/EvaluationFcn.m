function [ AUROC, AUPR, prec, tpr, fpr, Accuarcy,TPR,NMSE] = EvaluationFcn( C_hat,num_nodes,outer_mat)
%EVALUATIONFCN ָ�����
%   ���Ƚ�������Ϊ�Խ���Ϊ�㣬�����TPR��׼ȷ�ʣ�
for i=1:num_nodes
    C_hat(i,i)=0; % �Խ��߻�Ϊ0
end
c_hat = reshape(C_hat,num_nodes*num_nodes,1);% ��C_hat���������ʽ
% ����ʵ��Cת��Ϊ������ʽ
C_real = outer_mat;
for i=1:num_nodes
    C_real(i,i)=0; % �Խ��߻�Ϊ0
end
c_real = reshape(C_real,num_nodes*num_nodes,1);
% ��ֵ����
NMSE = norm(c_hat -c_real)/norm(c_real);
epsi = 0.5;
for i=1:size(c_hat)
    if abs(c_hat(i)-1) < epsi
        c_hat(i) = 1;
    end
    if abs(c_hat(i)-0) < epsi
        c_hat(i) = 0;
    end 
end



% ����ROC���ߺ�PR����
[AUROC,AUPR,prec,tpr,fpr]=prec_rec(c_hat,c_real);


% ����Accuarcy��TPR
ACC = 0;
TP = 0;
TP_FN = 0;
for i=1:size(c_hat)
    if c_real(i) == c_hat(i)
        ACC = ACC + 1;
    end
    if ((c_real(i) == 1) && (c_hat(i) == 1))
        TP = TP + 1;  
        TP_FN = TP_FN + 1;
    end
    if ((c_real(i) == 1) && (c_hat(i) == 0))
        TP_FN = TP_FN + 1;
    end
end

Accuarcy = ACC/(num_nodes*num_nodes); %����׼ȷ��
TPR =TP / sum(c_real); %����TPR
% NMSE = norm(c_hat -c_real)/norm(c_real);
end

