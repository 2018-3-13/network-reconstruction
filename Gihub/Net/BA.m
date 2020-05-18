  tic
  m0 = 5;
  m = 3;
  NodesNum = 20;
  A = sparse(NodesNum, NodesNum);
  A(1:m0,1:m0) = round(rand(m0));
  A = tril(A); 
  A = A+A';
  A = A - diag(diag(A));
 for i=m0+1:NodesNum
     Degree = sum(A(1:i-1,1:i-1));
     for j=2:i-1
         Degree(j) = Degree(j) + Degree(j-1);
     end
     LinksNum = 0;
     while LinksNum<m
         link = fix(rand()*Degree(i-1)+1);
         for j=1:i-2
             if link<=Degree(1) && A(i,1)==0
                 A(i,1) = 1;     
                 A(1,i) = 1;
                 LinksNum = LinksNum+1;
             elseif link>Degree(j) && link<=Degree(j+1) && A(i,j+1)==0
                 A(i,j+1) = 1;       A(j+1,i) = 1;
                 LinksNum = LinksNum+1;
            end
        end
     end
 end
 save('BA-1.mat','A');
%  Degree = sum(A);
%  list = unique(Degree);
%  num = zeros(1,length(list));
%  for i=1:length(list)
%      num(i) = length(find(list(i)==Degree));
%  end
%  toc
%  loglog(list,num ./ sum(num),'.','markersize',20)
%  xlabel('k'),ylabel('P(k)')