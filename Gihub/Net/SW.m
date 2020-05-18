 NodesNum = 20;
K =2;
p = 0.2;
A = sparse(NodesNum, NodesNum);
 for i=1:K/2
     A = A + diag(ones(1, length(diag(A, i))), i);
 end
 A = sparse(A);
 for i=1:K/2
     A(i, (NodesNum-K/2+i):NodesNum) = 1;
 end
 A = A + A';
 for i=1:K/2
     P = rand(1,NodesNum);
     P = find(P <= p);
     for j=1:length(P)
         while true
             x = fix(rand()*NodesNum)+1;
             if A(x, P(j)) == 0
                 A(x, P(j)) = 1;
                 A(P(j), x) = 1;
                 break;
             end
         end
         
         if P(j) <= NodesNum-i
             A(NodesNum-i, P(j)) = 0;
             A(P(j), NodesNum-i) = 0;
         else
             A(P(j)-NodesNum+i, P(j)) = 0;
             A(P(j), P(j)-NodesNum+i) = 0;
         end
  
     end
 end
 save('SW-11.mat','A');
%  clear i j P x
%  B3 =A^3;
%  B33 = B3(1:NodesNum+1:end);
%  B2 = A^2;
%  B22 = B2(1:NodesNum+1:end);
%  c = B33./(B22.*(B22-1));
%  c(isnan(c)) = 0;
%  C = mean(c);
%  clear B33 B3 B2 B22
%  
%  Paths = graphallshortestpaths(tril(A), 'directed', false);
%  Paths = tril(Paths);
%  L = sum(sum(Paths)) / nchoosek(NodesNum, 2);