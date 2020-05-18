function x = BSSl0(A, y, p0, K_Bin, sigma_min, sigma_decrease_factor, mu, L)

% For given A x = y, we reconstruct binary signal x from given A and y using boxed re-weigthed smoothed l0 optimization.

% input data:
% (1) A is the sensing matrix.
% (2) y is the observation signal.
% (3) p0 is the prob that x_i being 0.

% input parameters:
% (1) sigma_min is a small positive number near 0 that sigma converses to.
% (2) L is the inner loop iteration number
% (3) mu is the gradient descent factor

% (c) Tianlin Liu 2016


%% 20180802,hkk
% 限制s中部分元素非0，即1
% 假设限制s(K_Bin:end)中的元素非0即1
% s(1:K_Bin-1)为实数


% Initialization
A_pinv = pinv(A);
[~,N] = size(A);
x = A_pinv*y;
K = (1-p0)*N;

sigma = 2*max(abs(x));

M = log(sigma_min/sigma) / log(sigma_decrease_factor);

iter = 0;
% Main Loop
while sigma > sigma_min
    k = 1 + iter / M * K;
    iter = iter + 1;
    
    for i=1:L        
        delta = ConstrainedDelta(x,sigma,k,K,K_Bin)/k;
        x = x - mu*delta;
        x = x - A_pinv*(A*x-y);   % Projection
            
    end
    
    sigma = sigma * sigma_decrease_factor;
end


    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta = ConstrainedDelta(s,sigma,k,K,K_Bin)

[N,~] = size(s);
w = ((k+1)/2 - ((k-1)/2)*sign(s).*sign(1-s));

if w == (k+1)/2,
    w = 1;
end

%% Add by hkk
C_Bin = ones(1,N);
for i=1:K_Bin-1
C_Bin(i) = 0;
end
C_NBin = 1 - C_Bin;

delta0 = (1 - K/N).*s.*(exp((-abs(s).^2)/(2*sigma^2))).*w + (K/N)*(s-1).*(exp((-abs(s-1).^2)/(2*sigma^2))).*w;
% delta1 = 2*s;
% delta = C_Bin*delta0 + C_NBin*delta1;
delta = C_Bin'.*delta0;
% delta = delta0;
