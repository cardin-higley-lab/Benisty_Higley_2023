function [M, tC] = RiemannianMean(tC)
Np=size(tC,3);

M  = mean(tC, 3);
epsv=zeros(200,1);
for ii = 1 : 200
    A = M ^ (1/2);      %-- A = C^(1/2)
    B = A ^ (-1);       %-- B = C^(-1/2)
    NN=0;
    S = zeros(size(M));
    for jj = 1 : Np
        C = tC(:,:,jj);
        del = logm(B * C * B);
        if any(~isreal(del(:)))
            error('Corr matrix is not invertible');
        end
        S = S + A * del * A;
        NN=NN+1;
    end
    S = S / NN;
    
    M = A * expm(B * S * B) * A;
    
    eps = norm(S, 'fro');
   disp(eps);
    if eps < .10%1e-6
        break;
    end
end
% if eps>1e-6
%     error('Corr matrix is not invertible');
%        
% end
end