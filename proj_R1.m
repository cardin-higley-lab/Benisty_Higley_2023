function mX = proj_R1(mCref,tC)
K  = size(tC, 3);
M  = size(tC,1);
MM = M * (M + 1) / 2;
mX = zeros(MM, K);

mW = sqrt(2) * ones(M) - (sqrt(2) - 1) * eye(M);
for kk = 1 : K
    Skk      = logm(mCref * tC(:,:,kk) * mCref) .* mW;
    if all(~isreal(Skk(:)))
        mX=nan;
        return;
    end
%         
    mX(:,kk) = Skk(triu(true(size(Skk))));
end
end