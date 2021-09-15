function [aff_mat, sigma] = CalcInitAff2D( data, params )
% Calculates the affinity between the columns of the data using euc distnace
% Output:
%     aff_mat - N-by-N symmetric matrix of non-negative affinities
%--------------------------------------------------------------------------

data = data.';
euc_dist = squareform(pdist(data));
nn_dist = sort(euc_dist.').';
params.knn = min(params.knn, size(nn_dist, 2));
sigma = median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));

aff_mat = exp(-euc_dist.^2/(2*sigma^2));


