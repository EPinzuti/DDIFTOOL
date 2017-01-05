

% Subselect kernel matrix from subset idx of data
function K = sker(K,idx)
K=K(1:size(idx,2),1:size(idx,2));
%K=K(size(idx,2),size(idx,2));