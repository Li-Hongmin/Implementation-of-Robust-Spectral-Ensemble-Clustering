
function labels = SC(S, k, kmRep)

    if ~exist('kmRep', 'var')
        kmRep = 5;
    end

    if ~issparse(S)
        S = sparse(S);
    end

    % D is degrees of Similarity
    n = size(S,1);
    D = sum(S, 2) + (1e-10);
    D = sqrt(1./D); % D^(-1/2)
    D = spdiags(D, 0, n, n);
    L = D * S * D;
    % L is normalize laplacian martix
    if size(L, 1) > 500
        OPTS.disp = 0;
        [V, ~] = eigs(L, k, 'LM', OPTS);
        V = real(V);
    else
        [eigVectors, eigValues] = eig(full(L));
        eigValues = diag(eigValues);
        [eigValues, idx] = sort(eigValues,'descend');
        nEigVec = eigVectors(:, idx(1:k));
        V = real(nEigVec);
    end
    sq_sum = sqrt(sum(V.*V, 2)) + 1e-20;
    U = V ./ repmat(sq_sum, 1, k);

    labels = kmeans(U, k, 'EmptyAction', 'drop', ...
        'Replicates', kmRep);
end