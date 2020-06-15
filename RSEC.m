%##############################################################################%
% Created Date: Saturday January 4th 2020                                      %
% Author: Li Hongmin (li.hongmin.xa@alumni.tsukuba.ac.jp)                      %
%##############################################################################%
% This function is for Robust Spectral Ensemble Clustering, which is described
% in Tao, Zhiqiang, et al. "Robust spectral ensemble clustering." Proceedings of the 25th ACM International on Conference on Information and Knowledge Management. ACM, 2016.

function label = RSEC(baseCls, k, lambda1, lambda2)

    % obtain co-assocition matrix
    [n, m] = size(baseCls);
    S = zeros(n, n);

    for i = 1:m
        tmpIDX = baseCls(:, i);
        S = S + double((repmat(tmpIDX, 1, n) - repmat(tmpIDX', n, 1)) == 0);
        % %         fprintf('%d/%d\n',i,numIter);
    end

    S = S * 1.0 / m;

    %% initialize

    [J, Z, E, Y1, Y2] = deal(zeros(n, n));
    H = zeros(n, k);
    mu = 10;
    p = 1.2;
    % p > 1

    t = 1;

    Dz = eye(n);
    oldobj = 1e200;
    obj = 1e100;
    while  oldobj > obj
        % update J
        tmp = Z + 1 / mu * Y2;
        if norm(tmp) >0 
            [u, s, v] = svt(tmp, 'lambda', lambda1 / mu,'k',k);
            J = u*s*v';
        end
        % update Z
        Dsqrt = sqrt(Dz);
        Z = sqrt(S * S' + eye(n)) * (S * S' + J - S' * E + 1/mu*(S'*Y1 -Y2 + Dsqrt\H*H'/Dsqrt));
        % update E
        Q = S - S*Z + Y1/mu;
        tmp = zeros(1,n);
        for j = 1:n
            nq = norm(Q(j,:));
            if nq > lambda2/mu
                tmp(j) = (nq -lambda2/mu)/nq;
            end
        end
        E = diag(tmp);
        
        % compute D and L
        Lz = (Z+ Z')/2 + H*H';
        Dz = diag(sum(Lz,2));
        
       % set H as the smallest k eigenvector of Lz
       opts.disp = 0;
       [H, ~] = eigs(Lz, k, 'sa', opts);
       
       % update multipliers
       tmp = S-S*Z-E;
       Y1 = Y1 + mu*(tmp);
       Y2 = Y2 + mu * (Z -J);
       
       oldobj = obj;
       obj = norm(tmp,'fro')/norm(S,'fro');
       % set mu
       mu = p* mu;
       t = t +1;   
       
       
    end
    
    label=SC(Z,k);
%     label=kmeans(H,k);

end
