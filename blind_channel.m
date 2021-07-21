function h = blind_channel(X)
    %BLIND_CHANNEL Estimate channel h from null subspace
    L = 1;
    [U,~,~] = svd(X);
    null_sub = U(:,3:end);
    null_her = null_sub';
    %mat_Her = blkdiag(null_her,null_her)
    
    h = null(null_her);
end

