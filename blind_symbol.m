function s = blind_symbol(X)
    %BLIND_SYMBOL Estimate symbol from SVD
    L = 1;
    [U,S,V] = svd(X);
    len = length(X);
    delay = 2;
    sig = delay+1;
    
    % noise subspace 3:end
    noise_sub = V(:,sig:end);
    [A,B] = size(noise_sub);
    
    % Create block matrix
    null_mat = zeros(len+delay,delay*B);
    for i = 1:delay
        null_mat(i:A+i-1,(i-1)*B+1:i*B) = noise_sub;
    end
    
    % Find null space of block matrix
    rank_mat = rank(null_mat);
    [len_null_mat,~] = size(null_mat);
    [U1,~,~] = svd(null_mat);
    
    s = U1(sig-2:end-1,end-1);
    
end

