function [w,y] = cma_try(X,mu, w_init)
    %% CMA Blind Beamformer
    %   Input: 
    %   X = data matrix
    %   mu = step size
    %   w_init = intitial beamformer
    %   Output:
    %   w = final beamformer
    %   y = sequence of estimated signal samples
    %%
    [M,N] = size(X);
    w(:,1) = w_init;
    for n = 1:size(X,2)-1
        y(n) = w(:,n)'*X(:,n);
        w(:,n+1) = w(:,n)-mu.*X(:,n).*(abs(y(n))^2-1)*conj(y(n));
        w(:,n+1) = w(:,n+1)/norm(w(:,n+1));
    end
    
end

