function [w,y] = cma_try2(X,mu, w_init)
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
    w = w_init;
    w_new = zeros(M,1);
    for i = 1:N
%         y = w'*X(:,i);
%         w_new = w-mu.*X(:,i).*(y*conj(y)-1)*conj(y);
        J(i) = (abs(w'*X(:,i))^2) - 1;
        w_new = w-mu*X(:,i)*X(:,i)'*w*J(i);
        w = w_new;
        %y_plot(i) = (sqrt(y*conj(y)));
        y_plot(i) = abs(w'*X(:,i));
        
        
    end
    
    y = w'*X;
    
    figure()
    plot(y_plot)
    xlabel('Sample number')
    ylabel('|y|')
    title('CMA convergence plot |y| vs sample number')
end


