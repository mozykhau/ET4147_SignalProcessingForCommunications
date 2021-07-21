function [DOAs] = music(X,d,Delta)
    %% MUSIC algortihm to define DOA
    %   Input: 
    %   X = data matrix
    %   d = number of sources
    %   Output:
    %   theta = angle vector
    
    %%  Covariance matrix vs noise subspace (page 71, example 3.4)
    [M,N] = size(X);
    Rx = (X*X')/N;

    
    [V, D] = eig(Rx);
    noiseSub = V(:, 1:M-d); % Noise subspace of R
    %% Music Algrotihm
    theta = -90:1:90; %Peak search
    a = zeros(M, length(theta));
    cost = zeros(length(theta), 1);
    inv_cost =  zeros(length(theta), 1);
     for i = 1:length(theta)
         a(:,i) = gen_a(M,Delta,deg2rad(theta(i)));
         cost_nom = (norm(noiseSub'*a(:,i))).^2;
         cost_denom = (norm(a(:,i))).^2;
         cost(i) = cost_nom/cost_denom;
         inv_cost(i) = 1./cost(i);
     end
     
     [costSort, orgInd] = sort(inv_cost, 'descend');
     DOAs = orgInd(1:d, 1)-91;
     [DOAs,~] = sort(DOAs);
     %% Plot
%      figure()
%      plot(theta,10*log10(inv_cost))
%      title('Direction of Arrival estimation via Music')
%      xlabel('Direction of arrival')
%      ylabel('Power [dB]')

    
end

