function [Freqs] = musicfreq(X,d)
    %% MUSIC algortihm to define normalised frequency
    %   
    %   Input: 
    %   X = data matrix
    %   d = number of sources
    %   Delta = spacing between microphones
    %   Output:
    %   theta = angle vector
    
    %%  Covariance mat            x vs noise subspace (page 71, example 3.4)
    [~,N] = size(X);
    X = X.';
    Rx = (X*X')/N;
    [~,~,V] = svd(Rx);
    %signalSub = V(:, 1:d); % Signal subspace of R
    noiseSub = V(:,d+1:end); % Noise subspace of R
    
    %% Music Algrotihm
    w_fft = 100;
    [h] = fft(noiseSub(:,1),w_fft);
    den   =  abs(h).^2;
    for n = 2:size(noiseSub,2)
        h = fft(noiseSub(:,n),w_fft);
        den = den + abs(h).^2;
    end
    pseudoSpect = 1./den; % This is the pseudospectrum
    [~, orgInd] = sort(pseudoSpect, 'descend');
    Freqs = (orgInd(1:d, 1)-1)/w_fft; 
    [Freqs, ~] = sort(Freqs);
    
     %% Plot
%      figure()
%      freqplot = linspace(0,1,w_fft);
%      plot(freqplot,10*log10(pseudoSpect));
%      title('Signal Spectrum via Music')
%      xlabel('Fundamental frequency x 2\pi')
%      ylabel('Power [dB]')

    
end

