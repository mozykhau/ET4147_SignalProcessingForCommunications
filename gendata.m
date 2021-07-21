function [X,A,S] = gendata(M,N,Delta,theta,f,SNR)
    % generate a data matrix X = AS+N
    %
    % Input:
    % M = number of antennas
    % N = number of samples
    % theta = array of directions in radians
    % Delta = wavelength spacing
    % SNR = signal to noise ratio [dB]
    %
    % Output:
    % X = data matrix
    
    
    d = length(theta);
    
    
    %% Create noise matrix
    var_noise = 1; %variance of noise
    Noise = sqrt((var_noise^2)/2)*(randn(M,N)+1i*randn(M,N));
    %% Create signal matrix
    var_s = 10^(SNR/10); % Variance of Signal
    S = complex_sin_signal(f(:,1),0,N);
    for i = 1:d
        Ps = sum(abs(S(i,:)).^2)/N;
        Pn = sum(abs(Noise(i,:)).^2)/N;
        S(i,:) = sqrt(var_s*Pn)*S(i,:);
    end

    
    %% Create Array response matrix

    for i=1:d
        A(:,i) = gen_a(M,Delta,theta(i));
    end
    %% Create data matrix
    X = A*S+Noise;
    %X = A*S; % no noise
end


