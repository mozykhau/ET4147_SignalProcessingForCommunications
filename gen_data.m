function X = gen_data(M,N,Delta,theta,SNR)
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
    var_s = 10^(SNR/10)*var_noise; %variance of signal
    S = sqrt((var_s^2)/2)*(randn(d,N)+1i*randn(1,N));
    
    %% Create Array response matrix

    for i=1:d
        A(:,i) = gen_a(M,Delta,theta(i));
    end
    
    %% Create data matrix
    X = A*S+Noise;
end

