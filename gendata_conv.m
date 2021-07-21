function x = gendata_conv(s,P,N,sigma)
    % GENDATA_CONV Construct 1/P sampled version of the output 
    % x(t) = h(t)*s(t)
    %
    % Input:
    % s = signal
    % P = sampling rate
    % N = number of samples
    % sigma = std of zero-mean Gaussian noise
    %
    % Output:
    % x = data matrix
    len = N*P;
    L = 1;
    
    % Generate noise
    n = randn(1,len)+1i*randn(1,len);
    n = sigma*n./std(n);
    
    % Generate channel
    tau = 0;
    h = channel_exam(L,P,tau);
    
    s2(1,:) = kron(s,[1;zeros(P-1,1)]);
    x = conv(h,s2);
    x = x(1:len)+n;
end



