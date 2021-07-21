function x = gen_data1(h,s,P,N)
    % Construct 1/P sampled version of the output x(t) = h(t)*s(t)
    %
    % Input:
    % h = channel
    % s = signal
    % P = sampling rate
    % N = number of samples
    %
    % Output:
    % x = data matrix
    s2(1,:) = kron(s,[1;zeros(P-1,1)]);
    x = conv(h,s2);
    
    
end

