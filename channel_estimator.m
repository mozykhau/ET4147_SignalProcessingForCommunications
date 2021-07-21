function h = channel_estimator(x,s,L)
    %CHANNEL_ESTIMATOR Estimate channel from signal, data and length of
    %channel
    % Input:
    % x = data vector
    % s = signal vector
    % L = length of channel in symbol periods
    % Output:
    % h = channel
    delay =2;
    s = s(delay:499);
    
    
    mu = 0.01; % step size
    [~,len_x] = size(x);
    P = 4;
    len_s = length(s);
    a = zeros(2*P,2*L);
    
    for k = 1:len_s-1
        e = x(:,k)-a.*s(k)
        a_new = a+mu.*e*conj(s(k));
        a = a_new;
    end
    h = a(:,end);
end

