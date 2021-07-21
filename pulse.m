function g = pulse(tau,L,P)
    %Generate a rate 1/P sampled version of the pulse g(t)
    %
    % tau = delay [0,1)         0.5
    % L = length of sampling    3
    % P = rate 1/P              1000
    l = L*P;
    for i = 1:l
        t = (i-1)/P-tau-1;
        if abs(t) < 1
            g(i,1) = 1-abs(t);
        else
            g(i,1) = 0;
        end    
    end
    
    %plot(g)
end
    

