function s = sourceqpsk(N)
    %SOURCEQPSK generate qpsk signal
    %   Detailed explanation goes here
    a = 1/sqrt(2);
    k = randi([0, 1], [N, 1]);
    m = randi([0, 1], [N, 1]);
    s = a*((-1).^k(:)+1i*(-1).^m(:));
end

