function s = source(N)
    %Construct QPSK source sequence fo length N
    %
    % s = source sequence
    % N = length
    m = zeros(N,1);
    k = zeros(N,1);
    s = zeros(N,1);
    
    k = randi([0, 1], [N, 1]);
    m = randi([0, 1], [N, 1]);

    s(:,1) = ((1i).^(2-m(:,1))).*(-1+2*k(:,1));

end

