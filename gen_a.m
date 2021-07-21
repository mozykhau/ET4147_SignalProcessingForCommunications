function response = gen_a(M, Delta, Theta)
    
    % generate the array response a(theta) of a uniform linear array
    
    % M = number of elements 
    % Delta = spacing  wavelengths
    % Theta = direction of arrival
    % slide 14 of lecture 1
    
    response = exp(1i*2*pi*Delta*sin(Theta)*((1:M)-1));
    response = response.';
end

