function [signal] = complex_sin_signal(f,phase,N)
    %COMPLEX_SIN_SIGNAL Summary of this function goes here
    %   Detailed explanation goes here
    %   f = normalised frequency
    
    ampin = 1; %amplitude
    k =1:N;
    signal = ampin*exp(1j*2*pi*f*k + phase);
end

