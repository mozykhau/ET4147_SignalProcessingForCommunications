function h = channel_exam(L,P,tau)
    % CHANNEL Generate channel for input t    
    r = length(tau);
    l = L*P;
    h = zeros(l,1); 
    i = 1:l;
    for j = 1:r       
        h(:) = h(:,1)+prygskok(L,P,tau(j));       
    end 
end

