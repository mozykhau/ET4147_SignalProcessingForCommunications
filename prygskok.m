function pu = prygskok(L,P,tau)
    % PRYGSKOK Generate coding scheme
    l = L*P;
    for i = 1:l
        t = (i-1)/P-tau;
        if t<1 && t>=0
            if t < 0.25 || (t>=0.5 && t<0.75)
                pu(i,1) = 1;   
            else
                pu(i,1) = -1;
            end    
        else
            pu(i,1) = 0;
        end        
    end
%     i = 1:l
%     figure()
%     plot(i/P,pu)

end

