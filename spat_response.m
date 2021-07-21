function y = spat_response(w,Delta,theta_range,M)
    % Generate the spatial response of a given beamformer w as a function 
    % of the direction of a source with array response a();
    % w = beamformer
    % Delta = spacing  wavelengths
    % theta_range = direction range
    % slide 20 lecture 1
    l = length(theta_range);
    for i=1:l
        a(:,i) = gen_a(M,Delta,theta_range(i));
    end
    %y(:) = abs(w'*a(:,1));
    y = abs(w'*a);
    plot(rad2deg(theta_range), y)
    xlabel('Angle');
    ylabel('Power');
    xlim([-90, 90])
    str = {sprintf("M = "+M),sprintf("Delta = "+Delta)};
    title(str);
end