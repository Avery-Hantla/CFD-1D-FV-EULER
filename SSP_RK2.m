%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SSP RK2 Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q_np1 = SSP_RK2(Q_n,dQdt,dX,sigma,gamma,QBC,islimon,order,aa,dAdX)
    % Find flow variables
    [~,u,~,~,c] = flowvariables(Q_n,gamma);
    
    % Find current iteration time step
    dt = (dX.*sigma)./(abs(u)+c);
    dt = min(dt);
    
    % Calculate Qstar and find the residual of Qstar
    Q_star = Q_n + dt.*dQdt;
    res_star = res(Q_star,dX,order,islimon,QBC,gamma,aa,dAdX);

    % Find Q for next time step
    Q_np1 = 0.5.*(Q_n+Q_star+dt.*res_star);
end