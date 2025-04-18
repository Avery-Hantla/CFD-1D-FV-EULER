%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Function to plot Q matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotQ(X,Q,linetype,ishold)
    % Plot the analytical solution
    subplot(3,1,1);
    plot(X,Q(1,:),linetype); if ishold == true; hold on; end
    xlabel('x');
    ylabel('rho');
    % ylim([0.4,1.1]);
    
    subplot(3,1,2);
    plot(X,Q(2,:),linetype); if ishold == true; hold on; end
    xlabel('x');
    ylabel('rho*u');
    % ylim([-0.,0.8]);
    
    subplot(3,1,3);
    plot(X,Q(3,:),linetype); if ishold == true; hold on; end
    xlabel('x');
    ylabel('E_t');
    % ylim([0.8,2.0]);
end