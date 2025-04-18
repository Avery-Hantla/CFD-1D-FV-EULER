%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Function to calculate the residual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = res(Qbar,dX,order,islimon,QBC,gamma,Abar,dAdX)
    [QL_iphalf, QR_iphalf] = reconstruction(Qbar,dX,order,islimon,QBC);
    [F_iphalf,F_imhalf] = riemann(gamma,QL_iphalf,QR_iphalf);

    [rho,u,E,P,~] = flowvariables(Qbar,gamma);

    % Compute G
    Gbari = [-rho.*u.*(1./Abar).*dAdX;
    -rho.*u.^2.*(1./Abar).*dAdX;
    -(u.*(E+P)./Abar).*dAdX];
    
    res = Gbari - (F_iphalf-F_imhalf)./dX;
end