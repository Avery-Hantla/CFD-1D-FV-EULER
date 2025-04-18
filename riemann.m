%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Function that calulcates the rosuvinov flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F_iphalf,F_imhalf] = riemann(gamma,QL_iphalf,QR_iphalf)
    % Shift Indexs
    % [QL_iphalf, QR_iphalf] = reconstruction(Qbar,dX,order,islimon,QBC);

    % Find value for left side
    [rhoL,uL,EL,PL,cL] = flowvariables(QL_iphalf,gamma); 
    FL_iphalf = [rhoL.*uL;(rhoL.*uL.^2)+PL;uL.*(EL+PL)];

    % Find values for right side
    [rhoR,uR,ER,PR,cR] = flowvariables(QR_iphalf,gamma);
    FR_iphalf = [rhoR.*uR;(rhoR.*uR.^2)+PR;uR.*(ER+PR)];

    % Find average u and c
    u = (uL+uR)./2;
    c = (cL+cR)./2;

    % Solve for F i plus half
    F = ((FL_iphalf+FR_iphalf)./2) - ((abs(u)+c)./2).*(QR_iphalf-QL_iphalf);
    F_iphalf = F(:,2:end);
    F_imhalf = F(:,1:end-1);
end
