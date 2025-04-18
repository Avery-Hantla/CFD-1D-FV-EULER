%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Function to reconstruct cells for 1st and 2nd order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QL_iphalf, QR_iphalf] = reconstruction(Qbar,dX,order,islimon,QBC)
    if order == 2
        % Shift Indexs
        Qbar_i = circshift(Qbar,0,2);
        Qbar_ip1 = circshift(Qbar,-1,2);% Qbar(:,2:end);
        Qbar_im1 = circshift(Qbar,1,2);% Qbar(:,1:end-1);

        % Find Slope
        if islimon == true
            S_i = minmod(Qbar_i,Qbar_im1,Qbar_ip1,dX);
        else 
            S_i = (Qbar_ip1-Qbar_im1)./(2*dX);
            S_1 = (Qbar_i(:,2)-Qbar_i(:,1))./dX;
            S_N = (Qbar_i(:,end)-Qbar_i(:,end-1))./dX;
            S_i = [S_1,S_i(:,2:end-1),S_N];
        end
        
        % Reconstructs Q i+1/2
        QL_iphalf = Qbar_i + S_i*(0.5*dX);
        
        %S_ip1 = circshift(S_i,-1);
        QR_iphalf =  Qbar_i - S_i.*(0.5.*dX);

        % Put boundary conditions in Q Reconstruction
        QL_iphalf = [QBC(:,1),QL_iphalf];
        QR_iphalf = [QR_iphalf,QBC(:,2)];
        
    elseif order == 1 
        QL_iphalf = [QBC(:,1),Qbar];

        QR_iphalf = [Qbar,QBC(:,2)];
    end
end