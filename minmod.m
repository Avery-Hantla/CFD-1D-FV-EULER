%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Function to limit the 2nd order reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S_i = minmod(Qbar_i,Qbar_im1,Qbar_ip1,dX)
    % Find the left and right slopes of Qi
    SL = (Qbar_i(:,2:end-1)-Qbar_i(:,1:end-2))./(dX);
    SR = (Qbar_i(:,3:end)-Qbar_i(:,2:end-1))./(dX);

    S_i = sign(SL).*min(abs(SL),abs(SR));
    S_i((SL.*SR)<=0) = 0;

    % Add 1 and N cells
    S1 = (Qbar_i(:,2)-Qbar_i(:,1))./dX;
    SN = (Qbar_i(:,end)-Qbar_i(:,end-1))./dX;

    S_i = [S1,S_i,SN];
end