%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           MATLAB 1D Euler Code
%                               Avery Hantla
%                              November, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% Inputs
Xbounds = [-4, 4];              % Grid Boundarys
num_points = [100,200,400,800]; % Number of points
sigma = 0.75;                   % CFL Number
gamma = 1.4;                    % Specific Heat Ratios
order = 2;                      % Desired Order of Error
N = 10000;
islimiteron = true;            % Use limiter? true/false
flow = 2;                       % Flow problem 1/2
isplot = false;                  % Plot during sim? true/false

for zdx = 1:length(num_points)
    clear Qbar Q % Clear variables between mesh simulations
    
    % Initilize the domain 
    dX = (Xbounds(2) - Xbounds(1))/(num_points(zdx)-1);
    X = Xbounds(1):dX:Xbounds(2);
    
    % Find analytical solution on mesh domain size
    [X_analytical,q_analytical] = ExactNozzle(Xbounds,flow,num_points(zdx));
    
    % Make a area vector at points
    aa = zeros(1,length(X));
    for jdx = 1:length(X)
       x_int = X(jdx);
       aa(jdx) = area(x_int);
    end

    % Make area vector at cells 
    aa_im1 = aa(1:end-1); %circshift(aa,1);
    aa_ip1 = aa(2:end); %circshift(aa,-1);
    Abar = (aa_ip1+aa_im1)./2;

    % Find slope of area across cell
    dAdX = (aa_ip1-aa_im1)./(dX); 
    
    % Guess initial conditions
    rho = q_analytical(1,end)*ones(1,length(X)-1); % Initial Density
    u = (q_analytical(2,end)/rho(1))*ones(1,length(X)-1); % Initial Velocity
    P = (q_analytical(3,end)-0.5*rho(1)*u(1)^2)*(gamma-1)*ones(1,length(X)-1);
    E = ((P./(gamma-1))+ 0.5.*rho.*u.^2);

    Qbar(:,:,1) = [rho;rho.*u;(P./(gamma-1))+0.5*rho.*u.^2]; % Q(:,i)

    % Specify the boundary conditions
    QBC = [q_analytical(:,1),q_analytical(:,end)];

    % Set initial conditiosn for while loop
    res(1) = 10; n = 1;
    while res(n) > (10^(-6)) && n < N
        % for ndx  = 1:4000
        % Reconstruct Cells
        [QL_iphalf, QR_iphalf] = reconstruction(Qbar(:,:,n),dX,order,islimiteron,QBC);
  
        % Calculate rossuvinov flux
        [F_iphalf,F_imhalf] = riemann(gamma,QL_iphalf,QR_iphalf);
        
        % Calculate current flow variables
        [rho,u,E,P,~] = flowvariables(Qbar(:,:,end),gamma);

        % Compute G
        Gbari = [-rho.*u.*(1./Abar).*dAdX;
        -rho.*u.^2.*(1./Abar).*dAdX;
        -(u.*(E+P)./Abar).*dAdX];
        
        dQdt = Gbari - (F_iphalf-F_imhalf)./dX;
        Qbar(:,:,n+1) = SSP_RK2(Qbar(:,:,n),dQdt,dX,sigma,gamma,QBC,islimiteron,order,Abar,dAdX);

        Q(:,:,n) = (QL_iphalf+QR_iphalf)/2;
        
        % Save residual and disp residual/figure
        res(n+1) = max(dQdt,[],'all');
        if mod(n,25) == 0 
            fprintf('Maximum Residual is: %d \n',res(n))
            if isplot == true
                plotQ(X,Q(:,:,end),'-',false)
                drawnow 
            end
        end
        n=n+1;
    end

    % Find the L2 Error and save residuals and Q
    Q_save{:,:,zdx} = Q(:,:,end);
    res_save{:,zdx} = res;
    EL2(zdx) = sqrt((sum((Q(1,:,end)-q_analytical(1,:)).^2))/(length(X)));
    
end % End multiple mesh sizes loop

% Plot the conserved variables
plotQ(X_analytical,q_analytical,'-.',true)
for zdx = 1:length(num_points)
    plotQ(Xbounds(1):(Xbounds(2) - Xbounds(1))/(num_points(zdx)-1):Xbounds(2),Q_save{:,:,zdx},'-',false)
end

% Add legend
leg{1} = 'Analytical Solution';
for jdx = 1:length(num_points)
    leg{jdx+1} = sprintf('%d Points',num_points(jdx));
end
legend(leg,'Location','southeast')

% Find order of error
for zdx = 1:length(num_points)-1
    P = (log(EL2(zdx)/EL2(zdx+1)))/log(2);
    fprintf('Order of error between the mesh with %d points and %d points is: %d \n',num_points(zdx),num_points(zdx+1),P)
end

%% Plotting for report
[rho_a,u_a,~,P_a,~] = flowvariables(q_analytical,gamma);
figure(2); plot(X_analytical,rho_a,'-.'); xlabel('x'); ylabel('Density, rho'); hold on
figure(3); plot(X_analytical,u_a,'-.'); xlabel('x'); ylabel('Velocity, u'); hold on
figure(4); plot(X_analytical,P_a,'-.'); xlabel('x'); ylabel('Pressure, p'); hold on
for jdx = 1:length(num_points)
    [rho,u,~,P,~] = flowvariables(Q_save{:,:,jdx},gamma);
    figure(2); plot(linspace(Xbounds(1),Xbounds(2),length(rho)),rho)
    figure(3); plot(linspace(Xbounds(1),Xbounds(2),length(rho)),u)
    figure(4); plot(linspace(Xbounds(1),Xbounds(2),length(rho)),P)

    figure(5); semilogy(1:length(res_save{jdx}(2:end)),res_save{jdx}(2:end)); hold on
end
save_names = {'rho.png','velocity.png','pressure.png'};
locations = {'southeast','northwest','southeast'};
for zdx = 2:4
    figure(zdx); legend(leg,'Location',locations{zdx-1}); %saveas(gcf,save_names{zdx-1});
end
figure(5); legend(leg{2:end},'Location','northeast'); xlabel('Iterations, n'); ylabel('Residual'); %saveas(gcf,'res.png')