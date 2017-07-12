%% Instruction
  % to run different cases, change the manufacturer only!
%% Info
% Grid refiner is for grid refinement analysis. 
% Right now, it needs to be aware of the geometry and the material. So it
% can call the manufacturer to get the MMS problem and solution
% The geometry and material also need to be passed to the coupler, so the
% coupler can keep passing the info on to the modules, because it's part of
% the problem description. 
% It needs to know the geometry and is responsible for generating the grid
% and pass the grid information to the coupler. 
function [order_phi,order_T]=converger(fbType,mocSrc,assumedSoln)
% clear;
nGrids=6;%10;%8;
refinementRatio=2;

% Geometry
Tau=10; 

% % Case configure options
% % fbType options: 'noFeedback','linear','squareRootPlus1'
% fbType='noFeedback'; 
% % mocSrc options: 'flat_source','linear_source'
% mocSrc='flat_source';
% % assumedSoln options: 'const_quadratic','sine_sine','sqrtPlus1_quadratic"
% assumedSoln='sqrtPlus1_quadratic'; 

error_phi0_n=zeros(nGrids,1);
error_T_n=zeros(nGrids,1);
gridMeshSize_n=ones(nGrids,1);
N=16; % angular discretization, fixed not refined. 
for iGrid=1:nGrids
  J=5*refinementRatio^iGrid;
  gridMeshSize_n(iGrid)=Tau/J;
  iGrid
  % Material
  field1='Sig_t_j';          value1=ones(J,1);
  field2='Sig_ss_j';         value2=ones(J,1)*0.5;
  field3='Sig_gamma_j';      value3=ones(J,1)*0.4;
  field4='Sig_f_j';          value4=ones(J,1)*0.1;
  field5='nuSig_f_j';        value5=ones(J,1)*0.2;
  field6='thermal_cond_k_j'; value6=ones(J,1);
  field7='kappaSig_f_j';     value7=ones(J,1)*0.1; % kappa=1.0;
  mat = struct(field1,value1,field2,value2,field3,value3,... 
    field4,value4,field5,value5,field6,value6,field7,value7);

  % Define function handles 
  if strcmp(mocSrc,'linear_source')
    % call the manufacturer to get MMS problem and solution
    [phi0_j_ana,psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n, ...
          T_j_ana,T_L,T_R,p_MMS_j]=...
          manufacturer_LS(J,N,Tau,mat,assumedSoln,fbType);

    % call the coupler to solve the above manufactured problem
    [phi0_j,T_j]=coupler_LS(J,N,Tau,mat,psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n, ...
                T_L,T_R,p_MMS_j,fbType, ...
              phi0_j_ana,T_j_ana);
  else
    [phi0_j_ana,psi_b1_n,psi_b2_n,Q_MMS_j_n, ...
          T_j_ana,T_L,T_R,p_MMS_j]=...
          manufacturer(J,N,Tau,mat,assumedSoln,fbType);

    % call the coupler to solve the above manufactured problem
    [phi0_j,T_j]=coupler(J,N,Tau,mat,psi_b1_n,psi_b2_n,Q_MMS_j_n, ...
                T_L,T_R,p_MMS_j,fbType);
  end

  
  % Calculate the error compared to manufactured solution
  error_phi0_n(iGrid)=norm(phi0_j-phi0_j_ana,2)/sqrt(J);
  error_T_n(iGrid)=norm(T_j-T_j_ana,2)/sqrt(J);
  
%   %% Plot the solution over grid refinements
%   x=linspace(0,10,J);
%   figure(67); clf; hold on;
%   if iGrid==1
%     % title('scalar flux');
%     xlabel('mesh size [cm]');
%     ylabel('scalar flux');
%   end
%   % Plot the solution
%   plot(x,phi0_j,'-*');
%   
%   figure(68); clf; hold on;
%   if iGrid==1
%     % title('temperature');
%     xlabel('mesh size [cm]');
%     ylabel('temperature');
%   end
%   plot(x,T_j,'-o');

end

% Calculate the order of accuracy
order_phi_nMinus1=ones(nGrids-1,1);
order_T_nMinus1=ones(nGrids-1,1);
for j=1:nGrids-1
  order_phi_nMinus1(j)=log(error_phi0_n(j)/error_phi0_n(j+1)) / ...
    log(gridMeshSize_n(j)/gridMeshSize_n(j+1));
  order_T_nMinus1(j)=log(error_T_n(j)/error_T_n(j+1)) / ...
    log(gridMeshSize_n(j)/gridMeshSize_n(j+1));
end

% %% Visualize the results
% orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];
% 
% scalarFluxErrorRMS_plot_handle=figure(11);
% loglog(gridMeshSize_n,error_phi0_n,'*');
% % title('scalar flux error convergence');
% xlabel('mesh size [cm]');
% ylabel('scalar flux error RMS');
% 
% hold on;
% orderGuess=round(order_phi_nMinus1(end));
% errorStt=error_phi0_n(end)*refinementRatio^(orderGuess*(nGrids-1));
% firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
% secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
% thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
% fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
% loglog(orderPlotGrid,firstOrder,'--');
% loglog(orderPlotGrid,secondOrder,'--');
% loglog(orderPlotGrid,thirdOrder,'--');
% loglog(orderPlotGrid,fourthOrder,'--');
% legend('scalar flux error','1st Order','2nd Order',...
%   '3rd Order','4th Order','location','best');
% hold off;
% 
% temperatureErrorRM_plot_handle=figure(12);
% loglog(gridMeshSize_n,error_T_n,'*');
% % title('temperature error convergence');
% xlabel('mesh size [cm]');
% ylabel('temperature error RMS');
% 
% hold on;
% orderGuess=round(order_T_nMinus1(end));
% errorStt=error_T_n(end)*refinementRatio^(orderGuess*(nGrids-1));
% firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
% secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
% thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
% fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
% loglog(orderPlotGrid,firstOrder,'--');
% loglog(orderPlotGrid,secondOrder,'--');
% loglog(orderPlotGrid,thirdOrder,'--');
% loglog(orderPlotGrid,fourthOrder,'--');
% legend('temperature error','1st Order','2nd Order',...
%   '3rd Order','4th Order','location','best');
% hold off;
% 
% % Plot the solution
% scalarFlux_plot_handle=figure(13);r
% plot(phi0_j,'-*');
% % title('scalar flux');
% xlabel('mesh size [cm]');
% ylabel('scalar flux');
% 
% temperature_plot_handle=figure(14);
% plot(T_j,'-o');
% % title('temperature');
% xlabel('mesh size [cm]');
% ylabel('temperature');
% 
% % Save the plots
% phi0_RMS_fn=char(strcat('fbType_',fbType,'_mocSrc_',mocSrc,'_soln_',assumedSoln,'_','phi0_RMS'));
% T_RMS_fn=char(strcat('fbType_',fbType,'_mocSrc_',mocSrc,'_soln_',assumedSoln,'_','T_RMS'));
% phi0_fn=char(strcat('fbType_',fbType,'_mocSrc_',mocSrc,'_soln_',assumedSoln,'_','phi0'));
% T_fn=char(strcat('fbType_',fbType,'_mocSrc_',mocSrc,'_soln_',assumedSoln,'_','T'));
% 
% savefig(scalarFluxErrorRMS_plot_handle,phi0_RMS_fn)
% savefig(temperatureErrorRM_plot_handle,T_RMS_fn)
% savefig(scalarFlux_plot_handle,phi0_fn)
% savefig(temperature_plot_handle,T_fn)
% Display the problem description

disp '=================';
% % display(fbType);
% % display(mocSrc)
% % display(assumedSoln);
% % Display the result
error_phi0_n
error_T_n
order_phi_nMinus1
order_T_nMinus1
display(char(strcat('fbType_',fbType,'_mocSrc_',mocSrc,'_soln_',assumedSoln)));
display(char(num2str(order_phi_nMinus1(end))));
display(char(num2str(order_T_nMinus1(end))));


order_phi=order_phi_nMinus1(end);
order_T=order_T_nMinus1(end);

% aa=0.0;
end