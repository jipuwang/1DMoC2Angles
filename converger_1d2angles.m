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
function [order_phi]=converger_1d2angles(assumedSoln)
% clear;
nGrids=2%6;%10;%8;
refinementRatio=2;

% Geometry
Tau=10; 

% Case configure options
if ~exist('assumedSoln','var')
  assumedSoln='sine_sine_sine';
end

error_phi0_n=zeros(nGrids,1);
gridMeshSize_n=ones(nGrids,1);
N=16; % angular discretization, fixed not refined. 
I=16;

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

  [phi0_j_ana,psi_b1_n_i,psi_b2_n_i,Q_MMS_j_n_i]=...
        manufacturer_1d2angles(J,N,I,Tau,mat,assumedSoln);
      
  [phi0_j]=OneDMoC_2Angles(J,N,I,Tau,mat,...
    psi_b1_n_i,psi_b2_n_i,Q_MMS_j_n_i);

  % Calculate the error compared to manufactured solution
  error_phi0_n(iGrid)=norm(phi0_j-phi0_j_ana,2)/sqrt(J);
  
end

% Calculate the order of accuracy
order_phi_nMinus1=ones(nGrids-1,1);
for j=1:nGrids-1
  order_phi_nMinus1(j)=log(error_phi0_n(j)/error_phi0_n(j+1)) / ...
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

% % Plot the solution
% scalarFlux_plot_handle=figure(13);
% plot(phi0_j,'-*');
% % title('scalar flux');
% xlabel('mesh size [cm]');
% ylabel('scalar flux');

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
display(assumedSoln);
% % Display the result
error_phi0_n
order_phi_nMinus1
display(char(strcat('soln_',assumedSoln)));
display(char(num2str(order_phi_nMinus1(end))));

order_phi=order_phi_nMinus1(end);

% aa=0.0;
end