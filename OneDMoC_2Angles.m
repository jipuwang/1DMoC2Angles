% 1D MoC Module
% Input: 
%   Geometry Tau
%   Spatial discretization J (or mesh size)
%   Polar discretization N
%   Azimuthal discretization I
%   Material: all cross sections
%   Boundary conditions
%   Distributed source, can be MMS
% Output: 
%   Cell-averaged scalar flux

function [phi0_j]=OneDMoC_2Angles(J,N,I,Tau,mat,...
           psi_b1_n_i,psi_b2_n_i,Q_MMS_j_n_i)

%   Input parameter
  if ~exist('Tau','var')
    Tau=10;
  end
  if ~exist('J','var')
    J=5*2;%*2%*2*2*2*2*2*2*2*2
  end
  if ~exist('N','var')
    N=16;
  end
  if ~exist('I','var')
    I=16;
  end
  if ~exist('mat','var')
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
  end
  if ~exist('psi_b1_n_i','var')
    psi_b1_n_i=ones(N,I)*1.0/(2*pi);
  end
  if ~exist('psi_b2_n_i','var')
    psi_b2_n_i=ones(N,I)*1.0/(2*pi);
  end
  if ~exist('Q_MMS_j_n_i','var')
    Q_MMS_j_n_i=ones(J,N,I)*0.3/(2*pi); % removed *2.0 (angular quantity)
  end
  
  % Material
  Sig_ss_j=mat.Sig_ss_j;
  nuSig_f_j=mat.nuSig_f_j;
  Sig_t_j=mat.Sig_t_j;
%   Sig_ss_j=ones(J,1)*0.5;
%   nuSig_f_j=ones(J,1)*0.2;
%   Sig_t_j=ones(J,1);

  Sig_t_inv_j=1./Sig_t_j;
  
  % Default variables, can be customized. 
  maxIterate=2000;
  epsilon_phi0=1e-12;
  delta=1E-13;
  [mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);
  [alpha_i,weight_i]=lgwt(I,0,2*pi);alpha_i=flipud(alpha_i);
    
  h_j=ones(J,1)*Tau/J;
  % N rays to trace, each angle has only 1 ray, no ray-spacing
  % n for each angle, and j for FSR region index
  segLen_j_n_i=zeros(J,N,I);
  for n=1:N
    for j=1:J
      segLen_j_n_i(j,n,:)=h_j(j)/abs(mu_n(n));
    end
  end
  
  phi0_j_old=ones(J,1);
  q_j_n_i=zeros(J,N,I);
  for iIterate=1:maxIterate
    for j=1:J
      for n=1:N
        for i=1:I
          q_j_n_i(j,n,i)=Sig_ss_j(j)*phi0_j_old(j)/(4*pi)+Q_MMS_j_n_i(j,n,i);
        end
      end
    end
    phi0_j_new=zeros(J,1);
    % ray tracing
    for n=1:N/2 % backward direction
      for i=1:I
        psi_in=psi_b2_n_i(n,i);
        for j=J:-1:1
          exp_temp=exp(-Sig_t_j(j)*segLen_j_n_i(j,n,i));
          psi_out=psi_in*exp_temp+q_j_n_i(j,n,i)*Sig_t_inv_j(j)*(1-exp_temp);
          psi_avg=q_j_n_i(j,n,i)*Sig_t_inv_j(j)+(psi_in-psi_out)/Sig_t_j(j)/segLen_j_n_i(j,n,i);
          phi0_j_new(j)=phi0_j_new(j)+weight_n(n)*weight_i(i)*psi_avg;
          psi_in=psi_out;
        end
      end
    end
    for n=N/2+1:N % forward direction
      for i=1:I
        psi_in=psi_b1_n_i(n,i);
        for j=1:J
          exp_temp=exp(-Sig_t_j(j)*segLen_j_n_i(j,n,i));
          psi_out=psi_in*exp_temp+q_j_n_i(j,n,i)*Sig_t_inv_j(j)*(1-exp_temp);
          psi_avg=q_j_n_i(j,n,i)*Sig_t_inv_j(j)+(psi_in-psi_out)/Sig_t_j(j)/segLen_j_n_i(j,n,i);
          phi0_j_new(j)=phi0_j_new(j)+weight_n(n)*weight_i(i)*psi_avg;
          psi_in=psi_out;
        end
      end
    end

    % test for convergence
%     error=norm(phi0_j_new-phi0_j_old);
    error=max(abs(phi0_j_new-phi0_j_old)./(phi0_j_new+delta));
    if error<epsilon_phi0
      break;
    end
    phi0_j_old=phi0_j_new;
  end  

  phi0_j=phi0_j_new;
  
  figure(11);
  plot(phi0_j,'*-');
  hold on;
  grid on;
  phi0_MMS =@(x) (sin(pi/(size(phi0_j,1)).*x)+1)*37.102114262431876;
  fplot(phi0_MMS,[0,size(phi0_j,1)],'bo-');
  legend('numerical','analytical');
%   hold off;
  iIterate
  
end
