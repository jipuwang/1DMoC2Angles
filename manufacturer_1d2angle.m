% A sine and sine solution generator including
    % Discretized analytical solution
    % Manufactured boundary conditions
    % Manufactured source
function [phi0_MMS_j,psi_b1_n,psi_b2_n,Q_MMS_j_n,...
          T_MMS_j,T_L,T_R,q_MMS_j]=...
          manufacturer(J,N,Tau,mat,assumedSoln,fbType)
  % input parameters
  if ~exist('J','var')
    J=5*2;%*2%*2*2*2*2*2*2*2*2
  end
  if ~exist('N','var')
    N=16;
  end
  if ~exist('Tau','var')
    Tau=10;
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
  if ~exist('assumedSoln','var')
    assumedSoln='sine_sine';
  end
  if ~exist('fbType','var')
    fbType='linear';
%     fbType='noFeedback';
%     fbType='sqareRoot';
  end
  % Material
  Sig_t_j=mat.Sig_t_j;
  Sig_ss_j=mat.Sig_ss_j;
  Sig_gamma_j=mat.Sig_gamma_j;
  Sig_f_j=mat.Sig_f_j;
  nuSig_f_j=mat.nuSig_f_j;
  kappaSig_f_j=mat.kappaSig_f_j;
  k_F=mat.thermal_cond_k_j(1);

  h=Tau/J;
  [mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);
  
  %% Manufactured Solutions for both fields
  % They need to be pre-defined here due to temperature dependence on the
  % xs. 
  % Options includes: sine_sine, const_cubic, sqrtPlus1_quadratic, etc.
  switch(assumedSoln)
    case('sine_sine')
      % Manufactured neutronics solution \psi(x,\mu)=sin(pi*x/Tau), 0<x<Tau
      psi_MMS =@(x) sin(pi*x/Tau);
      psi_MMS_Diff =@(x) pi/Tau*cos(pi*x/Tau);
      % Manufactured TH solution T(x)=sin(pi*x/Tau), 0<x<Tau
      T_MMS =@(x) sin(pi*x/Tau);
      T_MMS_xx =@(x) -(pi*pi/Tau/Tau)*sin(pi*x/Tau);
    case('const_cubic')
      % Manufactured neutronics solution \psi(x,\mu)=1.0, 0<x<Tau
      psi_MMS =@(x) 1.0+x*0.0;
      psi_MMS_Diff =@(x) x*0.0;
      % Manufactured TH solution T(x)=x.^3, 0<x<Tau
      T_MMS =@(x) x.^3;
      T_MMS_xx =@(x) x*6.0;
    case('sqrtPlus1_quadratic')
      % Manufactured neutronics solution \psi(x,\mu)=1.0, 0<x<Tau
      psi_MMS =@(x) sqrt(x+1);
      psi_MMS_Diff =@(x) 0.5./sqrt(x+1);
      % Manufactured TH solution T(x)=x.^2, 0<x<Tau
      T_MMS =@(x) x.^2;
      T_MMS_xx =@(x) 2.0+x*0.0;
  end
  
  %% XS update due to temperature feedback!
  % Change in capture is reflected in change in total. 
  switch fbType
    case 'noFeedback'
      Sig_gamma =@(x) mat.Sig_gamma_j(1)+0.0*x;
    case 'linear'
      % Assumes the original xs is homogeneous
      T0=50;
      gamma_coeff=0.004;
      Sig_gamma =@(x) mat.Sig_gamma_j(1)+gamma_coeff*(T_MMS(x)-T0);
    case 'squareRootPlus1'
      T0=50;
      Sig_gamma =@(x) mat.Sig_gamma_j(1)*sqrt((T0+1)./(T_MMS(x)+1));
  end
  
  Sig_ss =@(x) Sig_ss_j(1)+x*0;
  Sig_f =@(x) Sig_f_j(1)+x*0;
  nuSig_f =@(x) nuSig_f_j(1)+x*0;
  Sig_t =@(x) Sig_ss(x)+Sig_f(x)+Sig_gamma(x);
  
  phi0_MMS =@(x) 2*psi_MMS(x);
  % MMS source: mu_n * derivative(psi_MMS) +Sig_t* psi_MMS ...
  % -(Sig_ss+nuSig_f)*0.5*phi0_MMS;
  Q_MMS =@(x,mu) mu*psi_MMS_Diff(x) +Sig_t(x).*psi_MMS(x) ...
    -(Sig_ss(x)+nuSig_f(x))*0.5.*phi0_MMS(x);
  
  %% For MoC MMS solution and problem
  % Boundary condition and source
  % psi expression evaluated at x=0
  psi_b1_n=psi_MMS(0)*ones(N,1); % n=N/2+1:N % mu>0
  % psi expression evaluated at x=Tau
  psi_b2_n=psi_MMS(Tau)*ones(N,1); % n=1:N/2 % mu<0

  phi0_MMS_j=zeros(J,1);
  Q_MMS_j_n=zeros(J,N); % preallocate memory, avg'ed over tau_(j-1/2) and tau_(j+1/2)

  for j=1:J
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j(j)=1/h*integral(phi0_MMS,x_L,x_R);
    for n=1:N
      % g = @(c) (integral(@(x) (x.^2 + c*x + 1),0,1));
      Q_MMS_j_n(j,n)=integral(@(x) Q_MMS(x,mu_n(n)),x_L,x_R)/h;
    end % n
  end % j
  
  %% For TH MMS solution and problem
  % Boundary condition and source
  % Left boundary, T_MMS evaluated at x=0;
  T_L=T_MMS(0);
  % Right boundary, T_MMS evaluated at x=Tau;
  T_R=T_MMS(Tau);
  
  % Discretized MMS solution
  T_MMS_j=zeros(J,1);
  T_MMS_xx_j=zeros(J,1);

  % MMS source
  q_MMS_j=zeros(J,1);
  for j=1:J
    x_L=(j-1)*h;x_R=j*h;
    T_MMS_j(j)=1/h*integral(T_MMS,x_L,x_R);
    T_MMS_xx_j(j)=1/h*integral(T_MMS_xx,x_L,x_R);
    q_MMS_j(j)=k_F*T_MMS_xx_j(j)+kappaSig_f_j(j)*phi0_MMS_j(j);
  end  

end
