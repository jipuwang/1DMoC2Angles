% A sine and sine solution generator including
    % Discretized analytical solution
    % Manufactured boundary conditions
    % Manufactured source
function [phi0_MMS_j,psi_b1_n_i,psi_b2_n_i,Q_MMS_j_n_i]=...
          manufacturer_1d2angle(J,N,I,Tau,mat,assumedSoln)
  % input parameters
  if ~exist('J','var')
    J=5*2;%*2%*2*2*2*2*2*2*2*2
  end
  if ~exist('N','var')
    N=16;
  end
  if ~exist('I','var')
    I=16;
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
    assumedSoln='sine_sine_sine';
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
  [alpha_i,weight_i]=lgwt(I,0,2*pi);alpha_i=flipud(alpha_i);

  %% Manufactured Solutions for both fields
  % They need to be pre-defined here due to temperature dependence on the
  % xs. 
  % Options includes: sine_sine, const_cubic, sqrtPlus1_quadratic, etc.
  switch(assumedSoln)
    case('sine_sine_sine')
      % Manufactured neutronics solution \psi(x,\mu)=sin(pi*x/Tau), 0<x<Tau
      angleDep =@(mu,alpha) sin(mu+1).*sin(alpha/(2*pi));
      psi_MMS =@(x,mu,alpha) (sin(pi/Tau*x)+1).*angleDep(mu,alpha);
      psi_MMS_Diff =@(x,mu,alpha) pi/Tau*cos(pi/Tau*x).*sin(mu+1).*sin(alpha/(2*pi));
    case('const_cubic')
    case('sqrtPlus1_quadratic')
  end
  
  Sig_ss =@(x) Sig_ss_j(1)+x*0;
  Sig_gamma =@(x) mat.Sig_gamma_j(1)+0.0*x;
  Sig_f =@(x) Sig_f_j(1)+x*0;
  nuSig_f =@(x) nuSig_f_j(1)+x*0;
  Sig_t =@(x) Sig_ss(x)+Sig_f(x)+Sig_gamma(x);
  
  % phi0_MMS =@(x) integral2(@(mu,alpha) psi_MMS(x,mu,alpha), -1,1, 0,2*pi);
  % The above function handle can support fplot, and evaluation, but not
  % further integral.  So I have to do it by hand. 
  angleIntegral=integral2(angleDep, -1,1, 0,2*pi);
  phi0_MMS =@(x) (sin(pi/Tau*x)+1)*angleIntegral;
  % MMS source: mu_n * derivative(psi_MMS) +Sig_t* psi_MMS ...
  % -(Sig_ss+nuSig_f)*0.5*phi0_MMS;
  Q_MMS =@(x,mu,alpha) mu*psi_MMS_Diff(x,mu,alpha) ...
    +Sig_t(x).*psi_MMS(x,mu,alpha) ...
    -1/(4*pi)*Sig_ss(x).*phi0_MMS(x);
  % no fission source
%     -(Sig_ss(x)+nuSig_f(x))/(4*pi).*phi0_MMS(x);
  
  %% For MoC MMS solution and problem
  % Boundary condition and source
  psi_b1_n_i=zeros(N,I);
  % psi expression evaluated at x=0
  for n=1:N
    for i=1:I
      psi_b1_n_i(n,i)=psi_MMS(0,mu_n(n),alpha_i(i)); % n=N/2+1:N % mu>0
    end
  end
  % psi expression evaluated at x=Tau
  psi_b2_n_i=zeros(N,I);
  for n=1:N
    for i=1:I
      psi_b2_n_i(n,i)=psi_MMS(Tau,mu_n(n),alpha_i(i)); % n=N/2+1:N % mu>0
    end
  end
  
  phi0_MMS_j=zeros(J,1);
  Q_MMS_j_n_i=zeros(J,N,I);

  for j=1:J
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j(j)=1/h*integral(phi0_MMS,x_L,x_R);
    for n=1:N
      for i=1:I
        % g = @(c) (integral(@(x) (x.^2 + c*x + 1),0,1));
        Q_MMS_j_n_i(j,n,i)= ...
          integral(@(x) Q_MMS(x,mu_n(n),alpha_i(i)),x_L,x_R)/h;
      end % i
    end % n
  end % j
  
  figure(1);
  surf(psi_b1_n_i);
  figure(2);
  fsurf(@(mu,alpha) psi_MMS(0,mu,alpha),[-1 1 0 2*pi]);
  figure(3);
  surf(psi_b2_n_i);
  figure(4);
  fsurf(@(mu,alpha) psi_MMS(Tau,mu,alpha),[-1 1 0 2*pi]);
  figure(5);
  surf(0.5*(Q_MMS_j_n_i(5,:,:)+Q_MMS_j_n_i(6,:,:)));
  figure(6);
  fsurf(@(mu,alpha) Q_MMS(Tau/2,mu,alpha),[-1 1 0 2*pi]);
  
end
