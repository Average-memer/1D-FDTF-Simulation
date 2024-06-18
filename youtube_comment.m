% % INITIALIZE MATLAB 
close all; clc; clear all; 

%Units
seconds = 1;
hertz = 1/seconds;
megahertz = 1e+06 * hertz;
gigahertz = 1e+09 * hertz;
terahertz = 1e+12 * hertz;

meters = 1;
nanometers = 1e-09 * meters;

%%constants
c0 = 299792458 * meters/seconds;

% FDTD Parameters
Nz = 180;
dz = 1e-02;
STEPS = 2000;
Er = ones(1,Nz);
Ur = ones(1,Nz);
nzsrc = round(Nz/2); %source point
nbc=1;
nfreq = (1 * gigahertz)/(100*megahertz)+1;
fmax = 1 * gigahertz;

dt = dz/(2*c0);
tau = 0.5/fmax;
freq = linspace(0, fmax, nfreq);
t0 = 5*tau;
mHx = (dt*c0)./Ur;
mEy = (dt*c0)./Er; 
za = (0:Nz-1)*dz;
t = (0:STEPS-1)*dt;

%source, Total-field/scattered-field source formulation
A = sqrt(Er(nzsrc)/Ur(nzsrc));
st=(nbc*dz)/(2*c0)+ dt/2;
Esrc = exp(-((t-t0)/tau).^2);
Hsrc = A*exp(-((t-t0+st)/tau).^2);

%Init Fourier Transform
K = exp(-1i*2*pi*freq.*dt);
%Init Reflectance, Transitance, Source. 
%(Bounderies grid|REF grid|SRC grid|...inner grids points...|TRN grid|Bounderies grid)
REF = zeros(1,nfreq);
TRN = zeros(1,nfreq);
SRC = zeros(1,nfreq);

% COMPUTE DEFAULT GRID RESOLUTION
Hx = zeros(1,Nz);
Ey = zeros(1, Nz);
H3=0; H2=0; H1=0; E3=0; E2=0; E1=0;

figure('color','w'); 

% % MAIN FDTD LOOP 
for T = 1 : STEPS
% Update H from E (Perfect Boundary Conditions)
   for nz = 1 : Nz-1 
     Hx(nz) = Hx(nz) + mHx(nz)*(Ey(nz+1) - Ey(nz))/dz;
   end
     Hx(Nz) = Hx(Nz) + mHx(nz)*(E3 - Ey(Nz))/dz;
    
%Handle H‐Field Source
    Hx(nzsrc-1) = Hx(nzsrc-1) - mHx(nzsrc-1)*Esrc(T)/dz;

%PML,Record H at Boundary 
H3=H2; 
H2=H1; 
H1=Hx(1);

% Update E from H (Perfect Boundary Conditions) 
    Ey(1) = Ey(1) + mEy(nz)*(Hx(1) - H3)/dz;
   for nz = 2 : Nz
    Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1))/dz;
   end
 
%Handle E‐Field Source  
    Ey(nzsrc) = Ey(nzsrc) - mEy(nzsrc)*Hsrc(T)/dz;
   
%PML,Record E at Boundary 
E3=E2;
E2=E1; 
E1=Ey(Nz);
% Inject Soft Source 
 % Ey(nzsrc) = Ey(nzsrc) + g(T);
  
  
  % Update Fourier Transforms
  for nf = 1 : nfreq %loop over each frequency
    REF(nf) = REF(nf) + (K(nf)^T)*Ey(1);  %Fourier in one time (1*dt) of Ey
    TRN(nf) = TRN(nf) + (K(nf)^T)*Ey(Nz); %Fourier in one time (Nz*dt) of Ey
    SRC(nf) = SRC(nf) + (K(nf)^T)*Esrc(T);%Fourier in time (T) of Esrc 
  end    
  
if ~mod(T,1)
    %show field
 subplot(211);
plot(za,Ey,'-b');
hold off;
plot(za,Hx,'-r') 
 axis tight;
 xlim([dz Nz*dz]);
 ylim([-1 0.1])
 xlabel('z'); 
 title(['FDTD after ',num2str(T) , ' Iteration']); 

 % COMPUTE REFLECTANCE AND TRANSMITANCE, Temporary Normalization
Re = abs(REF./SRC).^2;
Tr = abs(TRN./SRC).^2;

 %show REF, TRN
 subplot(212);
 plot(freq,Re,'-b');
 hold off;
 plot(freq,Tr,'-r');
 plot(freq,Re+Tr,'-g');
 xlim([freq(1) freq(nfreq)]);
 ylim([-1 1]);
 xlabel('Frequency (Hz)'); 
 title('REFLECTANCE AND TRANSMITANCE'); 
 hold off;

 %draw graphic
 drawnow 
end
end