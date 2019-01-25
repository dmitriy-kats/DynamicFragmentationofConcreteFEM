%Initialization

lengthB=0.25;                    
strainrate=200;
%elem=1000;
nsteps=5E6 ;          
nstep_tolerance=8000;   
factor_timestep=0.1;

rho1=2400;      %density (kg/m3)
E=40E9;        %Young's Modulus Pa
TS1=3.5E6;      %Tensile Strength Pa
Gf1=34;         %Fracture Energy (N/m)

EE=E.*ones(1,elem+1);
EE(end)=NaN;
c=sqrt(E/rho1);     %wave speed
%Geometry
h=lengthB/elem;     %element size (m)
A=h;                %Area m^2 
dt=factor_timestep*lengthB/(elem)/c;  

SigmaC=wblrnd(TS1*2,5,1,elem+1)+TS1;  


deltaC=Gf1*2./SigmaC; %Critical opening size at each node
clear vars mL vL muL sigmaL


%% Node setup

nodes=1:1:(elem+1); %node numbers
coords = linspace(0,lengthB,(elem+1)); %coordinates of nodes
%Regular nodes
uu=zeros(1,elem+1); %Displacement              
vv=(strainrate*(coords-coords(end)/2));        %% david: eq. 14, Initial condition, coords(end) = lengthB
aa=zeros(1,elem+1); %acceleration


mmE=rho1*A*h*ones(1,elem+1); %element mass     
mm=mean([mmE(1:end-1);mmE(2:end)]);
mm=[mmE(1)/2 mm];                           %% david: mass of first  = mmE/2
mm(end)=mm(end)/2;                          %% david: mass of last = mmE/2



%Cohesive elements node
cohTracker=false(1,elem+1); 
uuc=zeros(1,elem+1); %Displacement Cohesive
%vvc=zeros(1,elem+1); %Velocity Cohesive          
vvc=(strainrate*(coords-coords(end)/2));          %% equivalent to vvc=zeros.
aac=zeros(1,elem+1); %Acceleration Cohesive
mmc=zeros(1,elem+1); %cohesive mass
%Cohesive opening
delta=zeros(1,elem+1);
deltamax=zeros(1,elem+1);
DD=zeros(1,elem+1);
TT=zeros(1,elem+1);
TTmax=zeros(1,elem+1);

AvgForce=zeros(1,elem+1);
NodeForce=zeros(1,elem+1);
NodeForceCoh=zeros(1,elem+1);


t=0:dt:nsteps*dt;

E_diss=zeros(1,nsteps+1); %Dissipiated energy by positive delta
E_rev=zeros(1,nsteps+1);%Reversible energy by positive delta
E_coh=zeros(1,nsteps+1);%Reversible energy by positive delta
W_rev=zeros(1,nsteps+1);%Reversible energy by positive delta
KE=zeros(1,nsteps+1); %Kinetic Energy
KEoriginal=sum(1/2.*vv.*mm.*vv); %Original Kinetic Energy
KE(1)=KEoriginal;
PE=zeros(1,nsteps+1); %Potential Energy
SE=zeros(1,nsteps+1); %Strain Energy
ExtWork=zeros(1,nsteps+1); %External Work

numFrag=zeros(1,nsteps+1);





