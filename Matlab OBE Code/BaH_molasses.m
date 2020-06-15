% clear all;
import bloch.*
import BaH.* %Small module that includes state generation, dipole transition matrix elements calculation and Zeeman effect calculation in Barium Hydride

warning('off','all');


%***********************PHYSICAL CONSTANTS********************************%

%Frequency units used in this calculation are GHz, and time units of ns.


hbar=1.054*10^(-34); %[Js]
k_b=1.381*10^(-23); % [J/K]
c=299792000; %[m/s]
eps_0=8.854*10^(-12); %[F/m]
a0=5.29*10^(-11); % [m]
q_e=1.602*10^(-19); % [C]
r_expval_APi=1.632*a0; % [m]
r_expval_BS=1.67*a0; % [m]
Bohr_mag=1.39962449/1000; %[GHz/G]
g_S=2.002319;



%*************************STATE GENERATION********************************%

[StX,StA,StB]=generateStates([1],1/2,1,[0]); %For reference, check +BaH/generateStates.m

nx=size(StX,1);
na=size(StA,1);
nb=size(StB,1);

n=nx+na+nb; %Total number of states

%*************************TRANSITION DIPOLES******************************%
disp('Transition matrix')

rabi_matrix=zeros(n,n,3);

ind=0;
% na=0;
for p=[-1,0,1] %Three basic polarizations: -1 - right circular, +1 - left circular, 0 - pi
    ind=ind+1;
    %X<->A
    for f=1:na     %f - excited electronic states
        for i=1:nx %i - ground electronic states
            rabi_matrix(i,nx+f,ind)=dipoleTransitionMatrixElement(StX(i,:),StA(f,:),p); %For reference, +BaH/dipoleTransitionMatrixElement.m
            rabi_matrix(nx+f,i,ind)=dipoleTransitionMatrixElement(StA(f,:),StX(i,:),p);            
        end
    end
    %X<->B
    for f=1:nb     %f - excited electronic states
        for i=1:nx %i - ground electronic states
            rabi_matrix(i,na+nx+f,ind)=dipoleTransitionMatrixElement(StX(i,:),StB(f,:),p); %For reference, +BaH/dipoleTransitionMatrixElement.m
            rabi_matrix(na+nx+f,i,ind)=dipoleTransitionMatrixElement(StB(f,:),StX(i,:),p);            
        end
    end
end


%*************************ZEEMAN EFFECT***********************************%

disp('Zeeman effect matrix')

zeeman_matrix=sym(zeros(n,n));

syms g_api g_xs_12 g_xs_32 g_bs g B real; %Symbolic variables for g-factor and magnetic field


%We calculate the Zeeman effect directly only for the ground state, which
%is in Hund's case (b) and for which we found appropriate formulas
for i=1:4     %f - final states
   for f=1:4 %i - initial states
       zeeman_matrix(i,f)=(-1)^(-1)*g*B*zeemanElement(StX(i,:),StX(f,:),-1)-(-1)^(1)*g*B*zeemanElement(StX(i,:),StX(f,:),1); %For reference, +BaH/zeemanElement.m
%         zeeman_matrix(i,i)=g_xs_12*B*StX(i,end); 
    end
end
for i=5:12     %f - final states
   for f=5:12 %i - initial states
       zeeman_matrix(i,f)=(-1)^(-1)*g*B*zeemanElement(StX(i,:),StX(f,:),-1)-(-1)^(1)*g*B*zeemanElement(StX(i,:),StX(f,:),1); %For reference, +BaH/zeemanElement.m
%         zeeman_matrix(i,i)=g_xs_32*B*StX(i,end); 
    end
end
for i=1:nb  %f - final states
   for f=1:nb %i - initial states
       zeeman_matrix(na+nx+i,na+nx+f)=(-1)^(-1)*g*B*zeemanElement(StB(i,:),StB(f,:),-1)-(-1)^(1)*g*B*zeemanElement(StB(i,:),StB(f,:),1); %For reference, +BaH/zeemanElement.m
%         zeeman_matrix(nx+na+i,nx+na+i)=g_bs*B*StB(i,end); 
    end
end
for i=1:na  %f - final states
   zeeman_matrix(nx+i,nx+i)=g_api*B*StA(i,end-1); 
end




%*************************BRANCHING RATIOS********************************%

BR=zeros(n);

%To calculate branching ratios, we first add squares of all transition dipoles
%(so transition strengths), and then divide appropriate transition
%strengths by that sum.

%Transition strengths
transition_strengths=zeros(n);
for i=1:n
    for f=1:n
        for p=1:3
            transition_strengths(i,f)=transition_strengths(i,f)+rabi_matrix(i,f,p)^2;
        end
    end
end

%Sums of transition strengths for a given initial state 'i'
for i=1:n
    sums=0;
    for f=1:n
        sums=sums+transition_strengths(i,f);
    end
    for f=1:n
        BR(i,f)=transition_strengths(i,f)/sums; %(rotational) branching ratio
    end
end

%Initial states don't decay, so we remove those terms. Otherwise, they
%would simply indicate fractional transition strengths.
for i=1:nx
    BR(i,:)=0;
end


%****************************DISSIPATOR***********************************%

disp('Dissipator')

L=Dissipator(n);

syms Ga Gb real;
assume(Ga,'positive')
assume(Gb,'positive')

%All excited states are assumed to have the same total decay rate G.
DR=zeros(1,nx);
for i=1:na
    DR=[DR,Ga];
end
for i=1:nb
    DR=[DR,Gb];
end

%We use branching ratios table to generate the whole dissipator, instead of
%adding decays one-by-one. That's also why we needed DR vector.
L.fromBranching(BR,DR);


%****************************HAMILTONIAN**********************************%

%Symbolic variables
syms w_0 w_2 w_3 w_a w_b real; %energies
syms w_L1 w_L2 w_L3 real; %light frequencies
syms W_L1 W_L2 real; %Rabi frequencies.
syms d_L1 d_L2 real; %detunings
syms D_0 D_2 D_a D_b real; %splittings
syms chi_1 chi_2 real; %phases
syms x v k_1 k_2 t real; %additional variables
syms D_L1 D_L2 D_L3 real; %Overall detunings
syms S_L1 S_L2 real; %Profile shifts


%Hamiltonian

disp('Hamiltonian')

H=Hamiltonian(n);

H.addEnergies([w_0,w_0+D_0,w_0+D_0,w_0+D_0,...        
    w_2,w_2,w_2,...
    w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,w_2+D_2,...    %w_3,w_3,w_3,w_3,w_3,...%     
    w_a,w_a+D_a,w_a+D_a,w_a+D_a,...
    w_b,w_b+D_b,w_b+D_b,w_b+D_b]);


%Couplings

%X(J=1/2)<->A
for i=1:4
    for f=nx+1:nx+na
        if rabi_matrix(i,f,2)~=0
            H.addPolyCoupling(i,f,W_L1*rabi_matrix(i,f,2)*(cos(k_1*x+chi_1/2+d_L1*t)-1i*sin(k_1*x+chi_1/2+d_L1*t)),w_L1+S_L1);
            H.addPolyCoupling(i,f,W_L1*rabi_matrix(i,f,2)*(cos(k_1*x-chi_1/2-d_L1*t)-1i*sin(k_1*x-chi_1/2-d_L1*t)),w_L1+S_L1);
            H.addPolyCoupling(i,f,W_L1*rabi_matrix(i,f,2)*(cos(-k_1*x+chi_1/2-d_L1*t)-1i*sin(-k_1*x+chi_1/2-d_L1*t)),w_L1-S_L1);
            H.addPolyCoupling(i,f,W_L1*rabi_matrix(i,f,2)*(cos(-k_1*x-chi_1/2+d_L1*t)-1i*sin(-k_1*x-chi_1/2+d_L1*t)),w_L1-S_L1);
            H.couplings(i,f,1)=w_L1;
        end
    end
end

pol_angle=0;

%X(J=3/2)<->B
for i=5:12
    for f=nx+na+1:n
        if rabi_matrix(i,f,1)~=0
            prefactor=(-1)^(-1)/sqrt(2)*sin(pol_angle);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(k_2*x+chi_2/2+d_L2*t)-1i*sin(k_2*x+chi_2/2+d_L2*t)),w_L2+S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(k_2*x-chi_2/2-d_L2*t)-1i*sin(k_2*x-chi_2/2-d_L2*t)),w_L2+S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(-k_2*x+chi_2/2-d_L2*t)-1i*sin(-k_2*x+chi_2/2-d_L2*t)),w_L2-S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(-k_2*x-chi_2/2+d_L2*t)-1i*sin(-k_2*x-chi_2/2+d_L2*t)),w_L2-S_L2);
%               H.addCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1),w_L2);              
              H.couplings(i,f,1)=w_L2;
        end
        if rabi_matrix(i,f,2)~=0
            prefactor=cos(pol_angle);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,2)*(cos(k_2*x+chi_2/2+d_L2*t)-1i*sin(k_2*x+chi_2/2+d_L2*t)),w_L2+S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,2)*(cos(k_2*x-chi_2/2-d_L2*t)-1i*sin(k_2*x-chi_2/2-d_L2*t)),w_L2+S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,2)*(cos(-k_2*x+chi_2/2-d_L2*t)-1i*sin(-k_2*x+chi_2/2-d_L2*t)),w_L2-S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,2)*(cos(-k_2*x-chi_2/2+d_L2*t)-1i*sin(-k_2*x-chi_2/2+d_L2*t)),w_L2-S_L2);
%               H.addCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,2),w_L2);              
              H.couplings(i,f,1)=w_L2;
        end
        if rabi_matrix(i,f,3)~=0
            prefactor=-(-1)^(1)/sqrt(2)*sin(pol_angle);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(k_2*x+chi_2/2+d_L2*t)-1i*sin(k_2*x+chi_2/2+d_L2*t)),w_L2+S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(k_2*x-chi_2/2-d_L2*t)-1i*sin(k_2*x-chi_2/2-d_L2*t)),w_L2+S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(-k_2*x+chi_2/2-d_L2*t)-1i*sin(-k_2*x+chi_2/2-d_L2*t)),w_L2-S_L2);
            H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(-k_2*x-chi_2/2+d_L2*t)-1i*sin(-k_2*x-chi_2/2+d_L2*t)),w_L2-S_L2);
%             H.addCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3),w_L2);
            H.couplings(i,f,1)=w_L2;
        end
    end
end

% for i=8:12
%     for f=nx+na+1:n
%         if rabi_matrix(i,f,1)~=0
%             prefactor=(-1)^(-1)/sqrt(2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(k_2*x+chi_2/2+d_L2*t)-1i*sin(k_2*x+chi_2/2+d_L2*t)),w_L2+S_L2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(k_2*x-chi_2/2-d_L2*t)-1i*sin(k_2*x-chi_2/2-d_L2*t)),w_L2+S_L2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(-k_2*x+chi_2/2-d_L2*t)-1i*sin(-k_2*x+chi_2/2-d_L2*t)),w_L2-S_L2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1)*(cos(-k_2*x-chi_2/2+d_L2*t)-1i*sin(-k_2*x-chi_2/2+d_L2*t)),w_L2-S_L2);
%               H.addCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,1),w_L3);              
%               H.couplings(i,f,1)=w_L3;
%         end
%         if rabi_matrix(i,f,3)~=0
%             prefactor=-(-1)^(1)/sqrt(2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(k_2*x+chi_2/2+d_L2*t)-1i*sin(k_2*x+chi_2/2+d_L2*t)),w_L2+S_L2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(k_2*x-chi_2/2-d_L2*t)-1i*sin(k_2*x-chi_2/2-d_L2*t)),w_L2+S_L2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(-k_2*x+chi_2/2-d_L2*t)-1i*sin(-k_2*x+chi_2/2-d_L2*t)),w_L2-S_L2);
% %             H.addPolyCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3)*(cos(-k_2*x-chi_2/2+d_L2*t)-1i*sin(-k_2*x-chi_2/2+d_L2*t)),w_L2-S_L2);
%             H.addCoupling(i,f,prefactor*W_L2*rabi_matrix(i,f,3),w_L3);
%             H.couplings(i,f,1)=w_L3;
%         end
%     end
% end


H.defineEnergyDetuning(w_0,w_a,D_L1,w_L1);
% H.defineStateDetuning(10,20,D_L3);

H.hamiltonian=H.hamiltonian+zeeman_matrix;

H.defineZero(w_0);
H.unitaryTransformation();
H.defineEnergyDetuning(w_2,w_b,D_L2,w_L2);
H.defineZero(w_2);
H.unitaryTransformation();



for i=1:n
    for f=1:n
        if i~=f
            H.transformed(i,f)=simplify(expand(H.transformed(i,f)),'Steps',1000);
        end
    end
end

H.takeGradient(x);

H.transformed=simplify(subs(H.transformed,x,v*t));
DH=simplify(subs(H.hamGradient,x,v*t));

disp(H.transformed)

%%
%**************************NUMERICAL VALUES*******************************%

%Splittings. Everything is in units of 2\pi*frequency
Del_0=-2*pi*2*10^(-3); 
Del_a=2*pi*2*10^(-3);
Del_2=2*pi*39*10^(-3);
Del_b=-2*pi*52*10^(-3);

%Decay rates
Gamma_APi=0.98772/136.5; %1/Lifetime of the excited state [GHz]
Gamma_BS=0.95312/125.1; %1/Lifetime of the excited state [GHz]

%k-values
k_api=2*pi/1060.7868; %1/nm
k_bs=2*pi/905.3197; %1/nm

%Detunings
Det_L1=0;
Det_L2=Del_2/2-Del_b; %Perhaps it should be (Del_2-Del_b)/2
Det_L3=0;%(Del_2-Del_b);
% det_L1=100*Gamma_APi;
% det_L2=det_L1*k_api*Gamma_BS/(k_bs*Gamma_APi);
det_L2=300*Gamma_BS;
det_L1=det_L2*k_bs*Gamma_APi/(k_api*Gamma_BS);
Shift_L1=-20*Gamma_APi;
Shift_L2=(det_L2/det_L1)*(k_api/k_bs)*(-Shift_L1);

%Magnetic field
g_fact=g_S*Bohr_mag;
g_eff_APi=-0.51*Bohr_mag;
g_eff_BS=2.54*Bohr_mag;
g_eff_XS_12=-0.65*Bohr_mag;
g_eff_XS_32=0.56*Bohr_mag;
B_field=2*pi*10; %[G]


%Average Rabi rates [GHz]
Rabi_L1=sqrt(3/2)*det_L1;
Rabi_L2=sqrt(3/2)*det_L2;
% Rabi_L2=20*Gamma_BS;

disp('Rabi rates')
disp('R_L1 [Gamma_APi]')
disp(vpa(Rabi_L1/Gamma_APi,5))
disp('R_L2 [Gamma_BSigma]')
disp(vpa(Rabi_L2/Gamma_BS,5))

disp('Detunings')
disp('d_L1 [Gamma_APi]')
disp(vpa(det_L1/Gamma_APi,5))
disp('d_L2 [Gamma_BSigma]')
disp(vpa(det_L2/Gamma_BS,5))

disp('Spacing')
disp('Delta_2+Delta_b [Gamma_BS]')
disp(vpa((Del_2-Del_b)/Gamma_BS,5))

%Phases
angle_1=45/180*pi;
angle_2=-45/180*pi;


t_equil=-15/Gamma_APi; %[ns]


%%
%*************************INITIAL CONDITIONS******************************%

IC=zeros(n);

%We simply assume that 1/4 of population is in J=1/2 and 3/4 in J=3/2. Both
%equally distributed among the Zeeman sublevels
for i=1:4
    IC(i,i)=1/16;
end
for i=5:12
    IC(i,i)=3/32;
end


%*************************MASTER EQUATION*********************************%

disp('Optical Bloch Equations')
Eq=BlochEqns(H,L);
%%
% Eq.eqnsRHS=subs(Eq.eqnsRHS,[D_0,D_2,D_a,D_b,D_L1,D_L2,D_L3,Ga,Gb,g_xs_12,g_xs_32,g_api,g_bs,k_1,k_2],...
%     [Del_0,Del_2,Del_a,Del_b,Det_L1,Det_L2,Det_L3,Gamma_APi,Gamma_BS,g_eff_XS_12,g_eff_XS_32,g_eff_APi,g_eff_BS,k_api,k_bs]);
% DH=subs(DH,[D_0,D_2,D_a,D_b,D_L1,D_L2,D_L3,Ga,Gb,g_xs_12,g_xs_32,g_api,g_bs,k_1,k_2],...
%     [Del_0,Del_2,Del_a,Del_b,Det_L1,Det_L2,Det_L3,Gamma_APi,Gamma_BS,g_eff_XS_12,g_eff_XS_32,g_eff_APi,g_eff_BS,k_api,k_bs]);
Eq.eqnsRHS=subs(Eq.eqnsRHS,[D_0,D_2,D_a,D_b,D_L1,D_L2,D_L3,Ga,Gb,g,g_api,k_1,k_2],...
    [Del_0,Del_2,Del_a,Del_b,Det_L1,Det_L2,Det_L3,Gamma_APi,Gamma_BS,g_fact,g_eff_APi,k_api,k_bs]);
DH=subs(DH,[D_0,D_2,D_a,D_b,D_L1,D_L2,D_L3,Ga,Gb,g,g_api,k_1,k_2],...
    [Del_0,Del_2,Del_a,Del_b,Det_L1,Det_L2,Det_L3,Gamma_APi,Gamma_BS,g_fact,g_eff_APi,k_api,k_bs]);


Eq.necessaryVariables();
disp(symvar(DH))

%%
%*************************VELOCITY PROFILE********************************%
import bloch.*
 DHs=vpa(simplify(subs(DH,[B,S_L1,S_L2,W_L1,W_L2,chi_1,chi_2,d_L1,d_L2],...
     [B_field,Shift_L1,Shift_L2,Rabi_L1,Rabi_L2,angle_1,angle_2,det_L1,det_L2]))); 
 Fs=[];
 step=10;
figure
u=0;
for velocity=-30:2:30

    tic
    u=u+1;

    Forces=[];

    if Shift_L1-k_api*velocity==0 || Shift_L2-k_bs*velocity==0
        t_end=2*pi*(3*200/det_L1);
    else
        t_end=2*2*pi*(1/(abs(Shift_L1-k_api*velocity))+1/(abs(Shift_L2-k_bs*velocity))+100/det_L1+100/det_L2);
    end
    
     
    Eq.evolve(t_equil,t_end,IC,[B_field,Shift_L1,Shift_L2,Rabi_L1,Rabi_L2,angle_1,angle_2,det_L1,det_L2,velocity]);
%     Eq.evolve(t_equil,t_end,IC,[B_field,Shift_L1,Rabi_L1,Rabi_L2,angle_1,det_L1,velocity]);
    i0=find(Eq.evTime(:)>0,1);
        
    av_e=Eq.evolution(nx+1,nx+1,i0:end);
    for ii=nx+2:n
        av_e=av_e+Eq.evolution(ii,ii,i0:end);
    end
    av_e=mean(av_e);
    
    for i=i0:step:length(Eq.evTime(:)) 
        gradH=double(subs(DHs,[t,v],[Eq.evTime(i),velocity]));
        force=-real(trace(squeeze(Eq.evolution(:,:,i))*gradH));
        Forces=[Forces,force];
    end

    tot_force=double(trapz(Eq.evTime(i0:step:i),Forces)/(Eq.evTime(i)-Eq.evTime(i0)));
    tot_force=tot_force/(k_api*Gamma_APi*0.5);
    fprintf('Velocity %.2f m/s\n',velocity);
    fprintf('p_ee=%.4f \n',av_e);
    fprintf('F=%.4f hkG/2 \n',tot_force);
    fprintf('Done %.2f%% \n',u/(31)*100);
    
%     Eq.plotEvolution();


    Fs=[Fs,tot_force];

    clf
    plot(-30:2:velocity,Fs)
    xlim([-30 30])
    ylim([-40 40])
    drawnow
    toc
%     
    
end
%%
figure 

plot(Eq.evTime(:),squeeze(Eq.evolution(13,13,:))+squeeze(Eq.evolution(14,14,:))+squeeze(Eq.evolution(15,15,:))+squeeze(Eq.evolution(16,16,:))) %APi
hold on
% plot(Eq.evTime(:),squeeze(Eq.evolution(17,17,:))+squeeze(Eq.evolution(18,18,:))+squeeze(Eq.evolution(19,19,:))+squeeze(Eq.evolution(20,20,:))) %BS

% hold on
% plot(Eq.evTime(:),squeeze(Eq.evolution(17,17,:))) %BS F=0
% % plot(Eq.evTime(:),squeeze(Eq.evolution(18,18,:))+squeeze(Eq.evolution(19,19,:))+squeeze(Eq.evolution(20,20,:))) %BS F=1

% plot(Eq.evTime(:),squeeze(Eq.evolution(5,5,:))+squeeze(Eq.evolution(6,6,:))+squeeze(Eq.evolution(7,7,:))); %XS J=3/2 F=1 
% plot(Eq.evTime(:),squeeze(Eq.evolution(8,8,:))+squeeze(Eq.evolution(9,9,:))+squeeze(Eq.evolution(10,10,:))+squeeze(Eq.evolution(11,11,:))+squeeze(Eq.evolution(12,12,:))); %XS J=3/2 F=2

% plot(Eq.evTime(:),squeeze(Eq.evolution(1,1,:))+squeeze(Eq.evolution(2,2,:))+squeeze(Eq.evolution(3,3,:))+squeeze(Eq.evolution(4,4,:))) %XS J=1/2
% hold on
% plot(Eq.evTime(:),squeeze(Eq.evolution(5,5,:))+squeeze(Eq.evolution(6,6,:))+squeeze(Eq.evolution(7,7,:))+squeeze(Eq.evolution(8,8,:))+squeeze(Eq.evolution(9,9,:))+squeeze(Eq.evolution(10,10,:))+squeeze(Eq.evolution(11,11,:))+squeeze(Eq.evolution(12,12,:))) %XS J=3/2
%%
% %%
% Eq.plotEvolution()
% % 
%%
figure
plot(Forces/(k_api*Gamma_APi*0.5))
%%
mean(Forces/(k_api*Gamma_APi*0.5))
%%
DHss=vpa(subs(DHs,[v],[-20]),4)
% Eqq=squeeze(Eq.evolution(:,:,30000))
% Res=Eqq*DHss
%%
figure
Om=DHs
fplot(imag(Om))
xlim([0,2000])
%%
vpa(simplify(real(Om),'Steps',200),4)

%%





x=(2^(1/2)*W_L1*cos(k_1*t*v)*sin(chi_1/2)*sin(S_L1*t)*sin(d_L1*t)*2i)/3 - (2^(1/2)*W_L1*sin(k_1*t*v)*sin(chi_1/2)*cos(S_L1*t)*sin(d_L1*t)*2i)/3;
x=simplify(x)
