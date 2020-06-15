clear;
import bloch.*

warning('off','all');


% Hamiltonian
syms w1 w2 wa we1 we2 wb W1 W2 d1 d2 D1 D2 D01 D02 Dx v k real;
syms chi1 chi2 t x real ;
n=4;
h=Hamiltonian(n);

h.addEnergies([wa,wb,we1,we2]);

k_f=1;

h.addPolyCoupling(1,3,-W1*(cos(x+chi1/2+d1*t)-1i*sin(x+chi1/2+d1*t)),w1+D01);
h.addPolyCoupling(1,3,-W1*(cos(x-chi1/2-d1*t)-1i*sin(x-chi1/2-d1*t)),w1+D01);
h.addPolyCoupling(1,3,-W1*(cos(-x+chi1/2-d1*t)-1i*sin(-x+chi1/2-d1*t)),w1-D01);
h.addPolyCoupling(1,3,-W1*(cos(-x-chi1/2+d1*t)-1i*sin(-x-chi1/2+d1*t)),w1-D01);


h.addPolyCoupling(2,4,-W2*(cos(k_f*x+chi2/2+d2*t)-1i*sin(k_f*x+chi2/2+d2*t)),w2+D02);
h.addPolyCoupling(2,4,-W2*(cos(k_f*x-chi2/2-d2*t)-1i*sin(k_f*x-chi2/2-d2*t)),w2+D02);
h.addPolyCoupling(2,4,-W2*(cos(-k_f*x+chi2/2-d2*t)-1i*sin(-k_f*x+chi2/2-d2*t)),w2-D02);
h.addPolyCoupling(2,4,-W2*(cos(-k_f*x-chi2/2+d2*t)-1i*sin(-k_f*x-chi2/2+d2*t)),w2-D02);

% % 


h.addPolyCoupling(1,3,-W1*(cos(x+3*chi1/2+3*d1*t)-1i*sin(x+3*chi1/2+3*d1*t)),w1+D01);
h.addPolyCoupling(1,3,-W1*(cos(x-3*chi1/2-3*d1*t)-1i*sin(x-3*chi1/2-3*d1*t)),w1+D01);
h.addPolyCoupling(1,3,-W1*(cos(-x+3*chi1/2-3*d1*t)-1i*sin(-x+3*chi1/2-3*d1*t)),w1-D01);
h.addPolyCoupling(1,3,-W1*(cos(-x-3*chi1/2+3*d1*t)-1i*sin(-x-3*chi1/2+3*d1*t)),w1-D01);

h.addPolyCoupling(2,4,-W2*(cos(k_f*x+3*chi2/2+3*d2*t)-1i*sin(k_f*x+3*chi2/2+3*d2*t)),w2+D02);
h.addPolyCoupling(2,4,-W2*(cos(k_f*x-3*chi2/2-3*d2*t)-1i*sin(k_f*x-3*chi2/2-3*d2*t)),w2+D02);
h.addPolyCoupling(2,4,-W2*(cos(-k_f*x+3*chi2/2-3*d2*t)-1i*sin(-k_f*x+3*chi2/2-3*d2*t)),w2-D02);
h.addPolyCoupling(2,4,-W2*(cos(-k_f*x-3*chi2/2+3*d2*t)-1i*sin(-k_f*x-3*chi2/2+3*d2*t)),w2-D02);




h.couplings(1,3,1)=w1;
h.couplings(2,4,1)=w2;
% 
h.defineEnergyDetuning(wa,we1,D1,w1);


h.defineZero(wa);
h.unitaryTransformation();

h.defineZero(wb);
h.unitaryTransformation();
h.defineEnergyDetuning(wb,we2,D2,w2);
H=h.transformed;

DH=sym('DH',[n n]);
for i=1:n
    for j=1:n
        DH(i,j)=simplify(diff(H(i,j),x));
    end
end


h.transformed=simplify(subs(h.transformed,x,v*t));
disp(h.transformed)

H=h.transformed;



DH=simplify(subs(DH,x,v*t));
disp(DH)

%%
import bloch.*
L=Dissipator(n);

syms G1 G2 real;
assume(G1,'positive')
assume(G2,'positive')

L.addDecay(3,1,2*G1/3);
L.addDecay(4,2,2*G2/3);
L.addDecay(3,2,G1/3);
L.addDecay(4,1,G2/3);


%%

IC=zeros(4);
IC(1,1)=1/2;
IC(2,2)=1/2;
%%
Gamma1=1;
Gamma2=1
g_f=Gamma2/Gamma1;

det1=100*Gamma1;
det2=det1/(g_f*k_f);

off1=-20;
off2=20;
R1=1.16*det1;
R2=1.16*det2;

angle1=(25/180)*pi;
angle2=(-25/180)*pi;
rad_factor1=0.21;
rad_factor2=0.195;
t_equil=-15;
%%
import bloch.*
eq=BlochEqns(h,L);
DH=subs(DH,[d1,chi1,W1,d2,chi2,W2,D1,D2,D01,D02],[det1,angle1,R1,det2,angle2,R2,0,0,off1,off2]);
eq.eqnsRHS=subs(eq.eqnsRHS,[G1,G2,d1,chi1,W1,d2,chi2,W2,D1,D2,D01,D02],[Gamma1,Gamma2,det1,angle1,R1,det2,angle2,R2,0,0,off1,off2]);
eq.necessaryVariables();
disp(symvar(DH))
return
%%


    j=0;
    Fs=[];
    step=10;
    figure
    for velocity=-70:0.35:0

        Forces=[];

        j=j+1;
        disp(j)
      
        if off1-velocity==0 || off2-velocity==0
            t_end=2*pi*(4*200/det1);
        else
            t_end=4*pi*(1/(abs(off1-velocity))+1/(abs(off2-velocity))+200/det1);
        end
        eq.evolve(t_equil,t_end,IC,[velocity]);
        i0=find(eq.evTime(:)>0,1);
        av_e=squeeze(eq.evolution(3,3,:)+eq.evolution(4,4,:));
        
       

        for i=i0:step:length(eq.evTime(:)) 

            gradH=double(subs(DH,[t,v],[eq.evTime(i),velocity]));

            force=-real(trace(squeeze(eq.evolution(:,:,i))*gradH));

            Forces=[Forces,force];
            
        end

       tot_force=2*double(real(trapz(eq.evTime(i0:step:i),Forces)/(eq.evTime(i)-eq.evTime(i0))));
       tot_e=double(real(trapz(eq.evTime(i0:step:i),av_e(i0:step:i)))/(eq.evTime(i)-eq.evTime(i0)))
       disp(tot_e)
       disp(tot_force)

        
       
       fprintf('Done %.2f%% \n',j/(201)*100);
        
       
    Fs=[Fs,tot_force];

    clf

    plot(-70:0.35:velocity,Fs)
    xlim([-70 0])
    ylim([-10 75])
    ylabel('Force [F_{rad}]')
    xlabel('Velocity [\Gamma_A/k_A]')

    drawnow

    end
    return

    %%
    Fsfl=fliplr(Fs);
    Fs1=[Fs(:)',-Fsfl(2:end)];
    figure
    plot(-70:0.35:70,Fs1)
    xlim([-70,70])
    %%
    Fs=Fs1;
    figure
    plot(-45:0.2:45,movmean(Fs,20))
    hold on
    plot(-45:1:45,zeros(71),'k--')
    xlim([-45 45])
    ylim([-120 120])
    ylabel('Force [F_{rad}]')
    xlabel('Velocity [\Gamma_A/k_A]')

%%
Fs=Fs1;
save('TCF_molasses_det35.mat','Fs');

    