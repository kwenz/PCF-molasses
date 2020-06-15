clear all;
import bloch.*

warning('off','all');


% Hamiltonian

syms w1 w2 wa we W1 W2 W0 d1 d2 D0 D1 D00 D d v k real;
syms chi chi1 chi2 t x y phi phi3 a b real;
syms W;
n=2;
h=Hamiltonian(n);

h.addEnergies([wa,we]);

h.addPolyCoupling(1,2,-W1*(cos(x+chi1/2+d1*t)-1i*sin(x+chi1/2+d1*t)),w1+D0);
h.addPolyCoupling(1,2,-W1*(cos(x-chi1/2-d1*t)-1i*sin(x-chi1/2-d1*t)),w1+D0);
h.addPolyCoupling(1,2,-W1*(cos(-x+chi1/2-d1*t)-1i*sin(-x+chi1/2-d1*t)),w1-D0);
h.addPolyCoupling(1,2,-W1*(cos(-x-chi1/2+d1*t)-1i*sin(-x-chi1/2+d1*t)),w1-D0);

% h.addPolyCoupling(1,2,-W1*(cos(x+chi2/2+d1*t)-1i*sin(x+chi2/2+d1*t)),w1+D1);
% h.addPolyCoupling(1,2,-W1*(cos(x-chi2/2-d1*t)-1i*sin(x-chi2/2-d1*t)),w1+D1);
% h.addPolyCoupling(1,2,-W1*(cos(-x+chi2/2-d1*t)-1i*sin(-x+chi2/2-d1*t)),w1-D1);
% h.addPolyCoupling(1,2,-W1*(cos(-x-chi2/2+d1*t)-1i*sin(-x-chi2/2+d1*t)),w1-D1);
% 
% X=simplify(W0*(cos(x+chi/2+d*t)-1i*sin(x+chi/2+d*t))+W0*(cos(x-chi/2-d*t)-1i*sin(x-chi/2-d*t))+...
%     W0*(cos(-x+chi/2-d*t)-1i*sin(-x+chi/2-d*t))+W0*(cos(-x-chi/2+d*t)-1i*sin(-x-chi/2+d*t)),'Steps',500);

% h.addPolyCoupling(1,2,-W,w1);
% h.addPolyCoupling(1,2,-W1,w1+D0);
% h.addPolyCoupling(1,2,-W1,w1-D0);
% h.addPolyCoupling(1,2,-W1,w1-D0);

% h.addPolyCoupling(1,2,W2*(cos(y-chi1/2)-1i*sin(y-chi1/2))*exp(1i*phi),w1-d2+D0-Dy);
% h.addPolyCoupling(1,2,W2*(cos(y+chi1/2)-1i*sin(y+chi1/2))*exp(1i*phi),w1+d2+D0-Dy);
% h.addPolyCoupling(1,2,W2*(cos(-y-chi1/2)-1i*sin(-y-chi1/2))*exp(1i*phi),w1+d2+D0-Dy);
% h.addPolyCoupling(1,2,W2*(cos(-y+chi1/2)-1i*sin(-y+chi1/2))*exp(1i*phi),w1-d2+D0-Dy);

% h.addPolyCoupling(1,2,W2*(cos(x-chi2/2)-1i*sin(x-chi2/2)),w1-d2-D0);
% h.addPolyCoupling(1,2,W2*(cos(x+chi2/2)-1i*sin(x+chi2/2)),w1+d2-D0);
% h.addPolyCoupling(1,2,W2*(cos(-x-chi2/2)-1i*sin(-x-chi2/2)),w1+d2-D0);
% h.addPolyCoupling(1,2,W2*(cos(-x+chi2/2)-1i*sin(-x+chi2/2)),w1-d2-D0);

% % % 
% h.addPolyCoupling(1,2,W1*(cos(x-3*chi1/2)-1i*sin(x-3*chi1/2)),w1-3*d1+D0+Dx);
% h.addPolyCoupling(1,2,W1*(cos(x+3*chi1/2)-1i*sin(x+3*chi1/2)),w1+3*d1+D0+Dx);
% h.addPolyCoupling(1,2,W1*(cos(-x-3*chi1/2)-1i*sin(-x-3*chi1/2)),w1+3*d1-D0-Dx);
% h.addPolyCoupling(1,2,W1*(cos(-x+3*chi1/2)-1i*sin(-x+3*chi1/2)),w1-3*d1-D0-Dx);
% 
% h.addPolyCoupling(1,2,W1*(cos(x-5*chi/2)-1i*sin(x-5*chi/2)),w1-5*d1);
% h.addPolyCoupling(1,2,W1*(cos(x+5*chi/2)-1i*sin(x+5*chi/2)),w1+5*d1);
% h.addPolyCoupling(1,2,W1*(cos(-x-5*chi/2)-1i*sin(-x-5*chi/2)),w1+5*d1);
% h.addPolyCoupling(1,2,W1*(cos(-x+5*chi/2)-1i*sin(-x+5*chi/2)),w1-5*d1);

% h.addPolyCoupling(1,2,W1*(cos(x-7*chi/2)+1i*sin(x-7*chi/2))*exp(-1i*phi),w1-7*d1);
% h.addPolyCoupling(1,2,W1*(cos(x+7*chi/2)+1i*sin(x+7*chi/2))*exp(-1i*phi),w1+7*d1);
% h.addPolyCoupling(1,2,W1*(cos(-x-7*chi/2)+1i*sin(-x-7*chi/2))*exp(-1i*phi),w1+7*d1);
% h.addPolyCoupling(1,2,W1*(cos(-x+7*chi/2)+1i*sin(-x+7*chi/2))*exp(-1i*phi),w1-7*d1);



h.couplings(1,2,1)=w1;
% 
h.defineEnergyDetuning(wa,we,D,w1);


h.defineZero(wa);
disp(h.hamiltonian)
h.unitaryTransformation();

H=simplify(h.transformed,'Steps',500);

DH=sym('DHx',[n n]);

for i=1:n
    for j=1:n
        DH(i,j)=simplify(diff(H(i,j),x));

    end
end
% for i=1:n
%     for j=1:n
%         if i~=j
%             DH(i,j)=1/2*simplify(diff(X,x));
%         end
% 
%     end
% end

h.transformed=simplify(subs(h.transformed,[x,D0],[v*t,0]),'Steps',200);
% disp(h.transformed)

H=h.transformed;


DH=simplify(subs(DH,[x,D0],[v*t,0]));
disp(DH)

% X=subs(X,x,v*t);

%%
import bloch.*
L=Dissipator(2);

syms G G2 real;
assume(G,'positive')
assume(G2,'positive')

L.addDecay(2,1,G);

%%

IC=zeros(2);
IC(1,1)=3/4;
IC(2,2)=1/4;
%%
Gamma=1;
det1=100*Gamma;
R1=sqrt(3/2)*det1;
angle=(45/180)*pi;
shift=0;
center_detuning=0;

t_equil=-15;
%%
import bloch.*
eq=BlochEqns(h,L);
%%
disp('Hamiltonian')
disp(h.transformed)
% fprintf('with W = %s \n\n',char(X))
disp('Gradient of the Hamiltonian')
disp(DH)
disp('Dissipator') 
disp(L.dissipator)
disp('Equations')
disp(eq.equationsVector)
% disp(eq.eqnsRHS)

DH=vpa(subs(DH,[d1,W1,chi1],[det1,R1,angle,]));

eq.eqnsRHS=vpa(subs(eq.eqnsRHS,[G,D,d1,W1,chi1],[Gamma,center_detuning,det1,R1,angle]));
eq.necessaryVariables();
% Wv=subs(2*h.transformed(1,2),[G,D,d1,W1,chi1],[Gamma,center_detuning,det1,R1,angle]);

disp(symvar(DH))
% disp(DH)


    %% Velocity profile
   
    k=0;
    Fs=[];
    step=10;
    figure

    for velocity=0.2:0.2:50
  
        k=k+1;
        Forces=[];
        Omegas=[];
        DHv=subs(DH,v,velocity);
      
        if velocity==0
            t_end=2*pi*(4*200/det1);
        else
            t_end=4*pi*(1/abs(velocity)+200/det1);
        end
        
        eq.evolve(t_equil,t_end,IC,[velocity]);
        
        i0=find(eq.evTime(:)>0,1);
        av_e=mean(eq.evolution(2,2,i0:end));
    
        for i=i0:step:length(eq.evTime(:))        
            gradH=subs(DHv,t,eq.evTime(i));
%             Omega_t=subs(Wv,[t,v],[eq.evTime(i),velocity]);
%             Omegas=[Omegas,Omega_t];
              force=-real(double(trace(squeeze(eq.evolution(:,:,i))*gradH)));
%             force=double(2*real(gradH(2,1)*squeeze(eq.evolution(1,2,i))));
      
            Forces=[Forces,force];
        end

       tot_force=2*double(real(trapz(eq.evTime(i0:step:i),Forces)/(eq.evTime(i)-eq.evTime(i0))));

       disp(av_e)
       disp(tot_force)

       
    Fs=[Fs,tot_force];

    clf
    plot(0.2:0.2:velocity,Fs)
    xlim([0 50])
    ylim([-5 70])
    drawnow
    
   

    end

return
%%
Fs_final=[fliplr(Fs),Fs(1),Fs];
%%
figure
plot(-50:0.2:50,movmean(Fs,20))
xlabel('Velocity [\Gamma/k]')
ylabel('Force [hbark \Gamma/2]')


%%
Force=Fs;
save('BCF_fv_Delta_pm20_chi45_omega122.mat','Force')

%% Rabi rate profile
  Ints=[];
    Fs=[];
    t_end=2*pi*(1/(abs(velocity))+200/det1);
    figure
    for Rabi=0:1:200

        disp(Rabi)
        Forces1=[];

        DH1=subs(DHx,[W1],[Rabi]);

        eq.evolve(t_equil,t_end,IC,[Rabi]);
        i0=find(eq.evTime(:)>0,1);

        av_e=mean(eq.evolution(2,2,i0:end));
%         pope(k,j)=av_e;
        for i=i0:50:length(eq.evTime(:)) %
            gradH=double(subs(DH1,[t],[eq.evTime(i)]));
            force=-(real(trace(squeeze(eq.evolution(:,:,i))*gradH)));
            Forces1=[Forces1,force];

        end

       tot_force=double(trapz(eq.evTime(i0:50:end),Forces1)/(eq.evTime(end)));
       
       tot_force=(tot_force*2);
%        disp(av_e)
%        disp(tot_force)

       fprintf('Done %.2f%% \n',Rabi/200*100);
        
       
    Fs=[Fs,tot_force];
    clf
    plot(0:1:Rabi,Fs)
    xlim([0 200])
    ylim([0 160])
    xlabel('\Omega [\Gamma]')
    drawnow

    end
    
%%
figure
plot(0:1:Rabi,Fs)
hold on 
plot(0:10:200,zeros(21),'k--')
    xlim([0 200])
    ylim([-10 130])
    xlabel('\Omega [\Gamma]')
    drawnow
    
    
%%
pop=zeros(31,21);
pope=zeros(31,21);
k=0;
Fs=[];
for R1=90:2:150
    k=k+1;
    j=0;
    disp(R1)
    
    for an=35:1:55
        angle=(an/180)*pi;
        DH1=subs(DHx,[chi1,W1],[angle,R1]);
        Forces1=[];
        j=j+1;
        disp(an)
      
        eq.evolve(t_equil,t_end,IC,[R1,angle]);
        i0=find(eq.evTime(:)>0,1);

        av_e=mean(eq.evolution(2,2,i0:end));
        pope(k,j)=av_e;
        for i=i0:50:length(eq.evTime(:)) %
            gradH=double(subs(DH1,[t],[eq.evTime(i)]));
            force=-2*real(trace(squeeze(eq.evolution(:,:,i))*gradH));
            Forces1=[Forces1,force];

        end

       tot_force=double(trapz(eq.evTime(i0:50:end),Forces1)/(eq.evTime(end)));
       disp(av_e)
       disp(tot_force)
%        fprintf('Done %.2f%% \n',j/(6*9)*100);
        
       
       
        pop(k,j)=tot_force;
        pope(k,j)=real(av_e);

    end
end
%%
% Multiple velocity profiles

Rabs=[108,133];
Ang=[28,21];
Fs=zeros(2,601);
Pe=zeros(2,601);
    figure
for j=1:2
    disp(j)
    R1=Rabs(j);
%     R1=110;
    angle=Ang(j);
%     angle=48;
    angle=(angle/180)*pi;
    
    k=0;
    for velocity=0:0.1:60
        k=k+1;
        com_det=-velocity;
        
        DH1=subs(DH,[chi,W1],[angle,R1]);
        Forces1=[];
    
      
        if velocity==0
            t_end=2*pi*(5*200/det1);
        else
            t_end=2*pi*(1/(abs(velocity))+200/det1);
        end
        eq.evolve(t_equil,t_end,IC,[com_det,R1,angle,velocity]);
        i0=find(eq.evTime(:)>0,1);

        av_e=mean(eq.evolution(2,2,i0:end));
    
        for i=i0:50:length(eq.evTime(:)) %
            gradH=double(subs(DH1,[t,v],[eq.evTime(i),velocity]));
            force=-2*real(trace(squeeze(eq.evolution(:,:,i))*gradH));
            Forces1=[Forces1,force];
        end

       tot_force=double(trapz(eq.evTime(i0:50:end),Forces1)/(eq.evTime(end)));
       


       fprintf('Done %.2f%% \n',(velocity+0.1)/(60.1)*100);
        
       
    Fs(j,k)=tot_force;
    Pe(j,k)=real(av_e);
    
    clf
    for x=1:j
        plot(0:0.1:60,Fs(x,:))
        hold on
    end
    xlim([0 60])
    ylim([0 130])
    drawnow

    end
end




%%

figure
surf(90:2:150,35:1:55,real(pope'),'EdgeColor','none','facecolor','interp')
colormap('jet')
ylabel('\chi [deg]')
xlabel('\Omega_0 [\Gamma]')
view(2)
c=colorbar;
c.Label.Interpreter = 'latex';
c.Label.String='F [$\hbar k \Gamma/2$]';
c.Label.FontSize=14;
title(['v = 1 \Gamma/k']);

drawnow


%     
%%
figure
corrs=[];
for v=0:1:59
    A=Pees(:,v+1);%(Pees1(1:end,v+1)-mean(Pees1(1:end,v+1)))./std(Pees1(1:end,v+1));
    B=Forces(:,v+1);%(Forces(1:end,v+1)-mean(Forces(1:end,v+1)))./std(Forces(1:end,v+1));
%     disp(A)
%     disp(B)
    pplt=corrcoef(A,B);
    pplt=pplt(1,2);
    corrs=[corrs,pplt];
%     plot(v,pplt,'o')
%     hold on
    
%     xlabel('Velocity [$\Gamma/k$]','Interpreter','latex')
%     ylabel('$\overline{\rho_{ee}}$','Interpreter','latex')
end
plot(0:1:59,corrs)
ylabel('Correlation')
xlabel('Velocity [\Gamma/k]')
ylim([-1.05 1.05])
    
% legend('1','6','11','16','21','26','31')
% legend('\chi=30, \Omega=104\Gamma','\chi=28, \Omega=108\Gamma','\chi=26, \Omega=113\Gamma','\chi=24.5, \Omega=118\Gamma','\chi=23.5, \Omega=121\Gamma','\chi=22.5, \Omega=125\Gamma','\chi=21, \Omega=133\Gamma')
% legend('\chi=28, \Omega=108\Gamma','\chi=24.5, \Omega=118\Gamma','\chi=21, \Omega=133\Gamma')
%%
dx=2;
dy=1;
figure 
surf(90:2:150,35:1:55,real(FSt'),'EdgeColor','none','facecolor','interp')
colormap('jet')
view(2)
hold on
scatter3(Rabs,Ang,[200,200,200,200,200,200,200],'ko')
text(Rabs+dx,Ang+dy,[200,200,200,200,200,200,200],[{'1'},{'2'},{'3'},{'4'},{'5'},{'6'},{'7'}],'Fontsize',10,'color','k')
c=colorbar;
c.Label.String='\rho_{ee}';
xlabel('\Omega [\Gamma]')
ylabel('\chi')
title('v = 1 \Gamma/k')
% hold on
% surf(90:2:150,15:1:35,real(FSt'),'EdgeColor','none','facecolor','interp')




%%
R1=55;
R2=133;
R3=105;
Forces1=[];
Forces2=[];
Forces3=[];
ReOmegas1=[];
ReOmegas2=[];
ReOmegas3=[];
ImOmegas1=[];
ImOmegas2=[];
ImOmegas3=[];
Om=H(1,2);
Om=subs(Om,[d1,v,chi],[det1,velocity,angle]);



DH1=subs(DH,[W1],[R1]);
Om1=subs(Om,W1,R1);

eq.evolve(t_equil,t_end,IC,[R1]);
i0=find(eq.evTime(:)>0,1);
av_e=mean(eq.evolution(2,2,i0:end));
for i=i0:50:length(eq.evTime(:)) %
    gradH=double(subs(DH1,[t],[eq.evTime(i)]));
    force=-(real(trace(squeeze(eq.evolution(:,:,i))*gradH)));
    Forces1=[Forces1,force];
end

for i=i0:5:length(eq.evTime(:))
    Omg=single(subs(Om1,t,eq.evTime(i)));
    ReOmegas1=[ReOmegas1,real(Omg)];
    ImOmegas1=[ImOmegas1,imag(Omg)];
end

tot_force=double(trapz(eq.evTime(i0:50:end),Forces1)/(eq.evTime(end)));
       
tot_force=(tot_force*2)

figure
plot(Forces1,'.')
xlabel('Time [1/\Gamma]')
ylabel('Force [2F_{rad}]')
title(['\Omega = ',num2str(R1),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])
eq.plotEvolution()


% ReOmegas1=ReOmegas1./max(ReOmegas1);
% ImOmegas1=ImOmegas1./max(ImOmegas1);
u1=squeeze(eq.evolution(1,2,i0:5:end)+eq.evolution(2,1,i0:5:end));
v1=1i*(squeeze(eq.evolution(1,2,i0:5:end)-eq.evolution(2,1,i0:5:end)));
w1=squeeze(eq.evolution(2,2,i0:5:end)-eq.evolution(1,1,i0:5:end));
k1=u1.^2+v1.^2+w1.^2;
% t1=squeeze(eq.evolution(2,2,i0:5:end)+eq.evolution(1,1,i0:5:end));
uo1=real(u1.*ReOmegas1');
vo1=real(v1.*ImOmegas1');
ko1=uo1.^2+vo1.^2+R1*w1.^2;


DH2=subs(DH,[W1],[R2]);
Om2=subs(Om,W1,R2);

eq.evolve(t_equil,t_end,IC,[R2]);
i0=find(eq.evTime(:)>0,1);
av_e=mean(eq.evolution(2,2,i0:end));
for i=i0:50:length(eq.evTime(:)) %
    gradH=double(subs(DH2,[t],[eq.evTime(i)]));
    force=-(real(trace(squeeze(eq.evolution(:,:,i))*gradH)));
    Forces2=[Forces2,force];
end

for i=i0:5:length(eq.evTime(:))
    Omg=single(subs(Om2,t,eq.evTime(i)));
    ReOmegas2=[ReOmegas2,real(Omg)];
    ImOmegas2=[ImOmegas2,imag(Omg)];
end
tot_force=double(trapz(eq.evTime(i0:50:end),Forces2)/(eq.evTime(end)));
       
tot_force=(tot_force*2)


figure
plot(Forces2,'.')
xlabel('Time [1/\Gamma]')
ylabel('Force [2F_{rad}]')
title(['\Omega = ',num2str(R2),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])
eq.plotEvolution()

% ReOmegas2=ReOmegas2./max(ReOmegas2);
% ImOmegas2=ImOmegas2./max(ImOmegas2);
u2=squeeze(eq.evolution(1,2,i0:5:end)+eq.evolution(2,1,i0:5:end));
v2=1i*(squeeze(eq.evolution(1,2,i0:5:end)-eq.evolution(2,1,i0:5:end)));
w2=squeeze(eq.evolution(2,2,i0:5:end)-eq.evolution(1,1,i0:5:end));
k2=u2.^2+v2.^2+w2.^2;
% t2=squeeze(eq.evolution(2,2,i0:5:end)+eq.evolution(1,1,i0:5:end));
uo2=real(u2.*ReOmegas2');
vo2=real(v2.*ImOmegas2');
ko2=uo2.^2+vo2.^2+R2.*w2.^2;



DH3=subs(DH,[W1],[R3]);
Om3=subs(Om,W1,R3);

eq.evolve(t_equil,t_end,IC,[R3]);
i0=find(eq.evTime(:)>0,1);
av_e=mean(eq.evolution(2,2,i0:end));
for i=i0:50:length(eq.evTime(:)) %
    gradH=double(subs(DH3,[t],[eq.evTime(i)]));
    force=-(real(trace(squeeze(eq.evolution(:,:,i))*gradH)));
    Forces3=[Forces3,force];
end

for i=i0:5:length(eq.evTime(:))
    Omg=single(subs(Om3,t,eq.evTime(i)));
    ReOmegas3=[ReOmegas3,real(Omg)];
    ImOmegas3=[ImOmegas3,imag(Omg)];
end

tot_force=double(trapz(eq.evTime(i0:50:end),Forces3)/(eq.evTime(end)));
       
tot_force=(tot_force*2)

figure
plot(Forces3,'.')
xlabel('Time [1/\Gamma]')
ylabel('Force [2F_{rad}]')
title(['\Omega = ',num2str(R3),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])
eq.plotEvolution() 


% ReOmegas3=ReOmegas3./max(ReOmegas3);
% ImOmegas3=ImOmegas3./max(ImOmegas3);
u3=squeeze(eq.evolution(1,2,i0:5:end)+eq.evolution(2,1,i0:5:end));
v3=1i*(squeeze(eq.evolution(1,2,i0:5:end)-eq.evolution(2,1,i0:5:end)));
w3=squeeze(eq.evolution(2,2,i0:5:end)-eq.evolution(1,1,i0:5:end));
k3=u3.^2+v3.^2+w3.^2;
% t3=squeeze(eq.evolution(2,2,i0:5:end)+eq.evolution(1,1,i0:5:end));

uo3=real(u3.*ReOmegas3');
vo3=real(v3.*ImOmegas3');
ko3=uo3.^2+vo3.^2+R3^2.*w3.^2;
disp('done')

%%
uo1=uo1./sqrt(ko1);
vo1=vo1./sqrt(ko1);
wo1=w1./sqrt(ko1);
uo2=uo2./sqrt(ko2);
vo2=vo2./sqrt(ko2);
wo2=w2./sqrt(ko2);
uo3=uo3./sqrt(ko3);
vo3=vo3./sqrt(ko3);
wo3=R3.*w3./sqrt(ko3);


u1=u1./sqrt(k1);
v1=v1./sqrt(k1);
w1=w1./sqrt(k1);
u2=u2./sqrt(k2);
v2=v2./sqrt(k2);
w2=w2./sqrt(k2);
u3=u3./sqrt(k3);
v3=v3./sqrt(k3);
w3=w3./sqrt(k3);

%%
figure
plot(ReOmegas1(20000:20500))
hold on
plot(ImOmegas1(20000:20500))

%%
figure
scatter3(u1(20000:21000),v1(20000:21000),w1(20000:21000),'.')
hold on
[x y z] = sphere(20);
h = surf(x, y, z); 
set(h, 'FaceAlpha', 0.01,'EdgeColor',[0.8 0.8 0.8])

xlabel('u')
ylabel('v')
zlabel('w')
title(['\Omega = ',num2str(R1),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])

figure
scatter3(u2(20000:21000),v2(20000:21000),w2(20000:21000),'.')
hold on
[x y z] = sphere(20);
h = surf(x, y, z); 
set(h, 'FaceAlpha', 0.01,'EdgeColor',[0.8 0.8 0.8])

xlabel('u')
ylabel('v')
zlabel('w')
title(['\Omega = ',num2str(R2),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])


figure
scatter3(u3(20000:21000),v3(20000:21000),w3(20000:21000),'.')
hold on
[x y z] = sphere(20);
h = surf(x, y, z); 
set(h, 'FaceAlpha', 0.01,'EdgeColor',[0.8 0.8 0.8])

xlabel('u')
ylabel('v')
zlabel('w')
title(['\Omega = ',num2str(R3),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])
%%
figure
scatter3(uo1(20000:21000),vo1(20000:21000),wo1(20000:21000),'.')
hold on
[x y z] = sphere(20);
h = surf(x, y, z); 
set(h, 'FaceAlpha', 0.01,'EdgeColor',[0.8 0.8 0.8])

xlabel('u')
ylabel('v')
zlabel('w')
title(['\Omega = ',num2str(R1),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])

figure
scatter3(uo2(20000:21000),vo2(20000:21000),wo2(20000:21000),'.')
hold on
[x y z] = sphere(20);
h = surf(x, y, z); 
set(h, 'FaceAlpha', 0.01,'EdgeColor',[0.8 0.8 0.8])

xlabel('u')
ylabel('v')
zlabel('w')
title(['\Omega = ',num2str(R2),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])


figure
scatter3(uo3(20000:21000),vo3(20000:21000),wo3(20000:21000),'.')
hold on
[x y z] = sphere(20);
h = surf(x, y, z); 
set(h, 'FaceAlpha', 0.01,'EdgeColor',[0.8 0.8 0.8])

xlabel('u')
ylabel('v')
zlabel('w')
title(['\Omega = ',num2str(R3),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])
%%
figure
scatter3(1:11001,uo3(10000:21000),vo3(10000:21000),'.')
% hold on
% [x y z] = sphere(20);
% h = surf(x, y, z); 
% set(h, 'FaceAlpha', 0.01,'EdgeColor',[0.8 0.8 0.8])

xlabel('u')
ylabel('v')
zlabel('w')
title(['\Omega = ',num2str(R2),' \Gamma, \delta = 100 \Gamma, \chi = 21 deg'])
%%


