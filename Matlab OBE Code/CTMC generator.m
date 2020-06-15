% clear;
IF=load('BaH_Intensity_force_2CF.mat');
IF=IF.Intensity_force;


Int_profile=zeros(33);
Force_profile=zeros(33);
i=1;

for x=-0.4:0.025:0.4
    j=1;
    for y=-0.4:0.025:0.4
        Int_profile(i,j)=intensity(x,y);
        Force_profile(i,j)=interp1(IF(1,:),IF(2,:),intensity(x,y),'cubic');
        j=j+1;
    end
    i=i+1;
end

dlmwrite('Force_profile.txt',Force_profile)

%%
figure 
surf(-0.4:0.025:0.4,-0.4:0.025:0.4,Int_profile,'EdgeColor','none','facecolor','interp')
colormap('jet')
view(2)
hold on
xlim([-0.4 0.4])
ylim([-0.4 0.4])
%%
figure 
surf(-0.4:0.025:0.4,-0.4:0.025:0.4,Force_profile,'EdgeColor','none','facecolor','interp')
colormap('jet')
view(2)
hold on
xlim([-0.4 0.4])
ylim([-0.4 0.4])
% %%
% figure
% plot(IF(1,:),IF(2,:))




function w=intensity(x,y)

I0=1.1*45;
s=1.13/2;

w=I0*exp(-(x^2+y^2)/(2*s^2));
end