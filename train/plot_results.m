clear; close all;

%%%%%%%%%load data%%%%%%%%%
test = 1; % Default to test

load loss.out; load energy_train.out;load virial_train.out;load force_train.out;
energy=energy_train;
virial=virial_train;
force=force_train;
loss = loss(:,3:7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

for i = 1: size(virial)
    if(virial(i,2)==0)
        virial(i,1)=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1)
loglog(loss(:,:),'-','linewidth',1); hold on;
title('(a)','fontsize',15)
xlabel('Generation','fontsize',18,'interpreter','latex');
ylabel('Loss functions','fontsize',18,'interpreter','latex');
set(gca,'FontName','Times New Roman','fontsize',15,'ticklength',get(gca,'ticklength')*2);
legend('L1','L2','Energy','Force','Virial','Location','SouthWest');
axis tight

subplot(2,2,2)
plot(energy(:,2),energy(:,1),'.','markersize',20,'color',[0,0.45,0.74]); hold on;
x=get(gca,'xlim');
plot(x,x,'linewidth',2,'color', [0.85,0.33,0.1]);
title('(b)','fontsize',15)
xlabel('DFT energy (eV/atom)','fontsize',18,'interpreter','latex');
ylabel('NEP energy (eV/atom)','fontsize',18,'interpreter','latex');
set(gca,'FontName','Times New Roman','fontsize',15,'ticklength',get(gca,'ticklength')*2)
axis([-8.8 -4.2 -8.8 -4.2 ])
%axis tight

subplot(2,2,3)
plot(force(:,4:6),force(:,1:3),'c.','markersize',20,'color' ,[0,0.45,0.74]); hold on;
x=get(gca,'xlim');
plot(x,x,'linewidth',2,'color', [0.85,0.33,0.1]);
title('(c)','fontsize',15)
xlabel('DFT force (eV/\AA)','fontsize',18,'interpreter','latex');
ylabel('NEP force (eV/\AA)','fontsize',18,'interpreter','latex');
set(gca,'FontName','Times New Roman','fontsize',15,'ticklength',get(gca,'ticklength')*2);
axis([-30 30 -30 30])
%axis tight;

subplot(2,2,4)
plot(virial(:,2),virial(:,1),'.','markersize',20 ,'color' ,[0,0.45,0.74]); hold on;
x=get(gca,'xlim')*1.2;
plot(x,x,'linewidth',2,'color', [0.85,0.33,0.1]);
title('(d)','fontsize',15)
xlabel('DFT virial (eV/atom)','fontsize',18,'interpreter','latex');
ylabel('NEP virial (eV/atom)','fontsize',18,'interpreter','latex');
set(gca,'FontName','Times New Roman','fontsize',15,'ticklength',get(gca,'ticklength')*2);
axis([-6 11 -6 11])

set(gcf,'Position',[100 100 800 800]);
