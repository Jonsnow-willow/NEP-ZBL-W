figure
nep = [0, 0.019068558478153547, 0.05907028731268223, 0.09350702893908938, 0.11113330991944471, 0.11623021519991258, 0.1111391463490679, 0.09352263786248778, 0.05909397539638859, 0.019088964597844482, 0];
dft = [0, 0.02, 0.065, 0.102, 0.117, 0.113, 0.099, 0.084, 0.057, 0.021, 0.0];
MV2B = [0, 0.018, 0.057, 0.089, 0.106, 0.108, 0.105, 0.088, 0.058, 0.019, 0];
MV4S = [0, 0.018, 0.058, 0.089, 0.105, 0.109, 0.104, 0.088, 0.058, 0.020 0];

plot(linspace(0, 1, length(nep)),nep,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(dft)),dft,'s-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(MV2B)),MV2B,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(MV4S)),MV4S,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
xlabel('Reaction coordinate','fontsize',16,'interpreter','latex');
ylabel('Stacking fault Energy (eV/\AA$^2$)','fontsize',16,'interpreter','latex');

h = legend('NEP','DFT','MV2B', 'MV4S');
ylim([0 0.13])
set(h,'fontsize',13,'interpreter','latex')
