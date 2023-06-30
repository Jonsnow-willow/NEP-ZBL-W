figure
nep = [0, 0.015918666757837423, 0.05121355553937171, 0.08368180022410947, 0.1002255177080304, 0.10438729596588193, 0.10023050338010915, 0.08369681635953567, 0.05123516494484517, 0.015935983111235386, -5.931235975112001e-07];
dft = [0, 0.016, 0.051, 0.082, 0.101, 0.112, 0.100, 0.080, 0.049, 0.014, 0];
MV2B = [0, 0.015, 0.044, 0.070, 0.086, 0.093, 0.086, 0.069, 0.043, 0.013, 0];
MV4S = [0, 0.015, 0.044, 0.072, 0.088, 0.095, 0.088, 0.071, 0.043, 0.013 0];
GAP = [0.0, 0.017, 0.060 0.10 0.12 0.13 0.12 0.10 0.060 0.017 0.00042];


plot(linspace(0, 1, length(nep)),nep,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(dft)),dft,'s-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(MV2B)),MV2B,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(MV4S)),MV4S,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
xlabel('Reaction coordinate','fontsize',16,'interpreter','latex');
ylabel('Stacking fault Energy (eV/\AA$^2$)','fontsize',16,'interpreter','latex');

h = legend('NEP','DFT','MV2B', 'MV4S');
ylim([0 0.13])
set(h,'fontsize',13,'interpreter','latex')
