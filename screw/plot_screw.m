nep = [0. 0.00834422 0.03092597 0.0559019  0.0670981  0.06859237 0.06711217 0.05607725 0.03117593 0.00850427 0.00017849];
MV4S = [0,0,0,0.00272027,0.018451292,0.042325469,0.055881881,0.059555449,0.055453378,0.035437329,0.008974482,0,0,0,0];
MV2B = [0,0,0.0,0.001211684,0.014145402,0.033933558,0.047929706,0.051603274,0.045556091,0.027274916,0.007043813,0,0,0,0];
DFT_csanyi = [0,0.048115872,0.10106564,0.071874498,0.03106564];
DFT_Rodney = [0,0.007497994,0.030386776,0.057786872,0.079858771,0.087961804,0.078850907,0.056397047,0.028762638,0.00496068,0];
DFT_stoller = [0,0.005376344,0.023318889,0.046451613,0.068284385,0.078709677,0.067839833,0.045806452,0.021808698,0.003020382,0];

plot(linspace(0, 1, length(nep)),nep,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(MV4S)),MV4S,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(MV2B)),MV2B,'o-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(DFT_stoller)),DFT_stoller,'s-' ,'linewidth',2,MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(DFT_Rodney)),DFT_Rodney,'s-' ,'linewidth',2, MarkerFaceColor='auto'),hold on
plot(linspace(0, 1, length(DFT_csanyi)),DFT_csanyi,'s-' ,'linewidth',2, MarkerFaceColor='auto'),hold on

xlabel('Reaction coordinate','fontsize',16,'interpreter','latex');
ylabel('Energy per unit length (eV/b)','fontsize',16,'interpreter','latex');
h = legend('NEP','MV4S','MV2B','DFT\_stoller','DFT\_Rodney','DFT\_csanyi');
ylim([0 0.13])
set(h,'fontsize',12,'interpreter','latex')