clear, clc

load('Z1_Save_ASE_v0_V2');

MarkerSize = 8;
LineWidth = 2;
Fontsize = 14;

%%
figure
plot(SNR_dBs,ASE_Per_CSI_No_TrSq,'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
plot(SNR_dBs,ASE_Per_CSI_TrSq,'-r^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(SNR_dBs,ASE_Est_CSI_No_TrSq,'-mp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(SNR_dBs,ASE_Est_CSI_TrSq,'-ko','LineWidth',LineWidth,'MarkerSize',MarkerSize);
axis([-35 -5 0 55])
set(gca,'xtick',(-35:3:-5));
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel('SNR [dB]','Fontsize',18); ylabel('Average Spectral Efficiency [bit/s/Hz]','Fontsize',18)
set(gcf, 'position', [600 250 700 500]);
set(gca,'position',[0.082 0.12 0.90 0.86])
axis normal;
h1 = legend('Perfect CSI, no triple squint','Perfect CSI, triple squint',...
    'Estimated CSI, no triple squint','Estimated CSI, triple squint','Location','southeast');
set(h1,'Fontsize',Fontsize);
h_small_11 = axes('position',[0.607714 0.648333 0.245 0.175]);
axis(h_small_11);
plot([-8,-5],ASE_Per_CSI_No_TrSq(end-1:end),'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
plot([-8,-5],ASE_Est_CSI_No_TrSq(end-1:end),'-mp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
axis([-5.2 -4.9 53.78 53.9])
set(gca,'xtick',[-5.2 -5.1 -5 -4.9]);
set(gca,'xticklabels',{'-5.2','-5.1','-5','-4.9'}, 'linewidth',1.0);
annotation('line',[0.98143 0.77],[0.96 0.824],'LineWidth',1.0,'LineStyle','-');
h_small_12 = axes('position',[0.720238 0.41466667 0.245 0.175]);
axis(h_small_12);
plot([-8,-5],ASE_Per_CSI_TrSq(end-1:end),'-r^','LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
plot([-8,-5],ASE_Est_CSI_TrSq(end-1:end),'-ko','LineWidth',LineWidth,'MarkerSize',MarkerSize);
axis([-5.2 -4.9 51.18 51.3])
set(gca,'xtick',[-5.2 -5.1 -5 -4.9]);
set(gca,'xticklabels',{'-5.2','-5.1','-5','-4.9'}, 'linewidth',1.0);
annotation('line',[0.98143 0.8843],[0.92 0.588],'LineWidth',1.0,'LineStyle','-');
