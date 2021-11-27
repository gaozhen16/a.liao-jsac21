clear, clc

load('Z1_Save_Delay_CRB_v0_BD_V1_1');
RMSE_Delay_set_00 = sqrt(MSE_Delay_set_0(:,1));
RMSE_Delay_set_01 = sqrt(MSE_Delay_set_1(:,1));
RCRB_Delay_set_0 = sqrt(CRB_Delay_set(:,1));
load('Z1_Save_Delay_CRB_v0_BD_V2_1');
RMSE_Delay_set_10 = sqrt(MSE_Delay_set_0(:,2));
RMSE_Delay_set_11 = sqrt(MSE_Delay_set_1(:,2));
RCRB_Delay_set_1 = sqrt(CRB_Delay_set(:,2));

MarkerSize = 8;
LineWidth = 2;
Fontsize = 14;

%%
figure
semilogy(SNR_dBs,RMSE_Delay_set_01,'--b^','LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
semilogy(SNR_dBs,RMSE_Delay_set_00,'--ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RCRB_Delay_set_0,'--k','LineWidth',LineWidth+0.2);
semilogy(SNR_dBs,RMSE_Delay_set_11,'-b^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Delay_set_10,'-ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RCRB_Delay_set_1,'-k','LineWidth',LineWidth+0.2);
semilogy(SNR_dBs,RMSE_Delay_set_00,'--ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Delay_set_01,'--b^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Delay_set_10,'-ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Delay_set_11,'-b^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
axis([-120 20 3e-8 1e3])
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel('SNR [dB]','Fontsize',18); ylabel('RMSE','Fontsize',18)
set(gcf, 'position', [600 250 700 500]);
set(gca,'position',[0.098 0.12 0.883 0.87])
axis normal;
h1 = legend('Proposed scheme, triple squint','Conventional scheme, no triple squint','CRLB, SNR = -20 dB',...
    'Proposed scheme, triple squint','Conventional scheme, no triple squint','CRLB, SNR = 20 dB','Location','northeast');
set(h1,'Fontsize',Fontsize); % legend('boxoff')
text(-106.3848, 0.000316536,'SNR = -20 dB','color','k','FontSize',15,'Fontname','Times New Roman')
annotation('arrow',[0.53 0.357],[0.536 0.452],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
annotation('ellipse',[0.525 0.524 0.03043 0.1],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
text(-63.701, 6.0968e-07,'SNR = 20 dB','color','k','FontSize',15,'Fontname','Times New Roman')
annotation('arrow',[0.72 0.616],[0.284 0.226],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
annotation('ellipse',[0.715 0.276 0.03043 0.082],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
