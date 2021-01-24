%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wrong %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %% part b) Plot the expected maximum story drift (SDR) vs. ground motion intensity
% % Existing Building
% max_drift_EB=zeros(length(EB_NCL),NS);
% for i=1:length(EB_SDRs)
%     for j=1:size(EB_SDRs{i},1)
%         for k=1:NS
%             if max_drift_EB(j,k) < EB_SDRs{i}(j,k)
%                 max_drift_EB(j,k)=EB_SDRs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(EB_NCL,max_drift_EB(:,i),'linewidth',1.5)            
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('SDR_{max}')
% title('Existing Building')
% % Retrofitted building
% max_drift_RB=zeros(length(RB_NCL),NS);
% for i=1:length(RB_SDRs)
%     for j=1:size(RB_SDRs{i},1)
%         for k=1:NS
%             if max_drift_RB(j,k) < RB_SDRs{i}(j,k)
%                 max_drift_RB(j,k)=RB_SDRs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(RB_NCL,max_drift_RB(:,i),'linewidth',1.5)            
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('SDR_{max}')
% title('Retrofitted Building')
% 
% %% part c) plot expected peak floor acceleration (PFA) vs. ground motion intensity
% % Existing Building
% max_PFA_EB=zeros(length(EB_NCL),NS);
% for i=1:length(EB_PFAs)
%     for j=1:size(EB_PFAs{i},1)
%         for k=1:NS
%             if max_PFA_EB(j,k) < EB_PFAs{i}(j,k)
%                 max_PFA_EB(j,k)=EB_PFAs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(EB_NCL,max_PFA_EB(:,i),'linewidth',1.5)        %%%  EB_NCL == SaT1 ?????    
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('PFA')
% title('Existing Building')
% % Retrofitted building
% max_PFA_RB=zeros(length(RB_NCL),NS);
% for i=1:length(RB_PFAs)
%     for j=1:size(RB_PFAs{i},1)
%         for k=1:NS
%             if max_PFA_RB(j,k) < RB_PFAs{i}(j,k)
%                 max_PFA_RB(j,k)=RB_PFAs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(RB_NCL,max_PFA_RB(:,i),'linewidth',1.5)            
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('PFA')
% title('Retrofitted Building')
% 
% %% part d) probability of exceedance P(SDR_max > sdr_max) vs. maximum story drift
% sdr=linspace(0,0.02,1000);
% % Existing Building
% for i=1:length(EB_SDRs)
%     for j=1:size(EB_SDRs{i},1)
%             PFA_011(i,:)=EB_SDRs{i}(1,:);
%             if j==5
%                 PFA_055(i,:)=EB_SDRs{i}(5,:);
%             else if j==9
%                     PFA_099(i,:)=EB_SDRs{i}(9,:);
%                 end
%             end
%     end
% end
% for i=1:NS
%     % SaT1=0.11g
%     Mu_011(i)=mean(PFA_011(:,i));
%     Sigma_011(i)=std(PFA_011(:,i));
%     sigma_ln011(i)=sqrt(log((Sigma_011(i)/Mu_011(i))^2+1));
%     u_ln011(i)=log(Mu_011(i))-1/2*sigma_ln011(i)^2;
%     P_011(i,:)=1-logncdf(sdr,u_ln011(i),sigma_ln011(i));
%     % SaT1=0.55g
%     Mu_055(i)=mean(PFA_055(:,i));
%     Sigma_055(i)=std(PFA_055(:,i));
%     sigma_ln055(i)=sqrt(log((Sigma_055(i)/Mu_055(i))^2+1));
%     u_ln055(i)=log(Mu_055(i))-1/2*sigma_ln055(i)^2;
%     P_055(i,:)=1-logncdf(sdr,u_ln055(i),sigma_ln055(i));
%     % SaT1=0.99g
%     Mu_099(i)=mean(PFA_099(:,i));
%     Sigma_099(i)=std(PFA_099(:,i));
%     sigma_ln099(i)=sqrt(log((Sigma_099(i)/Mu_099(i))^2+1));
%     u_ln099(i)=log(Mu_099(i))-1/2*sigma_ln099(i)^2;
%     P_099(i,:)=1-logncdf(sdr,u_ln099(i),sigma_ln099(i));
% end
% 
% for i=1:NS
%     figure
%     plot(sdr,P_011(i,:),sdr,P_055(i,:),sdr,P_099(i,:))
%     box on; grid on;
%     xlabel('sdr_{max}')
%     ylabel('Probability of exceedance')
%     title(['Existing Building ',num2str(i),'^{st} Floor'])
%     legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
%     if i==2
%         xlim([0 0.004])
%     else if i==3
%             xlim([0 0.002])
%         end
%     end
% end
% % Retrofitted building
% for i=1:length(RB_SDRs)
%     for j=1:size(RB_SDRs{i},1)
%             PFA_011(i,:)=RB_SDRs{i}(1,:);
%             if j==5
%                 PFA_055(i,:)=RB_SDRs{i}(5,:);
%             else if j==9
%                     PFA_099(i,:)=RB_SDRs{i}(9,:);
%                 end
%             end
%     end
% end
% for i=1:NS
%     % SaT1=0.11g
%     Mu_011(i)=mean(PFA_011(:,i));
%     Sigma_011(i)=std(PFA_011(:,i));
%     sigma_ln011(i)=sqrt(log((Sigma_011(i)/Mu_011(i))^2+1));
%     u_ln011(i)=log(Mu_011(i))-1/2*sigma_ln011(i)^2;
%     P_011(i,:)=1-logncdf(sdr,u_ln011(i),sigma_ln011(i));
%     % SaT1=0.55g
%     Mu_055(i)=mean(PFA_055(:,i));
%     Sigma_055(i)=std(PFA_055(:,i));
%     sigma_ln055(i)=sqrt(log((Sigma_055(i)/Mu_055(i))^2+1));
%     u_ln055(i)=log(Mu_055(i))-1/2*sigma_ln055(i)^2;
%     P_055(i,:)=1-logncdf(sdr,u_ln055(i),sigma_ln055(i));
%     % SaT1=0.99g
%     Mu_099(i)=mean(PFA_099(:,i));
%     Sigma_099(i)=std(PFA_099(:,i));
%     sigma_ln099(i)=sqrt(log((Sigma_099(i)/Mu_099(i))^2+1));
%     u_ln099(i)=log(Mu_099(i))-1/2*sigma_ln099(i)^2;
%     P_099(i,:)=1-logncdf(sdr,u_ln099(i),sigma_ln099(i));
% end
% 
% for i=1:NS
%     figure
%     plot(sdr,P_011(i,:),sdr,P_055(i,:),sdr,P_099(i,:))
%     box on; grid on;
%     xlabel('sdr_{max}')
%     ylabel('Probability of exceedance')
%     title(['Retrofitted building ',num2str(i),'^{st} Floor'])
%     legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
%     if i==2
%         xlim([0 0.004])
%     else if i==3
%             xlim([0 0.002])
%         end
%     end
% end
% 
% %% part e) probability of exceedance P(PFA > pfa) vs. peak floor acceleration
% pfa=linspace(0,1,1000);
% % Existing Building
% for i=1:length(EB_PFAs)
%     for j=1:size(EB_PFAs{i},1)
%             PFA_011(i,:)=EB_PFAs{i}(1,:);
%             if j==5
%                 PFA_055(i,:)=EB_PFAs{i}(5,:);
%             else if j==9
%                     PFA_099(i,:)=EB_PFAs{i}(9,:);
%                 end
%             end
%     end
% end
% for i=1:NS
%     % SaT1=0.11g
%     Mu_011(i)=mean(PFA_011(:,i));
%     Sigma_011(i)=std(PFA_011(:,i));
%     sigma_ln011(i)=sqrt(log((Sigma_011(i)/Mu_011(i))^2+1));
%     u_ln011(i)=log(Mu_011(i))-1/2*sigma_ln011(i)^2;
%     P_011(i,:)=1-logncdf(pfa,u_ln011(i),sigma_ln011(i));
%     % SaT1=0.55g
%     Mu_055(i)=mean(PFA_055(:,i));
%     Sigma_055(i)=std(PFA_055(:,i));
%     sigma_ln055(i)=sqrt(log((Sigma_055(i)/Mu_055(i))^2+1));
%     u_ln055(i)=log(Mu_055(i))-1/2*sigma_ln055(i)^2;
%     P_055(i,:)=1-logncdf(pfa,u_ln055(i),sigma_ln055(i));
%     % SaT1=0.99g
%     Mu_099(i)=mean(PFA_099(:,i));
%     Sigma_099(i)=std(PFA_099(:,i));
%     sigma_ln099(i)=sqrt(log((Sigma_099(i)/Mu_099(i))^2+1));
%     u_ln099(i)=log(Mu_099(i))-1/2*sigma_ln099(i)^2;
%     P_099(i,:)=1-logncdf(pfa,u_ln099(i),sigma_ln099(i));
% end
% 
% for i=1:NS
%     figure
%     plot(pfa,P_011(i,:),pfa,P_055(i,:),pfa,P_099(i,:))
%     box on; grid on;
%     xlabel('pfa')
%     ylabel('Probability of exceedance')
%     title(['Existing Building ',num2str(i),'^{st} Floor'])
%     legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
% end
% 
% % Retrofitted building
% for i=1:length(RB_PFAs)
%     for j=1:size(RB_PFAs{i},1)
%             PFA_011(i,:)=RB_PFAs{i}(1,:);
%             if j==5
%                 PFA_055(i,:)=RB_PFAs{i}(5,:);
%             else if j==9
%                     PFA_099(i,:)=RB_PFAs{i}(9,:);
%                 end
%             end
%     end
% end
% for i=1:NS
%     % SaT1=0.11g
%     Mu_011(i)=mean(PFA_011(:,i));
%     Sigma_011(i)=std(PFA_011(:,i));
%     sigma_ln011(i)=sqrt(log((Sigma_011(i)/Mu_011(i))^2+1));
%     u_ln011(i)=log(Mu_011(i))-1/2*sigma_ln011(i)^2;
%     P_011(i,:)=1-logncdf(pfa,u_ln011(i),sigma_ln011(i));
%     % SaT1=0.55g
%     Mu_055(i)=mean(PFA_055(:,i));
%     Sigma_055(i)=std(PFA_055(:,i));
%     sigma_ln055(i)=sqrt(log((Sigma_055(i)/Mu_055(i))^2+1));
%     u_ln055(i)=log(Mu_055(i))-1/2*sigma_ln055(i)^2;
%     P_055(i,:)=1-logncdf(pfa,u_ln055(i),sigma_ln055(i));
%     % SaT1=0.99g
%     Mu_099(i)=mean(PFA_099(:,i));
%     Sigma_099(i)=std(PFA_099(:,i));
%     sigma_ln099(i)=sqrt(log((Sigma_099(i)/Mu_099(i))^2+1));
%     u_ln099(i)=log(Mu_099(i))-1/2*sigma_ln099(i)^2;
%     P_099(i,:)=1-logncdf(pfa,u_ln099(i),sigma_ln099(i));
% end
% 
% for i=1:NS
%     figure
%     plot(pfa,P_011(i,:),pfa,P_055(i,:),pfa,P_099(i,:))
%     box on; grid on;
%     xlabel('pfa')
%     ylabel('Probability of exceedance')
%     title(['Retrofitted Building ',num2str(i),'^{st} Floor'])
%     legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
% end



















% %% part b) Plot the expected maximum story drift (SDR) vs. ground motion intensity
% % Existing Building
% max_drift_EB=zeros(length(EB_NCL),NS);
% for i=1:length(EB_NCL)
%     for j=1:size(EB_SDRs{i},1)
%         for k=1:NS
%             if max_drift_EB(i,k) < EB_SDRs{i}(j,k)
%                 max_drift_EB(i,k)=EB_SDRs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(EB_NCL,max_drift_EB(:,i),'linewidth',1.5)            
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('SDR_{max}')
% title('Existing Building')
% % Retrofitted building
% max_drift_RB=zeros(length(RB_NCL),NS);
% for i=1:length(RB_NCL)
%     for j=1:size(RB_SDRs{i},1)
%         for k=1:NS
%             if max_drift_RB(i,k) < RB_SDRs{i}(j,k)
%                 max_drift_RB(i,k)=RB_SDRs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(RB_NCL,max_drift_RB(:,i),'linewidth',1.5)            
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('SDR_{max}')
% title('Retrofitted Building')
% 
% % mean(max_drift_EB(:,2))
% % mean(max_drift_RB(:,2))
% 
% %% part c) plot expected peak floor acceleration (PFA) vs. ground motion intensity
% % Existing Building
% max_PFA_EB=zeros(length(EB_NCL),NS);
% for i=1:length(EB_NCL)
%     for j=1:size(EB_PFAs{i},1)
%         for k=1:NS
%             if max_PFA_EB(i,k) < EB_PFAs{i}(j,k)
%                 max_PFA_EB(i,k)=EB_PFAs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(EB_NCL,max_PFA_EB(:,i),'linewidth',1.5)        %%%  EB_NCL == SaT1 ?????    
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('PFA')
% title('Existing Building')
% % Retrofitted building
% max_PFA_RB=zeros(length(RB_NCL),NS);
% for i=1:length(RB_NCL)
%     for j=1:size(RB_PFAs{i},1)
%         for k=1:NS
%             if max_PFA_RB(i,k) < RB_PFAs{i}(j,k)
%                 max_PFA_RB(i,k)=RB_PFAs{i}(j,k);
%             end
%         end
%     end  
% end
% figure
% for i=1:NS
%     plot(RB_NCL,max_PFA_RB(:,i),'linewidth',1.5)            
%     hold on;
% end
% box on; grid on;
% legend('1','2','3') %%%% ??????
% xlabel('SaT1')
% ylabel('PFA')
% title('Retrofitted Building')
% 
% part d) probability of exceedance P(SDR_max > sdr_max) vs. maximum story drift

