% CEE244_Final_Project_Part b to part o ;Weiyang_Bai

clear all; close all; clc
load ExistingBuildingPSDAData
load RetrofittedBuildingPSDAData

% EB = ExistingBuilding; RB = RetrofittedBuilding
% Story Drift Ratio
EB_SDRs=ExistingBuildingPSDAData.SDRs;
RB_SDRs=RetrofittedBuildingPSDAData.SDRs;

% Peak floor acceleration
EB_PFAs=ExistingBuildingPSDAData.PFAs;
RB_PFAs=RetrofittedBuildingPSDAData.PFAs;

% incremental dynamic analyses (IDAs) (Ground Motion)
EB_IDASas=ExistingBuildingPSDAData.IDASas;
RB_IDASas=RetrofittedBuildingPSDAData.IDASas;

% Collapse intensity associated
EB_CSas=ExistingBuildingPSDAData.CollapseSas;
RB_CSas=RetrofittedBuildingPSDAData.CollapseSas;

% Number of Sas (intensity levels)
EB_NSas=ExistingBuildingPSDAData.NumberOfSasRun;
RB_NSas=RetrofittedBuildingPSDAData.NumberOfSasRun;

% intensity levels to be used for non-collapse losses
EB_NCL=ExistingBuildingPSDAData.imsForNonCollapseLosses;
RB_NCL=RetrofittedBuildingPSDAData.imsForNonCollapseLosses;

%% extract data 
NS=3; % number of stories
%% Existing Building
% SDR data (story drift ratio)
EB_SDR_STORE={{}}; % the first layer of 3 cell is story 1, 2 & 3; Second layer is ground motion starting from 0.11g 
for i=1:NS
    for j=1:length(EB_NCL) 
            temp=[];
        for k=1:length(EB_SDRs)
            if j<=EB_NSas(k)
                temp=[temp,EB_SDRs{k}(j,i)];
            end
        end
        EB_SDR_STORE{i}{j}=temp;
    end
end
% Mean value of each ground motion
for i=1:NS
    for j=1:length(EB_NCL)
    mean_EB_SDR{i}(j)=mean(EB_SDR_STORE{i}{j});
    end
end
% Stander diviation of each ground motion
for i=1:NS
    for j=1:length(EB_NCL)
    std_EB_SDR{i}(j)=std(EB_SDR_STORE{i}{j});
    end
end

% PFA data (peak floor acceleration)
EB_PFA_STORE={{}}; % the first layer of 3 cell is story 1, 2 & 3; Second layer is ground motion starting from 0.11g 
for i=1:NS
    for j=1:length(EB_NCL) 
            temp=[];
        for k=1:length(EB_PFAs)
            if j<=EB_NSas(k)
                temp=[temp,EB_PFAs{k}(j,i)];
            end               
        end
        EB_PFA_STORE{i}{j}=temp;
    end
end
% Mean value of each ground motion
for i=1:NS
    for j=1:length(EB_NCL)
    mean_EB_PFA{i}(j)=mean(EB_PFA_STORE{i}{j});
    end
end
% Stander diviation of each ground motion
for i=1:NS
    for j=1:length(EB_NCL)
    std_EB_PFA{i}(j)=std(EB_PFA_STORE{i}{j});
    end
end


%% Retrofitted Building
RB_SDR_STORE={{}}; % same as EB_SDR_STORE={{}}
for i=1:NS
    for j=1:length(RB_NCL) 
            temp=[];
        for k=1:length(RB_SDRs)
            if j<=RB_NSas(k)
                temp=[temp,RB_SDRs{k}(j,i)];
            end               
        end
        RB_SDR_STORE{i}{j}=temp;
    end
end
% Mean value of each ground motion
for i=1:NS
    for j=1:length(RB_NCL)
    mean_RB_SDR{i}(j)=mean(RB_SDR_STORE{i}{j});
    end
end
% Stander diviation of each ground motion
for i=1:NS
    for j=1:length(RB_NCL)
    std_RB_SDR{i}(j)=std(RB_SDR_STORE{i}{j});
    end
end
% PFA data (peak floor acceleration)
RB_PFA_STORE={{}}; % the first layer of 3 cell is story 1, 2 & 3; Second layer is ground motion starting from 0.11g 
for i=1:NS
    for j=1:length(RB_NCL) 
            temp=[];
        for k=1:length(RB_PFAs)
            if j<=RB_NSas(k)
                temp=[temp,RB_PFAs{k}(j,i)];
            end               
        end
        RB_PFA_STORE{i}{j}=temp;
    end
end
% Mean value of each ground motion
for i=1:NS
    for j=1:length(RB_NCL)
    mean_RB_PFA{i}(j)=mean(RB_PFA_STORE{i}{j});
    end
end
% Stander diviation of each ground motion
for i=1:NS
    for j=1:length(RB_NCL)
    std_RB_PFA{i}(j)=std(RB_PFA_STORE{i}{j});
    end
end

%% part b)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the expected maximum story drift (SDR) vs. ground motion intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Existing Building
figure
for i=1:NS
    plot(EB_NCL,mean_EB_SDR{i},'linewidth',1.5)            
    hold on;
end
box on; grid on;
legend('1^{st}','2^{nd}','3^{rd}','location','northwest') 
xlabel('SaT1')
ylabel('SDR_{max}')
title('SDR of Existing Building')
cd figure
saveas (gcf,'P_b1.jpg')
cd ..

% Retrofitted Building
figure
for i=1:NS
    plot(RB_NCL,mean_RB_SDR{i},'linewidth',1.5)            
    hold on;
end
box on; grid on;
legend('1^{st}','2^{nd}','3^{rd}','location','northwest') 
xlabel('SaT1')
ylabel('SDR_{max}')
title('SDR of Retrofitted Building')
cd figure
saveas (gcf,'P_b2.jpg')
cd ..

%% part c)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot expected peak floor acceleration (PFA) vs. ground motion intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Existing Building
figure
% subplot (2,1,1)
for i=1:NS
    plot(EB_NCL,mean_EB_PFA{i},'linewidth',1.5)            
    hold on;
end
box on; grid on;
legend('1^{st}','2^{nd}','3^{rd}','location','northwest') 
xlabel('SaT1')
ylabel('PFA')
title('PFA of Existing Building')
cd figure
saveas (gcf,'P_c1.jpg')
cd ..


% Retrofitted Building
% subplot (2,1,2)
figure
for i=1:NS
    plot(RB_NCL,mean_RB_PFA{i},'linewidth',1.5)            
    hold on;
end
box on; grid on;
legend('1^{st}','2^{nd}','3^{rd}','location','northwest') 
xlabel('SaT1')
ylabel('PFA')
title('PFA of Retrofitted Building')
cd figure
saveas (gcf,'P_c2.jpg')
cd ..

%% part d)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probability of exceedance P(SDR_max > sdr_max) vs. maximum story drift %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Existing Building
sdr=linspace(0,0.2,1000);
for i=1:NS
    for j=1:length(EB_NCL)
        sigma_ln_EB_SDR{i}(j)=sqrt(log((std_EB_SDR{i}(j)/mean_EB_SDR{i}(j))^2+1));
        u_ln_EB_SDR{i}(j)=log(mean_EB_SDR{i}(j))-1/2*sigma_ln_EB_SDR{i}(j)^2;
        P_EB_SDR{i}{j}=1-logncdf(sdr,u_ln_EB_SDR{i}(j),sigma_ln_EB_SDR{i}(j));
    end
end
figure
for i=1:NS   
    subplot(3,1,i)
    plot(sdr,P_EB_SDR{i}{1},sdr,P_EB_SDR{i}{5},sdr,P_EB_SDR{i}{9})
    box on; grid on;
    xlabel('sdr_{max}')
    ylabel('Probability of exceedance')
    title(['Probability of Exceedance SDR-EB ',num2str(i),'^{st} Floor'])
    legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
    xlim([0 0.02])
end        
cd figure
saveas(gcf,sprintf('P_d1.jpg',i))
cd ..
% Retrofitted Building
for i=1:NS
    for j=1:length(RB_NCL)
        sigma_ln_RB_SDR{i}(j)=sqrt(log((std_RB_SDR{i}(j)/mean_RB_SDR{i}(j))^2+1));
        u_ln_RB_SDR{i}(j)=log(mean_RB_SDR{i}(j))-1/2*sigma_ln_RB_SDR{i}(j)^2;
        P_RB_SDR{i}{j}=1-logncdf(sdr,u_ln_RB_SDR{i}(j),sigma_ln_RB_SDR{i}(j));
    end
end
figure
for i=1:NS
    subplot(3,1,i)
    plot(sdr,P_RB_SDR{i}{1},sdr,P_RB_SDR{i}{5},sdr,P_RB_SDR{i}{9})
    box on; grid on;
    xlabel('sdr_{max}')
    ylabel('Probability of exceedance')
    title(['Probability of Exceedance SDR-RB ',num2str(i),'^{st} Floor'])
    legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
    xlim([0 0.02])
end        
cd figure
saveas(gcf,sprintf('P_d2.jpg',i))
cd ..
%% part e)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probability of exceedance P(PFA > pfa) vs. peak floor acceleration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Existing Building
pfa=linspace(0,10,1000);
for i=1:NS
    for j=1:length(EB_NCL)
        sigma_ln_EB_PFA{i}(j)=sqrt(log((std_EB_PFA{i}(j)/mean_EB_PFA{i}(j))^2+1));
        u_ln_EB_PFA{i}(j)=log(mean_EB_PFA{i}(j))-1/2*sigma_ln_EB_PFA{i}(j)^2;
        P_EB_PFA{i}{j}=1-logncdf(pfa,u_ln_EB_PFA{i}(j),sigma_ln_EB_PFA{i}(j));
    end
end
figure
for i=1:NS
    subplot(3,1,i)
    plot(pfa,P_EB_PFA{i}{1},pfa,P_EB_PFA{i}{5},pfa,P_EB_PFA{i}{9})
    box on; grid on;
    xlabel('pfa')
    ylabel('Probability of exceedance')
    title(['Probability of Exceedance PFA-EB ',num2str(i),'^{st} Floor'])
    legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
    xlim ([0 1])
end      
cd figure
saveas(gcf,sprintf('P_e1.jpg',i))
cd ..

% Retrofitted Building
for i=1:NS
    for j=1:length(RB_NCL)
        sigma_ln_RB_PFA{i}(j)=sqrt(log((std_RB_PFA{i}(j)/mean_RB_PFA{i}(j))^2+1));
        u_ln_RB_PFA{i}(j)=log(mean_RB_PFA{i}(j))-1/2*sigma_ln_RB_PFA{i}(j)^2;
        P_RB_PFA{i}{j}=1-logncdf(pfa,u_ln_RB_PFA{i}(j),sigma_ln_RB_PFA{i}(j));
    end
end
figure
for i=1:NS
    subplot(3,1,i)
    plot(pfa,P_RB_PFA{i}{1},pfa,P_RB_PFA{i}{5},pfa,P_RB_PFA{i}{9})
    box on; grid on;
    xlabel('pfa')
    ylabel('Probability of exceedance')
    title(['Probability of Exceedance PFA-RB ',num2str(i),'^{st} Floor'])
    legend('Sa_{T1} = 0.11g','Sa_{T1} = 0.55g','Sa_{T1} = 0.99g')
    xlim ([0 1])
end 
cd figure
saveas(gcf,sprintf('P_e2.jpg',i))
cd ..
%% part f)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% collapse fragility function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load PSHA.mat
% Existing Building
% Empirical & lognormal function
sort_EB_CSas=sort(EB_CSas); % sort the SaT1 (ground motion)
for i=1:length(EB_CSas)
    Rank_EB(i)=i/length(EB_CSas);
end
mean_EB_CSas=mean(EB_CSas);
std_EB_CSas=std(EB_CSas);
sigma_ln_EB_CSas=sqrt(log((std_EB_CSas/mean_EB_CSas)^2+1));
u_ln_EB_CSas=log(mean_EB_CSas)-1/2*sigma_ln_EB_CSas^2;
P_EB_Collapse=logncdf(SaT1,u_ln_EB_CSas,sigma_ln_EB_CSas);
% a. Probability of collapse at the MCE
% building A & B 
P_EB_MCE_A=interp1(SaT1,P_EB_Collapse,MCE_A);
P_EB_MCE_B=interp1(SaT1,P_EB_Collapse,MCE_B);
figure
plot (sort_EB_CSas,Rank_EB,'ro',SaT1,P_EB_Collapse,'linewidth',1.5),hold on;
box on;grid on;
plot(MCE_A,P_EB_MCE_A,'g^','linewidth',1.5);
text(MCE_A+0.2,P_EB_MCE_A,['P(C|IM = Sa_{T1}(',num2str(MCE_A),')',' = ',num2str(P_EB_MCE_A)])
xlabel('Sa_{T1}(g)')
ylabel('Probability of Collapse')
legend('Empirical','Lognormal','Sa_{T1}(MCE)','location','northwest')
title('MCE-A of EB Building')
cd figure
saveas(gcf,'P_f1.jpg')
cd ..

figure
plot (sort_EB_CSas,Rank_EB,'ro',SaT1,P_EB_Collapse,'linewidth',1.5),hold on;
box on;grid on;
plot(MCE_B,P_EB_MCE_B,'g^','linewidth',1.5);
text(MCE_B+0.2,P_EB_MCE_B,['P(C|IM = Sa_{T1}(',num2str(MCE_B),')',' = ',num2str(P_EB_MCE_B)])
xlabel('Sa_{T1}(g)')
ylabel('Probability of Collapse')
legend('Empirical','Lognormal','Sa_{T1}(MCE)','location','northwest')
title('MCE-B of EB Building')
cd figure
saveas(gcf,'P_f2.jpg')
cd ..

% b. The median collapse intensity
P_EB_median=0.5;
SaT1_EB_median=interp1(P_EB_Collapse,SaT1,P_EB_median);
figure
plot (sort_EB_CSas,Rank_EB,'ro',SaT1,P_EB_Collapse,'linewidth',1.5),hold on;
box on;grid on;
plot(SaT1_EB_median,P_EB_median,'ks','linewidth',1.5);
text(SaT1_EB_median+0.2,P_EB_median,['Sa_{T1}(median) = ',num2str(SaT1_EB_median)])
xlabel('Sa_{T1}(g)')
ylabel('Probability of Collapse')
legend('Empirical','Lognormal','Sa_{T1}(median)','location','northwest')
title('Median Collapse Intensity of EB Building')
cd figure
saveas(gcf,'P_f3.jpg')
cd ..
% c. CMR the collaps margin ratio
CMR_EB_A=SaT1_EB_median/MCE_A;
CMR_EB_B=SaT1_EB_median/MCE_B;
% d. The mean annual frequency of collapse
for i =1:length(P_PSHA_A)
    if i==1
        d_lambda_A(i)=P_PSHA_A(i)/dIM;
        d_lambda_B(i)=P_PSHA_B(i)/dIM;
    else if i==length(P_PSHA_A)
        d_lambda_A(i)=(P_PSHA_A(i)-P_PSHA_A(i-1))/dIM;
        d_lambda_B(i)=(P_PSHA_B(i)-P_PSHA_B(i-1))/dIM;
        else
        d_lambda_A(i)=(P_PSHA_A(i+1)-P_PSHA_A(i-1))/(dIM*2);
        d_lambda_B(i)=(P_PSHA_B(i+1)-P_PSHA_B(i-1))/(dIM*2);
        end
    end
end
for i=1:length(P_EB_Collapse)    
    lambdaC_EB_A(i)=P_EB_Collapse(i)*abs(d_lambda_A(i));   
    lambdaC_EB_B(i)=P_EB_Collapse(i)*abs(d_lambda_B(i));   
end
MAF_EB_A=sum(lambdaC_EB_A)*dIM; % Annual frequency of collapse for EB A
MAF_EB_B=sum(lambdaC_EB_B)*dIM; % Annual frequency of collapse for EB A
disp (['part f_d']);
disp (['Annual frequency of collapse for EB A is ',num2str(MAF_EB_A)]);
disp (['Annual frequency of collapse for EB B is ',num2str(MAF_EB_B)]);
figure
plot (SaT1,lambdaC_EB_A,'linewidth',1.5);box on;grid on;
xlabel('Sa_{T1}(g)')
ylabel('P(C|IM=x)x d\lambda/dIM')
title('EB-A Disaggregation of Collapse Risk')
cd figure
saveas(gcf,'P_f4.jpg')
cd ..

figure
plot (SaT1,lambdaC_EB_B,'linewidth',1.5);box on;grid on;
xlabel('Sa_{T1}(g)')
ylabel('P(C|IM=x)x d\lambda/dIM')
title('EB-B Disaggregation of Collapse Risk')
cd figure
saveas(gcf,'P_f5.jpg')
cd ..
% Retrofitted Building
sort_RB_CSas=sort(RB_CSas); % sort the SaT1 (ground motion)
for i=1:length(RB_CSas)
    Rank_RB(i)=i/length(RB_CSas);
end
mean_RB_CSas=mean(RB_CSas);
std_RB_CSas=std(RB_CSas);
sigma_ln_RB_CSas=sqrt(log((std_RB_CSas/mean_RB_CSas)^2+1));
u_ln_RB_CSas=log(mean_RB_CSas)-1/2*sigma_ln_RB_CSas^2;
P_RB_Collapse=logncdf(SaT1,u_ln_RB_CSas,sigma_ln_RB_CSas);
% a. Probability of collapse at the MCE
% building A & B 
P_RB_MCE_A=interp1(SaT1,P_RB_Collapse,MCE_A);
P_RB_MCE_B=interp1(SaT1,P_RB_Collapse,MCE_B);
figure
plot (sort_RB_CSas,Rank_RB,'ro',SaT1,P_RB_Collapse,'linewidth',1.5),hold on;
box on;grid on;
plot(MCE_A,P_RB_MCE_A,'g^','linewidth',1.5);
text(MCE_A+0.2,P_RB_MCE_A,['P(C|IM = Sa_{T1}(',num2str(MCE_A),')',' = ',num2str(P_RB_MCE_A)])
xlabel('Sa_{T1}(g)')
ylabel('Probability of Collapse')
legend('Empirical','Lognormal','Sa_{T1}(MCE)','location','northwest')
title('MCE-A of RB Building')
cd figure
saveas(gcf,'P_f6.jpg')
cd ..
figure
plot (sort_RB_CSas,Rank_RB,'ro',SaT1,P_RB_Collapse,'linewidth',1.5),hold on;
box on;grid on;
plot(MCE_B,P_RB_MCE_B,'g^','linewidth',1.5);
text(MCE_B+0.2,P_RB_MCE_B,['P(C|IM = Sa_{T1}(',num2str(MCE_B),')',' = ',num2str(P_RB_MCE_B)])
xlabel('Sa_{T1}(g)')
ylabel('Probability of Collapse')
legend('Empirical','Lognormal','Sa_{T1}(MCE)','location','northwest')
title('MCE-B of RB Building')
cd figure
saveas(gcf,'P_f7.jpg')
cd ..
% b. The median collapse intensity
P_RB_median=0.5;
SaT1_RB_median=interp1(P_RB_Collapse,SaT1,P_RB_median);
figure
plot (sort_RB_CSas,Rank_RB,'ro',SaT1,P_RB_Collapse,'linewidth',1.5),hold on;
box on;grid on;
plot(SaT1_RB_median,P_RB_median,'ks','linewidth',1.5);
text(SaT1_RB_median+0.2,P_RB_median,['Sa_{T1}(median) = ',num2str(SaT1_RB_median)])
xlabel('Sa_{T1}(g)')
ylabel('Probability of Collapse')
legend('Empirical','Lognormal','Sa_{T1}(median)','location','northwest')
title('Median Collapse Intensity of RB Building')
cd figure
saveas(gcf,'P_f8.jpg')
cd ..
% c. CMR the collaps margin ratio
CMR_RB_A=SaT1_RB_median/MCE_A;
CMR_RB_B=SaT1_RB_median/MCE_B;
% d. The mean annual frequency of collapse
for i =1:length(P_PSHA_A)
    if i==1
        d_lambda_A(i)=P_PSHA_A(i)/dIM;
        d_lambda_B(i)=P_PSHA_B(i)/dIM;
    else if i==length(P_PSHA_A)
        d_lambda_A(i)=(P_PSHA_A(i)-P_PSHA_A(i-1))/dIM;
        d_lambda_B(i)=(P_PSHA_B(i)-P_PSHA_B(i-1))/dIM;
        else
        d_lambda_A(i)=(P_PSHA_A(i+1)-P_PSHA_A(i-1))/(dIM*2);
        d_lambda_B(i)=(P_PSHA_B(i+1)-P_PSHA_B(i-1))/(dIM*2);
        end
    end
    
end
for i=1:length(P_EB_Collapse)    
    lambdaC_RB_A(i)=P_RB_Collapse(i)*abs(d_lambda_A(i));   
    lambdaC_RB_B(i)=P_RB_Collapse(i)*abs(d_lambda_B(i));   
end
MAF_RB_A=sum(lambdaC_RB_A)*dIM; % Annual frequency of collapse for EB A
MAF_RB_B=sum(lambdaC_RB_B)*dIM; % Annual frequency of collapse for EB A
disp (['Annual frequency of collapse for RB A is ',num2str(MAF_RB_A)]);
disp (['Annual frequency of collapse for RB B is ',num2str(MAF_RB_B)]);
figure
plot (SaT1,lambdaC_RB_A,'linewidth',1.5);box on;grid on;
xlabel('Sa_{T1}(g)')
ylabel('P(C|IM=x)x d\lambda/dIM')
title('RB-A Disaggregation of Collapse Risk')
cd figure
saveas(gcf,'P_f9.jpg')
cd ..
figure
plot (SaT1,lambdaC_RB_B,'linewidth',1.5);box on;grid on;
xlabel('Sa_{T1}(g)')
ylabel('P(C|IM=x)x d\lambda/dIM')
title('RB-B Disaggregation of Collapse Risk')
cd figure
saveas(gcf,'P_f10.jpg')
cd ..

%% part g)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Damage fragility functions for structural and non-structural and equapment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structural components
sdr=linspace(0,0.2,1000);
SC_DS1=log(0.005);
SC_DS2=log(0.015);
SC_DS3=log(0.02);
SC_DS4=log(0.05);
SC_std_DS=0.3;
P_SC_DS1=logncdf(sdr,SC_DS1,SC_std_DS);
P_SC_DS2=logncdf(sdr,SC_DS2,SC_std_DS);
P_SC_DS3=logncdf(sdr,SC_DS3,SC_std_DS);
P_SC_DS4=logncdf(sdr,SC_DS4,SC_std_DS);
figure
plot(sdr,P_SC_DS1,sdr,P_SC_DS2,sdr,P_SC_DS3,sdr,P_SC_DS4,'linewidth',1.5);box on; grid on;
xlabel('SDR')
ylabel('P(DS|SDR)')
legend('DS1','DS2','DS3','DS4','location','southeast')
title ('Damage Fragility Function of Structural Components')
cd figure
saveas(gcf,'P_g1.jpg')
cd ..
% non-structural components
NSC_DS1=log(0.001);
NSC_DS2=log(0.005);
NSC_DS3=log(0.015);
NSC_DS4=log(0.025);
NSC_std_DS=0.4;
P_NSC_DS1=logncdf(sdr,NSC_DS1,NSC_std_DS);
P_NSC_DS2=logncdf(sdr,NSC_DS2,NSC_std_DS);
P_NSC_DS3=logncdf(sdr,NSC_DS3,NSC_std_DS);
P_NSC_DS4=logncdf(sdr,NSC_DS4,NSC_std_DS);
figure
plot(sdr,P_NSC_DS1,sdr,P_NSC_DS2,sdr,P_NSC_DS3,sdr,P_NSC_DS4,'linewidth',1.5);box on; grid on;
xlabel('SDR')
ylabel('P(DS|SDR)')
legend('DS1','DS2','DS3','DS4','location','southeast')
title ('Damage Fragility Function of Non-Structural Components')
cd figure
saveas(gcf,'P_g2.jpg')
cd ..
% equipment components
pfa=linspace(0,10,1000);
EC_DS1=log(0.1);
EC_DS2=log(0.5);
EC_DS3=log(1.2);
EC_DS4=log(2.0);
EC_std_DS=0.6;
P_EC_DS1=logncdf(pfa,EC_DS1,EC_std_DS);
P_EC_DS2=logncdf(pfa,EC_DS2,EC_std_DS);
P_EC_DS3=logncdf(pfa,EC_DS3,EC_std_DS);
P_EC_DS4=logncdf(pfa,EC_DS4,EC_std_DS);
figure
plot(pfa,P_EC_DS1,pfa,P_EC_DS2,pfa,P_EC_DS3,pfa,P_EC_DS4,'linewidth',1.5);box on; grid on;
xlabel('PFA')
ylabel('P(DS|PFA)')
legend('DS1','DS2','DS3','DS4','location','southeast')
title ('Damage Fragility Function of Equipment Components')
cd figure
saveas(gcf,'P_g3.jpg')
cd ..

%% part h)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  the expected non-collapse direct(downtime not included) losses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LT_A=1500000;
LT_B=6000000;
% loss due to structural
E_LST_A=[0.08 0.25 0.7 1]*LT_A/3*0.25; % [DS1 DS2 DS3 DS4]; building A at each floor
E_LST_B=[0.08 0.25 0.7 1]*LT_B/3*0.25; % [DS1 DS2 DS3 DS4]; building B at each floor
for i=1:length(P_SC_DS1)
    PSC_DS0(i)=1-P_SC_DS1(i);
    PSC_DS1(i)=P_SC_DS1(i)-P_SC_DS2(i);
    PSC_DS2(i)=P_SC_DS2(i)-P_SC_DS3(i);
    PSC_DS3(i)=P_SC_DS3(i)-P_SC_DS4(i);
    PSC_DS4(i)=P_SC_DS4(i);
end
PSC_DS=[PSC_DS0;PSC_DS1;PSC_DS2;PSC_DS3;PSC_DS4];

L_SC_DS_A=PSC_DS1*E_LST_A(1)+PSC_DS2*E_LST_A(2)+PSC_DS3*E_LST_A(3)+PSC_DS4*E_LST_A(4); % Expected Loss for interior Partion; building A
L_SC_DS_B=PSC_DS1*E_LST_B(1)+PSC_DS2*E_LST_B(2)+PSC_DS3*E_LST_B(3)+PSC_DS4*E_LST_B(4); % Expected Loss for interior Partion; building B
% Existing Building
E_SaT1_SC_EB_A={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building A
E_SaT1_SC_EB_B={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building B
for i=1:NS
    for j=1:length(EB_NCL)
        for k=1:length(P_EB_SDR{i}{j})                                
            if k==1
                E_SaT1_SC_EB_A_temp{i}{j}(k)=abs(L_SC_DS_A(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k))/(0.0002));
                E_SaT1_SC_EB_B_temp{i}{j}(k)=abs(L_SC_DS_B(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k))/(0.0002));
            elseif k == length(P_RB_SDR{i}{j})
               E_SaT1_SC_EB_A_temp{i}{j}(k)=abs(L_SC_DS_A(k)*(P_EB_SDR{i}{j}(k)-P_EB_SDR{i}{j}(k-1))/(0.0002));
               E_SaT1_SC_EB_B_temp{i}{j}(k)=abs(L_SC_DS_B(k)*(P_EB_SDR{i}{j}(k)-P_EB_SDR{i}{j}(k-1))/(0.0002));
            else
               E_SaT1_SC_EB_A_temp{i}{j}(k)=abs(L_SC_DS_A(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k-1))/(2*0.0002));
               E_SaT1_SC_EB_B_temp{i}{j}(k)=abs(L_SC_DS_B(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k-1))/(2*0.0002));
            end
        end
        E_SaT1_SC_EB_A{i}(j)=sum(E_SaT1_SC_EB_A_temp{i}{j})*(0.0002);
        E_SaT1_SC_EB_B{i}(j)=sum(E_SaT1_SC_EB_B_temp{i}{j})*(0.0002);
    end
end
P_EB_NC=1-interp1(SaT1,P_EB_Collapse,EB_NCL'); % probability of non-collapse on each SaT1 level (0.11g 0.22g...)
NC_Loss_SC_EB_A=(E_SaT1_SC_EB_A{1}+E_SaT1_SC_EB_A{2}+E_SaT1_SC_EB_A{3}).*P_EB_NC; % E[Lst|SaT1,NC]*P(NC|SaT1)
NC_Loss_SC_EB_B=(E_SaT1_SC_EB_B{1}+E_SaT1_SC_EB_B{2}+E_SaT1_SC_EB_B{3}).*P_EB_NC;

% Retrofitted Building
E_SaT1_SC_RB_A={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building A
E_SaT1_SC_RB_B={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building B
for i=1:NS
    for j=1:length(RB_NCL)
       for k=1:length(P_RB_SDR{i}{j})
           if k==1
                E_SaT1_SC_RB_A_temp{i}{j}(k)=abs(L_SC_DS_A(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k))/(0.0002));
                E_SaT1_SC_RB_B_temp{i}{j}(k)=abs(L_SC_DS_B(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k))/(0.0002));
            elseif k == length(P_RB_SDR{i}{j})
               E_SaT1_SC_RB_A_temp{i}{j}(k)=abs(L_SC_DS_A(k)*(P_RB_SDR{i}{j}(k)-P_RB_SDR{i}{j}(k-1))/(0.0002));
               E_SaT1_SC_RB_B_temp{i}{j}(k)=abs(L_SC_DS_B(k)*(P_RB_SDR{i}{j}(k)-P_RB_SDR{i}{j}(k-1))/(0.0002));
            else
               E_SaT1_SC_RB_A_temp{i}{j}(k)=abs(L_SC_DS_A(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k-1))/(2*0.0002));
               E_SaT1_SC_RB_B_temp{i}{j}(k)=abs(L_SC_DS_B(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k-1))/(2*0.0002));
           end
        end
        E_SaT1_SC_RB_A{i}(j)=sum(E_SaT1_SC_RB_A_temp{i}{j})*(0.0002);
        E_SaT1_SC_RB_B{i}(j)=sum(E_SaT1_SC_RB_B_temp{i}{j})*(0.0002);
    end
end
P_RB_NC=1-interp1(SaT1,P_RB_Collapse,RB_NCL'); % probability of non-collapse on each SaT1 level (0.11g 0.22g...)
NC_Loss_SC_RB_A=(E_SaT1_SC_RB_A{1}+E_SaT1_SC_RB_A{2}+E_SaT1_SC_RB_A{3}).*P_RB_NC; % E[Lst|SaT1,NC]*P(NC|SaT1)
NC_Loss_SC_RB_B=(E_SaT1_SC_RB_B{1}+E_SaT1_SC_RB_B{2}+E_SaT1_SC_RB_B{3}).*P_RB_NC;

% loss due to non structural
E_LST_A=[0.12 0.2 0.6 1]*LT_A/3*0.5; % [DS1 DS2 DS3 DS4]; building A at each floor
E_LST_B=[0.12 0.2 0.6 1]*LT_B/3*0.5; % [DS1 DS2 DS3 DS4]; building B at each floor
for i=1:length(P_NSC_DS1)
    PNSC_DS0(i)=1-P_NSC_DS1(i);
    PNSC_DS1(i)=P_NSC_DS1(i)-P_NSC_DS2(i);
    PNSC_DS2(i)=P_NSC_DS2(i)-P_NSC_DS3(i);
    PNSC_DS3(i)=P_NSC_DS3(i)-P_NSC_DS4(i);
    PNSC_DS4(i)=P_NSC_DS4(i);
end
PNSC_DS=[PNSC_DS0;PNSC_DS1;PNSC_DS2;PNSC_DS3;PNSC_DS4];

L_NSC_DS_A=PNSC_DS1*E_LST_A(1)+PNSC_DS2*E_LST_A(2)+PNSC_DS3*E_LST_A(3)+PNSC_DS4*E_LST_A(4); % Expected Loss for interior Partion; building A
L_NSC_DS_B=PNSC_DS1*E_LST_B(1)+PNSC_DS2*E_LST_B(2)+PNSC_DS3*E_LST_B(3)+PNSC_DS4*E_LST_B(4); % Expected Loss for interior Partion; building B
% Existing Building
E_SaT1_NSC_EB_A={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building A
E_SaT1_NSC_EB_B={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building B
for i=1:NS
    for j=1:length(EB_NCL)
        for k=1:length(P_EB_SDR{i}{j})                                
            if k==1
                E_SaT1_NSC_EB_A_temp{i}{j}(k)=abs(L_NSC_DS_A(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k))/(0.0002));
                E_SaT1_NSC_EB_B_temp{i}{j}(k)=abs(L_NSC_DS_B(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k))/(0.0002));
            elseif k == length(P_RB_SDR{i}{j})
               E_SaT1_NSC_EB_A_temp{i}{j}(k)=abs(L_NSC_DS_A(k)*(P_EB_SDR{i}{j}(k)-P_EB_SDR{i}{j}(k-1))/(0.0002));
               E_SaT1_NSC_EB_B_temp{i}{j}(k)=abs(L_NSC_DS_B(k)*(P_EB_SDR{i}{j}(k)-P_EB_SDR{i}{j}(k-1))/(0.0002));
            else
               E_SaT1_NSC_EB_A_temp{i}{j}(k)=abs(L_NSC_DS_A(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k-1))/(2*0.0002));
               E_SaT1_NSC_EB_B_temp{i}{j}(k)=abs(L_NSC_DS_B(k)*(P_EB_SDR{i}{j}(k+1)-P_EB_SDR{i}{j}(k-1))/(2*0.0002));
            end
        end
        E_SaT1_NSC_EB_A{i}(j)=sum(E_SaT1_NSC_EB_A_temp{i}{j})*(0.0002);
        E_SaT1_NSC_EB_B{i}(j)=sum(E_SaT1_NSC_EB_B_temp{i}{j})*(0.0002);
    end
end
P_EB_NC=1-interp1(SaT1,P_EB_Collapse,EB_NCL'); % probability of non-collapse on each SaT1 level (0.11g 0.22g...)
NC_Loss_NSC_EB_A=(E_SaT1_NSC_EB_A{1}+E_SaT1_NSC_EB_A{2}+E_SaT1_NSC_EB_A{3}).*P_EB_NC; % E[Lst|SaT1,NC]*P(NC|SaT1)
NC_Loss_NSC_EB_B=(E_SaT1_NSC_EB_B{1}+E_SaT1_NSC_EB_B{2}+E_SaT1_NSC_EB_B{3}).*P_EB_NC;

% Retrofitted Building
E_SaT1_NSC_RB_A={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building A
E_SaT1_NSC_RB_B={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building B
for i=1:NS
    for j=1:length(RB_NCL)
       for k=1:length(P_RB_SDR{i}{j})
           if k==1
                E_SaT1_NSC_RB_A_temp{i}{j}(k)=abs(L_NSC_DS_A(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k))/(0.0002));
                E_SaT1_NSC_RB_B_temp{i}{j}(k)=abs(L_NSC_DS_B(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k))/(0.0002));
            elseif k == length(P_RB_SDR{i}{j})
               E_SaT1_NSC_RB_A_temp{i}{j}(k)=abs(L_NSC_DS_A(k)*(P_RB_SDR{i}{j}(k)-P_RB_SDR{i}{j}(k-1))/(0.0002));
               E_SaT1_NSC_RB_B_temp{i}{j}(k)=abs(L_NSC_DS_B(k)*(P_RB_SDR{i}{j}(k)-P_RB_SDR{i}{j}(k-1))/(0.0002));
            else
               E_SaT1_NSC_RB_A_temp{i}{j}(k)=abs(L_NSC_DS_A(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k-1))/(2*0.0002));
               E_SaT1_NSC_RB_B_temp{i}{j}(k)=abs(L_NSC_DS_B(k)*(P_RB_SDR{i}{j}(k+1)-P_RB_SDR{i}{j}(k-1))/(2*0.0002));
           end
        end
        E_SaT1_NSC_RB_A{i}(j)=sum(E_SaT1_NSC_RB_A_temp{i}{j})*(0.0002);
        E_SaT1_NSC_RB_B{i}(j)=sum(E_SaT1_NSC_RB_B_temp{i}{j})*(0.0002);
    end
end
P_RB_NC=1-interp1(SaT1,P_RB_Collapse,RB_NCL'); % probability of non-collapse on each SaT1 level (0.11g 0.22g...)
NC_Loss_NSC_RB_A=(E_SaT1_NSC_RB_A{1}+E_SaT1_NSC_RB_A{2}+E_SaT1_NSC_RB_A{3}).*P_RB_NC; % E[Lst|SaT1,NC]*P(NC|SaT1)
NC_Loss_NSC_RB_B=(E_SaT1_NSC_RB_B{1}+E_SaT1_NSC_RB_B{2}+E_SaT1_NSC_RB_B{3}).*P_RB_NC;

% Loss due to equipment
E_LST_A=[0.25 0.4 0.8 1]*LT_A/3*0.25; % [DS1 DS2 DS3 DS4]; building A at each floor
E_LST_B=[0.25 0.4 0.8 1]*LT_B/3*0.25; % [DS1 DS2 DS3 DS4]; building B at each floor
for i=1:length(P_EC_DS1)
    PEC_DS0(i)=1-P_EC_DS1(i);
    PEC_DS1(i)=P_EC_DS1(i)-P_EC_DS2(i);
    PEC_DS2(i)=P_EC_DS2(i)-P_EC_DS3(i);
    PEC_DS3(i)=P_EC_DS3(i)-P_EC_DS4(i);
    PEC_DS4(i)=P_EC_DS4(i);
end
PSC_DS=[PEC_DS0;PEC_DS1;PEC_DS2;PEC_DS3;PEC_DS4];

L_EC_DS_A=PEC_DS1*E_LST_A(1)+PEC_DS2*E_LST_A(2)+PEC_DS3*E_LST_A(3)+PEC_DS4*E_LST_A(4); % Expected Loss for interior Partion; building A
L_EC_DS_B=PEC_DS1*E_LST_B(1)+PEC_DS2*E_LST_B(2)+PEC_DS3*E_LST_B(3)+PEC_DS4*E_LST_B(4); % Expected Loss for interior Partion; building B
% Existing Building
E_SaT1_EC_EB_A={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building A
E_SaT1_EC_EB_B={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building B
for i=1:NS
    for j=1:length(EB_NCL)
        for k=1:length(P_EB_PFA{i}{j})                                
            if k==1
                E_SaT1_EC_EB_A_temp{i}{j}(k)=abs(L_EC_DS_A(k)*(P_EB_PFA{i}{j}(k+1)-P_EB_PFA{i}{j}(k))/(0.01));
                E_SaT1_EC_EB_B_temp{i}{j}(k)=abs(L_EC_DS_B(k)*(P_EB_PFA{i}{j}(k+1)-P_EB_PFA{i}{j}(k))/(0.01));
            elseif k == length(P_RB_SDR{i}{j})
               E_SaT1_EC_EB_A_temp{i}{j}(k)=abs(L_EC_DS_A(k)*(P_EB_PFA{i}{j}(k)-P_EB_PFA{i}{j}(k-1))/(0.01));
               E_SaT1_EC_EB_B_temp{i}{j}(k)=abs(L_EC_DS_B(k)*(P_EB_PFA{i}{j}(k)-P_EB_PFA{i}{j}(k-1))/(0.01));
            else
               E_SaT1_EC_EB_A_temp{i}{j}(k)=abs(L_EC_DS_A(k)*(P_EB_PFA{i}{j}(k+1)-P_EB_PFA{i}{j}(k-1))/(2*0.01));
               E_SaT1_EC_EB_B_temp{i}{j}(k)=abs(L_EC_DS_B(k)*(P_EB_PFA{i}{j}(k+1)-P_EB_PFA{i}{j}(k-1))/(2*0.01));
            end
        end
        E_SaT1_EC_EB_A{i}(j)=sum(E_SaT1_EC_EB_A_temp{i}{j})*(0.01);
        E_SaT1_EC_EB_B{i}(j)=sum(E_SaT1_EC_EB_B_temp{i}{j})*(0.01);
    end
end
P_EB_NC=1-interp1(SaT1,P_EB_Collapse,EB_NCL'); % probability of non-collapse on each SaT1 level (0.11g 0.22g...)
NC_Loss_EC_EB_A=(E_SaT1_EC_EB_A{1}+E_SaT1_EC_EB_A{2}+E_SaT1_EC_EB_A{3}).*P_EB_NC; % E[Lst|SaT1,NC]*P(NC|SaT1)
NC_Loss_EC_EB_B=(E_SaT1_EC_EB_B{1}+E_SaT1_EC_EB_B{2}+E_SaT1_EC_EB_B{3}).*P_EB_NC;

% Retrofitted Building
E_SaT1_EC_RB_A={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building A
E_SaT1_EC_RB_B={}; % E[Li|SaT1]; Expected Loss conditioned on each ground motion intensity; building B
for i=1:NS
    for j=1:length(RB_NCL)
       for k=1:length(P_RB_PFA{i}{j})
           if k==1
                E_SaT1_EC_RB_A_temp{i}{j}(k)=abs(L_EC_DS_A(k)*(P_RB_PFA{i}{j}(k+1)-P_RB_PFA{i}{j}(k))/(0.01));
                E_SaT1_EC_RB_B_temp{i}{j}(k)=abs(L_EC_DS_B(k)*(P_RB_PFA{i}{j}(k+1)-P_RB_PFA{i}{j}(k))/(0.01));
            elseif k == length(P_RB_SDR{i}{j})
               E_SaT1_EC_RB_A_temp{i}{j}(k)=abs(L_EC_DS_A(k)*(P_RB_PFA{i}{j}(k)-P_RB_PFA{i}{j}(k-1))/(0.01));
               E_SaT1_EC_RB_B_temp{i}{j}(k)=abs(L_EC_DS_B(k)*(P_RB_PFA{i}{j}(k)-P_RB_PFA{i}{j}(k-1))/(0.01));
            else
               E_SaT1_EC_RB_A_temp{i}{j}(k)=abs(L_EC_DS_A(k)*(P_RB_PFA{i}{j}(k+1)-P_RB_PFA{i}{j}(k-1))/(2*0.01));
               E_SaT1_EC_RB_B_temp{i}{j}(k)=abs(L_EC_DS_B(k)*(P_RB_PFA{i}{j}(k+1)-P_RB_PFA{i}{j}(k-1))/(2*0.01));
           end
        end
        E_SaT1_EC_RB_A{i}(j)=sum(E_SaT1_EC_RB_A_temp{i}{j})*(0.01);
        E_SaT1_EC_RB_B{i}(j)=sum(E_SaT1_EC_RB_B_temp{i}{j})*(0.01);
    end
end
P_RB_NC=1-interp1(SaT1,P_RB_Collapse,RB_NCL'); % probability of non-collapse on each SaT1 level (0.11g 0.22g...)
NC_Loss_EC_RB_A=(E_SaT1_EC_RB_A{1}+E_SaT1_EC_RB_A{2}+E_SaT1_EC_RB_A{3}).*P_RB_NC; % E[Lst|SaT1,NC]*P(NC|SaT1)
NC_Loss_EC_RB_B=(E_SaT1_EC_RB_B{1}+E_SaT1_EC_RB_B{2}+E_SaT1_EC_RB_B{3}).*P_RB_NC;

% Total non-collapse cost
NC_TotalLoss_EB_A=NC_Loss_SC_EB_A+NC_Loss_NSC_EB_A+NC_Loss_EC_EB_A;
NC_TotalLoss_EB_B=NC_Loss_SC_EB_B+NC_Loss_NSC_EB_B+NC_Loss_EC_EB_B;
NC_TotalLoss_RB_A=NC_Loss_SC_RB_A+NC_Loss_NSC_RB_A+NC_Loss_EC_RB_A;
NC_TotalLoss_RB_B=NC_Loss_SC_RB_B+NC_Loss_NSC_RB_B+NC_Loss_EC_RB_B;

figure 
plot(EB_NCL,NC_Loss_SC_EB_A,EB_NCL,NC_Loss_NSC_EB_A,EB_NCL,NC_Loss_EC_EB_A,EB_NCL,NC_TotalLoss_EB_A,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment')
title ('Existing Building (A)')
cd figure
saveas(gcf,'P_h1.jpg')
cd ..

figure
plot(EB_NCL,NC_Loss_SC_EB_B,EB_NCL,NC_Loss_NSC_EB_B,EB_NCL,NC_Loss_EC_EB_B,EB_NCL,NC_TotalLoss_EB_B,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total')
title ('Existing Building (B)')
cd figure
saveas(gcf,'P_h2.jpg')
cd ..

figure 
plot(RB_NCL,NC_Loss_SC_RB_A,RB_NCL,NC_Loss_NSC_RB_A,RB_NCL,NC_Loss_EC_RB_A,RB_NCL,NC_TotalLoss_RB_A,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment')
title ('Retrofitted Building (A)')
cd figure
saveas(gcf,'P_h3.jpg')
cd ..

figure
plot(RB_NCL,NC_Loss_SC_RB_B,RB_NCL,NC_Loss_NSC_RB_B,RB_NCL,NC_Loss_EC_RB_B,RB_NCL,NC_TotalLoss_RB_B,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total')
title ('Retrofitted Building (B)')
cd figure
saveas(gcf,'P_h4.jpg')
cd ..
%% part i)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Total expected collapse losses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LT_A=1500000;
LT_B=6000000;
P_EB_C=interp1(SaT1,P_EB_Collapse,EB_NCL'); % probability of collapse in each ground motion (ex 0.11g 0.22g...)
P_RB_C=interp1(SaT1,P_RB_Collapse,RB_NCL'); % probability of collapse in each ground motion (ex 0.11g 0.22g...)
Total_Loss_C_EB_A=P_EB_C*LT_A;
Total_Loss_C_EB_B=P_EB_C*LT_B;
Total_Loss_C_RB_A=P_RB_C*LT_A;
Total_Loss_C_RB_B=P_RB_C*LT_B;
figure
plot(EB_NCL,Total_Loss_C_EB_A,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Collapse Loss ($)')
title ('Existing Building (A)')
cd figure
saveas(gcf,'P_i1.jpg')
cd ..

figure
plot(EB_NCL,Total_Loss_C_EB_B,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Collapse Loss ($)')
title ('Existing Building (B)')
cd figure
saveas(gcf,'P_i2.jpg')
cd ..

figure
plot(RB_NCL,Total_Loss_C_RB_A,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Collapse Loss ($)')
title ('Retrofitted Building (A)')
cd figure
saveas(gcf,'P_i3.jpg')
cd ..

figure
plot(RB_NCL,Total_Loss_C_RB_B,'linewidth',1.5); box on; grid on;
xlabel ('SaT1 (g)')
ylabel ('Expected Collapse Loss ($)')
title ('Retrofitted Building (B)')
cd figure
saveas(gcf,'P_i4.jpg')
cd ..

%% part j)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% total probability theorem of the expected direct loss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotalLoss_EB_A=Total_Loss_C_EB_A+NC_TotalLoss_EB_A;
TotalLoss_EB_B=Total_Loss_C_EB_B+NC_TotalLoss_EB_B;
TotalLoss_RB_A=Total_Loss_C_RB_A+NC_TotalLoss_RB_A;
TotalLoss_RB_B=Total_Loss_C_RB_B+NC_TotalLoss_RB_B;

figure 
plot(EB_NCL,NC_Loss_SC_EB_A,EB_NCL,NC_Loss_NSC_EB_A,EB_NCL,NC_Loss_EC_EB_A,EB_NCL,NC_TotalLoss_EB_A,EB_NCL,Total_Loss_C_EB_A,EB_NCL,TotalLoss_EB_A,'linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','location','northwest')
title ('Existing Building (A)')
cd figure
saveas(gcf,'P_j1.jpg')
cd ..

figure
plot(EB_NCL,NC_Loss_SC_EB_B,EB_NCL,NC_Loss_NSC_EB_B,EB_NCL,NC_Loss_EC_EB_B,EB_NCL,NC_TotalLoss_EB_B,EB_NCL,Total_Loss_C_EB_B,EB_NCL,TotalLoss_EB_B,'linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','location','northwest')
title ('Existing Building (B)')
cd figure
saveas(gcf,'P_j2.jpg')
cd ..

figure 
plot(RB_NCL,NC_Loss_SC_RB_A,RB_NCL,NC_Loss_NSC_RB_A,RB_NCL,NC_Loss_EC_RB_A,RB_NCL,NC_TotalLoss_RB_A,RB_NCL,Total_Loss_C_RB_A,RB_NCL,TotalLoss_RB_A,'linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','location','northwest')
title ('Retrofitted Building (A)')
cd figure
saveas(gcf,'P_j3.jpg')
cd ..

figure
plot(RB_NCL,NC_Loss_SC_RB_B,RB_NCL,NC_Loss_NSC_RB_B,RB_NCL,NC_Loss_EC_RB_B,RB_NCL,NC_TotalLoss_RB_B,RB_NCL,Total_Loss_C_RB_B,RB_NCL,TotalLoss_RB_B,'linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','location','northwest')
title ('Retrofitted Building (B)')
cd figure
saveas(gcf,'P_j4.jpg')
cd ..

%% part k)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the expected annual direct losses and total direct losses over a 10-year period%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Existing Building
lambda_EB_A=interp1(SaT1,abs(d_lambda_A),EB_NCL'); % d_lambda_A is the slope of the building A Hazard Curve; Computed in part f
lambda_EB_B=interp1(SaT1,abs(d_lambda_B),EB_NCL'); % d_lambda_B is the slope of the building B Hazard Curve; Computed in part f
Annual_direct_losses_EB_A= TotalLoss_EB_A.*lambda_EB_A*0.11;
Annual_direct_losses_EB_B= TotalLoss_EB_B.*lambda_EB_B*0.11;

% Retrofitted Building
lambda_RB_A=interp1(SaT1,abs(d_lambda_A),RB_NCL'); % d_lambda_A is the slope of the building A Hazard Curve; Computed in part f
lambda_RB_B=interp1(SaT1,abs(d_lambda_B),RB_NCL'); % d_lambda_B is the slope of the building B Hazard Curve; Computed in part f
Annual_direct_losses_RB_A= TotalLoss_RB_A.*lambda_RB_A*0.11;
Annual_direct_losses_RB_B= TotalLoss_RB_B.*lambda_RB_B*0.11;

sum_Annual_direct_losses_EB_A=sum(Annual_direct_losses_EB_A);
sum_Annual_direct_losses_EB_B=sum(Annual_direct_losses_EB_B);
sum_Annual_direct_losses_RB_A=sum(Annual_direct_losses_RB_A);
sum_Annual_direct_losses_RB_B=sum(Annual_direct_losses_RB_B);
display('Part k')
display (['Annual direct losses of EB(A) is ',num2str(sum_Annual_direct_losses_EB_A/1000),'k']);
display (['Annual direct losses of EB(B) is ',num2str(sum_Annual_direct_losses_EB_B/1000),'k']);
display (['Annual direct losses of RB(A) is ',num2str(sum_Annual_direct_losses_RB_A/1000),'k']);
display (['Annual direct losses of RB(B) is ',num2str(sum_Annual_direct_losses_RB_B/1000),'k']);
display (['Total direct losses over a 10-year period of EB(A) is ',num2str(sum_Annual_direct_losses_EB_A/100),'k']);
display (['Total direct losses over a 10-year period of EB(B) is ',num2str(sum_Annual_direct_losses_EB_B/100),'k']);
display (['Total direct losses over a 10-year period of RB(A) is ',num2str(sum_Annual_direct_losses_RB_A/100),'k']);
display (['Total direct losses over a 10-year period of RB(B) is ',num2str(sum_Annual_direct_losses_RB_B/100),'k']);


%% part l)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Indirect Loss of Non-collapse & collapse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_SD=[30 240 480];   % time for structural damage. [DS2 DS3 DS4]
Time_NSD=[10 120 240];  % time for non-structural damage. [DS2 DS3 DS4]
Time_ED=[0 30 120];   % time for equipment damage. [DS2 DS3 DS4]
Time_Collapse=730;      % time for collapse damage. [DS2 DS3 DS4]

% Existing Building
for case_number=1:3
    switch case_number
        case 1
            log_DS=[SC_DS1, SC_DS2, SC_DS3, SC_DS4];
            Time=Time_SD;
            std_DS=SC_std_DS;
            u_ln_EB=u_ln_EB_SDR;
            sigma_ln_EB=sigma_ln_EB_SDR;
        case 2
            log_DS=[NSC_DS1, NSC_DS2, NSC_DS3, NSC_DS4];
            Time=Time_NSD;
            std_DS=NSC_std_DS;
            u_ln_EB=u_ln_EB_SDR;
            sigma_ln_EB=sigma_ln_EB_SDR;
        case 3
            log_DS=[EC_DS1, EC_DS2, EC_DS3, EC_DS4];
            Time=Time_ED;
            std_DS=EC_std_DS;
            u_ln_EB=u_ln_EB_PFA;
            sigma_ln_EB=sigma_ln_EB_PFA;
    end
    num_simulations=1000;
    repair_time_temp_EB=0;
    for Sa_level=1:length(EB_NCL)
        for Simulations=1:num_simulations
            for story=1:NS
                u1=rand();  % create a random probability for SDR
                u1_SDR=logninv(u1, u_ln_EB{story}(Sa_level),sigma_ln_EB{story}(Sa_level));
                P_DS=logncdf(u1_SDR,log_DS,std_DS);
                u2=rand();  % create a random value for design DS
                if u2<P_DS(4)
                    repair_time_temp_EB=Time(3);
                elseif u2<P_DS(3)
                    repair_time_temp_EB=Time(2);
                elseif u2<P_DS(2)
                    repair_time_temp_EB=Time(1);    
                end
                repair_time_temp2_EB(story)=repair_time_temp_EB;
            end
            repair_time_temp3_EB(Simulations)=max(repair_time_temp2_EB);
        end
        repair_time_EB(case_number,Sa_level)=sum(repair_time_temp3_EB)/num_simulations; % SC, NSC & EC repair time at each ground motion
    end
end   
max_repair_time_EB=max(repair_time_EB); % maximum repair time(SC, NSC & EC) for each ground motion

% Retrofitted Building
for case_number=1:3
    switch case_number
        case 1
            log_DS=[SC_DS1, SC_DS2, SC_DS3, SC_DS4];
            Time=Time_SD;
            std_DS=SC_std_DS;
            u_ln_RB=u_ln_RB_SDR;
            sigma_ln_RB=sigma_ln_RB_SDR;
        case 2
            log_DS=[NSC_DS1, NSC_DS2, NSC_DS3, NSC_DS4];
            Time=Time_NSD;
            std_DS=NSC_std_DS;
            u_ln_RB=u_ln_RB_SDR;
            sigma_ln_RB=sigma_ln_RB_SDR;
        case 3
            log_DS=[EC_DS1, EC_DS2, EC_DS3, EC_DS4];
            Time=Time_ED;
            std_DS=EC_std_DS;
            u_ln_RB=u_ln_RB_PFA;
            sigma_ln_RB=sigma_ln_RB_PFA;
    end
    num_simulations=1000;
    repair_time_temp_RB=0;
    for Sa_level=1:length(RB_NCL)
        for Simulations=1:num_simulations
            for story=1:NS
                u1=rand();  % create a random probability for SDR
                u1_SDR=logninv(u1, u_ln_RB{story}(Sa_level),sigma_ln_RB{story}(Sa_level));
                P_DS=logncdf(u1_SDR,log_DS,std_DS);
                u2=rand();  % create a random value for design DS
                if u2<P_DS(4)
                    repair_time_temp_RB=Time(3);
                elseif u2<P_DS(3)
                    repair_time_temp_RB=Time(2);
                elseif u2<P_DS(2)
                    repair_time_temp_RB=Time(1);    
                end
                repair_time_temp2_RB(story)=repair_time_temp_RB;
            end
            repair_time_temp3_RB(Simulations)=max(repair_time_temp2_RB);
        end
        repair_time_RB(case_number,Sa_level)=sum(repair_time_temp3_RB)/num_simulations; % SC, NSC & EC repair time at each ground motion
    end
end   
max_repair_time_RB=max(repair_time_RB); % maximum repair time(SC, NSC & EC) for each ground motion

% Indirect Down time Loss NC
down_time_loss_NC_EB_A=max_repair_time_EB*400.*P_EB_NC;
down_time_loss_NC_EB_B=max_repair_time_EB*5000.*P_EB_NC;
down_time_loss_NC_RB_A=max_repair_time_RB*400.*P_RB_NC;
down_time_loss_NC_RB_B=max_repair_time_RB*5000.*P_RB_NC;

% Indirect Down time Loss C
down_time_loss_C_EB_A=730*400.*P_EB_C;
down_time_loss_C_EB_B=730*5000.*P_EB_C;
down_time_loss_C_RB_A=730*400.*P_RB_C;
down_time_loss_C_RB_B=730*5000.*P_RB_C;

% Indirect Down time Total Loss
down_time_totalloss_EB_A=down_time_loss_NC_EB_A+down_time_loss_C_EB_A;
down_time_totalloss_EB_B=down_time_loss_NC_EB_B+down_time_loss_C_EB_B;
down_time_totalloss_RB_A=down_time_loss_NC_RB_A+down_time_loss_C_RB_A;
down_time_totalloss_RB_B=down_time_loss_NC_RB_B+down_time_loss_C_RB_B;

figure 
plot(EB_NCL,NC_Loss_SC_EB_A,EB_NCL,NC_Loss_NSC_EB_A,EB_NCL,NC_Loss_EC_EB_A,EB_NCL,NC_TotalLoss_EB_A,EB_NCL,Total_Loss_C_EB_A,EB_NCL,TotalLoss_EB_A,EB_NCL,down_time_totalloss_EB_A,'--','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','location','northwest')
title ('Existing Building (A)')
cd figure
saveas(gcf,'P_l1.jpg')
cd ..

figure
plot(EB_NCL,NC_Loss_SC_EB_B,EB_NCL,NC_Loss_NSC_EB_B,EB_NCL,NC_Loss_EC_EB_B,EB_NCL,NC_TotalLoss_EB_B,EB_NCL,Total_Loss_C_EB_B,EB_NCL,TotalLoss_EB_B,EB_NCL,down_time_totalloss_EB_B,'--','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','location','northwest')
title ('Existing Building (B)')
cd figure
saveas(gcf,'P_l2.jpg')
cd ..

figure 
plot(RB_NCL,NC_Loss_SC_RB_A,RB_NCL,NC_Loss_NSC_RB_A,RB_NCL,NC_Loss_EC_RB_A,RB_NCL,NC_TotalLoss_RB_A,RB_NCL,Total_Loss_C_RB_A,RB_NCL,TotalLoss_RB_A,RB_NCL,down_time_totalloss_RB_A,'--','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','location','northwest')
title ('Retrofitted Building (A)')
cd figure
saveas(gcf,'P_l3.jpg')
cd ..

figure
plot(RB_NCL,NC_Loss_SC_RB_B,RB_NCL,NC_Loss_NSC_RB_B,RB_NCL,NC_Loss_EC_RB_B,RB_NCL,NC_TotalLoss_RB_B,RB_NCL,Total_Loss_C_RB_B,RB_NCL,TotalLoss_RB_B,RB_NCL,down_time_totalloss_RB_B,'--','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','location','northwest')
title ('Retrofitted Building (B)')
cd figure
saveas(gcf,'P_l4.jpg')
cd ..
%% part m)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the expected loss(direct and indirect) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
plot(EB_NCL,NC_Loss_SC_EB_A,EB_NCL,NC_Loss_NSC_EB_A,EB_NCL,NC_Loss_EC_EB_A,EB_NCL,NC_TotalLoss_EB_A,EB_NCL,Total_Loss_C_EB_A,EB_NCL,TotalLoss_EB_A,EB_NCL,down_time_totalloss_EB_A,'--',EB_NCL,TotalLoss_EB_A+down_time_totalloss_EB_A,'-.','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','Indirect + Direct','location','northwest')
title ('Existing Building (A)')
cd figure
saveas(gcf,'P_m1.jpg')
cd ..

figure
plot(EB_NCL,NC_Loss_SC_EB_B,EB_NCL,NC_Loss_NSC_EB_B,EB_NCL,NC_Loss_EC_EB_B,EB_NCL,NC_TotalLoss_EB_B,EB_NCL,Total_Loss_C_EB_B,EB_NCL,TotalLoss_EB_B,EB_NCL,down_time_totalloss_EB_B,'--',EB_NCL,TotalLoss_EB_B+down_time_totalloss_EB_B,'-.','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','Indirect + Direct','location','northwest')
title ('Existing Building (B)')
cd figure
saveas(gcf,'P_m2.jpg')
cd ..

figure 
plot(RB_NCL,NC_Loss_SC_RB_A,RB_NCL,NC_Loss_NSC_RB_A,RB_NCL,NC_Loss_EC_RB_A,RB_NCL,NC_TotalLoss_RB_A,RB_NCL,Total_Loss_C_RB_A,RB_NCL,TotalLoss_RB_A,RB_NCL,down_time_totalloss_RB_A,'--',RB_NCL,TotalLoss_RB_A+down_time_totalloss_RB_A,'-.','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','Indirect + Direct','location','northwest')
title ('Retrofitted Building (A)')
cd figure
saveas(gcf,'P_m3.jpg')
cd ..

figure
plot(RB_NCL,NC_Loss_SC_RB_B,RB_NCL,NC_Loss_NSC_RB_B,RB_NCL,NC_Loss_EC_RB_B,RB_NCL,NC_TotalLoss_RB_B,RB_NCL,Total_Loss_C_RB_B,RB_NCL,TotalLoss_RB_B,RB_NCL,down_time_totalloss_RB_B,'--',RB_NCL,TotalLoss_RB_B+down_time_totalloss_RB_B,'-.','linewidth',1.5); box on; grid on
xlabel ('SaT1 (g)')
ylabel ('Expected Non-Collapse Loss ($)')
legend('Structural','Non-Structural','Equipment','Total NC','Collpase','Total Direct Loss','Total Indirect Loss','Indirect + Direct','location','northwest')
title ('Retrofitted Building (B)')
cd figure
saveas(gcf,'P_m4.jpg')
cd ..
%% part n)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the expected loss(direct and indirect) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Existing Building
lambda_EB_A=interp1(SaT1,abs(d_lambda_A),EB_NCL'); % d_lambda_A is the slope of the building A Hazard Curve; Computed in part f
lambda_EB_B=interp1(SaT1,abs(d_lambda_B),EB_NCL'); % d_lambda_B is the slope of the building B Hazard Curve; Computed in part f
Annual_direct_losses_EB_A= (TotalLoss_EB_A+down_time_totalloss_EB_A).*lambda_EB_A*0.11;
Annual_direct_losses_EB_B= (TotalLoss_EB_B+down_time_totalloss_EB_B).*lambda_EB_B*0.11;

% Retrofitted Building
lambda_RB_A=interp1(SaT1,abs(d_lambda_A),RB_NCL'); % d_lambda_A is the slope of the building A Hazard Curve; Computed in part f
lambda_RB_B=interp1(SaT1,abs(d_lambda_B),RB_NCL'); % d_lambda_B is the slope of the building B Hazard Curve; Computed in part f
Annual_direct_losses_RB_A= (TotalLoss_RB_A+down_time_totalloss_RB_A).*lambda_RB_A*0.11;
Annual_direct_losses_RB_B= (TotalLoss_RB_B+down_time_totalloss_RB_B).*lambda_RB_B*0.11;

sum_Annual_direct_losses_EB_A=sum(Annual_direct_losses_EB_A);
sum_Annual_direct_losses_EB_B=sum(Annual_direct_losses_EB_B);
sum_Annual_direct_losses_RB_A=sum(Annual_direct_losses_RB_A);
sum_Annual_direct_losses_RB_B=sum(Annual_direct_losses_RB_B);
display('Part n')
display (['Annual direct & indirect losses of EB(A) is ',num2str(sum_Annual_direct_losses_EB_A/1000),'k']);
display (['Annual direct & indirect losses of EB(B) is ',num2str(sum_Annual_direct_losses_EB_B/1000),'k']);
display (['Annual direct & indirect losses of RB(A) is ',num2str(sum_Annual_direct_losses_RB_A/1000),'k']);
display (['Annual direct & indirect losses of RB(B) is ',num2str(sum_Annual_direct_losses_RB_B/1000),'k']);
display (['Total direct & indirect losses over a 10-year period of EB(A) is ',num2str(sum_Annual_direct_losses_EB_A/100),'k']);
display (['Total direct & indirect losses over a 10-year period of EB(B) is ',num2str(sum_Annual_direct_losses_EB_B/100),'k']);
display (['Total direct & indirect losses over a 10-year period of RB(A) is ',num2str(sum_Annual_direct_losses_RB_A/100),'k']);
display (['Total direct & indirect losses over a 10-year period of RB(B) is ',num2str(sum_Annual_direct_losses_RB_B/100),'k']);

%% part o)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Perform a deaggregation on magnitude and distance at Sa_MCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input data
m_min=6;
m_max=8;
L=120;  % km
L_A=20; % km
L_B=120-75; % km
D_A=10; % km
D_B=50; % km
dIM=0.01;
SaT1=0:dIM:4.5; % set up pga
n=length(SaT1); % number of terms
m=linspace(m_min,m_max,10); % Magnitude of the earthquake
b=1;
lamda_m_min=0.01;
r=linspace(10,110,10);        % Distance from the site
load ('PSHA');

% analysis
% Characterizing Distribution of Earthquake Magnitudes
% lamda_m=10^(a-b*m); % lamda_m is number of earthquake with magnitude > m per year
FM=(1-10.^(-b*(m-m_min)))/(1-10.^(-b*(m_max-m_min)));   %  m_min < m < m_max
fM=b*log(10)*10.^(-b*(m-m_min))/(1-10.^(-b*(m_max-m_min))); 
% probability P(M=m);
for i=1:length(FM)
    if i==1
        P_m(i)=FM(i);
    else
        P_m(i)=FM(i)-FM(i-1);
    end
end

% Characterizing Distribution of Distances to Future Earthquakes
% for building A
% compute FR
for i=1:length(r)
    if r(i) <= D_A
        FR_A(i)=0;
    else if r(i) <= sqrt(L_A^2 + D_A^2) && D_A <= r(i)
            FR_A(i)=2*(sqrt(r(i)^2-D_A^2))/L;
        else if r(i) >= sqrt(L_A^2 + D_A^2) && sqrt((L-L_A)^2 + D_A^2) >= r(i)
                FR_A(i)=(L_A+(sqrt(r(i)^2-D_A^2)))/L;
            else
                FR_A(i)=1;
            end
        end
    end
end
% compute fR
for i=1:length(r) 
    if r(i) <= sqrt(L_A^2 + D_A^2) && D_A <= r(i)
        fR_A(i)=2*(r(i)^2-D_A^2)^(-1/2)*r(i)/L;
    else
        fR_A(i)=0;
    end
end
% probability P(R=r);
for i=1:length(FR_A)
    if i==1
        P_r_A(i)=FR_A(i);
    else
        P_r_A(i)=FR_A(i)-FR_A(i-1);
    end
end

for i=1:length(r)
    if r(i) <= D_B
        FR_B(i)=0;
    else if r(i) <= sqrt(L_B^2 + D_B^2) && D_B <= r(i)
            FR_B(i)=2*(sqrt(r(i)^2-D_B^2))/L;
        else if r(i) >= sqrt(L_B^2 + D_B^2) && sqrt((L-L_B)^2 + D_B^2) >= r(i)
                FR_B(i)=(L_B+(sqrt(r(i)^2-D_B^2)))/L;
            else
                FR_B(i)=1;
            end
        end
    end
end
% compute fR
for i=1:length(r) 
    if r(i) <= sqrt(L_B^2 + D_B^2) && D_B <= r(i)
        fR_B(i)=2*(r(i)^2-D_B^2)^(-1/2)*r(i)/L;
    else
        fR_B(i)=0;
    end
end
% probability P(R=r);
for i=1:length(FR_B)
    if i==1
        P_r_B(i)=FR_B(i);
    else
        P_r_B(i)=FR_B(i)-FR_B(i-1);
    end
end


lambda_MCE_A = interp1(SaT1,P_PSHA_A,MCE_A);
lambda_MCE_B = interp1(SaT1,P_PSHA_B,MCE_B);
sigmaln_SaT1=0.6;
P_PSHA_A=zeros(1,length(SaT1));
P_PSHA_B=zeros(1,length(SaT1));
for i=1:length(r)
    for j=1:length(r)
        uln_SaT1=0.420+0.073*(m(i)-6.75)+(-0.59+0.0293*(m(i)-4.5))*log(r(j))-0.007*(r(j)-1)+0.5646;
        P_SaT1_A=1-normcdf(log(MCE_A),uln_SaT1,sigmaln_SaT1); % probability of pga > x
        P_SaT1_B=1-normcdf(log(MCE_B),uln_SaT1,sigmaln_SaT1); % probability of pga > x
        % combine Probability SaT1>x at all magnitude and distance
        P_MR_A(i,j)=(lamda_m_min*P_SaT1_A*P_m(i)*P_r_A(j))/lambda_MCE_A;
        P_MR_B(i,j)=(lamda_m_min*P_SaT1_B*P_m(i)*P_r_B(j))/lambda_MCE_B;
    end
end

figure
bar3(P_MR_A,0.5);xlabel('Distance (Site-to-Fault) (km)');ylabel('Magnitude');zlabel('Probability of Exceedence')
title('Deaggregation at Building A with MCE-A');grid on; 
set(gca,'xticklabel',r);set(gca,'yticklabel',m-0.5)
cd figure
saveas(gcf,'P_o1.jpg')
cd ..
figure
bar3(P_MR_B,0.5);xlabel('Distance (Site-to-Fault) (km)');ylabel('Magnitude');zlabel('Probability of Exceedence')
title('Deaggregation at Building B with MCE-B');grid on; 
set(gca,'xticklabel',r);set(gca,'yticklabel',m-0.5)
cd figure
saveas(gcf,'P_o2.jpg')
cd ..