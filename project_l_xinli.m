%% part l
clc;close all;clear all;
cd data
load P_exceedPFA_RB.mat 
load P_exceedPFA_EB.mat 
load P_exceedSDR_RB.mat 
load P_exceedSDR_EB.mat 
load fragCurve_ST.mat
load fragCurve_NST.mat
load fragCurve_EQT.mat
numSimulations = 1e2;
dIM = 0.11;
EBim = 0.11:dIM:2.42;
RBim = 0.11:dIM:3.3;
IM_col = 0.1:0.01:5; % increment of collapse fragility
Iloss_siteA = 0.4; %in thousands
Iloss_siteB = 5; %in thousands
RT_col = 730;
%% Existing building
repairTime = zeros(1,length(EBim));
repairTime_col = zeros(1,length(EBim));
for indexIM = 1:length(EBim)
    indexIM
    repairTime_IM = zeros(3,3); % row: sotry; col: component
    for indexComponent = 1:3
        % define component parameters
        switch indexComponent
        case 1 % 1: Structual Component
            P_exceedEDP = P_exceedSDR_EB{indexIM};
            fragCurve = fragCurve_ST;
            dEDP = 0.0002; edp = 0:dEDP:0.1; %sdr
        case 2 % 2: Non-Structual Component
            P_exceedEDP = P_exceedSDR_EB{indexIM};
            fragCurve = fragCurve_NST;
            dEDP = 0.0002; edp = 0:dEDP:0.1; %sdr
        case 3 % 3: Equipment Component
            P_exceedEDP = P_exceedPFA_EB{indexIM};
            fragCurve = fragCurve_EQT;
            dEDP = 0.01; edp = 0:dEDP:10; %PFA
        end
     	
        for indexStory = 1:3
            thisP_exceedEDP = P_exceedEDP(indexStory,:);
            
            % total num of simulations
            for indexSimu = 1:numSimulations
                u1 = rand();
                u2 = rand();
                tempIndex1 = find(thisP_exceedEDP<=u1,1);
                if length(tempIndex1) == 0 % if generate u1 is tooo small
                   tempIndex1 = length(edp);
                end
                thisEDP = interp1(thisP_exceedEDP(tempIndex1-1:tempIndex1),edp(tempIndex1-1:tempIndex1),u1);
                for i=1:4
                    thisDSCurve(i) = interp1(edp,fragCurve(i,:),thisEDP);
                end
                if u2 <= thisDSCurve(4) % DS4
                    DS_name = 'DS4';
                elseif u2<= thisDSCurve(3) %DS3
                    DS_name = 'DS3';
                elseif u2<= thisDSCurve(2) %DS2
                    DS_name = 'DS2';
                else %u2<= thisDSCurve(1) %DS1
                    DS_name = 'DS1';
                end
                thisRepairT(indexSimu) = findRepairTime(indexComponent,DS_name);  
            end
            repairTime_IM(indexStory,indexComponent) = sum(thisRepairT)/numSimulations;
        end
    end
    repairTime(indexIM) = max(max(repairTime_IM));
    P_C_IM(indexIM) = interp1(IM_col,colFragCurve_EB,EB_im(indexIM));
    repairTime_col(indexIM) = RT_col * P_NC_IM(indexIM);
end
repairTime_EB = repairTime;
ILoss_NC_EB_siteA = repairTime * Iloss_siteA; %in k;
ILoss_NC_EB_siteB = repairTime * Iloss_siteB;   %in k;
ILoss_C_EB_siteA = repairTime_col * Iloss_siteA; %in k;
ILoss_C_EB_siteB = repairTime_col * Iloss_siteB; %in k;

save ILoss_NC_EB_siteA.mat ILoss_NC_EB_siteA
save ILoss_NC_EB_siteB.mat ILoss_NC_EB_siteB
save ILoss_C_EB_siteA.mat ILoss_C_EB_siteA
save ILoss_C_EB_siteB.mat ILoss_C_EB_siteB
%% Retrofitted building
repairTime = zeros(1,length(RBim));
repairTime_col = zeros(1,length(EBim));
for indexIM = 1:length(RBim)
    indexIM
    repairTime_IM = zeros(3,3); % row: sotry; col: component
    for indexComponent = 1:3
        % define component parameters
        switch indexComponent
        case 1 % 1: Structual Component
            P_exceedEDP = P_exceedSDR_RB{indexIM};
            fragCurve = fragCurve_ST;
            dEDP = 0.0002; edp = 0:dEDP:0.1; %sdr
        case 2 % 2: Non-Structual Component
            P_exceedEDP = P_exceedSDR_RB{indexIM};
            fragCurve = fragCurve_NST;
            dEDP = 0.0002; edp = 0:dEDP:0.1; %sdr
        case 3 % 3: Equipment Component
            P_exceedEDP = P_exceedPFA_RB{indexIM};
            fragCurve = fragCurve_EQT;
            dEDP = 0.01; edp = 0:dEDP:10; %PFA
        end
     	
        for indexStory = 1:3
            thisP_exceedEDP = P_exceedEDP(indexStory,:);
            % total num of simulations
            for indexSimu = 1:numSimulations
                u1 = rand();
                u2 = rand();
                tempIndex1 = find(thisP_exceedEDP<=u1,1);
                if length(tempIndex1) == 0 % if generate u1 is tooo small
                   tempIndex1 = length(edp);
                end
                thisEDP = interp1(thisP_exceedEDP(tempIndex1-1:tempIndex1),edp(tempIndex1-1:tempIndex1),u1);
                for i=1:4
                    thisDSCurve(i) = interp1(edp,fragCurve(i,:),thisEDP);
                end
                if u2 <= thisDSCurve(4) % DS4
                    DS_name = 'DS4';
                elseif u2<= thisDSCurve(3) %DS3
                    DS_name = 'DS3';
                elseif u2<= thisDSCurve(2) %DS2
                    DS_name = 'DS2';
                else %u2<= thisDSCurve(1) %DS1
                    DS_name = 'DS1';
                end
                thisRepairT(indexSimu) = findRepairTime(indexComponent,DS_name);  
            end
            repairTime_IM(indexStory,indexComponent) = sum(thisRepairT)/numSimulations;
        end
    end
    repairTime(indexIM) = max(max(repairTime_IM));
    P_C_IM(indexIM) = interp1(IM_col,colFragCurve_RB,RB_im(indexIM));
    repairTime_col(indexIM) = RT_col * P_NC_IM(indexIM);
end
repairTime_RB = repairTime;
ILoss_NC_RB_siteA = repairTime * Iloss_siteA; %in k;
ILoss_NC_RB_siteB = repairTime * Iloss_siteB;   %in k;
ILoss_C_RB_siteA = repairTime_col * Iloss_siteA; %in k;
ILoss_C_RB_siteB = repairTime_col * Iloss_siteB; %in k;
%% plot


