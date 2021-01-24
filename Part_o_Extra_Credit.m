% CEE_244_Project_Part_o_Extra_Credit
% Weiyang Bai
clc; clear all; close all
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

%%
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
%%
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
