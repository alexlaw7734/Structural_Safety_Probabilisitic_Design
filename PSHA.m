% CEE_244_PSHA (probabilistic seismic hazard analysis)
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
m=linspace(m_min,m_max,n); % Magnitude of the earthquake
b=1;
lamda_m_min=0.01;
r=linspace(9,102,n);        % Distance from the site


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

% figure
% plot(m,FM); box on; grid on;
% title('FM')
% figure
% plot(m,fM); box on; grid on;
% title('fM')


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

% figure
% plot(r,FR_A); box on; grid on;
% title('FR_A')
% figure
% plot(r,fR_A); box on; grid on;
% title('fR_A')

% for building B
% compute FR
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

% figure
% plot(r,FR_A); box on; grid on;
% title('FR_B')
% figure
% plot(r,fR_A); box on; grid on;
% title('fR_B')

% Predicting Ground Motion Intensity
% dIM=0.01;
% SaT1=0:dIM:4.5; % set up pga
sigmaln_SaT1=0.6;
P_PSHA_A=zeros(1,length(SaT1));
P_PSHA_B=zeros(1,length(SaT1));
for k=1:length(SaT1)
    for i=1:length(SaT1)
        for j=1:length(SaT1)
            uln_SaT1=0.420+0.073*(m(i)-6.75)+(-0.59+0.0293*(m(i)-4.5))*log(r(j))-0.007*(r(j)-1)+0.5646;
            P_SaT1=1-normcdf(log(SaT1(k)),uln_SaT1,sigmaln_SaT1); % probability of pga > x
            % combine Probability SaT1>x at all magnitude and distance
            P_PSHA_A(k)=P_PSHA_A(k)+lamda_m_min*P_SaT1*P_m(i)*P_r_A(j);
            P_PSHA_B(k)=P_PSHA_B(k)+lamda_m_min*P_SaT1*P_m(i)*P_r_B(j);
        end
    end
end
% 2% probability of exceedance in 50 years
lamda=-log(1-0.02)/50;
returnperiod=1/lamda;
%
MCE_A=interp1(P_PSHA_A,SaT1,lamda);
MCE_B=interp1(P_PSHA_B,SaT1,lamda);
%%
figure
semilogy(SaT1,P_PSHA_A,SaT1,P_PSHA_B,'--','linewidth',2);hold on;
plot(MCE_A,lamda,'o',MCE_B,lamda,'g*')
box on;grid on
title('Hazard Curve')
xlabel('SaT1')
ylabel('Probability of exceedance')
legend ('Building A','Building B',['2% Probability of exceedance, MCE-A is ',num2str(MCE_A)],['2% Probability of exceedance, MCE-B is ',num2str(MCE_A)])
saveas(gcf,'Hazard Curve.jpg')

% 2% probability of exceedance in 50 years
lamda=-log(1-0.02)/50;
returnperiod=1/lamda;
%
MCE_A=interp1(P_PSHA_A,SaT1,lamda);
MCE_B=interp1(P_PSHA_B,SaT1,lamda);
%%
save ('PSHA','SaT1', 'P_PSHA_A', 'P_PSHA_B','MCE_A','MCE_B','dIM')