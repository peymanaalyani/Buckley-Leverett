
% This is the Matlab script base on Buckley-Levertt
% related to advanced Enhanced Oil Recovery course /
% Written by Peyman Alyani / Spring 2022 / At Amirkabir University
% of Technology / Faculty of Petroleum Engineering
% written Based on formulation in book Enhanced Oil Recovery Willhite

clear all
clc

global no nw Siw Sor
% set the Variables
uo=5; % oil viscosity
uw=[1 50]; % waters viscosity for case 1 & case 2
Sw=0:0.001:1; % water saturation
dSw=1/(length(Sw)-1);
Siw=0.15; % water saturation
Sor=0.35; % residual oil saturation
Krw0=0.6; % water relative permeability
Kro0=0.8; % oil relative permeability
nw=2; % for water
no=4; % for oil


for i=1:length(uw)
    M(i)=Krw0*uw(i)/Kro0/uo; % Mobility ratio
end


% Normalized water and oil saturation
Snw=Sw; % Normalized water saturation
Sno=Snw; % Normalized oil saturation
Snw=(Sw-Siw)./(1-Siw-Sor); % calculte Normalized water saturation

for i=1:length(Snw)
   if Snw(i)<=0
       Snw(i)=eps;
   elseif Snw(i)>=1
       Snw(i)=1-eps;
   end
end

Sno=1-Snw; % calculate normalized oil saturation

% Relative Permeability Curve

Krw=Krw0*Snw.^nw;
Kro=Kro0*Sno.^no;

% Fractional flow
fw=ones(length(M),length(Krw));
fo=fw;
for i=1:length(uw)
    fw(i,:)=1./(1+Kro.*uo./(Krw.*uw(i))); % water fractional flow
    fo=1-fw; % oil fractional flow
end


% slope of fractional flow
dfds=fw;
deltafs=fw;
for i=1:length(uw)
    dfds(i,:)=((fw(i,:).^2)./M(i)).*(((1-Snw)).^no)./(Snw.^nw).*(no./(1-Snw)+nw./(Snw)); % derivative of water fractional dlow
    deltafs(i,:)=fw(i,:)./(Snw);
end

% find_shock=dfds-deltafs; % find shock function
Snw_shock=ones(1,length(M));
fw_shock=Snw_shock;
Sw_shock=fw_shock;
dfds_shock=fw_shock;
deltafs_shock=fw_shock;
t_BT=fw_shock;
ER_BT_Snw=fw_shock;

for i=1:length(uw)
    
    
    Snw_shock(1,i)=fzero('find_shock_Com',0.7,[],M(i)); % find the shock normalized water saturation
    fw_shock(1,i)=1/(1+((1-Snw_shock(1,i))^no/(Snw_shock(1,i)^nw))/M(i)); % find the shock water fraction flow
    Sw_shock(1,i)=Snw_shock(1,i)*(1-Siw-Sor)+Siw; % find the shock water saturation
    dfds_shock(1,i)=((fw_shock(1,i)^2)/M(i))*(((1-Snw_shock(1,i))^no)/(Snw_shock(1,i))^nw)*(no/(1-Snw_shock(1,i))+nw/(Snw_shock(1,i))); % derivative of water fractional dlow
    deltafs_shock(1,i)=fw_shock(1,i)/(Snw_shock(1,i)); % slope of fractional flow
    

    t_BT(1,i)=1/deltafs_shock(1,i);
    ER_BT_Snw(1,i)=Snw_shock(1,i)-(fw_shock(1,i)-1)/dfds_shock(1,i); % find ER at water B.T.
    ER_BT_Sw(1,i)=ER_BT_Snw(1,i)*(1-Siw-Sor)+Siw; 
    
    i_Sw_shock(i)=find(Sw>(Sw_shock(1,i))&Sw<(Sw_shock(1,i)+dSw));
    
    if Sw_shock(1,i)>=1-Sor-dSw
        Snw_shock(1,i)=1;
        Sw_shock(1,i)=1-Sor;
        fw_shock(1,i)=1;
        dfds_shock(1,i)=1;
        deltafs_shock(1,i)=1;
        i_Sw_shock(i)=find(Sw>(1-Sor-dSw) & Sw<(1-Sor+dSw));
    end
end

% Calculate the Recovery Factor

ER_WF_Snw=fw;
Wcut_WF_Snw=fw;
Ocut_WF_Snw=fw;
T_WF_Snw=fw;

V_Snw=fw;

for i=1:length(uw) 
    for j=1:length(Snw)
        if j<i_Sw_shock(i)
            V_Snw(i,j)=deltafs_shock(1,i);
            Wcut_WF_Snw(i,j)=0;
        else
            V_Snw(i,j)=dfds(i,j);
            Wcut_WF_Snw(i,j)=fw(i,j);
        end
    end
    
    if Snw(1,i)==1
       for j=1:length(Snw)
           if j<i_Sw_shock
               V_Snw(i,j)=1;
               Wcut_WF_Snw(i,j)=0;
           else
               Wcut_WF_Snw(i,j)=1;
           end
       end
    end

end

Ocut_WF_Snw=1-Wcut_WF_Snw;

for i=1:length(uw)
    T_WF_Snw(i,:)=1./V_Snw(i,:);
    n_ER_WF_Snw(i)=T_WF_Snw(i,1)/(i_Sw_shock(i)-1);
    ER_WF_Snw(i,1:i_Sw_shock(i))=0:n_ER_WF_Snw(i):T_WF_Snw(i,1);
end

for i=1:length(uw)
    for k=1:length(Snw)
        if k<i_Sw_shock(i)
            ER_WF_Snw(i,k)=ER_WF_Snw(i,k);
            T_WF_Snw(i,k)=ER_WF_Snw(i,k);
        else
            ER_WF_Snw(i,k)=Snw(k)-(fw(i,k)-1)/dfds(i,k);
        end
    end
end

ER_WF_Sw=ER_WF_Snw;
ER_WF_Sw=ER_WF_Snw*(1-Sor);

%  plot the relative permeability curve

figure(1) 
plot(Sw(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),Krw(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),'b',...
    Sw(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),Kro(Siw*(length(Sw)-1):(1-Sor)*(length(Sw)-1)),'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Relative Permeability','fontsize',16)
legend({'Krw','Kro'},'fontsize',12)

figure(2) 
plot(Sw,fw(1,:),'b',Sw,fw(2,:),'r','linewidth',2)
axis([0 1 0 1])
xlabel('Sw, Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)
hold on

for i=1:length(uw)
    S_shock(i,:)=[Siw,Sw_shock(i)];
    F_shock(i,:)=[0,fw_shock(i)];
    S_shock_end(i,:)=[Siw,ER_BT_Sw(i)];
    F_shock_end=[0,1];
end

plot(S_shock(1,:),F_shock(1,:),'b-',S_shock(2,:),F_shock(2,:),'r-',...
    S_shock_end(1,:),F_shock_end,'b--',S_shock_end(2,:),F_shock_end,'r--')
hold off

legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))]},'location','best','fontsize',12)

figure(3) 
plot(Snw,fw(1,:),'b',Snw,fw(2,:),'r','linewidth',2)
axis([0 1 0 1])
xlabel('Snw, Normalized Water Saturation','fontsize',16)
ylabel('Water Fractional Flow','fontsize',16)

hold on

for i=1:length(uw)
    Sn_shock(i,:)=[0,Snw_shock(i)];
    F_shock(i,:)=[0,fw_shock(i)];
    Sn_shock_end(i,:)=[0,ER_BT_Snw(i)];
    F_shock_end=[0,1];
end

plot(Sn_shock(1,:),F_shock(1,:),'b-',Sn_shock(2,:),F_shock(2,:),'r-',...
    Sn_shock_end(1,:),F_shock_end,'b--',Sn_shock_end(2,:),F_shock_end,'r--','linewidth',1)
hold off

legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))]},'location','best','fontsize',12)

figure(4) 
plot(T_WF_Snw(1,:),Wcut_WF_Snw(1,:),'b',T_WF_Snw(2,:),Wcut_WF_Snw(2,:),'r','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Water Cut','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))]},'location','best','fontsize',12)

figure(5) 
plot(T_WF_Snw(1,:),Ocut_WF_Snw(1,:),'b',T_WF_Snw(2,:),Ocut_WF_Snw(2,:),'r','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Oil Cut','fontsize',16)
legend({['M^0=',num2str(M(1))],['M^0=',num2str(M(2))]},'location','best','fontsize',12)

figure(6) 
plot(T_WF_Snw(1,:),ER_WF_Snw(1,:),'b',T_WF_Snw(2,:),ER_WF_Snw(2,:),'r','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Recovery Factor, Normalized','fontsize',16)
legend({'WatVis=50 cP','WatVis=1 cP'},'location','NorthEast','fontsize',11)

figure(7) 
plot(T_WF_Snw(1,:),ER_WF_Sw(1,:),'b',T_WF_Snw(2,:),ER_WF_Sw(2,:),'r','linewidth',2)
axis([0 5 0 1])
xlabel('Dimensionless Time, MPV','fontsize',16)
ylabel('Recovery Factor','fontsize',16)
legend({'WatVis=50 cP','WatVis=1 cP'},'location','NorthEast','fontsize',11)
figure(8) 
t=0.1;%td=0.1 pv
X_Snw=t.*V_Snw;
plot(X_Snw(1,:),Sw,'b',X_Snw(2,:),Sw,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance, T_D=0.1 pv','fontsize',16)
ylabel('Water Saturation','fontsize',16)
legend({'WatVis=50 cP','WatVis=1 cP'},'location','NorthEast','fontsize',11)

figure(9) 
t=0.3;%td=0.3 pv
X_Snw=t.*V_Snw;
plot(X_Snw(1,:),Sw,'b',X_Snw(2,:),Sw,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance, T_D=0.3 pv','fontsize',16)
ylabel('Water Saturation','fontsize',16)
legend({'WatVis=50 cP','WatVis=1 cP'},'location','NorthEast','fontsize',11)

figure(10) 
t=0.7;%td=0.7 pv
X_Snw=t.*V_Snw;
plot(X_Snw(1,:),Sw,'b',X_Snw(2,:),Sw,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance, T_D=0.7 pv','fontsize',16)
ylabel('Water Saturation','fontsize',16)
legend({'WatVis=50 cP','WatVis=1 cP'},'location','NorthEast','fontsize',11)

figure(11) 
t=1;%td=1 pv
X_Snw=t.*V_Snw;
plot(X_Snw(1,:),Sw,'b',X_Snw(2,:),Sw,'r','linewidth',2)
axis([0 1 0 1])
xlabel('Dimensionless Distance, T_D=1 pv','fontsize',16)
ylabel('Water Saturation','fontsize',16)
legend({'WatVis=50 cP','WatVis=1 cP'},'location','NorthEast','fontsize',11)