clc
clear
gamma=1.4; M1=2.75; P1=12.11*(10^3); T1=216.5; rho1=P1/(T1*287); [~, ~, P, ~, ~] = flowisentropic(gamma, M1); P10=P1/P;
X=[1.1000	1.2000	1.3000	1.4000	1.5000	1.6000	1.7000	1.8000	1.9000	2.0000	2.2000	2.4000	2.6000	2.8000	3.0000	3.2000	3.4000	3.6000	3.8000	4.0000	5.0000	6.0000	7.0000	8.0000	9.0000	10.0000];
Y=[1.5152	3.9442	6.6621	9.4272	12.1127	14.6515	17.0119	19.1833	21.1675	22.9735	26.1028	28.6814	30.8137	32.5875	34.0734	35.3275	36.3934	37.3059	38.0922	38.7739	41.1177	42.4398	43.2546	43.7908	44.1619	44.4290];
c=polyfit(X,Y,6);
theta_max=poly2sym(c);
clear c Y X
fplot(theta_max, [ 1.1, 10]);
AA=zeros(1, 3500);
P40_max=0;
i=1;
for theta1=0.1 : 0.5 : double(subs(theta_max,M1))- 0.5
[M2,beta1,~,~] = OSW(M1,theta1);
    M1n=M1*sind(beta1);
    [~, ~, ~, ~, ~, P0, ~] = flownormalshock(gamma, M1n);
    P20=P10*P0;
     for theta2=.1: 0.5: double(subs(theta_max,M2))-0.5   
         [M3,beta2,~,~] = OSW(M2,theta2);
         M2n=M2*sind(beta2);
         [~, ~, ~, ~, ~, P0, ~] = flownormalshock(gamma, M2n);
         P30=P20*P0;
        if M3>1.000
            [~, ~, ~, ~, ~, P0, ~] = flownormalshock(gamma, M3);
            P40=P30*P0;
        else
            P40=P30;
        end
        AA(i)=P40;
        temp=max(AA);
        if temp> P40_max
            P40_max=temp;
            theta1_opt=theta1;
            theta2_opt=theta2;
        end
        i=i+1;
    end   
end

[M2,beta1,~,~] = OSW(M1,theta1_opt);
 M1n=M1*sind(beta1);
[~, T, P, rho, ~, P0, ~] = flownormalshock(gamma, M1n);
P2=P*P1;
T2=T*T1;
rho2=rho*rho1;
P20=P10*P0;
[M3,beta2,~,~] = OSW(M2,theta2_opt);
 M2n=M2*sind(beta2);
[~, T, P, rho, ~, P0, ~] = flownormalshock(gamma, M2n);
P3=P*P2;
T3=T*T2;
rho3=rho*rho2;
P30=P20*P0;
[~, T, P, rho, ~, P0, ~] = flownormalshock(gamma, M3);
M4= shock_wave(M3,gamma);
P4=P*P3;
T4=T*T3;
rho4=rho*rho3;
P40=P30*P0;

A_in=0.0001:0.0005:0.0625; 
for i=1:length(A_in)
    A_inD = A_in(i);
M_inD=M4;
[mach, T, P, rho, area] = flowisentropic(gamma, M_inD);
T0_inD=T4/T; 
A_star=(A_inD/area);
A_outD = 0.0625;
area2 = (A_outD/A_star);
[M_outD, T2, P2, rho2, area2] = flowisentropic(gamma, area2,'sub');
T_outD(i)=T2*T0_inD;
M_outD2(i)=M_outD;
P_outD(i)=P40*P2;
P0_outD(i)=P40;
T0_outD(i)=T0_inD;
M_inH=M_outD2(i); T_tot1=T0_outD(i); 
P01=P0_outD(i); P1=P_outD(i); T1=T_outD(i);
[mach_inH, T, P, rho, velocity, T0, P0] = flowrayleigh(gamma, M_inH);
T_star_tot=T_tot1/T0;
P0_star=P01/P0;
Pstar=P1/P;
Tstar=T1/T;
if (T_star_tot <= 2200)
        M_outH1=1;
        T0_outH=T_star_tot;
        [mach_outH, T2, P2, rho2, velocity2, T02, P02] = flowrayleigh(gamma, M_outH1);          
elseif(T_star_tot>2200) 
        T0_outH=2200;
        T_ratio=T0_outH/T_star_tot;
        [mach_outH, T2, P2, rho2, velocity2, T02, P02] = flowrayleigh(gamma, T_ratio,'totaltsub');
end
M_outH(i)=mach_outH;
T0_outH2(i)=T0_outH;
T_outH(i)=T2*Tstar;
P_outH(i)=P2*Pstar;
P0_outH(i)=P02*P0_star;
p0_in_N(i)=P0_outH(i)*0.95;
Mi=M_outH(i); Pe=1.01*(10^5); 
A=0.0625; T0=T0_outH2(i); 
P_ratio_N_exit=Pe/p0_in_N(i);
[mach, T, P, rho, area] = flowisentropic(gamma, P_ratio_N_exit, 'pres');
Te=T*T0; 
if Mi==1
    A_throat=A;
    Ae=area*A_throat;
elseif Mi<1
    Ae=A;
    A_throat=Ae/area;
end
Me_N(i)=mach;
Te_N(i)=Te;
Pe_N(i)=Pe;
A_star_N(i)=A_throat;
Ae_N(i)=Ae;
rho_e(i)=rho;
P1=12.11*(10^3); M1=2.75; T1=216.5;
mdot(i)=rho4*A_inD*M4*sqrt(gamma*287*T4);
V1=M1*sqrt(gamma*287*T1);
V2=Me_N(i)*sqrt(gamma*287*Te_N(i));
Thrust(i)=mdot(i)*(V2-V1)+(Pe_N(i)-P1)*A_in(i);
end 
plot(A_in,Thrust);
xlabel("Inlet Area (m^2)");
ylabel("Thrust (N)")
[Max_Thrust, Max_thrust_Index]=max(Thrust);
A_opt=A_in(Max_thrust_Index);
mdot_opt=mdot(Max_thrust_Index);
M_out_Diffuser=M_outD2(Max_thrust_Index);
P0_out_Diffuser=P0_outD(Max_thrust_Index);
P_out_Diffuser=P_outD(Max_thrust_Index);
T_out_Diffuser=T_outD(Max_thrust_Index);
T0_out_Diffuser=T0_outD(Max_thrust_Index);
M_out_Combustor=M_outH(Max_thrust_Index);
P0_out_Combustor=P0_outH(Max_thrust_Index);
T0_out_Combustor=T0_outH2(Max_thrust_Index);
P_out_Combustor=P_outH(Max_thrust_Index);
T_out_Combustor=T_outH(Max_thrust_Index);
P0_in_Nozzle=p0_in_N(Max_thrust_Index);
T_out_Nozzle=Te_N(Max_thrust_Index);
P_out_Nozzle=Pe_N(Max_thrust_Index);
M_out_Nozzle=Me_N(Max_thrust_Index);
Ae_Nozzle=Ae_N(Max_thrust_Index);
Ath_Nozzle=A_star_N(Max_thrust_Index);
%Finding wedge length
Wedge2=A_opt/tand(beta2-theta2_opt);
Hyp=A_opt/sind(beta2-theta2_opt);
angle1=beta2-beta1+theta1_opt;
angle2=beta1-theta1_opt;
Wedge1=Hyp/sind(angle2)*sind(angle1);