% This script post-process the exported DRW model numerical results
% and compare (validate) the numerical results against Snyder-Lumley
% experimental data.

clc;
close all;
clear all; 

%===========User inputs==========%
Velocity_at_inlet=2;      % in [m/s]
U_ave=Velocity_at_inlet;

%Value of time scale constant coefficient from RANS FLUENT
CL=0.01;     % in [-] (NOTE! CT=2*CL) 
%=======End User inputs==========%

%===Loading DRW Data (Extracted .txt data from RANS FLUENT via UDF)===%
n_data_drw=10;      %Number of exported data columns
file='FlowTracers_TidalChannel_GITaylor_Calibration_500try_CL0p01.txt';
fid=fopen(file,'r');

[data,count]=fscanf(fid,'%f',[n_data_drw inf]);
fclose(fid);

%Particle coordinate
X_drw_all(:,1)=data(1,:);  
Y_drw_all(:,1)=abs(data(2,:)-50);  %Inlet center is at (0,50,0) 
Z_drw_all(:,1)=data(3,:);

%Particle mean velocity
v_drw_all(:,1)=data(4,:);

%Particle fluctuation velocity
u_prim_drw_all(:,1)=data(5,:);
v_prim_drw_all(:,1)=data(6,:);
w_prim_drw_all(:,1)=data(7,:);
 
%Particle ID and time scales
P_ID_all(:,1)=data(9,:)+1;  %"+1" is to start from Particle ID from 1

T_drw_all(:,1)=data(8,:); 
TL_drw_all(:,1)=data(10,:);
number_of_tries=max(P_ID_all);
%===End of loading DRW Data===%

%==========Sorting loaded DRW model data=============%

%Finding the minimum number of data for each realization (stop). This is
%required since DRW is a stochastic model and path of different particle
%realization is different!

for m=1 : number_of_tries
    [column]=find(P_ID_all==m);
    stop_mat(m,1)=length(column);
end
stop=min(stop_mat);

for m=1 : number_of_tries
    
    %Isolating the values of EACH realization for DRW model.
    index=find(P_ID_all==m);
    X_drw(:,1)=X_drw_all(index);
    Y_drw(:,1)=Y_drw_all(index);
    Z_drw(:,1)=Z_drw_all(index);
    
    v_drw(:,1)=v_drw_all(index);
    
    u_prim_drw(:,1)=u_prim_drw_all(index);     
    v_prim_drw(:,1)=v_prim_drw_all(index);
    w_prim_drw(:,1)=w_prim_drw_all(index);
    
    T_drw(:,1)=T_drw_all(index);  
    
    %Saving the sorted positions of each realization
    Xp(1:stop,m) = X_drw(1:stop,1);
    Yp(1:stop,m) = Y_drw(1:stop,1);
    Zp(1:stop,m) = Z_drw(1:stop,1);
    
    %Calculating the displacement of each realization in each direction
    dx(1:stop,m) = X_drw(1:stop,1) - 0;
    dy(1:stop,m) = Y_drw(1:stop,1) - (U_ave*T_drw(1:stop,1));
    dz(1:stop,m) = Z_drw(1:stop,1) - 0;
    
    %Saving the sorted vel. fluctuation of each realization
    v_prim(1:stop,m) = v_prim_drw(1:stop,:);
    u_prim(1:stop,m) = u_prim_drw(1:stop,:);
    w_prim(1:stop,m) = w_prim_drw(1:stop,:);

    %Saving the sorted time steps of each realization
    Tp(1:stop,m) = T_drw(1:stop,:);
    
     if m~=number_of_tries
         clear X_drw Y_drw Z_drw...
               v_drw u_prim_drw v_prim_drw w_prim_drw...
               T_drw
     end
   
end 
%=======End Sorting DRW model data for post-processing=============%

%===Result(1) - Plot of all particle realizations dispersion===%
figure(1)
plot(Yp,Zp,'.')
title(['Realizationzs path - C_L = ',num2str(CL)])
xlabel('Y [m]')
ylabel('Z [m]')
xlim([min(Yp(1,:)) min(Yp(stop,:))])
ylim([min(Zp(stop,:)) max(Zp(stop,:))])
grid on
%===End of Result (1)===%

%===Evaluation of the mean Lagrangian time scale (TL) predicted by DRW===%
count=1;
for n=1 : length(u_prim_drw_all)-1
    if  abs(u_prim_drw_all(n+1)-u_prim_drw_all(n))>10^-3 &&...
        TL_drw_all(n+1,1)~=0
        
        TL(count,1)=TL_drw_all(n+1,1);
        count=count+1;
        
    end
end
I_DRW=mean(TL)
%====End of TL calculation============================% 
 
clear X_drw_all Y_drw_all Z_drw_all V_drw_all 
clear u_prim_drw_all v_prim_drw_all w_prim_drw_all
clear T_drw_all

%===Statistical analysis on sorted DRW data===%
for j=1 : stop

    dy_rms(j,1)=std(dy(j,:),1);
    dx_rms(j,1)=std(dx(j,:),1);
    dz_rms(j,1)=std(dz(j,:),1);

    u_prime_rms(j,1)=std(u_prim(j,:),1);
    v_prime_rms(j,1)=std(v_prim(j,:),1);
    w_prime_rms(j,1)=std(w_prim(j,:),1);
    
    
    Tp_ave(j,1)=mean(Tp(j,:));
     
end
%===End of statistical analysis on sorted DRW data===%

%=======Autocorrelation function and it's integral==========%
count2=1;
for t=2
    count1=1;
        for tau=1 : stop-t
            for j=1:number_of_tries
        
                R_num(1,j) = (w_prim(t,j)*w_prim(t+tau,j));
                R_denum1(1,j)=(w_prim(t,j)^2);   
                R_denum2(1,j)=(w_prim(t+tau,j)^2);
             
            end

            R_tau(tau,1)=mean(R_num)/(sqrt(mean(R_denum1))*sqrt(mean(R_denum2)));
            kesi(tau,1)=Tp_ave(t+tau)-Tp_ave(t);
            clear R_num;clear R_denum;clear R_denum1;clear R_denum2;
        end
    %Limiting integration to the point were autocorrelation functio
    %becomes zero for the first time. The value of this integral will
    %more accurate, limiting will not be required, by increasing number
    %of particles.
    I_Taylor(count2,1)=trapz(R_tau(1:70,1))*mean(diff(Tp_ave(1:70,1)))
    
    figure(2)
    plot(kesi,R_tau,'r*-')
    title(['Velocity Autocorrelation Function - C_L = ',num2str(CL)],'FontSize',14,'FontWeight','bold','Color','k')
    xlabel('\zeta','FontSize',14,'FontWeight','bold','Color','k')
    ylabel('R_{\zeta}','FontSize',14,'FontWeight','bold','Color','k')
    grid on
    ylim([min(R_tau) 1])
    count2=count2+1;
    hold on
end
%==============%


%===Comparison of DRW data against G.I. taylor dispersion theory===%

%Evaluation of Root Mean Square of particle dispersion in three different
%directions based on the G.I. Taylor's disperssion theory.
LHS_u=(sqrt(2*I_Taylor*Tp_ave)).*u_prime_rms;
LHS_v=(sqrt(2*I_Taylor*Tp_ave)).*v_prime_rms;
LHS_w=(sqrt(2*I_Taylor*Tp_ave)).*w_prime_rms;

figure(3)
plot(Tp_ave,dx_rms,'.')
hold on
plot(Tp_ave,LHS_u(1:stop,1),'r.')
hold on
title(['RMS dispersion in X-Dir - C_L = ',num2str(CL)])
xlabel('t [sec]')
ylabel('RMS(X)')
hleg1 = legend('DRW prediction','G.I Taylor theory');
set(hleg1,'Location','NorthWest')
grid on

figure(4)
plot(Tp_ave,dy_rms,'.')
hold on
plot(Tp_ave,LHS_v(1:stop,1),'r.')
hold on
title(['RMS dispersion in Y-Dir - C_L = ',num2str(CL)])
xlabel('t [sec]')
ylabel('RMS(Y)')
hleg1 = legend('DRW prediction','G.I Taylor theory');
set(hleg1,'Location','NorthWest')
grid on

figure(5)
plot(Tp_ave,dz_rms,'.')
hold on
plot(Tp_ave,LHS_w(1:stop,1),'r.')
hold on
title(['RMS dispersion in Z-Dir - C_L = ',num2str(CL)])
xlabel('t [sec]')
ylabel('RMS(Z)')
hleg1 = legend('DRW prediction','G.I Taylor theory');
set(hleg1,'Location','NorthWest')
grid on
%===End of comparison of DRW data against G.I. taylor dispersion theory===%
