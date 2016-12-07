% This script post-process the exported DRW model numerical results
% and compare (validate) the numerical results against Snyder-Lumley
% experimental data.

clc;
close all;
clear all; 

%===========User inputs==========%
Velocity_at_inlet=6.5;      % in [m/s]
U_ave=Velocity_at_inlet;

%Value of time scale constant coefficient from RANS FLUENT
CL=0.15;     % in [-] (NOTE! CT=2*CL) 
%=======End User inputs==========%

%===Loading DRW Data (Extracted .txt data from RANS FLUENT via UDF)===%
n_data_drw=10;      %Number of exported data columns
file='HollowGlass_Particle_SL_Exp_validation_500try_CT0p3.txt';
fid=fopen(file,'r');

[data,count]=fscanf(fid,'%f',[n_data_drw inf]);
fclose(fid);

%Particle coordinate
X_drw_all(:,1)=data(1,:)-0.5588;   %Inlet center is at (0.5588,0,0)
Y_drw_all(:,1)=data(2,:); 
Z_drw_all(:,1)=data(3,:);

%Particle mean velocity
u_drw_all(:,1)=data(4,:);

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
    
    u_drw(:,1)=u_drw_all(index);
    
    u_prim_drw(:,1)=u_prim_drw_all(index);     
    v_prim_drw(:,1)=v_prim_drw_all(index);
    w_prim_drw(:,1)=w_prim_drw_all(index);
    
    T_drw(:,1)=T_drw_all(index);  
    
    %Saving the sorted positions of each realization
    Xp(1:stop,m) = X_drw(1:stop,1);
    Yp(1:stop,m) = Y_drw(1:stop,1);
    Zp(1:stop,m) = Z_drw(1:stop,1);
    
    %Calculating the displacement of each realization in each direction
    dx(1:stop,m) = X_drw(1:stop,1) - (U_ave*T_drw(1:stop,1));
    dy(1:stop,m) = Y_drw(1:stop,1) - 0;
    dz(1:stop,m) = Z_drw(1:stop,1) - 0;
    
    %Saving the sorted vel. fluctuation of each realization
    v_prim(1:stop,m) = v_prim_drw(1:stop,:);
    u_prim(1:stop,m) = u_prim_drw(1:stop,:);
    w_prim(1:stop,m) = w_prim_drw(1:stop,:);

    %Saving the sorted time steps of each realization
    Tp(1:stop,m) = T_drw(1:stop,:);
    
     if m~=number_of_tries
         clear X_drw Y_drw Z_drw...
               u_drw u_prim_drw v_prim_drw w_prim_drw...
               T_drw
     end
   
end 
%=======End Sorting DRW model data for post-processing=============%

%===Result(1) - Plot of all particle realizations dispersion===%
figure(1)
plot(Zp,Xp,'.')
title(['Realizationzs path - C_L = ',num2str(CL)])
xlabel('Z [m]')
ylabel('X [m]')
xlim([-0.06 0.06])
grid on
axis([min(Zp(stop,:)) max(Zp(stop,:)) 0 min(Xp(stop,:))])
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
    
    I_Taylor(count2,1)=trapz(R_tau)*mean(diff(Tp_ave))
    
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


%Plot of numerical (DRW) vs. Experimental (SL) lateral dispersion (z-axis)

figure(3)
%Numerical (DRW) data
plot(Tp_ave.*10^3,(dz_rms.^2)*10^4,'r.')
title(['Mean square of hollow glass particle lateral dispersion - C_L = ',num2str(CL)],'FontSize',14,'FontWeight','bold','Color','k')
hold on

%Experimental (SL) data extracted from SL original paper
plot(412.407862, 4.517544,'kO','markers',9)
hold on
plot(339.434889, 3.421053,'kO','markers',9)
hold on
plot(275.307125, 2.631579,'kO','markers',9)
hold on
plot(218.918919, 1.907895,'kO','markers',9)
hold on
plot(169.164619, 1.469298,'kO','markers',9)
hold on
plot(126.044226, 0.986842,'kO','markers',9)
hold on
plot(95.085995,	0.635965,'kO','markers',9)
hold on
plot(56.388206,	0.285088,'kO','markers',9)
hold on
plot(25.429975,	0.131579,'kO','markers',9)
hold on

title(['Particle Mean Square Lateral Dispersion - C_L = ',num2str(CL)])
xlabel('t [ms]')
ylabel('Mean Square Z [cm2]')
hleg1 = legend('DRW prediction','SL Experiment');
set(hleg1,'Location','NorthWest')
xlim([0 500])
ylim([0 7])
grid on
