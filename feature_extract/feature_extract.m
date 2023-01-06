clear;clc;close all;
%load data
data=csvread('a_pulse.csv');
%select dataset
t0=100000;
t1=400000;
time=data(:,1);
U=data(:,7);
%voltage to resistance
R0=U.*1000./(4.78-U);
N=length(R0);
t=0:0.001:(N-1)*0.001;
lev=6;
wf='sym8';
%denoise
s_s=wden(R0,'heursure','s','one',lev,wf);
%base
[c,l]=wavedec(R0,11,'sym8');
base=wrcoef('a',c,l,'sym8',11);
%denosie and debase
s_ss=s_s-base;
%noise
err=R0-s_ss;
%normalization
s=s_ss./base;

%plot
figure(1)
subplot(6,1,1)
plot(t,R0);xlabel('t/s');ylabel('origin');
%axis([0,3,2.4,3.4]);
subplot(6,1,2)
plot(t,s_s);xlabel('t/s');ylabel('denoise');
subplot(6,1,3)
plot(t,s_ss);xlabel('t/s');ylabel('debase');
subplot(6,1,4)
plot(t,err);xlabel('t/s');ylabel('noise');
subplot(6,1,5)
plot(t,base);xlabel('t/s');ylabel('base') ;
subplot(6,1,6)
plot(t,s);xlabel('t/s');ylabel('normalization');

s=s(t0:t1,1);
t=t(1,t0+1:t1+1);
%find p peak
pulse=find(s>max(s)/7);
L=length(pulse);
ke=1;j=1;
for i=1:length(pulse)-1
    if((pulse(i+1)-pulse(i))<2)
        ke=ke+1;
    else
        max1(j)=find(s==max(s(pulse(i-ke+1:i))));
        j=j+1;
        ke=1;
    end
end
max1(j)=find(s==max(s(pulse(L-ke+1:L))));

%location to time
time=zeros(1,length(max1)-1);
for i=1:length(time)
    time(i)=(max1(i+1)-max1(i))*0.001+(t0+1)*0.001;
end
time;

%find outlier
odd=find(abs(time-mean(time))>0.3);
if length(odd)>1
    for i=1:length(odd)
        time_odd=time(odd);
        isB=ismember(time,time_odd);
        time_quodd = time(~isB);
        max1_odd=max1(odd);
        isB1=ismember(max1,max1_odd);
        max1_quodd = max1(~isB);
    end
else
    time_quodd=time;        
end

%average beating time and heart rate
T=mean(time_quodd);
rp=60/T;
%1st derivation
ds=diff(s);
dt=t(1:end-1);
%2nd derivation
dds=diff(s,2);
ddt=t(1:end-2);
%plot
figure (2)
subplot(3,1,1)
plot(t,s);xlabel('t/s');ylabel('origin');
subplot(3,1,2)
plot(dt,ds);xlabel('t/s');ylabel('1st derivation');
subplot(3,1,3)
plot(ddt,dds);xlabel('t/s');ylabel('2nd derivation');


figure(3)
plot(t,s);xlabel('t/s');ylabel('Intensity (a.u.)');
hold on
%find starting point of every cardiac cycle
l=length(max1);
loc1=ones(l);
for i=1:length(max1)
    yp(i)=s(max1(i));
    xp(i)=find(s==yp(i));
    tp(i)=0.001*(max1(i)-1)+(t0+1)*0.001; 
    if max1(i)-250<1
        continue
    else
        yd_1(i)=min(s(max1(i)-250:max1(i)));
    end
    if yd_1(i)==s(max1(i)-250)        
        xd_1(i)=find(s==yd_1(i));
        kk=0;
        for ei=xd_1(i):xp(i)
            if ds(ei)<0
                kk=kk+1;
            else
                kk=kk;
            end
        end
        if kk>10
            y_1=s(xd_1(i):xp(i));
            y_2=yp(i)-y_1;
            [pk_1,loc_1]=findpeaks(y_2);
            loc1(i)=max(loc_1);
            xd(i)=xd_1(i)+loc1(i);
        else
            xd(i)=xd_1(i);
        end
        yd(i)=s(xd(i));          
        td(i)=0.001*(xd(i)-1)+(t0+1)*0.001;
    else
        xd_1(i)=find(s==yd_1(i));
        kk=0;
        for ei=xd_1(i):xp(i)
            if ds(ei)<0
                kk=kk+1;
            else
                kk=kk;
            end
        end
        if kk>10
            y_1=s(xd_1(i):xp(i));
            y_2=yp(i)-y_1;
            [pk_1,loc_1]=findpeaks(y_2);
            loc1(i)=max(loc_1);
            xd(i)=xd_1(i)+loc1(i);
        else
            xd(i)=xd_1(i);
        end
        yd(i)=s(xd(i));          
        td(i)=0.001*(xd(i)-1)+(t0+1)*0.001;
    end
end

%find base for every cardiac cycle
for i=1:length(max1)-3
    if xd(i)==0||xd(i+1)==0
        continue
    else
        base_k(i)=(yd(i+1)-yd(i))/(xd(i+1)-xd(i));
        base_b(i)=yd(i)-base_k(i)*xd(i);
    	h0(i)=yp(i)-base_k(i)*xp(i)-base_b(i);
        y_min(i)=min(s(xd(i):xd(i+1)));
        x_min(i)=find(s==y_min(i));
        h1(i)=y_min(i)-base_k(i)*x_min(i)-base_b(i);
        if h1(i)/h0(i)<-0.5
            continue
        else
            xp_f(i)=xp(i);
            xd_f(i)=xd(i);
            xe_f(i)=xd(i+1);
            yp_f(i)=yp(i);
            tp_f(i)=tp(i);
            yd_f(i)=yd(i);
            td_f(i)=td(i);
            ye_f(i)=yd(i+1);
            te_f(i)=td(i+1);
        end
    end
     if td_f(i)==0||te_f(i)==0
        continue
    else
        plot(tp_f(i),yp_f(i),'rs');
        plot(td_f(i),yd_f(i),'bo');
        plot(te_f(i),ye_f(i),'rx');
    end
end
%delet abnormal data
tp_f = tp_f(:,any(tp_f,1));
td_f = td_f(:,any(td_f,1));
te_f = te_f(:,any(te_f,1));
yp_f = yp_f(:,any(yp_f,1));
yd_f = yd_f(:,any(yd_f,1));
ye_f = ye_f(:,any(ye_f,1));
xp_f = xp_f(:,any(xp_f,1));
xd_f = xd_f(:,any(xd_f,1));
xe_f = xe_f(:,any(xe_f,1));


%segmentation
l=length(te_f);
Sa=zeros(1024,2*l);
pulse_time_min=zeros(1,l);pulse_time_max=zeros(1,l);H_P=zeros(1,l);H_T=zeros(1,l);H_Dic=zeros(1,l);T_P_ratio=zeros(1,l);Dic_P_ratio=zeros(1,l);up_t=zeros(1,l);down_t=zeros(1,l);total_t=zeros(1,l);RWTT=zeros(1,l);PPT=zeros(1,l);LVET=zeros(1,l);up_10=zeros(1,l);down_10=zeros(1,l);up_25=zeros(1,l);down_25=zeros(1,l);up_33=zeros(1,l);down_33=zeros(1,l);up_50=zeros(1,l);down_50=zeros(1,l);down_up_ratio_10=zeros(1,l);down_up_ratio_25=zeros(1,l);down_up_ratio_33=zeros(1,l);down_up_ratio_50=zeros(1,l);width_10=zeros(1,l);width_30=zeros(1,l);width_50=zeros(1,l);width_70=zeros(1,l);up_k=zeros(1,l);up_max_k=zeros(1,l);up_max_t=zeros(1,l);DPT=zeros(1,l);
for i=1:l
    ke(i)=xe_f(i)-xd_f(i);
    ke(i)=fix(ke(i));
    if ke(i)>1024||ke(i)<600
        continue
    else
        base_k_f(i)=(ye_f(i)-yd_f(i))/(xe_f(i)-xd_f(i));
        base_b_f(i)=yd_f(i)-base_k_f(i)*xd_f(i);
        t=td_f(i):0.001:te_f(i);
        x=xd_f(i):xe_f(i);
        y_b=base_k_f(i)*x+base_b_f(i);
        y=s(xd_f(i):xe_f(i));
        y_db=y-y_b';            
        %y=(y-min(y))/(max(y)-min(y));
        Sa(1:ke(i)+1,2*i-1)=t(1:ke(i)+1)-T0;
        Sa(1:ke(i)+1,2*i)=y_db(1:ke(i)+1);
        pulse_time_min(i)=min(t(1:ke(i)+1)-T0);
        pulse_time_max(i)=max(t(1:ke(i)+1)-T0); 
        %extract feature
        %up wave
        up_y_db=y_db(1:xp_f(i)-xd_f(i));
        %dowmn wave
        down_y_db=y_db(xp_f(i)-xd_f(i):xe_f(i)-xd_f(i));
        %1st derivation
        d_y_db=diff(y_db);
        %2nd derivation
        dd_y_db=diff(y_db,2);        
              
        %Up_t
        up_t(i)=xp_f(i)-xd_f(i);
        %Down_t
        down_t(i)=xe_f(i)-xp_f(i);
        %total time
        total_t(i)=xe_f(i)-xd_f(i);
        %H_P
        H_P(i)=max(y_db);
        %up_k
        up_k(i)=H_P(i)/up_t(i);
        %up_max_k
        up_max_k(i)=max(d_y_db);
        %yp_max_t
        up_max_t(i)=find(up_max_k(i)==d_y_db);
        %DPT
        DPT(i)=up_t(i)-up_max_t(i);
        %H_10
        h_10=0.1*H_P(i);
        up_10(i)=length(find(up_y_db>h_10));
        %D_10
        down_10(i)=0;
        above_10=find(down_y_db>h_10);
        l_above=length(above_10);
        for j=1:l_above-1
            if above_10(j+1)-above_10(j)<2
                down_10(i)=down_10(i)+1;
            else
                down_10(i)=down_10(i);
                continue
            end
        end
        %H_25
        h_25=0.25*H_P(i);
        up_25(i)=length(find(up_y_db>h_25));
        %D_25
        down_25(i)=0;
        above_25=find(down_y_db>h_25);
        l_above=length(above_25);
        for j=1:l_above-1
            if above_25(j+1)-above_25(j)<2
                down_25(i)=down_25(i)+1;
            else
                down_25(i)=down_25(i);
                continue
            end
        end
        %H_33
        h_33=0.33*H_P(i);
        up_33(i)=length(find(up_y_db>h_33));
        %D_33
        down_33(i)=0;
        above_33=find(down_y_db>h_33);
        l_above=length(above_33);
        for j=1:l_above-1
            if above_33(j+1)-above_33(j)<2
                down_33(i)=down_33(i)+1;
            else
                down_33(i)=down_33(i);
                continue
            end
        end
        %Up_50
        h_50=0.5*H_P(i);
        up_50(i)=length(find(up_y_db>h_50));
        %Down_50
        down_50(i)=0;
        above_50=find(down_y_db>h_50);
        l_above=length(above_50);
        for j=1:l_above-1
            if above_50(j+1)-above_50(j)<2
                down_50(i)=down_50(i)+1;
            else
                down_50(i)=down_50(i);
                continue
            end
        end
        %Width_10
        width_10(i)=up_10(i)+down_10(i);                   
        %Width_130
        h_30=0.3*H_P(i);
        width_30(i)=length(find(y_db>h_30));
        %Width_50
        h_50=0.5*H_P(i);
        width_50(i)=length(find(y_db>h_50));
        %Width_70
        h_70=0.7*H_P(i);
        width_70(i)=length(find(y_db>h_70));
        %down_up_ratio_10 down_up_ratio_25 down_up_ratio_33  down_up_ratio_50
        down_up_ratio_10(i)=down_10(i)/up_10(i);
        down_up_ratio_25(i)=down_25(i)/up_25(i);
        down_up_ratio_33(i)=down_33(i)/up_33(i);
        down_up_ratio_50(i)=down_50(i)/up_50(i);
        %Reflected notch
        max_d_y_db(i)=max(d_y_db(up_t(i)+50:up_t(i)+245));
        if max_d_y_db(i)>0
            %loc_max_d_y_db(i)=find(d_y_db==max_d_y_db(i));
            for j=up_t(i)+50:up_t(i)+250
                if d_y_db(j)<0
                    continue
                else
                    width_W(i)=j;
                    break
                end
            end
            for j_j=width_W(i)+20:up_t(i)+300
                if d_y_db(j_j)>0
                    continue
                else
                    width_T(i)=j_j;
                    break
                end
            end
        else
            max_dd_y_db(i)=max(dd_y_db(up_t(i)+50:up_t(i)+250));
            dd_0(i)=find(dd_y_db==max_dd_y_db(i));
            max_d_y_db(i)=max(d_y_db(dd_0(i):dd_0(i)+100));
            width_T(i)=find(d_y_db==max_d_y_db(i));
        end
        %Dicrotic notch
        hhh=size(d_y_db);
        kkk=min(width_T(i)+355,hhh);
        max_d_y_db(i)=max(d_y_db(width_T(i)+50:kkk));
        if max_d_y_db(i)>0
            for jj=width_T(i)+50:width_T(i)+355
                if d_y_db(jj)<0
                    continue
                else
                    width_N(i)=jj;
                    break
                end
            end
        else
            continue
        end

        %Dicrotic wave
       kkk=min(width_N(i)+255,hhh);
       min_d_y_db(i)=min(d_y_db(width_N(i)+50:kkk));
       if min_d_y_db(i)<0
            for j=width_N(i)+50:width_N(i)+255
                if d_y_db(j)>0
                    continue
                else
                    width_Dic(i)=j;
                    break
                end
            end
       else
           continue
       end
    end
    %origin
    xD_f(i)=xd_f(i);
    tD_f(i)=0.001*xD_f(i)+t0*0.001;
    yD_f(i)=y_db(1);
    %P
    xP_f(i)=xp_f(i);
    tP_f(i)=0.001*xP_f(i)+t0*0.001;
    yP_f(i)=y_db(xP_f(i)-xD_f(i)+1);
    %T
    xT_f(i)=xd_f(i)+width_T(i)-1;
    tT_f(i)=0.001*xT_f(i)+t0*0.001;
    yT_f(i)=y_db(xT_f(i)-xD_f(i)+1);
    %V
    xN_f(i)=xd_f(i)+width_N(i)-1;
    tN_f(i)=0.001*xN_f(i)+t0*0.001;
    yN_f(i)=y_db(xN_f(i)-xD_f(i)+1);
    %D
    xDic_f(i)=xd_f(i)+width_Dic(i)-1;
    tDic_f(i)=0.001*xDic_f(i)+t0*0.001;
    yDic_f(i)=y_db(xDic_f(i)-xD_f(i)+1);
    
    %RWTT
    RWTT(i)=xT_f(i)-xP_f(i);
    %PPT
    PPT(i)=xDic_f(i)-xP_f(i);
    %LVET
    LVET(i)=xN_f(i)-xD_f(i);
    %T
    H_T(i)=yT_f(i);
    %T
    T_P_ratio(i)=H_T(i)/H_P(i);
    %H_D
    H_Dic(i)=yDic_f(i);
    %H_D/H_P
    Dic_P_ratio(i)=H_Dic(i)/H_P(i);
    
    if (mod(i,10)==0)
        figure(3+i)
        subplot(3,1,1)
        plot(t,y_db);xlabel('t/s');ylabel('Intensity (a.u.)');
        hold on
        plot(tD_f(i),yD_f(i),'r*');
        plot(tP_f(i),yP_f(i),'b+');
        plot(tT_f(i),yT_f(i),'rs');
        plot(tN_f(i),yN_f(i),'b*');
        plot(tDic_f(i),yDic_f(i),'rx');
        subplot(3,1,2)
        plot(td_f(i):0.001:te_f(i)-0.001,d_y_db);xlabel('t/s');ylabel('Intensity (a.u.)');
        subplot(3,1,3)
        plot(td_f(i):0.001:te_f(i)-0.002,dd_y_db);xlabel('t/s');ylabel('Intensity (a.u.)');
    end
end

%save
Sa = Sa(:,any(Sa,1));
dlmwrite('ls_left_1/bp_pulse/a.csv',Sa,'delimiter',',','precision',15)
% save feature
feature=[pulse_time_min;pulse_time_max;H_P;H_T;H_Dic;T_P_ratio;Dic_P_ratio;up_t;down_t;total_t;RWTT;PPT;LVET;up_10;down_10;up_25;down_25;up_33;down_33;up_50;down_50;down_up_ratio_10;down_up_ratio_25;down_up_ratio_33;down_up_ratio_50;width_10;width_30;width_50;width_70;up_k;up_max_k;up_max_t;DPT];
feature = feature(:,all(feature,1));
feature=feature';
label=["pulse_time_min","pulse_time_max","H_P","H_T","H_Dic","T_P_ratio","Dic_P_ratio","up_t","down_t","total_t","RWTT","PPT","LVET","up_10","down_10","up_25","down_25","up_33","down_33","up_50","down_50","down_up_ratio_10","down_up_ratio_25","down_up_ratio_33","down_up_ratio_50","width_10","width_30","width_50","width_70","up_k","up_max_k","up_max_t","DPT"];
label_feature=[label;feature];
dlmwrite('ls_left_1/bp_pulse/f_feature.csv',feature,'delimiter',',','precision',15)
