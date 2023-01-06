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
        Sa(1:ke(i)+1,2*i-1)=t(1:ke(i)+1);
        Sa(1:ke(i)+1,2*i)=y_db(1:ke(i)+1);
    end
end

%save segmentation data
Sa = Sa(:,any(Sa,1));
dlmwrite('a.csv',Sa,'delimiter',',','precision',15)
