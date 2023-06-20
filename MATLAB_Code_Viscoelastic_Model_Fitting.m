figure(1);
plot(bead_displacement,'o');
start=15;
time=[0;0.2;0.4;0.6;0.8;1;1.2;1.4;1.6;1.8;2.0;2.2;2.4;2.6;2.8;3.0];
d=[];
final_d=bead_displacement(500,1);
for i =1:5
    d = [d,bead_displacement(start+(i-1)*90:(i-1)*90+start+15,:)];
end
d1=d(:,1);
d2=d(:,2);
d3=d(:,3);
d4=d(:,4);
d5=d(:,5);
x1=[d1(16,1),d2(16,1),d3(16,1),d4(16,1),d5(16,1)];
y1=[d1(1,1),d2(1,1),d3(1,1),d4(1,1),d5(1,1),];
Deformation=x1-y1;
D1=d1-double(d1(1,1));
D2=d2-double(d2(1,1));
D3=d3-double(d3(1,1));
D4=d4-double(d4(1,1));
D5=d5-double(d5(1,1));
D =[D1,D2,D3,D4,D5];
r=[d1(1,1),d2(1,1),d3(1,1),d4(1,1),d5(1,1),final_d(1,1)];
residual=[];
 for i =2:6
    residual=[residual,r(1,i)-r(1,i-1)];
 end
figure(2);
subplot(2,2,1);
plot(bead_displacement);
title('Bead tracking');
xlim([0 500]);
subplot(2,2,2);
plot(d,"-- *");
title('Creep deformation'); 
subplot(2,2,3);
plot(Deformation,"-- *");
title('Deformation');
subplot(2,2,4);
plot (residual,"-- *");
title('Residual');
Fitting_Result=[];
for i =1:5
x = time;  
y = D(:,i);     
fun=@(a)(10*10^-9/(6*pi*1.4*10^-6)*(1/a(1)*(1-exp(-a(1)*x/a(2)))+1/a(3)*x)*10^6-y); 
a0 = [100 50 50];
A = lsqnonlin(fun,a0); 
Fitting_Result=[Fitting_Result;[A(1),A(2),A(3)]];
xx = 0:0.2:3;
yy = 10*10^-9/(6*pi*1.4*10^-6)*(1/Fitting_Result(i,1)*(1-exp(-Fitting_Result(i,1)*x/Fitting_Result(i,2)))+1/Fitting_Result(i,3)*x)*10^6;

figure(3)
subplot(1,5,i);
title('fitting result')
plot(xx,yy,'r-')
hold on
plot(x,y,'bp')
end
figure(4)
subplot(1,4,1);
plot(Fitting_Result(:,1),"-- *");
title('k1');
subplot(1,4,2);
plot(Fitting_Result(:,2),"-- *");
title('v1'); 
subplot(1,4,3);
plot(Fitting_Result(:,3),"-- *");
title('v2');
subplot(1,4,4);
plot(Fitting_Result(:,2)/Fitting_Result(:,1),"-- *");
title('time constant v1/k1');
figure('Units','centimeter','Position',[0 0 14 3.5]);
plot(bead_displacement);
xlim([0 500]);