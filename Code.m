clc; clear all;
% Δημιουργώ το γραμμικοποιημένο αρχικό σύστημα στον χώρο κατάστασης
m=4; J=0.05; r=0.3; g=9.81; c=0.07;
A=[ 0 0 0 1 0 0
	0 0 0 0 1 0
	0 0 0 0 0 1
	0 0 -g -c/m 0 0
	0 0 0 0 -c/m 0
	0 0 0 0 0 0];	
B=[ 0 0; 0 0; 0 0 ;1/m 0; 0 1/m; r/J 0];
C=[ 1 0 0 0 0 0
	0 1 0 0 0 0];	
D=[0 0 ; 0 0];
sys=ss(A,B,C,D);
%Οι πίνακες ελεγξιμότητας και παρατηρησιμότητας στο σύνεχες
e=rank(ctrb(sys));
if e==size(sys.A,1), fprintf('Continuous System is Controllable \n'), end
o=rank(obsv(sys));
if o==size(sys.A,1), fprintf('Continuous System is Observable \n'), end

eig(sys); %Βλέπουμε ότι οι ιδιοτιμές του συστήματος δεν ικανοποιούν το κριτήριο ευστάθειας
step(sys); % Βημαρική απόκριση του συνεχούς συστήματος

%Υπολογισμός συναρτήσεων μεταφοράς 
[nc1,dc1]=ss2tf(A,B,C,D,1);
[nc2,dc2]=ss2tf(A,B,C,D,2);
tfc11=tf(nc1(1,:),dc1);
tfc22=tf(nc2(2,:),dc2);
%Διαγράμματα Bode για την εύρεση συχνότητας εύρους ζώνης 
bode(tfc11); bode(tfc22);

%Μετασχηματίζουμε το συνεχές σύστημα σε διακριτό
Ts=0.1;
sysd=c2d(sys,Ts,'zoh');
eig(sysd); % Οι ιδιοτιμές του διακριτού
%Οι πίνακες ελεγξιμότητας και παρατηρησιμότητας του διακριτού
ed=rank(ctrb(sysd));
if ed==size(sysd.A,1), fprintf('Discrete System is Controllable \n'), end
od=rank(obsv(sysd));
if od==size(sysd.A,1), fprintf('Discrete System is Observable \n'), end

%Υλοποίηση του ελεγκτή Linear Quadratic Regulator
Q=eye(6);R=eye(2); %Επιλογή των πινάκων Q,R ώστε να ισχύουν οι ιδιότητες που αναφέραμε 
[K,~,E]=dlqr(sysd.A,sysd.B,Q,R);% Linear Quadratic Regulator για διακριτό σύστημα
%Δημιουργία ευσταθούς συστήματος
sysnew=ss(sysd.A-sysd.B*K,zeros(6,2),sysd.C-sysd.D*K,zeros(2,2),Ts);
eig(sysnew); % Ιδιοτιμές sysnew
% Διαγράμματα που δείχνουν ότι το σύστημα συγκλίνει στην κατάσταση ισορροπίας
x01=[17 78 pi/2 0 0 0];
x02=[87 9 0 5 0 1];
initial(sysnew,x01);
initial(sysnew,x02); 
%Το σύστημα μετά την τοποθέτηση της reference εισόδου 
sysn=ss(sysd.A-sysd.B*K,sysd.B,sysd.C,zeros(2,2),Ts);
%Υπολογιμός συναρτήσεων μεταφοράς για τα Dc Gain
[n1,d1]=ss2tf(sysn.A,sysn.B,sysn.C,sysn.D,1);
[n2,d2]=ss2tf(sysn.A,sysn.B,sysn.C,sysn.D,2);
tf11=tf(n1(1,:),d1,Ts);
tf22=tf(n2(2,:),d2,Ts);
%Κερδη προαντιστάθμισης
K1=1/dcgain(tf11);
K2=1/dcgain(tf22);
Pro=[K1 0 ; 0 K2]; % Προαντισταθμιτής
%Τελικό σύστημα
sysnn=ss(sysd.A-sysd.B*K,sysd.B*Pro,sysd.C,zeros(2,2),Ts);

f=0.01; t=0:Ts:180;%Χρόνος και συχνότητα του σήματος reference
%Κυκλική τροχιά 
xideal=2*sin(2*pi*f*t)+3;
yideal=2*sin(2*pi*f*t+pi/2)+3;
U=[xideal;yideal];
yr=lsim(sysnn,U,t);
plot(yr(:,1),yr(:,2),'r',xideal,yideal,'b--')
grid minor;
axis([-1 6 0 6]);
title('Cycle Trajectory');
legend('System','Reference Trajectory');

pause();
%Υλοποίηση τροχίας σε σχήμα απείρου 
xinf=5*cos(2*pi*f*t)+7;
yinf=2.5*sin(4*pi*f*t)+7;
Uinf=[xinf; yinf];
yrr=lsim(sysnn,Uinf,t);
plot(yrr(:,1),yrr(:,2),'r',xinf,yinf,'b--')
grid minor;
axis([-2 13 0 10]);
title('Infinity Trajectory');
legend('System','Reference Trajectory');

