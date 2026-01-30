function [Psingle_nor, Psingle_nand, Psingle_and, Psingle_or, Psingle_xor, Psingle_xnor,Psingle_sr] = trial_tot(D,Naverage,seed,A_input,input0)
  % D=0.03;
% Naverage=10;A_input=0.10;
% seed=1;
rng(seed,"twister");
%% measure
swtich_time_point=0.0;
Psingle_nand(1:Naverage)=0.0;
Psingle_nor(1:Naverage)=0.0;
Psingle_or(1:Naverage)=0.0;
Psingle_and(1:Naverage)=0.0;
Psingle_xor(1:Naverage)=0.0;
Psingle_xnor(1:Naverage)=0.0;
Psingle_sr(1:Naverage)=0.0;
flag(1:Naverage)=0;
% per_switch_time(1:Naverage)=0; switch_time(1:Naverage)=0;
stat_num=0;
xm(1:Naverage)=0;
stat_success_NAND(1:Naverage)=0.0;
stat_success_NOR(1:Naverage)=0.0;
stat_success_OR(1:Naverage)=0.0;
stat_success_AND(1:Naverage)=0.0;
stat_success_XOR(1:Naverage)=0.0;
stat_success_XNOR(1:Naverage)=0.0;
stat_success_SR(1:Naverage)=0.0;
% bb(1:Naverage)=0.0;
%% 
I1 = input0 + A_input*(-2);
I2 = input0 + A_input*0;
I3 = input0 + A_input*2;
xm_values(1) = find_xm_only(I1, 1.5, 0.001);
xm_values(2) = find_xm_only(I2, 1.5, 0.001);
xm_values(3) = find_xm_only(I3, 1.5, 0.001);
%%
step=0.01; Tinterval=1000; T0=100; Ttot=20*Tinterval; Nt=round(Ttot/step)+1; 
t=0:step:Ttot;

x(1:Naverage)=-1.0;
prex=x;
input1(1:Naverage)=0.0;
input2(1:Naverage)=0.0;
% input0(1:Naverage)=0.52;

% output(1:Naverage,1:Nt)=0.0;
% input(1:Naverage,1:Nt)=0.0;
%% time loop
for n  = 1:Nt

    A=rand(2,Naverage);
    if mod(n-1,round(Tinterval/step))==0 %% 
         % input1=(round(A(1,:))*2-1)*0.1;
         % input2=(round(A(2,:))*2-1)*0.1;
         input1=(round(A(1,:))-0.5);
         input2=(round(A(2,:))-0.5);
         
        swtich_time_point = t(n);
        prex=x;

        Ilogic=input1+input2;
        xm(Ilogic==-1) = xm_values(1);
        xm(Ilogic==0) = xm_values(2);
        xm(Ilogic==1) = xm_values(3);
        flag((x>xm))=1;
        flag((x<xm))=-1;
        
        inputorig=(input1*A_input*2+input2*A_input*2+input0);
        c=0.8;
        Iinput = 4*2*c* (1-inputorig.^2) .* inputorig.^2 -c;
        % Iinput=a^2*inputorig.*(1-inputorig).*(1-a*inputorig.*(1-inputorig))-c;
    end

    ran_noise = rand(2,Naverage);
    WhiteNoise=sqrt(-4.0*D*step.*log(ran_noise(1,:))).*cos(2.0*pi.*ran_noise(2,:));
     % WhiteNoise=D*sqrt(-step*log(ran_noise(1,:))).*cos(2.0*pi.*ran_noise(2,:));

    dx = x + step*(4*x-20*x.^3+Iinput)+WhiteNoise;  
    x = dx;

    if t(n)-swtich_time_point>=T0
          stat_num=stat_num+1;
          stat_success_NOR(round(Ilogic)==-1 & x>xm)=stat_success_NOR(round(Ilogic)==-1& x>xm)+1;
          stat_success_NOR(round(Ilogic)==0 & x<xm)=stat_success_NOR(round(Ilogic)==0 & x<xm)+1;
          stat_success_NOR(round(Ilogic)==1 & x<xm)=stat_success_NOR(round(Ilogic)==1 & x<xm)+1;
           
         stat_success_AND(round(Ilogic)==-1 & x<xm)=stat_success_AND(round(Ilogic)==-1 & x<xm)+1;
         stat_success_AND(round(Ilogic)==0 & x<xm)=stat_success_AND(round(Ilogic)==0 & x<xm)+1;
         stat_success_AND(round(Ilogic)==1 & x>xm)=stat_success_AND(round(Ilogic)==1 & x>xm)+1;
    
          stat_success_OR(round(Ilogic)==-1 & x<xm)=stat_success_OR(round(Ilogic)==-1 & x<xm)+1;
          stat_success_OR(round(Ilogic)==0 & x>xm)=stat_success_OR(round(Ilogic)==0 & x>xm)+1;
          stat_success_OR(round(Ilogic)==1 & x>xm)=stat_success_OR(round(Ilogic)==1 & x>xm)+1;
          
         stat_success_NAND(round(Ilogic)==-1 & x>xm)=stat_success_NAND(round(Ilogic)==-1 & x>xm)+1;
         stat_success_NAND(round(Ilogic)==0 & x>xm)=stat_success_NAND(round(Ilogic)==0 & x>xm)+1;
         stat_success_NAND(round(Ilogic)==1 & x<xm)=stat_success_NAND(round(Ilogic)==1 & x<xm)+1;
    
         stat_success_XOR(round(Ilogic)==-1 & x<xm)=stat_success_XOR(round(Ilogic)==-1 & x<xm)+1;
         stat_success_XOR(round(Ilogic)==0 & x>xm)=stat_success_XOR(round(Ilogic)==0 & x>xm)+1;
         stat_success_XOR(round(Ilogic)==1 & x<xm)=stat_success_XOR(round(Ilogic)==1 & x<xm)+1;
              
         stat_success_XNOR(round(Ilogic)==-1 & x>xm)=stat_success_XNOR(round(Ilogic)==-1 & x>xm)+1;
         stat_success_XNOR(round(Ilogic)==0 & x<xm)=stat_success_XNOR(round(Ilogic)==0 & x<xm)+1;
         stat_success_XNOR(round(Ilogic)==1 & x>xm)=stat_success_XNOR(round(Ilogic)==1 & x>xm)+1;
     
         stat_success_SR(round(Ilogic)==-1 & x<xm)=stat_success_SR(round(Ilogic)==-1 & x<xm)+1;
         stat_success_SR(round(Ilogic)==0 & prex<xm & x<xm)=stat_success_SR(round(Ilogic)==0 & prex<xm & x<xm)+1;
         stat_success_SR(round(Ilogic)==0 & prex>xm & x>xm)=stat_success_SR(round(Ilogic)==0 & prex>xm & x>xm)+1;
         stat_success_SR(round(Ilogic)==1 & x>xm)=stat_success_SR(round(Ilogic)==1 & x>xm)+1;
    end  
end

for iii = 1:Naverage

     if stat_success_NOR(iii)==stat_num, Psingle_nor(iii)=1;end
     if stat_success_NAND(iii)==stat_num, Psingle_nand(iii)=1;end
     if stat_success_AND(iii)==stat_num, Psingle_and(iii)=1;end
     if stat_success_OR(iii)==stat_num, Psingle_or(iii)=1;end
     if stat_success_XOR(iii)==stat_num, Psingle_xor(iii)=1;end
     if stat_success_XNOR(iii)==stat_num, Psingle_xnor(iii)=1;end
     if stat_success_SR(iii)==stat_num, Psingle_sr(iii)=1;end
    
end