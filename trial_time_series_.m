clear;clc
Naverage=100;seed=1;Nflag=99;Ainput=0.15;
D=0.03; 
input0=0.495;
rng(seed,"twister");

a=4.0;
c=0.5;
%%
step=0.01; Tinterval=1000; T0=100; Ttot=7*Tinterval; Nt=round(Ttot/step)+1; 
t=0:step:Ttot;

x=-1.0;
input1=0.0;
input2=0.0;
output(1:Nt)=0.0;
input(1:Nt)=0.0;
k=0;
%% time loop
for n  = 1:Nt
    % A=rand(2,Naverage);
    if mod(n-1,round(Tinterval/step))==0 %% 
        k=k+1;
        % input1=(round(A(1,Nflag))-0.5);
        % input2=(round(A(2,Nflag))-0.5);
            % swtich_time_point = t(n);
            % flag((x>0.0))=1;
            % flag((x<-0.0))=-1;
        % Ilogic=input1*Ainput*2+input2*Ainput*2;
        if k == 1 
            Ilogic=0.5;
        elseif k == 2
            Ilogic=0;
        elseif k == 3
            Ilogic=-0.5;
        elseif k == 4 
            Ilogic=0;
        elseif k == 5
            Ilogic=0.5;
        elseif k == 6 
            Ilogic=-0.5;
        elseif k == 7 
            Ilogic=0.5;
        end
      
        inputorig=(Ilogic+input0);
        Iinput=6.4 * (1-inputorig.^2) .* inputorig.^2 -0.8;
    end
    
    ran_noise = rand(2,Naverage);
    WhiteNoise=sqrt(-4.0*D*step.*log(ran_noise(1,Nflag))).*cos(2.0*pi.*ran_noise(2,Nflag));

    dx = x + step*(4*x-20*x.^3+Iinput)+WhiteNoise;    
    x = dx;
   
    output(n)=dx;
    input(n)=Ilogic;
    % output(:,n)=dx(1,:);
    % input(:,n)=Iinput(1,:);
end