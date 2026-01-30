clear;clc
rng(1);
%%
Naverage=1000; 
Dmin=0.001;Dmax=0.1;
D = Dmin * ((Dmax/Dmin).^(1/95)) .^(0:95);
input0_=-1:0.02:0;
A_input=0.25;
%%
tic
for j=1:length(input0_)
    input0=input0_(j);
       parfor d=1:length(D)
           [P_single_nor(j,d,:), P_single_nand(j,d,:), P_single_and(j,d,:),P_single_or(j,d,:),...
               P_single_xor(j,d,:), P_single_xnor(j,d,:),P_single_sr(j,d,:)] = trial_tot(D(d),Naverage,j,A_input,input0);
        end
        Pand_d_nor(j,:,:) = P_single_nor(j,:,:);
        Pand_d_nand(j,:,:) = P_single_nand(j,:,:);
        Pand_d_and(j,:,:) = P_single_and(j,:,:);
        Pand_d_or(j,:,:) = P_single_or(j,:,:);
        Pand_d_xor(j,:,:) = P_single_xor(j,:,:);
        Pand_d_xnor(j,:,:) = P_single_xnor(j,:,:);
        Pand_d_sr(j,:,:) = P_single_sr(j,:,:);
end
Pand_nor = mean(Pand_d_nor, 3);
Pand_nand = mean(Pand_d_nand, 3);
Pand_and = mean(Pand_d_and, 3);
Pand_or = mean(Pand_d_or, 3);
Pand_xor = mean(Pand_d_xor, 3);
Pand_xnor = mean(Pand_d_xnor, 3);
Pand_sr = mean(Pand_d_sr, 3);