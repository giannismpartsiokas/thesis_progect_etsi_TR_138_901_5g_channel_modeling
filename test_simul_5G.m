Mt = 2;
Mr = 2;
tiers = 1;
PRB_max = 132;
% subcarriers_per_PRB = 12;
N_subcarriers = PRB_max;
rrm_type = 0;
SNR_req = 9.6;
max_users_per_BS = 17;
PRB_per_MS = 1;
Pt_BS_max = 10^(0.5);
load('H_generated.mat');

input_matrix = [Mt Mr tiers N_subcarriers rrm_type SNR_req max_users_per_BS PRB_per_MS];
    
MC_repeats = 1e+2;
Pt = zeros(1,MC_repeats);
Pt_central_cell = zeros(1,MC_repeats);
SINR_user = zeros(1,MC_repeats);
reject = zeros(1,MC_repeats);
outage = 0;

for index = 1:1:MC_repeats
 index_rnd = randsrc(1,1,1:1e+7);
 disp(index);
 [Pt(1,index),Pt_central_cell(1,index),SINR_user(1,index),reject(1,index)] = simul_5g(input_matrix,H_final(index_rnd,:,:));
 disp(mean(Pt_central_cell(1,1:index)));
 disp([min(SINR_user(1,1:index)) mean(SINR_user(1,1:index)) max(SINR_user(1,1:index))]);
 disp(mean(reject(1,1:index)));
 if(Pt_central_cell(1,index)>Pt_BS_max)
     outage = outage + 1
 end
end


P_out = outage/index;

c = Pt_central_cell;
color = 'r';
[p,w,bins] = pdfplot(c,color,flag);

figure;

c = SINR_user;
color = 'b';
[p,w,bins] = pdfplot(c,color,flag);

figure;
ecdf(c);
save results_Mt_2_Mr_2_MSs_per_BS_17_PRB_per_MS_1;