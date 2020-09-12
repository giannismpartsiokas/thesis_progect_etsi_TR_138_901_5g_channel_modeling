Mt = 2;
Mr = 2;
tiers = 2;
PRB_max = 132;
% subcarriers_per_PRB = 12;
N_subcarriers = PRB_max;
rrm_type = 1;
SNR_req = 16.4;
PRB_per_MS = 15;
Pt_BS_max = 10^(0.5);
AMC = 1;

subcarrier_spacing = 60000;
% Throughput_per_Cell = max_users_per_BS*PRB_per_MS*12*subcarrier_spacing;
% Total_Throughput = Throughput_per_Cell*7;

input_matrix = [Mt Mr tiers N_subcarriers rrm_type SNR_req PRB_per_MS AMC];
    
MC_repeats = 1e+2;
Pt = zeros(1,MC_repeats);
Pt_central_cell = zeros(1,MC_repeats);
users = zeros(1,MC_repeats);
users_central_cell = zeros(1,MC_repeats);
throughput_total = zeros(1,MC_repeats);
throughput_central_cell = zeros(1,MC_repeats);
SINR_user = zeros(1,MC_repeats);
reject = zeros(1,MC_repeats);
outage = 0;

load H_generated.mat;

for index = 1:1:MC_repeats
 disp(index);
 [Pt(1,index),Pt_central_cell(1,index),SINR_user(1,index),users(1,index),users_central_cell(1,index),throughput_total(1,index),throughput_central_cell(1,index)] = simul_5gusers_distri(input_matrix,H_final);
 disp([mean(Pt(1,1:index)) mean(Pt_central_cell(1,1:index))]);
 disp([min(SINR_user(1,1:index)) mean(SINR_user(1,1:index)) max(SINR_user(1,1:index))]);
 disp([mean(users(1,1:index)) mean(users_central_cell(1,1:index))]);
end

% P_out = outage/index;
% 
% c = Pt_central_cell;
% color = 'r';
% [p,w,bins] = pdfplot(c,color,flag);
% 
% figure;
% 
% ecdf(c);
% grid on;

save results_Mt_2_Mr_2_PRB_per_MS_15_ADAPTIVE;