function [Pt,Pt_central_cell,SINR_user,users,users_central_cell,throughput_total,throughput_central_cell] = simul_5gusers_distri(input_matrix,H_final)

Mt = input_matrix(1);
Mr = input_matrix(2);
tiers = input_matrix(3);
N_subcarriers = input_matrix(4);
rrm_type = input_matrix(5);
SNR_req = input_matrix(6);
PRB_per_MS = input_matrix(7);
AMC = input_matrix(8);

index_h = randsrc(1,1,1:1e+7);
%%%%  UMa Scenario %%%%%

h_BS = 25;
h_UT = 1.5;
fc = 28;
R = 500;
PL_max = 160;
Pt_max_per_PRB = 0.0240;
Pt_BS_max = 10^(0.5);

No = -174 + 10*log10(100*(1e+6));
No = 10^(No/10)/1000;
SNR_req  = 10^(SNR_req/10);
G_m = 4;
G_b = 18;

sectors = 3;
xo = 0; yo = 0; zo = 25;                                                              
number_of_cells = 1;                                 
for j = 1:1:tiers
    number_of_cells = number_of_cells + j*6;
end

Rc = R*sqrt((3*number_of_cells*sqrt(3))/(2*pi));     
Rc_cell = R*sqrt((3*sqrt(3))/(2*pi));
position = cell_array_design(xo,yo,R,tiers,0); 
position(:,3) = zo;

x_min = -(1 + 2*tiers)*R; 
x_max =  (1 + 2*tiers)*R;                           
y_min = -(1 + 2*tiers)*R*(sqrt(3)/2);               
y_max =  (1 + 2*tiers)*R*(sqrt(3)/2); 
total_sectors = number_of_cells*sectors;

flag = 0;

%%%%%%%%%%%%%%%% initialization of matrices %%%%%%%%%%%%%%%%%%%%

subcarriers_of_BSs = zeros(number_of_cells,N_subcarriers);
subcarriers_of_users = zeros(1,N_subcarriers);
users  = zeros(1,total_sectors);
users_BS  = zeros(1,number_of_cells);
total_loss_user = zeros(1,total_sectors);
total_loss_user_own_sector = zeros(1,1);
sector_user = zeros(1,1);
base_station_user = zeros(1,1);
sector_of_BS_user = zeros(1,1);
t_user = zeros(1,Mt,N_subcarriers);
order_of_subcarriers_of_users = zeros(1,1);
H = zeros(1,total_sectors,Mr,Mt,N_subcarriers);
P_users = zeros(1,N_subcarriers);
modul_type = zeros(1,N_subcarriers);
Pt_BS = zeros(size(position,1),1);
Pt_central_cell = 0;

j  = 0;


while(flag==0)
   users_BS_previous = users_BS;
   [x_pos,y_pos,z_pos,PL,base_station,sector,angle]  = users_distribute(fc,x_min,x_max,y_min,y_max,Rc,position,PL_max,h_UT);
   users_BS(1,base_station) = users_BS(1,base_station) + 1;
   j = j + 1;
   reject = 0;
   [A,theta_BS] = attenuation_omni(x_pos,y_pos,position,number_of_cells);
   temp_loss = zeros(1,total_sectors);
   temp_loss(1,:) = PL(1,:) + A(1,:) - G_m - G_b;
   users(1,3*(base_station-1) + sector) = users(1,3*(base_station-1) + sector) + 1; 
   sector_user(sum(users),1) = 3*(base_station-1) + sector;
   base_station_user(sum(users),1) = base_station;
   sector_of_BS_user(sum(users),1) = sector;
   total_loss_user(sum(users),:) = temp_loss;
   total_loss_user_own_sector(1,sum(users)) = total_loss_user(sum(users),sector_user(sum(users),1));
   d2D = sqrt((x_pos-position(base_station,1))^2 + (y_pos-position(base_station,2))^2);
   for n = 1:1:total_sectors
       for PRB_index = 1:1:N_subcarriers
           index_h = index_h + 1;
           if(index_h==1e+7)
              index_h = randsrc(1,1,1:1e+7);
           end
           H(sum(users),n,1:Mr,1:Mt,PRB_index) = H_final(index_h,1:Mr,1:Mt);
       end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   final_temp_index = zeros(1,PRB_per_MS);
   subcarriers_of_BSs_previous = subcarriers_of_BSs;
   t_user_previous = t_user;
   temp_P = zeros(1,N_subcarriers);
   temp_modul_type = zeros(1,N_subcarriers);
   reject_flag = 0;
   pool_of_subcarriers = find(subcarriers_of_BSs(base_station_user(sum(users),1),1:N_subcarriers)==0); 
   if(length(pool_of_subcarriers)<PRB_per_MS)
       flag = 1;
       reject_flag = 1;
   else
       for sub_index = 1:1:PRB_per_MS
            temp_gain = 0;
            pool_of_subcarriers = find(subcarriers_of_BSs(base_station_user(sum(users),1),1:N_subcarriers)==0); 
            if(rrm_type==1)
                for channel_index = 1:1:length(pool_of_subcarriers)
                     current_h = zeros(Mr,Mt);
                     current_h(1:Mr,1:Mt) = H(sum(users),sector_user(sum(users),1),1:Mr,1:Mt,pool_of_subcarriers(channel_index));
                     [v,d] = eig(current_h'*current_h);
                     temp_t = zeros(Mt,1);
                     temp_t(1:Mt,1) = v(:,diag(d)==max(diag(d)));
                     if((temp_t'*(current_h'*current_h)*temp_t)>temp_gain)
                         temp_gain = temp_t'*(current_h'*current_h)*temp_t;
                         temp_index = pool_of_subcarriers(channel_index);
                         t_user(sum(users),1:Mt,temp_index) = temp_t(1:Mt,1);
                         final_temp_index(1,sub_index) = temp_index;
                     end
                end
                temp_P(1,temp_index) = SNR_req*No*10^(total_loss_user_own_sector(1,sum(users))/10);
                temp_P(1,temp_index) = temp_P(1,temp_index)/abs(temp_gain);
                temp_modul_type(1,temp_index) = 1;
                if(AMC==1)
                    if(temp_P(1,temp_index)<Pt_max_per_PRB)
                        SNR_req = 16.4;
                        temp_P(1,temp_index) = SNR_req*No*10^(total_loss_user_own_sector(1,sum(users))/10);
                        temp_P(1,temp_index) = temp_P(1,temp_index)/abs(temp_gain);
                        temp_modul_type(1,temp_index) = 4;
                        if(temp_P(1,temp_index)<Pt_max_per_PRB)
                            SNR_req = 22.7;
                            temp_P(1,temp_index) = SNR_req*No*10^(total_loss_user_own_sector(1,sum(users))/10);
                            temp_P(1,temp_index) = temp_P(1,temp_index)/abs(temp_gain);
                            temp_modul_type(1,temp_index) = 6;
                            if(temp_P(1,temp_index)>Pt_max_per_PRB)
                                SNR_req = 16.4;
                                temp_P(1,temp_index) = SNR_req*No*10^(total_loss_user_own_sector(1,sum(users))/10);
                                temp_P(1,temp_index) = temp_P(1,temp_index)/abs(temp_gain);
                                temp_modul_type(1,temp_index) = 4;
                            end
                        else
                            SNR_req = 9.6;
                            temp_P(1,temp_index) = SNR_req*No*10^(total_loss_user_own_sector(1,sum(users))/10);
                            temp_P(1,temp_index) = temp_P(1,temp_index)/abs(temp_gain);
                            temp_modul_type(1,temp_index) = 1;
                        end
                    end
                end
                subcarriers_of_BSs(base_station_user(sum(users),1),temp_index) = 1;
            else
                channel_index = randsrc(1,1,1:length(pool_of_subcarriers));
                temp_index = pool_of_subcarriers(channel_index);
                current_h = zeros(Mr,Mt);
                current_h(1:Mr,1:Mt) = H(sum(users),sector_user(sum(users),1),1:Mr,1:Mt,temp_index);
                [v,d] = eig(current_h'*current_h);
                temp_t = zeros(Mt,1);
                temp_t(1:Mt,1) = v(:,diag(d)==max(diag(d)));
                temp_gain = temp_t'*(current_h'*current_h)*temp_t;
                t_user(sum(users),1:Mt,temp_index) = temp_t(1:Mt,1);
                temp_P(1,temp_index) = SNR_req*No*10^(total_loss_user_own_sector(1,sum(users))/10);
                temp_P(1,temp_index) = temp_P(1,temp_index)/abs(temp_gain);
                final_temp_index(1,sub_index) = temp_index;
                subcarriers_of_BSs(base_station_user(sum(users),1),temp_index) = 1;
            end
            if(temp_P(1,temp_index)>Pt_max_per_PRB)
                reject_flag = 1;
            end
       end
   end
   if(reject_flag>0)
       reject = reject + 1;
       sector_user(sum(users)) = [];
       base_station_user(sum(users)) = [];
       sector_of_BS_user(sum(users)) = [];
       total_loss_user(sum(users),:) = [];
       total_loss_user_own_sector(sum(users)) = [];
       users(1,3*(base_station-1) + sector) = users(1,3*(base_station-1) + sector) - 1;  
       subcarriers_of_BSs  = subcarriers_of_BSs_previous;
       t_user = t_user_previous;
       users_BS = users_BS_previous;
   else
       order_of_subcarriers_of_users(sum(users),1:PRB_per_MS) = final_temp_index;
       subcarriers_of_users(sum(users),final_temp_index) = 1;
       P_users(sum(users),:) = temp_P;
       modul_type(sum(users),:) = temp_modul_type;
       if(base_station==1)
            Pt_central_cell = Pt_central_cell + real(sum(P_users(sum(users),:)));
       end
       Pt_BS(base_station,1) = Pt_BS(base_station,1) + real(sum(P_users(sum(users),:)));
       if(Pt_BS(base_station,1)>Pt_BS_max)
           flag = 1;
       end
   end
%    disp(j);
%    disp(Pt_BS.');
%    disp(sum(users));
end

Pt = real(sum(Pt_BS));
index_user = 1;
SINR_user = 0;

if(sum(users)>0)
    sub_index = 1;
    temp_index = find(subcarriers_of_users(index_user,:)==1);
    temp_index = temp_index(sub_index);
    co_channel_MSs = find(subcarriers_of_users(:,temp_index)==1);
    h = zeros(Mr,Mt);
    h(1:Mr,1:Mt) = H(index_user,sector_user(index_user,1),1:Mr,1:Mt,temp_index);
    temp_t(1:Mt,1) = t_user(index_user,1:Mt,temp_index);
    temp_gain = abs(temp_t'*(h'*h)*temp_t)^2*P_users(index_user,temp_index);
    temp_gain = temp_gain/10^(total_loss_user_own_sector(1,index_user)/10);
    interference = No*abs(temp_t'*(h'*h)*temp_t);
    temp_t_user = zeros(Mt,1);
    for index_MS = 1:1:length(co_channel_MSs)
        if(co_channel_MSs(index_MS)~=index_user)
            temp_t_user(1:Mt,1) = t_user(co_channel_MSs(index_MS),1:Mt,temp_index);
            h_int(1:Mr,1:Mt) = H(co_channel_MSs(index_MS),sector_user(index_user,1),1:Mr,1:Mt,temp_index);
            interference = interference + abs(temp_t'*(h'*h_int)*temp_t_user)^2*(P_users(co_channel_MSs(index_MS),temp_index)/10^(total_loss_user(co_channel_MSs(index_MS),sector_user(index_user,1))/10));
        end
    end

    temp_SINR_user = temp_gain/interference;
    %     if(temp_SINR_user>SINR_user)
            SINR_user = temp_SINR_user;
    %     end
    % disp([Pt Pt_central_cell]);
end


users_central_cell = sum(users(1,1:sectors));
users = sum(users);

subcarrier_spacing = 60000;
% Throughput_per_Cell = max_users_per_BS*PRB_per_MS*12*subcarrier_spacing;
% Total_Throughput = Throughput_per_Cell*7;
throughput_total = 0;
throughput_central_cell = 0;

if(AMC==1)
    for index_users = 1:1:users;
        for index_PRB = 1:1:N_subcarriers
            if(modul_type(index_users,index_PRB)~=0)
                throughput_total = throughput_total + 2^(modul_type(index_users,index_PRB))*subcarrier_spacing*12*PRB_per_MS;
                if(base_station_user(index_users,1)==1)
                   throughput_central_cell = throughput_central_cell + 2^(modul_type(index_users,index_PRB))*subcarrier_spacing*12*PRB_per_MS;
                end
            end
        end
    end
end

% pause;