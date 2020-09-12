function h_final = channel_3gpp_5G(Mt,Mr,fc,h_UT,d2D)

N_clusters = 20;
N_rays_per_cluster = 20;
h = zeros(Mr,Mt,N_clusters);
h_final = zeros(Mr,Mt);

offset_angles = zeros(1,N_rays_per_cluster);
offset_angles(1,1) = 0.0447; 
offset_angles(1,2) = -0.0447;
offset_angles(1,3) = 0.1413;
offset_angles(1,4) = -0.1413;
offset_angles(1,5) = 0.2492;
offset_angles(1,6) = -0.2492;
offset_angles(1,7) = 0.3715;
offset_angles(1,8) = -0.3715;
offset_angles(1,9) = 0.5129;
offset_angles(1,10) = -0.5129;
offset_angles(1,11) = 0.6797;
offset_angles(1,12) = -0.6797;
offset_angles(1,13) = 0.8844;
offset_angles(1,14) = -0.8844;
offset_angles(1,15) = 1.1481;
offset_angles(1,16) = -1.1481;
offset_angles(1,17) = 1.5195;
offset_angles(1,18) = -1.5195;
offset_angles(1,19) = 2.1551;
offset_angles(1,20) = -2.1551;

mean_DS_LOS = -6.955 - 0.0963*log10(fc); 
std_DS_LOS = 0.66; 
mean_DS_NLOS = -6.28 - 0.204*log10(fc);
std_DS_NLOS = 0.39; 

DS = 10^(mean_DS_NLOS + std_DS_NLOS*randn);

mean_ASD_LOS  = 1.06 + 0.1114*log10(fc); 
std_ASD_LOS  = 0.28; 
mean_ASD_NLOS = 1.5 - 0.1144*log10(fc);
std_ASD_NLOS = 0.28;

ASD = 10^(mean_ASD_NLOS+std_ASD_NLOS*randn);
ASD = min(ASD,104);

mean_ASA_LOS = 1.81; 
std_ASA_LOS = 0.2;
mean_ASA_NLOS = 2.08 - 0.27*log10(fc);
std_ASA_NLOS = 0.11;

ASA = 10^(mean_ASA_NLOS + std_ASA_NLOS*randn);
ASA = min(ASA,104);

mean_ZSA_LOS  = 0.95;
std_ZSA_LOS = 0.16;
mean_ZSA_NLOS = -0.3236*log10(fc) + 1.512;
std_ZSA_NLOS = 0.16;

ZSA = 10^(mean_ZSA_NLOS + std_ZSA_NLOS*randn);
ZSA = min(ZSA,52); 

mean_ZSD_NLOS = max(-0.5, -2.1*(d2D/1000) -0.01*(h_UT - 1.5)+0.9);
std_ZSD_NLOS  = 0.49;

ZSD = 10^(mean_ZSD_NLOS + std_ZSD_NLOS*randn);
ZSD = min(ZSD,52);

a_fc = 0.208*log10(fc)- 0.782;
b_fc = 25;
c_fc = -0.13*log10(fc )+ 2.03;
e_fc = 7.66*log10(fc)  - 5.96;

m_offset_ZOD = e_fc - 10^(a_fc*log10(max(b_fc,d2D)) + c_fc -0.07*(h_UT-1.5));

C_ASA = 15;
C_ASD = 3;
C_ZSA = 7;

mlg_ZSD = mean_ZSD_NLOS;

r_t_LOS = 2.5;
r_t_NLOS = 2.3;
m_xpr = 7;
s_xpr = 3;
C_DS = max(0.25,6.5622-3.4084*log10(fc));

t_n = zeros(1,N_clusters);
P_n = zeros(1,N_clusters);
phi_AOA = zeros(N_clusters,N_rays_per_cluster);
phi_AOD = zeros(N_clusters,N_rays_per_cluster);
theta_ZOA = zeros(N_clusters,N_rays_per_cluster);
theta_ZOD = zeros(N_clusters,N_rays_per_cluster);

C_f = ones(1,N_clusters);
C_f(1,4) = 0.779;
C_f(1,5) = 0.860;
C_f(1,8) = 1.018;
C_f(1,10) = 1.090;
C_f(1,11) = 1.123;
C_f(1,12) = 1.146;
C_f(1,14) = 1.190;
C_f(1,15) = 1.211;
C_f(1,16) = 1.226;
C_f(1,19) = 1.273;
C_f(1,20) = 1.289;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_theta = ones(1,N_clusters);
C_theta(1,8) = 0.889;
C_theta(1,10) = 0.957;
C_theta(1,11) = 1.031;
C_theta(1,12) = 1.104;
C_theta(1,15) = 1.1088;
C_theta(1,19) = 1.184;
C_theta(1,20) = 1.178;

%%%% Generate Cluster Delays %%%%%%%%%%
for index_clusters = 1:1:N_clusters
    X_n = rand;
    t_n(1,index_clusters) = -r_t_NLOS*DS*log(X_n);
end

t_n = sort(t_n-min(t_n),'ascend');

%%%% Generate Cluster Delays %%%%%%%%%%
for index_clusters = 1:1:N_clusters
    z_n = 3*randn;
    P_n(1,index_clusters) = exp(-(t_n(1,index_clusters)*(r_t_NLOS-1))/r_t_NLOS/DS)*10^(-z_n/10);
end

P_n = P_n/sum(P_n);

%%%% Generate arrival angles %%%%
for index_clusters = 1:1:N_clusters
    temp_phi_AOA = 2*(ASA/1.4)*sqrt(-log(P_n(1,index_clusters)/max(P_n)));
    temp_phi_AOA = temp_phi_AOA/C_f(1,index_clusters);
    X_n = randsrc(1,1,[-1,1]);
    Y_n = randn*(ASA/7);
    temp_phi_AOA = X_n*temp_phi_AOA + Y_n;
%      + AOA_user(user)
    for index_subpaths = 1:1:N_rays_per_cluster
        phi_AOA(index_clusters,index_subpaths) = temp_phi_AOA + offset_angles(1,index_subpaths)*C_ASA;
    end
end

%%%% Generate departures angles %%%%s

for index_clusters = 1:1:N_clusters
    temp_phi_AOD = 2*(ASD/1.4)*sqrt(-log(P_n(1,index_clusters)/max(P_n)));
    temp_phi_AOD = temp_phi_AOD/C_f(1,index_clusters);
    X_n = randsrc(1,1,[-1,1]);
    Y_n = randn*(ASD/7);
    temp_phi_AOD = X_n*temp_phi_AOD + Y_n;
%      + AOD_user(user)
    for index_subpaths = 1:1:N_rays_per_cluster
        phi_AOD(index_clusters,index_subpaths) = temp_phi_AOD + offset_angles(1,index_subpaths)*C_ASD;
    end
end


%%%% Generate ZOA angles %%%%

for index_clusters = 1:1:N_clusters
    temp_theta_ZOA = -2*ZSA*log(P_n(1,index_clusters)/max(P_n));
    temp_theta_ZOA = temp_theta_ZOA/C_theta(1,index_clusters);
    X_n = randsrc(1,1,[-1,1]);
    Y_n = randn*(ZSA/7);
    temp_theta_ZOA = X_n*temp_theta_ZOA + Y_n;
%      + ZOA_user(user)
    for index_subpaths = 1:1:N_rays_per_cluster
        theta_ZOA(index_clusters,index_subpaths) = temp_theta_ZOA + offset_angles(1,index_subpaths)*C_ZSA;
    end
end

%%%% Generate ZOD angles %%%%

for index_clusters = 1:1:N_clusters
    temp_theta_ZOD = -2*ZSD*log(P_n(1,index_clusters)/max(P_n));
    temp_theta_ZOD = temp_theta_ZOD/C_theta(1,index_clusters);
    X_n = randsrc(1,1,[-1,1]);
    Y_n = randn*(ZSD/7);
    temp_theta_ZOD = X_n*temp_theta_ZOD + Y_n + m_offset_ZOD;
    theta_output = temp_theta_ZOD;
%      + ZOD_user(user)
    for index_subpaths = 1:1:N_rays_per_cluster
        theta_ZOD(index_clusters,index_subpaths) = temp_theta_ZOD + offset_angles(1,index_subpaths)*(3/8)*10^mlg_ZSD;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_theta_theta = zeros(N_clusters,N_rays_per_cluster);
F_theta_phi = zeros(N_clusters,N_rays_per_cluster);
F_phi_theta = zeros(N_clusters,N_rays_per_cluster);
F_phi_phi = zeros(N_clusters,N_rays_per_cluster);

for index_clusters = 1:1:N_clusters
    for index_subpaths = 1:1:N_rays_per_cluster
        F_theta_theta(index_clusters,index_subpaths) = -pi + rand*2*pi;
        F_theta_phi(index_clusters,index_subpaths) = -pi + rand*2*pi;
        F_phi_theta(index_clusters,index_subpaths) = -pi + rand*2*pi;
        F_phi_phi(index_clusters,index_subpaths) = -pi + rand*2*pi;
    end
end

angle_matrix = zeros(2,2);

for u = 1:1:Mr
    for s = 1:1:Mt
        for n = 3:1:N_clusters
            for m = 1:1:N_rays_per_cluster
               x = m_xpr + s_xpr*randn;
               k = 10^(x/10);
               angle_matrix(1,1) = exp(1i*F_theta_theta(n,m));
               angle_matrix(1,2) = sqrt(k^(-1))*exp(1i*F_theta_phi(n,m));
               angle_matrix(2,1) = sqrt(k^(-1))*exp(1i*F_phi_theta(n,m));
               angle_matrix(2,2) = exp(1i*F_theta_theta(n,m)); 
               Fr = ones(2,1);
               Ft = ones(2,1);
               h(u,s,n) = h(u,s,n) + sqrt(P_n(1,n)/N_rays_per_cluster)*Fr.'*angle_matrix*Ft;
            end
%             if(u==1)&&(s==1)
% %               disp(abs(h(u,s,n)));
% %             disp(sqrt(P_n(1,n)/N_rays_per_cluster));
%             disp(Fr.'*angle_matrix*Ft);
%             end
        end
    end
end

%%%%%%%%%%%%%% calculation of channel for the strongest clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for u = 1:1:Mr
    for s = 1:1:Mt
        for n = 1:1:2
            for i = 1:1:3 
                if(i==1)
                    m_set = [1 2 3 4 5 6 7 8 19 20];
                    delay_offset = 0;
                elseif(i==2)
                    m_set = [9 10 11 12 17 18];
                    delay_offset = 1.28*C_DS;
                else
                    m_set = [13 14 15 16];
                    delay_offset = 2.56*C_DS;
                end
                for m_index = 1:1:length(m_set)
                   m = m_set(m_index); 
                   x = m_xpr + s_xpr*randn;
                   k = 10^(x/10);
                   angle_matrix(1,1) = exp(1i*F_theta_theta(n,m));
                   angle_matrix(1,2) = sqrt(k^(-1))*exp(1i*F_theta_phi(n,m));
                   angle_matrix(2,1) = sqrt(k^(-1))*exp(1i*F_phi_theta(n,m));
                   angle_matrix(2,2) = exp(1i*F_theta_theta(n,m)); 
                   Fr = ones(2,1);
                   Ft = ones(2,1);
                   h(u,s,n) = h(u,s,n) + sqrt(P_n(1,n)/N_rays_per_cluster)*Fr.'*angle_matrix*Ft;
                end
%                disp(abs(h(u,s,n)));
            end
        end
    end
end

for u = 1:1:Mr
    for s = 1:1:Mt
        for n = 1:1:N_clusters
           h_final(u,s) = h_final(u,s) + h(u,s,n);
        end
    end
end
