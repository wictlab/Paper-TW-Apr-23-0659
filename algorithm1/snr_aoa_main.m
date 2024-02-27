close all;clear all;clc;
timer = 0;
Nt = 1;
Nr = 64;
Ny = 12;
Nx = 12;
L = 3;
d = 1;
H = zeros(Nr, Nt);
Gaoa=300;
search_area_aoa=(-90:180/Gaoa:90-180/Gaoa);
A_aoa=zeros(Nr,Gaoa);
grid_aoa=zeros(1,Gaoa);
for i=1:Gaoa
    grid_aoa(:,i)=sin(search_area_aoa(i)*pi/180);
    A_aoa(:,i)=exp(-1i*2*pi*(0:Nr-1)'*sin(search_area_aoa(i)*pi/180));
end
Gaod=300;
search_area_aod=(-90:180/Gaod:90-180/Gaod);
A_aod=zeros(Nt,Gaod);
grid_aod=zeros(1,Gaod);
for i=1:Gaoa
    grid_aod(:,i)=sin(search_area_aod(i)*pi/180);
end

SNR_list = -5:2.5:7.5;
sample_num = 110;
nmse_result = zeros(numel(SNR_list),1);
mean_error_mmse=zeros(numel(SNR_list),1);
mean_error_esprit=zeros(numel(SNR_list),1);
mean_error_omp=zeros(numel(SNR_list),1);
alpha_true = zeros(numel(SNR_list),sample_num,L);
phi_t_true = zeros(numel(SNR_list),sample_num,L);
phi_r_true = zeros(numel(SNR_list),sample_num,L);
L_result = zeros(numel(SNR_list),sample_num);
theta_result = zeros(numel(SNR_list),sample_num,2,10);
z_result = zeros(numel(SNR_list),sample_num,10);
for snr_ii = 1:numel(SNR_list)
    snr = SNR_list(snr_ii);
    noise = sqrt(10^(-snr/10)/2);
    for sample_ii = 1:sample_num
        H = zeros(Nr,Nt);
        H_vec = zeros(Nr*Nt,1);
        alpha = zeros(L,1);
        alpha(1) = exp(1i*2*pi*rand(1));
        alpha(2:L) = (normrnd(0, 0.1, L-1, 1) + 1i*normrnd(0, 0.1, L-1, 1)) / sqrt(2);
        while (find(abs(alpha)<0.01))
            alpha(2:L) = (normrnd(0, 0.1, L-1, 1) + 1i*normrnd(0, 0.1, L-1, 1)) / sqrt(2);
        end
        alpha = sort(alpha, 'descend');
        
        aod_taps=(rand(1,L)-0.5)*2*90;
        aoa_taps=(rand(1,L)-0.5)*2*90;
alphatrue=aoa_taps; 
%         phi_t = 2*rand(L,1)-1;%virtual AoD -1到1
%         phi_r = 2*rand(L,1)-1;%virtual AoA

        for l = 1:L
            at = exp(-1i*2*pi*[0:Nt-1]'*d*sin(aod_taps(l)*pi/180));
            ar = exp(-1i*2*pi*[0:Nr-1]'*d*sin(aoa_taps(l)*pi/180));
            H_vec = H_vec + alpha(l)*kron(ar,at);
            H = H + alpha(l)*(ar*at');
        end
        X = 1/sqrt(Nt)*exp(-1i*2*pi*rand(Nt,Nx));
        Y = (H*X + noise*(normrnd(0, 1, Nr, Nx) + 1i*normrnd(0, 1, Nr, Nx)));
        %Y_vec = (H_vec*X_vec + noise*(normrnd(0, 1, Nr*Nt, Nx) + 1i*normrnd(0, 1, Nr*Nt, Nx)));
        Y_vec=vec(Y);
        Rth = noise*sqrt(Ny*Nr);
        for i=1:Gaod
            A_aod(:,i)=exp(-1i*2*pi*(0:Nt-1)'*sin(search_area_aod(i)*pi/180));
            AM(:,i)=X'*A_aod(:,i);
        end
        dict=kron(AM,A_aoa);
        [theta_es,z_es,err]=proposed(dict,Y,X,Nx,Nt,Nr,Ny,Rth,Y_vec,Gaoa,Gaod,grid_aoa,grid_aod);

        mean_error_esprit(snr_ii) = mean_error_esprit(snr_ii)+ESPRIT(alphatrue,Ny,Nr,snr);

        A = exp(-1i*2*pi*[0:Nr-1]'*sin(alphatrue*pi/180));          % 流型矩阵
        noise1 =sqrt(1/2)*(randn(Nr, Ny)+1j*randn(Nr, Ny));
        Vj=sqrt((   10^(snr/10)   )/2);
        S=Vj*(randn(L, Ny)+1j*randn(L, Ny));
        Htrue=A*S;
        X1=Htrue+noise1;
        H_LS=X1;
        I=ones(Nr,1);
        I=diag(I);
        SNR = 10.^(snr/10);                                       % 转换
        noise_var_sqrt = sqrt(1./SNR);                               % 噪声方差 SNR=S/N，噪声功率为方差
        sigma_2 = abs(noise_var_sqrt).^2;
        H_mmse=H_LS*H_LS'*pinv(H_LS*H_LS'+I.*sigma_2)*H_LS;
        mean_error_mmse(snr_ii) = mean_error_mmse(snr_ii)+norm(H_mmse-Htrue,'fro') / norm(Htrue,'fro');
       
        H_es = zeros(Nr, Nt);
        at = zeros(Nt,1);
        ar = zeros(Nr,1);
        for l = 1:numel(z_es)
            at = exp(-1i*2*pi*[0:Nt-1]'*theta_es(1,l));
            ar = exp(-1i*2*pi*[0:Nr-1]'*theta_es(2,l));
            H_es = H_es + z_es(l)*ar*at';
        end
        h_OMP50=OMP(A_aoa,X1(:,1),L);
        H_OMP50=zeros(Nr,1);
        for h1=1:Gaoa
            if h_OMP50(h1)~=0
                H_OMP50=H_OMP50+h_OMP50(h1)*A_aoa(:,h1);
            end
        end
        mean_error_omp(snr_ii) = mean_error_omp(snr_ii)+norm(Htrue(:,1)-H_OMP50,'fro')^2/norm(Htrue(:,1),'fro')^2;
        %nmse_sample = sum(sum(abs(H-H_es).^2))/sum(sum(abs(H).^2));
        nmse_sample =norm(H-H_es,'fro')^2/norm(H,'fro')^2;
        disp(['snr=' num2str(SNR_list(snr_ii)) ' sample_ii=' num2str(sample_ii) ' nmse=' num2str(nmse_sample) ' err=' num2str(err)]);
        
       %% spectral efficiency
%         [U_perfectCSI,S,V_perfectCSI] = svd(H);
%         P_perfectCSI = V_perfectCSI(:,1:N_RF);
%         Q_perfectCSI = U_perfectCSI(:,1:N_RF);
%         R_perfectCSI = 0.5*noise*noise*(Q_perfectCSI'*Q_perfectCSI);
%         SE_perfectCSI = log2(det(eye(N_RF)+(1/N_RF)*inv(R_perfectCSI)*((Q_perfectCSI'*H*P_perfectCSI)*(Q_perfectCSI'*H*P_perfectCSI)')));
%         
%         [U_estCSI_SURE,S,V_estCSI_SURE] = svd(H_es);
%         P_estCSI_SURE = V_estCSI_SURE(:,1:N_RF);
%         Q_estCSI_SURE = U_estCSI_SURE(:,1:N_RF);
%         R_estCSI_SURE = 0.5*noise*noise*(Q_estCSI_SURE'*Q_estCSI_SURE);
%         SE_estCSI_SURE = log2(det(eye(N_RF)+(1/N_RF)*inv(R_estCSI_SURE)*((Q_perfectCSI'*H*P_estCSI_SURE)*(Q_perfectCSI'*H*P_estCSI_SURE)')));
        
%         alpha_true(snr_ii, sample_ii,:) = alpha;
%         phi_t_true(snr_ii, sample_ii,:) = phi_t;
%         phi_r_true(snr_ii, sample_ii,:) = phi_r;
        nmse_result(snr_ii) = nmse_result(snr_ii)+nmse_sample;
        L_sample = numel(z_es);
        L_result(snr_ii, sample_ii) = L_sample;
        theta_result(snr_ii, sample_ii, :, 1:L_sample) = theta_es(:,:);
        z_result(snr_ii, sample_ii, 1:L_sample) = z_es(:);
    end
end
nmse_result = nmse_result./sample_num;
mean_error_esprit=mean_error_esprit./sample_num;
mean_error_mmse=mean_error_mmse./sample_num;
mean_error_omp=mean_error_omp./sample_num;
figure;
semilogy(SNR_list,nmse_result,'-*','LineWidth',2);
 hold on;
semilogy(SNR_list,mean_error_esprit,'-*','LineWidth',2);
 hold on;
  semilogy(SNR_list,mean_error_mmse,'-*','LineWidth',2);
 hold on;
   semilogy(SNR_list,mean_error_omp,'-*','LineWidth',2);
 hold on;
 legend('proposed','ESPRIT','MMSE','OMP');
ylabel('NMSE');
xlabel('SNR in dB');
grid on 