close all;clear all;clc;
timer = 0;
Nt = 2;
Nr = 16;
Ny = 6;
Nx = 6;
L_list = [2:1:6];
d = 1;
gaoa=[10 80 300];
gaod=[10 80 300];
snr = 5;
sample_num = 100;
nmse_result = zeros(numel(L_list),size(gaoa,2));
mean_error_mmse=zeros(numel(L_list),size(gaoa,2));
mean_error_omp=zeros(numel(L_list),size(gaoa,2));
for g_num=1:size(gaoa,2)
    Gaoa=gaoa(g_num);
    Gaod=gaod(g_num);
    search_area_aoa=(-90:180/Gaoa:90-180/Gaoa);
    search_area_aod=(-90:180/Gaod:90-180/Gaod);
    A_aoa=zeros(Nr,Gaoa);
    grid_aoa=zeros(1,Gaoa);
    for i=1:Gaoa
        grid_aoa(:,i)=sin(search_area_aoa(i)*pi/180);
        A_aoa(:,i)=exp(-1i*2*pi*(0:Nr-1)'*sin(search_area_aoa(i)*pi/180));
    end
    A_aod=zeros(Nt,Gaod);
    grid_aod=zeros(1,Gaod);
    for i=1:Gaoa
        grid_aod(:,i)=sin(search_area_aod(i)*pi/180);
    end
    for L_ii = 1:numel(L_list)
        L = L_list(L_ii)
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
            for l = 1:L
                at = exp(-1i*2*pi*[0:Nt-1]'*d*sin(aod_taps(l)*pi/180));
                ar = exp(-1i*2*pi*[0:Nr-1]'*d*sin(aoa_taps(l)*pi/180));
                H_vec = H_vec + alpha(l)*kron(ar,at);
                H = H + alpha(l)*(ar*at');
            end
            X = (1/sqrt(Nt)*exp(-1i*2*pi*rand(Nt,Nx)));
            Y = (H*X + noise*(normrnd(0, 1, Nr, Nx) + 1i*normrnd(0, 1, Nr, Nx)));
            temp=kron(X.'*(at').',ar);
            Y_vec=vec(Y);

            Rth = noise*sqrt(Ny*Nr);
            for i=1:Gaod
                A_aod(:,i)=exp(-1i*2*pi*(0:Nt-1)'*sin(search_area_aod(i)*pi/180));
            end
            AM=X.'*((A_aod').');
            dict=kron(AM,A_aoa);
            [theta_es,z_es,err]=proposed(dict,Y,X,Nx,Nt,Nr,Ny,Rth,Y_vec,Gaoa,Gaod,grid_aoa,grid_aod,L);

            X1 = (1/sqrt(Nt)*exp(-1i*2*pi*rand(Nt,1)));
            Y1 = (H*X1 + noise*(normrnd(0, 1, Nr, 1) + 1i*normrnd(0, 1, Nr, 1)));
            AM1=X1.'*((A_aod').');
            dict1=kron(AM1,A_aoa);

            H_LS=Y1;
            I=ones(Nr,1);
            I=diag(I);
            SNR = 10.^(snr/10);                                       % 转换
            noise_var_sqrt = sqrt(1./SNR);                               % 噪声方差 SNR=S/N，噪声功率为方差
            sigma_2 = abs(noise_var_sqrt).^2;
            H_mmse=H_LS*H_LS'*pinv(H_LS*H_LS'+I.*sigma_2)*H_LS;
            H_mmse=H_mmse/X1;
            mean_error_mmse(L_ii,g_num) = mean_error_mmse(L_ii,g_num)+norm(H_mmse-H,'fro') / norm(H,'fro');

            H_es = zeros(Nr, Nt);
            at = zeros(Nt,1);
            ar = zeros(Nr,1);
            for l = 1:numel(z_es)
                at = exp(-1i*2*pi*[0:Nt-1]'*theta_es(1,l));
                ar = exp(-1i*2*pi*[0:Nr-1]'*theta_es(2,l));
                H_es = H_es + z_es(l)*ar*at';
            end
            h_OMP50=OMP(dict1,Y1,L);
            H_OMP50=zeros(Nr,1);
            for h1=1:Gaoa*Gaod
                if h_OMP50(h1)~=0
                    Pos_t= floor((h1-1)/Gaoa)+1;%存储这一列的序号 先t再r
                    Pos_r = mod(h1-1,Gaoa)+1;%存储这一列的序号
                    at = exp(-1i*2*pi*[0:Nt-1]'*grid_aod(Pos_t));
                    ar = exp(-1i*2*pi*[0:Nr-1]'*grid_aoa(Pos_r));
                    H_OMP50=H_OMP50+h_OMP50(h1)*ar*at';
                end
            end
            mean_error_omp(L_ii,g_num) = mean_error_omp(L_ii,g_num)+norm(H_OMP50-H,'fro')^2/norm(H,'fro')^2;
            nmse_sample =norm(H-H_es,'fro')^2/norm(H,'fro')^2;
            nmse_result(L_ii,g_num) = nmse_result(L_ii,g_num)+nmse_sample;
        end
    end
end
nmse_result = nmse_result./sample_num;
mean_error_mmse=mean_error_mmse./sample_num;
mean_error_omp=mean_error_omp./sample_num;
figure;
semilogy(L_list,nmse_result(:,1),'-*','LineWidth',2);
hold on;
semilogy(L_list,nmse_result(:,2),'-*','LineWidth',2);
hold on;
semilogy(L_list,nmse_result(:,3),'-*','LineWidth',2);
hold on;
semilogy(L_list,mean_error_mmse(:,1),'-*','LineWidth',2);
hold on;
semilogy(L_list,mean_error_omp(:,1),'-*','LineWidth',2);
hold on;
semilogy(L_list,mean_error_omp(:,2),'-*','LineWidth',2);
hold on;
semilogy(L_list,mean_error_omp(:,3),'-*','LineWidth',2);
hold on;
legend('proposed,10','proposed,80','proposed,300','MMSE','OMP,10','OMP,80','OMP,300');
ylabel('NMSE');
xlabel('SNR in dB');
grid on