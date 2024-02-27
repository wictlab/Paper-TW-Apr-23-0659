function[rmse_esprit]=ESPRIT(theta_true,T,N,snr)
%% 参数设定
c = 3e8;                                              % 光速
fc = 10e9;                                           % 载波频率
%N=32;
lambda = c/fc;                                        % 波长
d = lambda/2;                                         % 阵元间距，可设 2*d = lambda
twpi = 2.0*pi; 
k=2*pi/lambda;
derad = pi/180; 
theta = theta_true*derad;                              % 待估计角度
z = (0:d:(N-1)*d)';     % 阵元坐标分布

M = length(z);                                      % 阵元数
K = length(theta);                                    % 信源数
T = 15;                                              % 快拍数
%snr = 1;                                             % 信噪比

%% 阵列接收信号仿真模拟
A = exp(-1j*k*z*sin(theta));          % 流型矩阵
noise =sqrt(1/2)*(randn(N, T)+1j*randn(N, T));
Vj=sqrt((   10^(snr/10)   )/2);
S=Vj*(randn(K, T)+1j*randn(K, T));
X1=A*S+noise;
%  % unit random training symbols
%       s = 1;
%       
%       % uniform randomly generated values for the switches(signal)
%       p = randsrc(1, 1, [-1 1 1j -1j]);   
%       p = 1 / norm(p,2) * p;
%       % based on the input/output model from the reference paper
%       X1 = A*S + noise;
%   %X1=awgn(X,snr,'measured')

%% ESPRIT 算法
% 计算协方差矩阵
R = X1*X1'/T;
% 特征值分解并取得信号子空间
[U,D] = eig(R);                                       % 特征值分解
[D,I] = sort(diag(D));                                % 将特征值排序从小到大
U = fliplr(U(:, I));                                  % 对应特征矢量排序，fliplr 之后，较大特征值对应的特征矢量在前面
Us = U(:, 1:K);                                       % 信号子空间
% 角度估计
Ux = Us(1:M-1, :);
Uy = Us(2:M, :);

% 方法二：完全最小二乘法
Uxy = [Ux, Uy];
Uxy = Uxy'*Uxy;
[U,D] = eig(Uxy);
[D,I] = sort(diag(D));
F = fliplr(U(:,I));
F0 = F(1:K, K+1:K*2);                                 % F0是F的左上角部分
F1 = F(K+1:K*2, K+1:K*2);                             % F1是F的右下角部分
Psi = -F0/F1;

[T,Phi] = eig(Psi);
Theta = asin(-angle(diag(Phi))/pi)/derad;             % 估计角度
Theta = sort(Theta).';
disp('估计结果：');
disp(Theta);
[Theta,~]=sort(Theta, 'ascend');
rmse_esprit=norm(Theta-theta_true,'fro') / 2;
theta = Theta*derad; 
At=zeros(M,K);
for i=1:K
    At(:,i)=exp(-1j*k*z*sin(theta(i)));
end
theta_ls = pinv(At)*X1;%最小二乘解
H=At*theta_ls;
rmse_esprit=norm(H-A*S,'fro') / norm(A*S,'fro');