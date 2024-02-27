function [ theta ] = OMP( A1,y,t )
    A=A1;
    [M,N] = size(A);%传感矩阵A为M*N矩阵
    theta = zeros(N,1);%用来存储恢复的theta(列向量)
    At = zeros(M,t);%用来迭代过程中存储A被选择的列
    Pos_theta = zeros(1,t);%用来迭代过程中存储A被选择的列序号
    r_n = y;%初始化残差(residual)为y
    for ii=1:t%迭代t次，t为输入参数
        product = A'*r_n;%传感矩阵A各列与残差的内积
        [val,pos] = max(abs(product));%找到最大内积绝对值，即与残差最相关的列
        At(:,ii) = A(:,pos);%存储这一列
        Pos_theta(ii) = pos;%存储这一列的序号
        %y=At(:,1:ii)*theta，以下求theta的最小二乘解(Least Square)
        A(:,pos)=zeros(M,1);
        theta_ls = pinv(At(:,1:ii))*y;%最小二乘解
        %At(:,1:ii)*theta_ls是y在At(:,1:ii)列空间上的正交投影
        r_n = y - At(:,1:ii)*theta_ls;%更新残差        
    end
    theta(Pos_theta)=theta_ls;%恢复出的theta
end