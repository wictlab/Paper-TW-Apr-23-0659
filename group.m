run('C:\Users\jeffrey\Desktop\殷本全\殷本全\殷本全毕业论文\终期报告\code\Map.m');
%H=generate_mmwave_channel(16, 8, 4, 4);
Mr=16;%基站天线数
Mt=8;%用户天线数

f0=60;%GHz中心频率
nums_ue=36;%用户数量
nums_bs=3;%基站数
ax=32.4;%信号衰减指数
bx=2;%路径损耗指数
d0=1;%参考距离
Nx=6.8;%阴影衰落
c=299792458;%m/s 光速


%%%%%%%%%%%%%%%%%%%%%%%%%系统初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distance=zeros(nums_ue,nums_bs);%定义一个距离矩阵
UE= vertcat(u1, u2,u3);
for i=1:nums_bs
    bs_position=[BS(i,1),BS(i,2)];
    for j=1:nums_ue
        ue_position=[UE(j,1),UE(j,2)];
        distance(j,i)=sqrt((BS(i,1)-UE(j,1))^2+(BS(i,2)-UE(j,2))^2);
    end
end

Loss=zeros(nums_ue,nums_bs);%路径损耗系数矩阵
for i=1:nums_bs
    for j =1:nums_ue
           Loss(j,i)=ax+10*bx*log10(distance(j,i)/d0)+20*log10(4*pi*f0*d0*10^6/c)+Nx;
    end
end

H0 = cell(nums_ue, nums_bs);%定义一个cell数组放产生的信道矩阵H

for i=1:nums_ue
    for j=1:nums_bs
         %H0{i,j}=generate_mmwave_channel(Mr,Mt,randi([5,10]),randi([5,20]));
         H0{i,j}=generate_mmwave_channel(Mr,Mt,10,10);
    end
end

H=cell(nums_ue,nums_bs);%路径损耗系数修正后的矩阵
for i=1:nums_ue
    for j=1:nums_bs
        H{i,j}=H0{i,j}/Loss(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%用户分组第一阶段%%%%%%%%%%%%%%%%%%%%%%
%步骤一
%T联合干扰矩阵，i,j代表地i个小区内地第j用户的联合干扰矩阵
T=cell(nums_ue,nums_bs);
for i=1:nums_ue
    T{i,1}=[H{i,2}',H{i,3}'];
end
for i=1:nums_ue
    T{i,2}=[H{i,1}',H{i,3}'];
end
for i=1:nums_ue
    T{i,3}=[H{i,1}',H{i,2}'];
end
%步骤二
U=cell(nums_ue,nums_bs);
S=cell(nums_ue,nums_bs);%对角矩阵
V=cell(nums_ue,nums_bs);%右奇异矩阵

for i=1:nums_ue
    for j=1:nums_bs
    [U{i,j},S{i,j},V{i,j}]=svd(T{i,j});
    end
end

lamda=zeros(nums_ue,nums_bs);
for i=1:nums_ue
    for j=1:nums_bs
    lamda(i,j)=S{i,j}(Mt,Mt);
    end
end

theta2=zeros(nums_ue,nums_bs);%小区间干扰
for i=1:nums_ue
    for j=1:nums_bs 
    theta2(i,j)=lamda(i,j)*lamda(i,j);
    end
end

UE_information=zeros(nums_ue,2);
UE_num=zeros(nums_ue,1);
UE_Cellnum=zeros(nums_ue,1);
for  i=1:nums_ue
    UE_num(i,1)=i;
end

for i=1:nums_ue
    [min_value,min_index]=min(theta2(i,:));
    min_indices=find(theta2(i,:)==min_value);
    UE_Cellnum(i,1)=min_indices;
end
UE_information=[UE_num,UE_Cellnum];
%W=cell(nums_ue,1);%波束成形矩阵
%for i=1:nums_ue
%    W{i,1}=V{i,1}(:,Mr);
%end
%F=cell(nums,1);%预编码矩阵

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%用户分组第二阶段%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%小区一

nums_cell1=0;
for i=1:nums_ue
    if(UE_information(i,2)==1)
        nums_cell1=nums_cell1+1;
    end
end
cell1=cell(nums_cell1,8);
j=0;
for i=1:nums_ue
    if(UE_information(i,2)==1)
        j=j+1;
        cell1{j,1}=j;
        cell1{j,2}=H{i,1};%第二列是信道矩阵
        cell1{j,7}=theta2(i,1);%%第七列放一下小区间干扰
    end
end
for i=1:nums_cell1
    cell1{i,3}=(cell1{i,2}')*cell1{i,2};%第三列所示H'H
end
for i=1:nums_cell1
    [V,D]=eig(cell1{i,3});
    %disp('特征值矩阵 D:');
    %disp(D);
    cell1{i,4}=D;
    %第四列就是lamda

    %disp('特征向量矩阵 V:');
    %disp(V);
    cell1{i,5}=V;
%第五列是xi
end
for i=1:nums_cell1
 tao=zeros(Mt,1);
    for j=1:Mt
        tao=cell1{i,4}(j,j)*cell1{i,5}(j,:);
    end
    cell1{i,6}=tao;
end

row=zeros(nums_cell1,nums_cell1);%信道相关性矩阵
for i=1:nums_cell1
    for j=1:nums_cell1
        if i==j
            row(i,j)=0;
        end
        row(i,j)=norm(((cell1{i,6})'*cell1{i,6})'*((cell1{j,6})'*cell1{j,6}))/(norm(((cell1{i,6})'*cell1{i,6}))*norm(((cell1{j,6})'*cell1{j,6}))); 
    end
end