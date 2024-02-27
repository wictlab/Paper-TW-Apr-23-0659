run('C:\Users\jeffrey\Desktop\殷本全\殷本全\殷本全毕业论文\终期报告\code\Map.m');
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
%%%%%%%%%%%%%%%%%%%用户分组第一阶段%%%%%%%%%%%%%%%%%%%%%%%步骤一
%T联合干扰矩阵，i,j代表地i个小区内第j用户的联合干扰矩阵
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
%for i=1:nums_cell1
 %   [V1,D1]=eig(cell1{i,3});
%    cell1{i,4}=V1;
 %   cell1{i,5}=D1;
%end
B=5;%带宽总和
P=50;%每个基站传输功率
n=0.01;
disturb=0.1;

  result=[];a=0;
  for a_x=1:10
      a=0.1*a_x;%test 分组阈值
a=0.3;

nums_group=0;
group_information=cell(100,100);%第一列是组数，第二列是成员信息
%将第一个用户加入组1
nums_group=nums_group+1;
group_information{1,1}=nums_group;
group1=[];group1(1,1)=1;
group1_nums=1;
group_matrix=zeros(nums_cell1,1);%用来标记用户是否进行了分组
group_matrix(1,1)=1;
for j=2:nums_cell1
    if row(1,j)>a
        group1_nums=group1_nums+1;
        group1(group1_nums,1)=j;
        group_matrix(j,1)=1;
    end
end
group_information{nums_group,2}=group1;
%%第二阶段步骤3
%for jishu=1:nums_cell1
    while(any(group_matrix==0))
        countzeros=sum(group_matrix==0);
        nums_group=nums_group+1;
        [startindex,colindex]=find(group_matrix==0);
        group_new=[]; 
        group_new(1,1)=startindex(1,1);
        nums_group_new=1;
        group_matrix(startindex(1))=1;
        j=countzeros;
        [startindex_find,colindex_find]=find(group_matrix==0);%找到另外的还没分组的用户
        for i=1:j-1
            if row(startindex(1,1),startindex_find(i,1))>a
               nums_group_new=nums_group_new+1;
               group_new(nums_group_new,1)=startindex_find(i);
               group_matrix(startindex_find(i,1))=1;
            end
        end
    group_information{nums_group,2}=group_new;
    j=0;
    clear group_new;
    end
%end
for i=1:nums_group
    group_information{i,1}=i;
end
%%第二阶段步骤4
for i=1:nums_group
    nums_group_i=length(group_information{i,2});
    for j=1:nums_group_i
        nums_less_than_a=0;
        for k=j+1:nums_group_i
            if k==j
                continue
            end
            if row(group_information{i,2}(j,1),group_information{i,2}(k,1))<a
                nums_less_than_a=nums_less_than_a+1;
            end
        end
        if nums_less_than_a>1
            %添加新组
            nums_group=nums_group+1;
            group_new=[];
            group_new=group_information{i,2}(j);
            group_information{nums_group,1}=nums_group;
            group_information{nums_group,2}=group_new;
            nums_group_i= nums_group_i-1;%%回溯防止越界
            %从该分组中剔除该用户
            group_information{i,2}(j)=[];
        end
    end
end
disp(nums_group);
%%%%%%%%%%%%SINR计算%%%%%%%%%%
%%此阶段不考虑资源分配的问题
%%%每个小区 每个用户平均分配功率， P=5w
%%%每个小区 总带宽B,
%SINR=P/theta2+噪声
P_every_ue=P/nums_cell1;
B_every_group=B;%每个组复用带宽B，但是组内要分带宽

for i=1:nums_cell1%获得每个用户的SINR
    cell1{i,8}=P_every_ue/(cell1{ i,7}+n+1.3*(nums_group-1));
end%第八列放置SINR
R_sum=0;
for i=1:nums_group
    B_every_ue=B_every_group/length(group_information{i,2});
    for j=1:length(group_information{i,2})
        R=B_every_ue*log2(1+cell1{group_information{i,2}(j,1),8});%-(nums_group-1)*0.4;
        R_sum=R_sum+R;
    end
end

disp(R_sum);
 result(a_x,1)=R_sum;     
 end
