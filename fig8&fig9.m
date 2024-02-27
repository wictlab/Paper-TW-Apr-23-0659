rng(1) 
nums_cell1=0;
Mt1=Mt(nr);
Mr1=Mr(nr);
v1=v(nr);
for i=1:nums_ue
    if(UE_information(i,2)==3)
        nums_cell1=nums_cell1+1;
    end
end
cell1=cell(nums_cell1,9);%第九列放W
j=0;
for i=1:nums_ue
    if(UE_information(i,2)==3)
        j=j+1;
        numsrows=size(W{i,1},1);
        cell1{j,9}=W{i,1}(1:numsrows/2, :);%%%%%
        cell1{j,1}=j;
        cell1{j,2}=H{i,3};%第二列是信道矩阵
        cell1{j,7}=theta2(i,3);%%第七列放一下小区间干扰
    end
end
for i=1:nums_cell1
    cell1{i,3}=(cell1{i,2}')*cell1{i,2};%第三列所示H'H
end
for i=1:nums_cell1
    [V,D]=eig(cell1{i,3});
    cell1{i,4}=D;%第四列就是lamda
    cell1{i,5}=V;%第五列是xi
end
for i=1:nums_cell1
    tao=zeros(Mt1,1);
    for j=1:Mt1
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
B=100000000;%带宽总和
n0=-174;
n=10^(n0/10)/1000;
disturb=0.1;
D=40;
omega=2*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%小区三分组0.5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_3=[];
result_MI_3=[];
a=0.5;%test 分组阈值
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
    if row(1,j)<a
        group1_nums=group1_nums+1;
        group1(group1_nums,1)=j;
        group_matrix(j,1)=1;
    end
end
group_information{nums_group,2}=group1;
%%第二阶段步骤3
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
        if row(startindex(1,1),startindex_find(i,1))<a
            nums_group_new=nums_group_new+1;
            group_new(nums_group_new,1)=startindex_find(i);
            group_matrix(startindex_find(i,1))=1;
        end
    end
    group_information{nums_group,2}=group_new;
    j=0;
    clear group_new;
end
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
            if row(group_information{i,2}(j,1),group_information{i,2}(k,1))>a
                nums_less_than_a=nums_less_than_a+1;
            end
        end
        if nums_less_than_a>1 %添加新组
            nums_group=nums_group+1;
            group_new=[];
            group_new=group_information{i,2}(j);
            group_information{nums_group,1}=nums_group;
            group_information{nums_group,2}=group_new;
            nums_group_i= nums_group_i-1;%%回溯防止越界 %从该分组中剔除该用户
            group_information{i,2}(j)=[];
        end
    end
end
every_group=cell(nums_group,3);%第一列组内用户数量，第二列组的带宽,第三列用户成员
for i=1:nums_group
    every_group{i,1}=length(group_information{i,2});
    every_group{i,2}=B*(length(group_information{i,2})/nums_cell1);%组内复用带宽，组间带宽正交
    every_group{i,3}=group_information{i,2};
end
H_cell1=[];
trace1=0;
for i=1:nums_cell1
    H_cell1=[H_cell1,cell1{i,2}'];
end
H_cell1=(H_cell1')*H_cell1;
for i=1:length(H_cell1(1,:))
    for j=1:length(H_cell1(:,1))
        if i==j
            trace1=trace1+H_cell1(i,j);
        end
    end
end
b=[0.1:0.1:0.7];
com1=zeros(size(b,2),1);
pos1=zeros(size(b,2),1);
com3=zeros(size(b,2),1);
pos3=zeros(size(b,2),1);
com5=zeros(size(b,2),1);
pos5=zeros(size(b,2),1);
B =100000000;%hz 带宽
P = 12; % maximum transmit power 24dBm
Ucur = 0;%总用户数
Uset=zeros(nums_group,100);
for i=1:nums_group
    Ucur=Ucur+every_group{i,1};%每组用户个数
    for j=1:every_group{i,1}
        Uset(i,j)=every_group{i,3}(j,1);%每组用户编号组
    end
end

for i=1:size(b,2)
    Bo=(B*(1-b(i))/nums_group).*ones(nums_group,1);
Pc = (P*b(i)/Ucur).*ones(Ucur,1); 
Po = (P*(1-b(i))/Ucur).*ones(Ucur,1); 
Bc=(B*b(i)/nums_group).*ones(nums_group,1);
    [totcom,totpos,Pc1,Po1]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b(i),cell1,0);%ADMM
    [totcom_a,totpos_a,x,y]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b(i),cell1,1);%平均
    com1(i)=totcom;
    pos1(i)=totpos;
    com3(i)=totcom_a;
    pos3(i)=totpos_a;
    [totcom_b,totpos_b,x,y]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc1,Po1,B,nums_group,every_group,trace1,b(i),cell1,1);%平均
    com5(i)=totcom_b;
    pos5(i)=totpos_b;
end

com7=zeros(size(b,2),1);
pos7=zeros(size(b,2),1);
com8=zeros(size(b,2),1);
pos8=zeros(size(b,2),1);
com9=zeros(size(b,2),1);
pos9=zeros(size(b,2),1);
P = 5; 
for i=1:size(b,2)
    Bo=(B*(1-b(i))/nums_group).*ones(nums_group,1);
Pc = (P*b(i)/Ucur).*ones(Ucur,1);
Po = (P*(1-b(i))/Ucur).*ones(Ucur,1); 
Bc=(B*b(i)/nums_group).*ones(nums_group,1);
    [totcom,totpos,Pc1,Po1]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b(i),cell1,0);%ADMM
    [totcom_a,totpos_a,x,y]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b(i),cell1,1);%平均
    com7(i)=totcom;
    pos7(i)=totpos;
    com8(i)=totcom_a;
    pos8(i)=totpos_a;
    [totcom_b,totpos_b,x,y]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc1,Po1,B,nums_group,every_group,trace1,b(i),cell1,1);%平均
    com9(i)=totcom_b;
    pos9(i)=totpos_b;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%小区三分组0.35%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=12;
result_3=[];
result_MI_3=[];
a=0.35;%test 分组阈值
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
    if row(1,j)<a
        group1_nums=group1_nums+1;
        group1(group1_nums,1)=j;
        group_matrix(j,1)=1;
    end
end
group_information{nums_group,2}=group1;
%%第二阶段步骤3
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
        if row(startindex(1,1),startindex_find(i,1))<a
            nums_group_new=nums_group_new+1;
            group_new(nums_group_new,1)=startindex_find(i);
            group_matrix(startindex_find(i,1))=1;
        end
    end
    group_information{nums_group,2}=group_new;
    j=0;
    clear group_new;
end
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
            if row(group_information{i,2}(j,1),group_information{i,2}(k,1))>a
                nums_less_than_a=nums_less_than_a+1;
            end
        end
        if nums_less_than_a>1%添加新组
            nums_group=nums_group+1;
            group_new=[];
            group_new=group_information{i,2}(j);
            group_information{nums_group,1}=nums_group;
            group_information{nums_group,2}=group_new;
            nums_group_i= nums_group_i-1;%%回溯防止越界%从该分组中剔除该用户
            group_information{i,2}(j)=[];
        end
    end
end
P_every_ue=P/nums_cell1;
every_group=cell(nums_group,3);%第一列组内用户数量，第二列组的带宽,第三列用户成员
for i=1:nums_group
    every_group{i,1}=length(group_information{i,2});
    every_group{i,2}=B*(length(group_information{i,2})/nums_cell1);%组内复用带宽，组间带宽正交
    every_group{i,3}=group_information{i,2};
end
R_sum=0;
MI_sum=0;
H_cell1=[];
trace1=0;
for i=1:nums_cell1
    H_cell1=[H_cell1,cell1{i,2}'];
end
H_cell1=(H_cell1')*H_cell1;
for i=1:length(H_cell1(1,:))
    for j=1:length(H_cell1(:,1))
        if i==j
            trace1=trace1+H_cell1(i,j);
        end
    end
end
b=[0.1:0.1:0.7];
com2=zeros(size(b,2),1);
pos2=zeros(size(b,2),1);
com4=zeros(size(b,2),1);
pos4=zeros(size(b,2),1);
com6=zeros(size(b,2),1);
pos6=zeros(size(b,2),1);
Ucur = 0;%总用户数
Uset=zeros(nums_group,100);
for i=1:nums_group
    Ucur=Ucur+every_group{i,1};%每组用户个数
    for j=1:every_group{i,1}
        Uset(i,j)=every_group{i,3}(j,1);%每组用户编号组
    end
end

for i=1:size(b,2)
    Bo=(B*(1-b(i))/nums_group).*ones(nums_group,1);
Pc = (P*b(i)/Ucur).*ones(Ucur,1); 
Po = (P*(1-b(i))/Ucur).*ones(Ucur,1); 
Bc=(B*b(i)/nums_group).*ones(nums_group,1);
    [totcom,totpos,Pc1,Po1]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b(i),cell1,0);
    com2(i)=totcom;
    pos2(i)=totpos;
    [totcom_a,totpos_a,~,~]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b(i),cell1,1);
    com4(i)=totcom_a;
    pos4(i)=totpos_a;
    [totcom_b,totpos_b,~,~]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc1,Po1,B,nums_group,every_group,trace1,b(i),cell1,1);
    com6(i)=totcom_b;
    pos6(i)=totpos_b;
end

figure;
plot(com1,pos1,'o - k','linewidth',5)%所提算法
hold on
plot(com2,pos2,'^ -- k','linewidth',5)%无分组算法=每个用户单独一个组=a=1
hold on
plot(com3,pos3,'^ - b','linewidth',5)
hold on
plot(com4,pos4,'o -- b','linewidth',5)
hold on
plot(com5,pos5,'^ - r','linewidth',5)
hold on
plot(com6,pos6,'o -- r','linewidth',5)
hold on
legend('\fontname{Times new Roman}\fontsize{30}{Proposed joint algorithm,a=0.2}','\fontname{Times new Roman}\fontsize{30}{Proposed joint algorithm,a=0.7}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBPA, a=0.2}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBPA, a=0.7}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBA+proposed, a=0.2}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBA+proposed, a=0.7}');
xlabel('\fontname{Times new Roman}\fontsize{30}{total data rate/(bits/s)}');
ylabel('\fontname{Times new Roman}\fontsize{30}{total positioning estimation rate/(bits/s)}')
axis([-inf,inf -inf,inf]);
grid on;
set(gca,'Fontname','Times new roman','Fontsize',15)

figure;
plot(com1,pos1,'^ - k','linewidth',5)%所提算法
hold on
plot(com7,pos7,'o -- k','linewidth',5)%无分组算法=每个用户单独一个组=a=1
hold on
plot(com3,pos3,'^ - b','linewidth',5)
hold on
plot(com8,pos8,'o -- b','linewidth',5)
hold on
plot(com5,pos5,'^ - r','linewidth',5)
hold on
plot(com9,pos9,'o -- r','linewidth',5)
hold on
legend('\fontname{Times new Roman}\fontsize{30}{Proposed joint algorithm,P=12W}','\fontname{Times new Roman}\fontsize{30}{Proposed joint algorithm,P=5W}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBPA, ,P=12W}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBPA, P=5W}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBA+proposed, ,P=12W}','\fontname{Times new Roman}\fontsize{30}{Algorithm 2, EBA+proposed, P=5W}');
xlabel('\fontname{Times new Roman}\fontsize{30}{total data rate/(bits/s)}');
ylabel('\fontname{Times new Roman}\fontsize{30}{total positioning rate/(bits/s)}')
axis([-inf,inf -inf,inf]);
grid on;
set(gca,'Fontname','Times new roman','Fontsize',15)