function [idx,T1,T2,T3] = KNEG_cl(data,K)
%idx：Clustering label
%T1：Time taken to construct KNEG
%T2：Time taken to construct subclusters
%T3：Time taken to merge or split
%writen by Zexuan Fei
disp('dyknn start!!');
%data pre-processing
data(:,all(data==0, 1))=[];
n=size(data,1);
data=data_norm(data);%Normalization
k=floor(sqrt(n));
kf=0;
%% Construct KNEG
t1=clock;
[knn_neigh,knn_dist] = knnsearch(data,data,'k',k);
knn_neigh0=knn_neigh;%Backup result
knn_oppsite=ones(n,k);
for i=2:k
    for j=1:n
        neigh=knn_neigh(j,i);
        if neigh~=0
            loc=find(knn_neigh(neigh,:)==j);
            if isempty(loc) || loc>i
                knn_neigh(neigh,loc)=0;
                knn_dist(neigh,loc)=inf;
                knn_oppsite(j,i)=0;
            end
        end
    end
end

%% obtain peak density point
colordef black;
n_np_list=[];
is_np=[];
is_np2=[];
is_np_list=[];
np_list=[];
% last_flag=0;
for kk=1:k
    is_np=zeros(n,1);
    knn_oppsite2=~knn_oppsite;
    judge=sum(knn_neigh(:,1:kk).*knn_oppsite2(:,1:kk),2);
    is_np(judge==0)=1;
    is_np_list=[is_np_list,is_np];
    n_np=length(find(is_np==1));
    n_np_list=[n_np_list,n_np];
    if kk~=1
        if n_np_list(end-1)-n_np_list(end)<=n^0.5 %suspensive condition
            kf=kk;
            np_list=(find(is_np==1));
            break;
        end
    end
end
if n_np<K
    is_np=is_np2;
    np_list=(find(is_np==1));
    n_np=length(find(is_np==1));
end
t2=clock;
T1=etime(t2,t1);

%% Get subclusters(saved in group and np_group)
group=zeros(n,1);
group(is_np==1)=1:length(find(is_np==1));
count=0;
t3=clock;
for i=1:n_np    
    startnp=np_list(i);
    np_around=find(knn_neigh(startnp,1:kf)>0);
    for j=2:size(np_around,2)
        if is_np(knn_neigh(startnp,np_around(j)))==1
            group(knn_neigh(startnp,np_around(j)))=group(startnp);
            count=count+1;
        end
    end
end
%Initial partition
for i=1:n
    if group(i)>0
        continue;
    end
    next_is=i;
    next_i_loc=[];
    loc=find(knn_neigh(i,:)>0);
    for j=1:size(loc,2)
        if knn_oppsite(i,loc(j))==0
            next_i_loc=loc(j);
            break;
        end
    end
    if isempty(next_i_loc)
        next_i_loc=loc(2);
    end
    next_i=knn_neigh(i,next_i_loc);
    while group(next_i)==0 
        next_i0=next_i;
        next_is=[next_is,next_i];
        loc=find(knn_neigh(next_i,:)>0);
        for j=2:size(loc,2)
            if knn_oppsite(next_i,loc(j))==0
                next_i_loc=loc(j);
                break;
            end
        end
        if isempty(next_i_loc)
            next_i_loc=loc(2);
        end
        next_i=knn_neigh(next_i,next_i_loc);
        trace_back=0;
        while ~isempty(find(next_is==next_i,1))
            loc=find(knn_neigh(next_i,:)>0);
            flag=0;
            for j=2:size(loc,2)
                    if isempty(find(next_is==knn_neigh(next_i,loc(j)),1))
                        next_i_loc=loc(j);
                        flag=1;
                        break;
                    end
            end
            if flag==1 
                next_i=knn_neigh(next_i,next_i_loc);
                trace_back=0;
            else
                if size(next_is,2)>trace_back%backtrack
                    next_i=next_is(end-trace_back);
                    trace_back=trace_back+1;
                else
                    %can't backtrack
                    is_np(next_is(1))=1;
                    group(next_is(1))=length(find(is_np==1));
                    group(next_is)=group(next_is(1));
                    n_np=n_np+1;
                    np_list=[np_list;next_is(1)];
                    trace_back=0;
                    break;
                end

            end
        end
    end
    group(next_is)=group(next_i);
end

%weight np_group
np_list0=(find(is_np_list(:,kf-1)==1));

is_visited_np=[];
num_delete=0;
for i=1:n_np
    np_cur=np_list(i);
    if isempty(find(is_visited_np==group(np_cur),1))
        is_visited_np=[is_visited_np,group(np_cur)];
    else
        np_list(i)=0;
        num_delete=num_delete+1;
        is_np(np_cur)=0;
        continue;
    end
end
np_list(np_list==0)=[];
n_np=n_np-num_delete;
%formation of small clusters
max_pec=tabulate(group(:));
max_pec=max(max_pec(:,3));
maxsize=ceil(size(data,1)*max_pec/100)+1;
np_group_mem=zeros(maxsize,n_np);
np_group_pri=zeros(maxsize,n_np);
for i=1:n_np
    np_cur=np_list(i);
    mem=find(group==group(np_cur));
    [~,pri]=ismember(mem,knn_neigh0(np_cur,1:kf));
    fill=zeros(maxsize-size(mem,1),1);
    mem=[mem;fill];
    pri=[pri;fill];
    np_group_mem(:,i)=mem;
    np_group_pri(:,i)=pri;
end
t4=clock;
T2=etime(t4,t3);


%% Calculate the connection between np_group
t5=clock;
relation=zeros(n_np,n_np);
relation_count=zeros(n_np,n_np);
is_visited_np=[];
np_group_pri(np_group_pri==0)=kf+0.5;
for i=1:n_np
    if isempty(find(is_visited_np==group(np_list(i))))
        is_visited_np=[is_visited_np,group(np_list(i))];
    else
        continue;
    end
    np_around=np_group_mem(:,i);
    for j=1:find(np_around==0)-1
        np_around_around=knn_neigh(np_around(j),1:floor(2*kf-2));
        for m=1:size(np_around_around,2)
            if np_around_around(m)~=0
                [row,col]=find(np_group_mem==np_around_around(m));
                for n=1:size(row,1)
                    if(col(n)==i)
                        continue;
                    end                    
                    relation(i,col(n))=relation(i,col(n))+(2*kf+2-np_group_pri(row(n),col(n))-np_group_pri(j,i))/2;
                    relation_count(i,col(n))=relation_count(i,col(n))+1;
                end
            end
        end
    end
end
%% the correlation among density peak points
relation_norm=data_norm(relation);
in_out_rate=zeros(n_np,n_np);
for i=1:n_np
    x=data(np_list(i),:);
    x_row=find(relation(i,:)>0);
    if ~isempty(x_row)
        for j=1:size(x_row,2)
            y=data(np_list(x_row(j)),:);
            in_out_rate(i,x_row(j))=relation(x_row(j),i)/relation(i,x_row(j));
        end
    end
end
similarity=in_out_rate;
np_similar=zeros(n_np,n_np);
similarity(similarity>1)=0;
np_most_outline=zeros(n_np,n_np);
max_chubian=[];
relation2=relation;
same_cluster=0.8;
for i=1:n_np
    x=data(np_list(i),:);
    a=similarity(i,:);
    x_row=find((similarity(i,:)<1)&(similarity(i,:)>same_cluster));
    for j=1:size(x_row,2)
        y=data(np_list(x_row(j)),:);
        np_similar(i,x_row(j))=1;
    end
    [x_row_max,index_max]=max(relation(i,:));
    relation2(i,index_max)=0;
    x_num=length(find(group==group(np_list(i))));
    y_num=length(find(group==group(np_list(index_max))));
    x_chu=relation(i,index_max);
    y=data(np_list(index_max),:);   
    if x_row_max==0
        max_chubian=[max_chubian;i,i,1,x_num,x_num];
        y=x;
    else
        max_chubian=[max_chubian;i,index_max,in_out_rate(i,index_max),x_num,y_num];
    end
end

%% Merger or split
s=max_chubian(:,1);
t=max_chubian(:,2);
single_chu=max_chubian(~ismember(s,t),:);
single_chu(:,2)=[];
for i=1:size(single_chu,1)
    relation(single_chu(i,1),:)=relation(single_chu(i,1),:).*single_chu(i,2);
end
%construct weakly connected components
w=1./(max_chubian(:,3));
G = graph(s,t);
bins = conncomp(G);
finished=0;

if max(bins)<=K %split
    toall=max_chubian(:,end)+max_chubian(:,end-1);
    max_chubian(:,end-2)=max_chubian(:,end-2).*toall;
    np_relation=max_chubian(:,1:3);
    [np_relation,loc]=sortrows(np_relation,3);%relative membership degree
    for i=1:size(np_relation,1)
        if max(bins)==K
            finished=1;
            break;
        end
        G = rmedge(G,[s(loc(i))],[t(loc(i))]);
        bins = conncomp(G);
    end
end
%show result
np_cluster=1;
np_clusters=[];
np_datanum=[];
nps=cell(max(bins),1);
while ~isempty(find(bins==np_cluster,1))
    index=find(bins==np_cluster);
    for j=1:size(index,2)
        np_clusters=[np_clusters,max_chubian(index(j),1)];
    end
    np_clusters=unique(np_clusters);
    nps(np_cluster)={np_clusters};
    sum_datanum=0;
    for j=1:size(np_clusters,2)
        group(group==group(np_list(np_clusters(j))))=np_cluster;%出错了！
    end
    sum_datanum=length(find(group==np_cluster));
    np_datanum=[np_datanum;sum_datanum];
    np_clusters=[];
    np_cluster=np_cluster+1;
end
% merge
nps_relation=[];
if finished==0
    nps_relation=[];
    for i=1:size(nps,1)-1
        nps_i=nps{i};
        for j=i+1:size(nps,1)
            nps_j=nps{j};
            i2j=relation(nps_i,nps_j);
            j2i=relation(nps_j,nps_i);
            i2j_num=sum(i2j(:));
            j2i_num=sum(j2i(:));
            toall=i2j_num+j2i_num;
            if i2j_num<=j2i_num
                in_out_rate2=i2j_num/j2i_num;
            else
                in_out_rate2=j2i_num/i2j_num;
            end
            if toall~=0
                nps_relation=[nps_relation;[i,j,in_out_rate2,toall]];
            end
        end
    end
end
group_nps=1:size(nps,1);
if ~isempty(nps_relation) && finished==0
    nps_relation(:,3)=nps_relation(:,3).*nps_relation(:,4);
    nps_relation_sorted=sortrows(nps_relation,-3);
    is_incluster=zeros(size(nps,1),1);
    u=1:size(np_datanum,1);
    nums_sum=[u',np_datanum];
    for i=1:size(nps_relation_sorted,1)
        if ~isempty(nums_sum)
        end
        if is_incluster(nps_relation_sorted(i,1))==0 && is_incluster(nps_relation_sorted(i,2))==0
            a=group_nps(nps_relation_sorted(i,1));
            b=group_nps(nps_relation_sorted(i,2));
            group_nps(nps_relation_sorted(i,2))=group_nps(nps_relation_sorted(i,1));
            is_incluster(nps_relation_sorted(i,1))=1;
            is_incluster(nps_relation_sorted(i,2))=1;
            b_num=np_datanum(nps_relation_sorted(i,2));
            nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
            nums_sum(nums_sum(:,1)==b,:)=[];
        elseif is_incluster(nps_relation_sorted(i,1))==1 && is_incluster(nps_relation_sorted(i,2))==1
            if group_nps(nps_relation_sorted(i,1))~=group_nps(nps_relation_sorted(i,2))
                a=group_nps(nps_relation_sorted(i,1));
                b=group_nps(nps_relation_sorted(i,2));
                group_nps(group_nps==group_nps(nps_relation_sorted(i,2)))=group_nps(nps_relation_sorted(i,1));
                b_num=nums_sum(nums_sum(:,1)==b,2);
                nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
                nums_sum(nums_sum(:,1)==b,:)=[];
            end
        elseif  is_incluster(nps_relation_sorted(i,1))==0
            a=group_nps(nps_relation_sorted(i,2));
            b=group_nps(nps_relation_sorted(i,1));
            group_nps(nps_relation_sorted(i,1))=group_nps(nps_relation_sorted(i,2));
            is_incluster(nps_relation_sorted(i,1))=1;
            b_num=np_datanum(nps_relation_sorted(i,1));
            nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
            nums_sum(nums_sum(:,1)==b,:)=[];
        
        elseif  is_incluster(nps_relation_sorted(i,2))==0
            a=group_nps(nps_relation_sorted(i,1));
            b=group_nps(nps_relation_sorted(i,2));
            group_nps(nps_relation_sorted(i,2))=group_nps(nps_relation_sorted(i,1));
            is_incluster(nps_relation_sorted(i,2))=1;
            b_num=np_datanum(nps_relation_sorted(i,2));
            nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
            nums_sum(nums_sum(:,1)==b,:)=[];
        end
        if length(unique(group_nps))==K
            break;
        end
    end
else
end
 group_f=zeros(size(group,1),1);
 [group_nps_u]=unique(group_nps);
 np_cluster=1;
 for i=1:size(group_nps_u,2)
     group_nps_id=group_nps_u(i);
     onegroup=find(group_nps==group_nps_id);
     for j=1:size(onegroup,2)
         group_f(group==onegroup(j))=np_cluster;
     end
     np_cluster=np_cluster+1;
 end
 t6=clock;
 T3=etime(t6,t5);
idx=group_f;
disp('dyknn end!!');
end