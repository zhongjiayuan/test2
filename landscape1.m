close all
clc;
clear;

load weighted_delta_sd;
load weighted_delta_entropy;
load weighted_edges;
e=exp(1);

psize=size(weighted_delta_sd);
patients_num=[3,106,124,39,59,62,10,21];
%aver_CI=zeros(psize(1),psize(2));
yl=0;
for i=1:psize(1)
    tmp=zeros(psize(2),1);
    for l=1:psize(2)
        for s=1:patients_num(l)
            complex_index1(i,l,s)=weighted_delta_entropy(i,l,s)*(weighted_edges(i,l,s));  %*weighted_delta_sd(i,l,s);
            complex_index2(i,l,s)=weighted_delta_sd(i,l,s);
            complex_index3(i,l,s)=weighted_delta_entropy(i,l,s)*weighted_delta_sd(i,l,s)/weighted_edges(i,l,s);
            %aver_CI(i,l)=aver_CI(i,l)+complex_index3(i,l,s);
            tmp(l)=tmp(l)+complex_index3(i,l,s);
        end
        tmp(l)=tmp(l)/patients_num(l);
        %aver_CI(i,l)=aver_CI(i,l)/patients_num(l);
    end
    if sum(tmp)>0.001
        yl=yl+1;
        aver_CI(yl,:)=zscore(tmp);
    end
    miao(i,:)=tmp;
end

% minCI=min(min(aver_CI));
% maxCI=max(max(aver_CI));
% II_adjust=(aver_CI-minCI)/(maxCI-minCI);
II_adjust=aver_CI;
psize=size(II_adjust);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for t=1:psize(2)
%     if t==7
%         II_adjust(find(II(:,t)>3),t)=rand(1,size(II(find(II(:,t)>3),t),1))*3.5;
%         continue
%     end
%     II_adjust(find(II(:,t)>3),t)=rand(1,size(II(find(II(:,t)>3),t),1))*3;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
[I,idx]=sort(II_adjust(:,7),'descend');
for i=1:psize(1)
    II_new(i,:)=II_adjust(idx(i),:);
end
%psize=size(II_new);

clear Y Z T num
Y=pdist(II_new(:,7),@cofun2);
Z=linkage(Y);
[H,ter,perm] = dendrogram(Z, 0, 'orientation', 'left', ...
    'colorthreshold', 'default');
% tmp1=II';
% tmp2=tmp1(:);
% tmp2=zscore(tmp2);
% for i=1:msize(1)
%     for j=1:psize(2)
%         zs_II(i,j)=tmp2(psize(2)*(i-1)+j);
%     end
% end
clear II_final;
for i=1:psize(1)
    II_final(i,:)=II_new(perm(i),:);
end
II_draw=II_new(1:psize(1),:);
surf([1:psize(2)],[1:psize(1)],II_draw);
shading interp;