clear;
clc;
close all;


e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.01;

p=zeros(1,25);

p(1)=-0.5;
p(2)=-0.475;
p(3)=-0.45;
p(4)=-0.4;
p(5)=-0.375;
p(6)=-0.35;
p(7)=-0.325;
p(8)=-0.3;
p(9)=-0.275;
p(10)=-0.25;
p(11)=-0.225;
p(12)=-0.2;
p(13)=-0.175;
p(14)=-0.15;
p(15)=-0.125;
p(16)=-0.1;
p(17)=-0.08;
p(18)=-0.05;
p(19)=-0.02;
%p(20)=-0.0001;
p(20)=-0.001;
p(21)=0.02;
p(22)=0.05;
p(23)=0.08;
p(24)=0.1;
p(25)=0.15;







D= [-2 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
T=1;
TT=100;
kld1=zeros(1,24);

sample_num=2;
CC=zeros(8,4,sample_num);




CC1=zeros(8,4,20);
pp=-0.5:5/1900:-0.45;
%pp=-0.5:2/70:-0.3;

for ll=1:20
    qq(ll)=0.96^(1/abs(pp(ll)));
    E=[-2/5*qq(ll) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
    J=D*E*inv(D);
    for i=1:N-1;
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8;
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
        end
    end
    CC1(:,1,ll)=X(:,2000); 
end


pre_TC=reshape(CC1(:,1,:),8,20);
for s=1:TT
    for l=2:25
        q(l)=0.96^(1/abs(p(l)));
        E=[-2/5*q(l) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:8
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,2000);
        end
        TC=reshape(CC(:,l,:),8,sample_num);
        
        for i=1:8
             mu0=mean(pre_TC(i,1:20));
             sigma0=std(pre_TC(i,1:20));
             tumor(i)= normpdf(TC(i,1),mu0,sigma0);
             pvalue(i,1:20)=normpdf(pre_TC(i,1:20),mu0,sigma0);
             casecdf(i)=normcdf(TC(i,1),mu0,sigma0);
             
             %casecdf(i)=normpdf(TC(i,1),mu0,sigma0);
             if casecdf(i)<(normcdf(mu0-1*sigma0,mu0,sigma0))
                 casecdf(i)=1-casecdf(i);
             end
             concdf(i,1:20)=normcdf(pre_TC(i,1:20),mu0,sigma0);
             %concdf(i,1:10)=normpdf(pre_TC(i,1:10),mu0,sigma0);
        end
        %5
        [tmp_com_idx,index]=sort(tumor);
        for k=1:8
            mor(k)=casecdf(index(k));
            con(k)=mean(concdf(index(k),1:20));
        end
        kld1(l-1)= kld1(l-1)+0.5*(sum(con .* log(con./mor))+sum(mor.* log(mor./con)));
    end
        
        
         %mu1=mean(mor);
         %sigma1=std(mor);
         %mor= normpdf(mor,mu1,sigma1);
        
         %mu2=mean(con);
         %sigma2=std(con);
         %con= normpdf(con,mu2,sigma2);
        
        
        
        
        %mor=mor/sum(mor);
        %con=con/sum(con);
        
      
   
end

kld1=kld1/TT;

%%
%subplot(1,2,1);
save kld1.mat;
%plot(p(2:25),kld,'k-','LineWidth',2);
%title('KL');
%subplot(1,2,2);

plot(p(2:25),kld1,'k-','LineWidth',2);
title('KL-sym');




clear;
clc;
close all;


e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.05;

p=zeros(1,25);

p(1)=-0.5;
p(2)=-0.475;
p(3)=-0.45;
p(4)=-0.4;
p(5)=-0.375;
p(6)=-0.35;
p(7)=-0.325;
p(8)=-0.3;
p(9)=-0.275;
p(10)=-0.25;
p(11)=-0.225;
p(12)=-0.2;
p(13)=-0.175;
p(14)=-0.15;
p(15)=-0.125;
p(16)=-0.1;
p(17)=-0.08;
p(18)=-0.05;
p(19)=-0.02;
%p(20)=-0.0001;
p(20)=-0.001;
p(21)=0.02;
p(22)=0.05;
p(23)=0.08;
p(24)=0.1;
p(25)=0.15;







D= [-2 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
T=1;
TT=100;
kld2=zeros(1,24);

sample_num=2;
CC=zeros(8,4,sample_num);




CC1=zeros(8,4,20);
pp=-0.5:5/1900:-0.45;
%pp=-0.5:2/70:-0.3;

for ll=1:20
    qq(ll)=0.96^(1/abs(pp(ll)));
    E=[-2/5*qq(ll) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
    J=D*E*inv(D);
    for i=1:N-1;
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8;
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
        end
    end
    CC1(:,1,ll)=X(:,2000); 
end


pre_TC=reshape(CC1(:,1,:),8,20);
for s=1:TT
    for l=2:25
        q(l)=0.96^(1/abs(p(l)));
        E=[-2/5*q(l) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:8
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,2000);
        end
        TC=reshape(CC(:,l,:),8,sample_num);
        
        for i=1:8
             mu0=mean(pre_TC(i,1:20));
             sigma0=std(pre_TC(i,1:20));
             tumor(i)= normpdf(TC(i,1),mu0,sigma0);
             pvalue(i,1:20)=normpdf(pre_TC(i,1:20),mu0,sigma0);
             casecdf(i)=normcdf(TC(i,1),mu0,sigma0);
             
             %casecdf(i)=normpdf(TC(i,1),mu0,sigma0);
             if casecdf(i)<(normcdf(mu0-1*sigma0,mu0,sigma0))
                 casecdf(i)=1-casecdf(i);
             end
             concdf(i,1:20)=normcdf(pre_TC(i,1:20),mu0,sigma0);
             %concdf(i,1:10)=normpdf(pre_TC(i,1:10),mu0,sigma0);
        end
        %5
        [tmp_com_idx,index]=sort(tumor);
        for k=1:8
            mor(k)=casecdf(index(k));
            con(k)=mean(concdf(index(k),1:20));
        end
        kld2(l-1)= kld2(l-1)+0.5*(sum(con .* log(con./mor))+sum(mor.* log(mor./con)));
    end
        
        
         %mu1=mean(mor);
         %sigma1=std(mor);
         %mor= normpdf(mor,mu1,sigma1);
        
         %mu2=mean(con);
         %sigma2=std(con);
         %con= normpdf(con,mu2,sigma2);
        
        
        
        
        %mor=mor/sum(mor);
        %con=con/sum(con);
        
      
   
end

kld2=kld2/TT;

%%
%subplot(1,2,1);
save kld2.mat;
%plot(p(2:25),kld,'k-','LineWidth',2);
%title('KL');
%subplot(1,2,2);

plot(p(2:25),kld2,'k-','LineWidth',2);
title('KL-sym');




clear;
clc;
close all;


e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.1;

p=zeros(1,25);

p(1)=-0.5;
p(2)=-0.475;
p(3)=-0.45;
p(4)=-0.4;
p(5)=-0.375;
p(6)=-0.35;
p(7)=-0.325;
p(8)=-0.3;
p(9)=-0.275;
p(10)=-0.25;
p(11)=-0.225;
p(12)=-0.2;
p(13)=-0.175;
p(14)=-0.15;
p(15)=-0.125;
p(16)=-0.1;
p(17)=-0.08;
p(18)=-0.05;
p(19)=-0.02;
%p(20)=-0.0001;
p(20)=-0.001;
p(21)=0.02;
p(22)=0.05;
p(23)=0.08;
p(24)=0.1;
p(25)=0.15;







D= [-2 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
T=1;
TT=100;
kld3=zeros(1,24);

sample_num=2;
CC=zeros(8,4,sample_num);




CC1=zeros(8,4,20);
pp=-0.5:5/1900:-0.45;
%pp=-0.5:2/70:-0.3;

for ll=1:20
    qq(ll)=0.96^(1/abs(pp(ll)));
    E=[-2/5*qq(ll) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
    J=D*E*inv(D);
    for i=1:N-1;
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8;
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
        end
    end
    CC1(:,1,ll)=X(:,2000); 
end


pre_TC=reshape(CC1(:,1,:),8,20);
for s=1:TT
    for l=2:25
        q(l)=0.96^(1/abs(p(l)));
        E=[-2/5*q(l) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:8
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,2000);
        end
        TC=reshape(CC(:,l,:),8,sample_num);
        
        for i=1:8
             mu0=mean(pre_TC(i,1:20));
             sigma0=std(pre_TC(i,1:20));
             tumor(i)= normpdf(TC(i,1),mu0,sigma0);
             pvalue(i,1:20)=normpdf(pre_TC(i,1:20),mu0,sigma0);
             casecdf(i)=normcdf(TC(i,1),mu0,sigma0);
             
             %casecdf(i)=normpdf(TC(i,1),mu0,sigma0);
             if casecdf(i)<(normcdf(mu0-1*sigma0,mu0,sigma0))
                 casecdf(i)=1-casecdf(i);
             end
             concdf(i,1:20)=normcdf(pre_TC(i,1:20),mu0,sigma0);
             %concdf(i,1:10)=normpdf(pre_TC(i,1:10),mu0,sigma0);
        end
        %5
        [tmp_com_idx,index]=sort(tumor);
        for k=1:8
            mor(k)=casecdf(index(k));
            con(k)=mean(concdf(index(k),1:20));
        end
        kld3(l-1)= kld3(l-1)+0.5*(sum(con .* log(con./mor))+sum(mor.* log(mor./con)));
    end
        
        
         %mu1=mean(mor);
         %sigma1=std(mor);
         %mor= normpdf(mor,mu1,sigma1);
        
         %mu2=mean(con);
         %sigma2=std(con);
         %con= normpdf(con,mu2,sigma2);
        
        
        
        
        %mor=mor/sum(mor);
        %con=con/sum(con);
        
      
   
end

kld3=kld3/TT;

%%
%subplot(1,2,1);
save kld3.mat;
%plot(p(2:25),kld,'k-','LineWidth',2);
%title('KL');
%subplot(1,2,2);

plot(p(2:25),kld3,'k-','LineWidth',2);
title('KL-sym');







clear;
clc;
close all;


e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.2;

p=zeros(1,25);

p(1)=-0.5;
p(2)=-0.475;
p(3)=-0.45;
p(4)=-0.4;
p(5)=-0.375;
p(6)=-0.35;
p(7)=-0.325;
p(8)=-0.3;
p(9)=-0.275;
p(10)=-0.25;
p(11)=-0.225;
p(12)=-0.2;
p(13)=-0.175;
p(14)=-0.15;
p(15)=-0.125;
p(16)=-0.1;
p(17)=-0.08;
p(18)=-0.05;
p(19)=-0.02;
%p(20)=-0.0001;
p(20)=-0.001;
p(21)=0.02;
p(22)=0.05;
p(23)=0.08;
p(24)=0.1;
p(25)=0.15;







D= [-2 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
T=1;
TT=100;
kld4=zeros(1,24);

sample_num=2;
CC=zeros(8,4,sample_num);




CC1=zeros(8,4,20);
pp=-0.5:5/1900:-0.45;
%pp=-0.5:2/70:-0.3;

for ll=1:20
    qq(ll)=0.96^(1/abs(pp(ll)));
    E=[-2/5*qq(ll) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
    J=D*E*inv(D);
    for i=1:N-1;
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8;
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
        end
    end
    CC1(:,1,ll)=X(:,2000); 
end


pre_TC=reshape(CC1(:,1,:),8,20);
for s=1:TT
    for l=2:25
        q(l)=0.96^(1/abs(p(l)));
        E=[-2/5*q(l) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:8
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,2000);
        end
        TC=reshape(CC(:,l,:),8,sample_num);
        
        for i=1:8
             mu0=mean(pre_TC(i,1:20));
             sigma0=std(pre_TC(i,1:20));
             tumor(i)= normpdf(TC(i,1),mu0,sigma0);
             pvalue(i,1:20)=normpdf(pre_TC(i,1:20),mu0,sigma0);
             casecdf(i)=normcdf(TC(i,1),mu0,sigma0);
             
             %casecdf(i)=normpdf(TC(i,1),mu0,sigma0);
             if casecdf(i)<(normcdf(mu0-1*sigma0,mu0,sigma0))
                 casecdf(i)=1-casecdf(i);
             end
             concdf(i,1:20)=normcdf(pre_TC(i,1:20),mu0,sigma0);
             %concdf(i,1:10)=normpdf(pre_TC(i,1:10),mu0,sigma0);
        end
        %5
        [tmp_com_idx,index]=sort(tumor);
        for k=1:8
            mor(k)=casecdf(index(k));
            con(k)=mean(concdf(index(k),1:20));
        end
        kld4(l-1)= kld4(l-1)+0.5*(sum(con .* log(con./mor))+sum(mor.* log(mor./con)));
    end
        
        
         %mu1=mean(mor);
         %sigma1=std(mor);
         %mor= normpdf(mor,mu1,sigma1);
        
         %mu2=mean(con);
         %sigma2=std(con);
         %con= normpdf(con,mu2,sigma2);
        
        
        
        
        %mor=mor/sum(mor);
        %con=con/sum(con);
        
      
   
end

kld4=kld4/TT;

%%
%subplot(1,2,1);
save kld4.mat;
%plot(p(2:25),kld,'k-','LineWidth',2);
%title('KL');
%subplot(1,2,2);

plot(p(2:25),kld4,'k-','LineWidth',2);
title('KL-sym');






