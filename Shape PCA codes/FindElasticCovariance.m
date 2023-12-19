function U = FindElasticCovariance(mu,q)

close all;

reparamFlag = 1;
[a,b,n]=size(q);

muc=q_to_curve(mu);
figure(1), clf;
plot(muc(1,:),muc(2,:),'LineWidth',3)
axis equal off; view([1 90]);

for i=1:n
    i
    tmp = ElasticShootingVector(mu,q(:,:,i),reparamFlag);
    VV(i,1:b) = tmp(1,:);
    VV(i,b+1:2*b) = tmp(2,:);
end

K = cov(VV);
[U,S,V] = svd(K);

T = length(VV(1,:))/2;
s = sqrt(diag(S));

i=1;
figure(11); clf; hold on;
for t=-3:1:3
    u = t*s(i)*U(:,i)'/2;
    UU(1,1:T) = u(1:T);
    UU(2,1:T) = u(T+1:2*T);
    UUn=UU;
    
    q = ElasticShooting(mu,UUn);
   
    
    p = q_to_curve(q);
    
    
    if (t==0)
        plot(t*0.32+p(1,:),p(2,:),'r','LineWidth',3);
    else
        plot(t*0.32+p(1,:),p(2,:),'LineWidth',3);
    end
    axis equal off; view([1 90])
end

i = 2;
figure(12); clf; hold on;
for t=-3:1:3
    u = t*s(i)*U(:,i)'/2;
    UU(1,1:T) = u(1:T);
    UU(2,1:T) = u(T+1:2*T);
    UUn=UU;
    
    q = ElasticShooting(mu,UUn);   
    p = q_to_curve(q);
    

    if (t==0)
        plot(t*0.31+p(1,:),p(2,:),'r','LineWidth',3);
    else
        plot(t*0.31+p(1,:),p(2,:),'LineWidth',3);
    end
    axis equal off; view([1 90])
end

i = 3;
figure(13); clf; hold on;
for t=-3:1:3
    u = t*s(i)*U(:,i)'/2;
    UU(1,1:T) = u(1:T);
    UU(2,1:T) = u(T+1:2*T);
    UUn=UU;
    
    q = ElasticShooting(mu,UUn);   
    p = q_to_curve(q);
    
    if (t==0)
        plot(t*0.31+p(1,:),p(2,:),'r','LineWidth',3);
    else
        plot(t*0.31+p(1,:),p(2,:),'LineWidth',3);
    end
    axis equal off; view([1 90])
end

mun=mean(VV);
idx=1;
figure(17); clf; hold on;
while(idx<7)
        v = mun + zeros(1,2*T);
        for i=1:b
            % Gaussian
            c(i) = randn*s(i);
            v = v + c(i)*U(:,i)';
        end
        
        vn(1,1:T) = v(1:T);
        vn(2,1:T) = v(T+1:2*T);
        qsamp = ElasticShooting(mu,vn);
        psamp = q_to_curve(qsamp);

        z = plot(psamp(1,:)+.31*idx,psamp(2,:),'LineWidth',3);
        axis equal off; view([1 90]);
        idx=idx+1;
end