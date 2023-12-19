function [mu,q,E] = FindElasticMean(Data)

close all;
Niter=30;
n=size(Data,3);
figure(1); clf; 
for i=1:n
    X = ReSampleCurve(Data(:,:,i),100);
    q(:,:,i) = curve_to_q(X);
end

del = 0.5;
mu = sum(q,3)/n;
mu=mu/sqrt(InnerProd_Q(mu,mu));
mu=ProjectC(mu);

for iter =1:Niter
    
    vm = 0;    
    for i=1:n
        [iter i]
        v = ElasticShootingVector(mu,q(:,:,i),1);
        vm = vm + v;
        
    end
       vm = vm/n;
       
       E(iter) = norm(vm,'fro');
       mu = ElasticShooting(mu,del*vm);
       
       iter=iter+1;
       E
    
end

figure(101); 
plot(E);

figure(21); clf; 
p = q_to_curve(mu);
plot(p(1,:),p(2,:),'LineWidth',3);
axis equal;