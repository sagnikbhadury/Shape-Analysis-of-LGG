function d = GeodesicElasticClosed_fast(p1,p2)

% input p1 and p2 as 2xn matrices
% to turn off figures set figs=0

figs=0;

stp = 6;

if figs
figure(1); clf; hold on;
plot(p1(1,:),p1(2,:),'b','LineWidth',2);
plot(p2(1,:),p2(2,:),'r','LineWidth',2);
axis equal;
axis xy off;
end

p1 = ReSampleCurve(p1,100);
p2 = ReSampleCurve(p2,100);

q1 = curve_to_q(p1);
q2 = curve_to_q(p2);

[q2n,R] = Find_Rotation_and_Seed_unique_fast(q1,q2,1);
q2n = q2n/sqrt(InnerProd_Q(q2n,q2n));
q2n = ProjectC(q2n);
p2 = q_to_curve(q2);
p2n = R*p2;

d = acos(InnerProd_Q(q1,q2n));

if figs
alpha = geodesic_sphere_Full(q1,q2n,stp);
Path_Plot(alpha,p2n,10,'b',[73,6]);
axis xy; view([1 90]);
X1=q_to_curve(q1);
X2=q_to_curve(q2n);
figure(2),clf,hold on;
plot([X1(1,:) X1(1,1)],[X1(2,:) X1(2,1)],'-','LineWidth',2);
plot([X1(1,:) X1(1,1)],[X1(2,:) X1(2,1)],'*','LineWidth',2);
plot(X1(1,1),X1(2,1),'*r','LineWidth',10);
plot(X1(1,25),X1(2,25),'*g','LineWidth',10);
plot(X1(1,50),X1(2,50),'*k','LineWidth',10);
plot(X1(1,75),X1(2,75),'*c','LineWidth',10);
axis equal off; view([1 90])
figure(3),clf,hold on;
plot([X2(1,:) X2(1,1)],[X2(2,:) X2(2,1)],'-','LineWidth',2);
plot([X2(1,:) X2(1,1)],[X2(2,:) X2(2,1)],'*','LineWidth',2);
plot(X2(1,1),X2(2,1),'*r','LineWidth',10);
plot(X2(1,25),X2(2,25),'*g','LineWidth',10);
plot(X2(1,50),X2(2,50),'*k','LineWidth',10);
plot(X2(1,75),X2(2,75),'*c','LineWidth',10);
axis equal off; view([1 90])
end