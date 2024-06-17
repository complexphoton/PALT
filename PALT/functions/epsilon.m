function y = epsilon(x,na,Lft,Rt)
Passive=[Lft,[1;na;0],Rt];
d=Passive(1,:); d0=sum(Lft(1,:)); x = x+d0;
X = [0,cumsum(d)]; EPSILON=[Passive(2,:).^2];
y=ones(size(x));

for n=1:length(X)-1
y(x>X(n)&x<=X(n+1))=EPSILON(n);
end

