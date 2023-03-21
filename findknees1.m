% find knees
function [knees]=findknees1(non_fit)
x=non_fit(:,1);
y=non_fit(:,2);
angle=zeros(size(non_fit,1)-1,1);
for i=2:1:length(x)-1
    angle1=atan(abs((y(i-1)-y(i))/(x(i-1)-x(i))));
    angle2=atan(abs((x(i+1)-x(i))/(y(i+1)-y(i))));
    angle(i)=angle1+angle2;
end
[~,IX]=max(angle);
knees=IX(1);
end