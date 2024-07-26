function p = dumbSpline(x1,x2,y1,y2,m1,m2)


a = (m1 + m2 - 2*(y2-y1)/(x2-x1))/(x1-x2)^2;
b = (m2-m1)/(2*(x2-x1))-1.5*(x1+x2)*a;
c = m1 - 3*x1^2*a-2*x1*b;
d = y1- x1^3*a-x1^2*b-x1*c;

p =[a b c d];
end

