function Performances=calculateperformances(nonGround, trueResults)

erTable=crosstab(nonGround(:),trueResults(:));
a=erTable(1,1);
b=erTable(1,2);
c=erTable(2,1);
d=erTable(2,2);
ac=a+c; bd=b+d;
ab=a+b; cd=c+d;
e=a+b+c+d;
h=ac/e; i=bd/e;
f=ab/e; g=cd/e;
TI=b/ab;
TII=c/cd;
TE=(b+c)/e;
Pra=(a+d)/e;
Pre=f*h+g*i;
kappa=(Pra-Pre)/(1-Pre);

Performances.crossTable=erTable;
Performances.ab=ab;
Performances.f=f*100;
Performances.cd=cd;
Performances.g=g*100;
Performances.ac=ac;
Performances.bd=bd;
Performances.h=h*100;
Performances.i=i*100;
Performances.e=e;
Performances.TI=100*TI;
Performances.TII=100*TII;
Performances.TE=100*TE;
Performances.Pra=100*Pra;
Performances.Pre=100*Pre;
Performances.kappa=100*kappa;
