function [f1,f2] = package_loewe_3drug(x1,x2,x3,a,b,lam1,lam2,lam3,h1,h2,h3)
all3 = [x1,lam1,h1;x2,lam2,h2;x3,lam3,h3];
sort3 = table2array(sortrows(array2table(all3),3,'ascend')); %sort based on Hill coefficient
hmin = sort3(1,3);
hmax = sort3(end,3);
f1 = b+(a-b)./(1.+(sum((sort3(:,1)./sort3(:,2)).^(sort3(:,3)/hmin))).^hmin);
f2 = b+(a-b)./(1.+(sum((sort3(:,1)./sort3(:,2)).^(sort3(:,3)/hmax))).^hmax);
