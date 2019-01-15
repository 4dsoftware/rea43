function f = package_bliss_3drug(x1,x2,x3,a,b,lam1,lam2,lam3,h1,h2,h3)
f = b+(a-b)./(1.+(x1/lam1).^h1)./(1.+(x2/lam2).^h2)./(1.+(x3/lam3).^h3);