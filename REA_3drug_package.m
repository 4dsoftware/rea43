function [siai,bot] = REA_3drug_package(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Envelope Analysis (REA) for 3 drugs
%
% Di Du, Ph.D.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes
%
%Response envelope analysis (REA) is a tool to quantitatively determine 
%combination effects including synergy, additivity, and antagonism. 
%
%Inputs: 
%   data, the data formatted as the following:
%   column 1, concentrations of drug 1, in uM
%   column 2, concentrations of drug 2, in uM
%   column 3, concentrations of drug 3, in uM
%   column 4, corresponding survival rates in percentage

%Outputs:
%   siai is a vector that contains SI and AI.
%
%References:
%
%1. Du, D. et al, submited, 2017
%2. Wood, K. et al, PNAS, 2012
%--------------------------------------------------------------------------

ndose1 = length(unique(data(:,1))); %number of doses for drug 1
ndose2 = length(unique(data(:,2))); %number of doses for drug 2
ndose3 = length(unique(data(:,3))); %number of doses for drug 2
m1 = min(data(:,1));
m2 = min(data(:,2));
m3 = min(data(:,3));

%single drug response curve for drug 1
[dose1,id1] = sort(unique(data((data(:,2) == m2)&(data(:,3) == m3),1))); %concentrations of drug 1 
dose1 = dose1(2:end); %truncated concentrations of drug 1 to exclude 0
base1 = dose1(1)^2/dose1(3); 
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
surv1 = data((data(:,2) == 0)&(data(:,3) == 0),4)/100; %unaffected fraction of drug 1
surv1 = surv1(id1);  %sort surv1
surv1 = surv1(2:end); %truncated unaffected fraction of drug 1 to exclude 0

%single drug response curve for drug 2
[dose2,id2] = sort(unique(data((data(:,1) == m1)&(data(:,3) == m3),2))); %concentrations of drug 2 
dose2 = dose2(2:end); %truncated concentrations of drug 2 to exclude 0
base2 = dose2(1)^2/dose2(3);
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
surv2 = data((data(:,1) == 0)&(data(:,3) == 0),4)/100; %unaffected fraction of drug 2
surv2 = surv2(id2);  %sort surv2
surv2 = surv2(2:end); %truncated unaffected fraction of drug 2 to exclude 0

%single drug response curve for drug 3
[dose3,id3] = sort(unique(data((data(:,1) == m1)&(data(:,2) == m2),3))); %concentrations of drug 2 
dose3 = dose3(2:end); %truncated concentrations of drug 2 to exclude 0
base3 = dose3(1)^2/dose3(3);
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
surv3 = data((data(:,1) == 0)&(data(:,2) == 0),4)/100; %unaffected fraction of drug 2
surv3 = surv3(id3);  %sort surv2
surv3 = surv3(2:end); %truncated unaffected fraction of drug 2 to exclude 0

%combination data for drug 1 and 2
surv123 = zeros(length(dose1),length(dose2),length(dose3)); %unaffected fraction in a matrix
surv123_ar = zeros(length(dose1)*length(dose2)*length(dose3),4); %unaffected fraction in an array
q = 0; %index
for i = 1:length(dose1)
    for j = 1:length(dose2)
        for k = 1:length(dose3)
            q = q + 1;
            surv123(i,j,k) = data((data(:,1) == dose1(i) & data(:,2) == dose2(j) & ...
                data(:,3) == dose3(k)),4)/100; %unaffected fraction in a matrix
            surv123_ar(q,:) = [dose1(i) dose2(j) dose3(k) surv123(i,j,k)]; %unaffected fraction in an array
        end 
    end
end

%% 2 calculate Hill parameters for individual drugs
%first regression 
bot = min([surv1' surv2' surv3']); %initial guess of S0, assay background
c1 = package_envelope_hill([median(dose1) 2],[1;surv1],[base1;dose1],bot); %regression for drug 1
lam1 = c1(1); h1 = c1(2); %EC50 and Hill slope for drug 1
c2 = package_envelope_hill([median(dose2) 2],[1;surv2],[base2;dose2],bot); %regression for drug 2
lam2 = c2(1); h2 = c2(2); %EC50 and Hill slope for drug 2
c3 = package_envelope_hill([median(dose3) 2],[1;surv3],[base3;dose3],bot); %regression for drug 2
lam3 = c3(1); h3 = c3(2); %EC50 and Hill slope for drug 2

%second regression with adjusted S0, assay background
c = package_envelope_s0_3drug([lam1,h1,bot,lam2,h2,lam3,h3],surv1,dose1,surv2,dose2,surv3,dose3); %regression for S0
bot = c(3); %new S0, assay background
c1 = package_envelope_hill([median(dose1) 2],[1;surv1],[base1;dose1],bot); %regression for drug 1
lam1 = c1(1); h1 = (c1(2)); %EC50 and Hill slope for drug 1
c2 = package_envelope_hill([median(dose2) 2],[1;surv2],[base2;dose2],bot); %regression for drug 2
lam2 = c2(1); h2 = (c2(2)); %EC50 and Hill slope for drug 2
c3 = package_envelope_hill([median(dose3) 2],[1;surv3],[base3;dose3],bot); %regression for drug 2
lam3 = c3(1); h3 = c3(2); %EC50 and Hill slope for drug 2
%% 3 visualization
%initialize figure object
figure('position',[50 50 700 700]);

%aesthetics parameters
trim = 0.1; %trim the x scale and y scale this amount beyond the data points for 2-D representation
ft = 9; %font size

%define 3-d markers
asp1 = (log10(max(dose1))-log10(base1)+2*trim)/1.2; %define aspect ratio for axis 1
asp2 = (log10(max(dose2))-log10(base2)+2*trim)/1.2; %define aspect ratio for axis 2
asp3 = (log10(max(dose3))-log10(base3)+2*trim)/1.2; %define aspect ratio for axis 2

alpha = 0.7; %transparency
mksz = 12; %redefine marker size for this panel

%label the data points
ff = zeros(1,length(dose1)*length(dose2)*length(dose3)); %synergy labels
gg = zeros(1,length(dose1)*length(dose2)*length(dose3)); %antagonism labels
for k = 1:(length(dose1)*length(dose2)*length(dose3))
    [f1,f2] = package_loewe_3drug(surv123_ar(k,1),surv123_ar(k,2),surv123_ar(k,3),...
        1,bot,lam1,lam2,lam3,h1,h2,h3);
    f3 = package_bliss_3drug(surv123_ar(k,1),surv123_ar(k,2),surv123_ar(k,3),...
        1,bot,lam1,lam2,lam3,h1,h2,h3);
    if surv123_ar(k,4) < min([f1 f2 f3])
        ff(k) = 1;
    end
    if surv123_ar(k,4) > max([f1 f2 f3])
        gg(k) = 1;
    end
end

%find the largest island of synergy
ff = reshape(ff,length(dose1),length(dose2),length(dose3)); %reshape synergy label array to tensor
CC = bwconncomp(ff,26); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,id] = max(numPixels);
if isempty(id) == 0
    reg = CC.PixelIdxList{id}'; %data list corresponding to the largest island of antagonism
else
    reg = [];
end

%find the largest island of antagonism
gg = reshape(gg,length(dose1),length(dose2),length(dose3)); %reshape antagonism label array to tensor
CC = bwconncomp(gg,26); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,id] = max(numPixels);
if isempty(id) == 0
    reg2 = CC.PixelIdxList{id}'; %data list corresponding to the largest island of antagonism
else
    reg2 = [];
end

%examine if additivity region is larger than the other two
adv = 1 - gg - ff; %additivity region
adv = reshape(adv,length(dose1),length(dose2),length(dose3)); %reshape additivity label array to tensor
CC = bwconncomp(adv,6); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,id] = max(numPixels);
if isempty(id) == 0
    reg3 = CC.PixelIdxList{id}'; %data list corresponding to the largest island of antagonism
else
    reg3 = [];
end

%visualize the island
plot3(log10(surv123_ar(reg,1)),log10(surv123_ar(reg,2)),log10(surv123_ar(reg,3)),...
            'r.','markersize',mksz);
hold on
plot3(log10(surv123_ar(reg2,1)),log10(surv123_ar(reg2,2)),log10(surv123_ar(reg2,3)),...
            'b.','markersize',mksz);
hold off
box on
xlim([log10(min(dose1))-trim log10(max(dose1))+trim]);
ylim([log10(min(dose2))-trim log10(max(dose2))+trim]);
zlim([log10(min(dose3))-trim log10(max(dose3))+trim]);
xlabel('drug1','fontsize',ft); %x axis label
ylabel('drug2','fontsize',ft); %y axis label
zlabel('drug3','fontsize',ft); %z axis label
set(gca,'fontsize',ft);

%added 6/29/2018, linear scale correction
d1p = [dose1(1)^2/dose1(2);dose1];
d2p = [dose2(1)^2/dose2(2);dose2];
d3p = [dose3(1)^2/dose3(2);dose3];

nn = 0;
area1 = 0; %volume instead
if isempty(reg) == 0
dif = zeros(1,length(reg));
for k = reg
    nn = nn + 1;
    [f1,f2] = package_loewe_3drug(surv123_ar(k,1),surv123_ar(k,2),surv123_ar(k,3),...
        1,bot,lam1,lam2,lam3,h1,h2,h3);
    f3 = package_bliss_3drug(surv123_ar(k,1),surv123_ar(k,2),surv123_ar(k,3),...
        1,bot,lam1,lam2,lam3,h1,h2,h3);
    flow = min([f1 f2 f3]);
    id1 = find(d1p == surv123_ar(k,1));
    id2 = find(d2p == surv123_ar(k,2));
    id3 = find(d3p == surv123_ar(k,3));
    dif(nn) = abs(surv123_ar(k,4) - flow)*log(d1p(id1)/d1p(id1-1))*...
        log(d2p(id2)/d2p(id2-1))*log(d3p(id3)/d3p(id3-1));
    area1 = area1 + log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1))*log(d3p(id3)/d3p(id3-1));
end
else 
dif = 0;   
end

%difference of the measured data from the envelope for antagonism
nn = 0;
area2 = 0;
if isempty(reg2) == 0
dif2= zeros(1,length(reg2));
for k = reg2
    nn = nn + 1;
    [f1,f2] = package_loewe_3drug(surv123_ar(k,1),surv123_ar(k,2),surv123_ar(k,3),...
        1,bot,lam1,lam2,lam3,h1,h2,h3);
    f3 = package_bliss_3drug(surv123_ar(k,1),surv123_ar(k,2),surv123_ar(k,3),...
        1,bot,lam1,lam2,lam3,h1,h2,h3);
    flow = max([f1 f2 f3]);
    id1 = find(d1p == surv123_ar(k,1));
    id2 = find(d2p == surv123_ar(k,2));
    id3 = find(d3p == surv123_ar(k,3));
    dif2(nn) = abs(surv123_ar(k,4) - flow) * log(d1p(id1)/d1p(id1-1))*...
        log(d2p(id2)/d2p(id2-1))*log(d3p(id3)/d3p(id3-1));
    area2 = area2 + log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1))*log(d3p(id3)/d3p(id3-1));
end
else 
dif2 = 0;  
end

%difference of the measured data from the envelope for antagonism
nn = 0;
area3 = 0;
if isempty(reg3) == 0
for k = reg3
    nn = nn + 1;
    id1 = find(d1p == surv123_ar(k,1));
    id2 = find(d2p == surv123_ar(k,2));
    id3 = find(d3p == surv123_ar(k,3));
    area3 = area3 + log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1))*log(d3p(id3)/d3p(id3-1));
end
end

if area3 >= max(area1,area2)
    dif = 0;
    dif2 = 0;
end

tarea = 0;
for i = 2:length(d1p)
    for j = 2:length(d2p)
        for k = 2:length(d3p)
            tarea = tarea + log(d1p(i)/d1p(i-1))*log(d2p(j)/d2p(j-1))*log(d3p(k)/d3p(k-1));
        end
    end
end
si = sum(dif)/tarea; %SI
ai = sum(dif2)/tarea; %AI

%view(150,20);
siai = [si ai]; %final output
