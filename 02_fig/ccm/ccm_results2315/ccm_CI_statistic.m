close all; clear ;clc;

A=csvread('GR_xmap_NPI.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('GR_xmap_NPI_rho_statistics.csv',data_a);
%% 2
close all; clear ;clc;
A=csvread('NPI_xmap_GR.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('NPI_xmap_GR_rho_statistics.csv',data_a);
%% 3
close all; clear ;clc;
A=csvread('GR_xmap_SI.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('GR_xmap_SI_rho_statistics.csv',data_a);
%% 4
close all; clear ;clc;
A=csvread('SI_xmap_GR.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('SI_xmap_GR_rho_statistics.csv',data_a);
%% 5
close all; clear ;clc;
A=csvread('GR_xmap_Vero.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('GR_xmap_Vero_rho_statistics.csv',data_a);
%% 6
close all; clear ;clc;
A=csvread('Vero_xmap_GR.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('Vero_xmap_GR_rho_statistics.csv',data_a);
%% 7
close all; clear ;clc;
A=csvread('SI_xmap_NPI.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('SI_xmap_NPI_rho_statistics.csv',data_a);
%% 8
close all; clear ;clc;
A=csvread('NPI_xmap_SI.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('NPI_xmap_SI_rho_statistics.csv',data_a);
%% 9
close all; clear ;clc;
A=csvread('SI_xmap_Vero.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=mean(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('SI_xmap_Vero_rho_statistics.csv',data_a);
%% 10
close all; clear ;clc;
A=csvread('Vero_xmap_SI.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('Vero_xmap_SI_rho_statistics.csv',data_a);
%% 11
close all; clear ;clc;
A=csvread('TR_xmap_SI.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('TR_xmap_SI_rho_statistics.csv',data_a);
%% 12
close all; clear ;clc;
A=csvread('SI_xmap_TR.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('SI_xmap_TR_rho_statistics.csv',data_a);
%% 13
close all; clear ;clc;
A=csvread('TR_xmap_GR.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('TR_xmap_GR_rho_statistics.csv',data_a);
%% 14
close all; clear ;clc;
A=csvread('GR_xmap_TR.csv');
a=size(A,1)/100;
data_a=zeros(a,4);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
for jj=1:a
    data_a(jj,1)=A((jj-1)*100+1,1);
    data_a(jj,3)=median(A((jj-1)*100+1:(jj)*100,2));
    data_a(jj,2)=prctile((A((jj-1)*100+1:(jj)*100,2)),97.5);
    data_a(jj,4)=prctile((A((jj-1)*100+1:(jj)*100,2)),2.5);
end
csvwrite('GR_xmap_TR_rho_statistics.csv',data_a);




