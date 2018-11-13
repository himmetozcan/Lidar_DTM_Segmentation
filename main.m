clear
close all
clc

%% load data

load utah1;
% load utah5;
% load utah7;
% load utah9;

I = dsm;
%% preprocessing: small objects

th = 1; % meters
ws = 5;
Is = imopen(I, strel('disk', ws));
Os = (I - Is) > th;

%% method
tic
Xkj = votesegdsm( Is );

% [Xkj, kernelMap, Xe] = votesegdsm( Is, 'kernelsigma', 5, 'votemaxwinsize', 10 );
% figure,imagesc(kernelMap);axis image
% figure,imagesc(Is);axis image
% hold on,plot(Xe(:,1),Xe(:,2),'k^','MarkerFaceColor','r')
toc
%% Results

Objs = Xkj + 2 * (Os > 0);
Objs = Objs > 0;

trueResults=(dsm-dtm_NFS)>1;

tmp=2*Objs+trueResults;

r=ones(size(dsm)); g=ones(size(dsm)); b=ones(size(dsm));
g(tmp==3)=1; r(tmp==3)=.5; b(tmp==3)=0;
r(tmp==2)=1; g(tmp==2)=0; b(tmp==2)=0;
b(tmp==1)=1; r(tmp==1)=0; g(tmp==1)=0;
im=[];
im(:,:,1)=r; im(:,:,2)=g; im(:,:,3)=b;

figure,imshow(im,[])

disp(' ')
disp('Calculating performance ... ')
performances=calculateperformances(Objs(:), trueResults(:));

disp(['Kappa(%): ', num2str(performances.kappa)])
disp(['Total Error(%): ', num2str(performances.TE), ', Type-I Error(%): ', num2str(performances.TI),', Type-II Error(%): ', num2str(performances.TII)])

