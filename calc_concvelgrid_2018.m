% get ice concentration data on velocity spatial grid
% (output field has spatial dimensions of velocity, 
%  but temporal dimension of concentration)
%
% Paul Holland 22/1/15
%

close all
clear all

% ++++++++++++++++++++
% START MAIN POLE LOOP
% ++++++++++++++++++++
pole=1; %Antarctica

slat=-70;
cmer=0;

%disp('amsre/fowler: takes 5:30 (Antarctic) or 8:00 (Arctic)')
%load(strcat('amsre_2018_l1_pole',num2str(pole)));
%load(strcat('fowler_2018_l1_pole',num2str(pole)));
%outfname=strcat('amsre_fowlergrid_mean_l1_pole',num2str(pole));
%outfname=strcat('amsre_fowlergrid_2018_nanmean_l1_pole',num2str(pole));


load(strcat('CDR_2018_l1_pole',num2str(pole)));
load(strcat('kimura_2018_l1_pole',num2str(pole)));


%outfname=strcat('amsre_kimuragrid_mean_l1_pole',num2str(pole));
outfname=strcat('CDR_kimuragrid_2018_nanmean_l1_pole',num2str(pole));

%disp('bootstrap/kimura: takes ???? (Antarctic) or 2:00 (Arctic)')
%load(strcat('bootstrap_l1_pole',num2str(pole)));
%load(strcat('kimura_l1_pole',num2str(pole)));
%%outfname=strcat('bootstrap_kimuragrid_mean_l1_pole',num2str(pole));
%outfname=strcat('bootstrap_kimuragrid_nanmean_l1_pole',num2str(pole));

%disp('nasateam/fowler: takes 0:40 (Antarctic) or 1:10 (Arctic)')
%load(strcat('nasateam_l1_pole',num2str(pole)));
%load(strcat('fowler_l1_pole',num2str(pole)));
%outfname=strcat('nasateam_fowlergrid_mean_l1_pole',num2str(pole));
%outfname=strcat('nasateam_fowlergrid_nanmean_l1_pole',num2str(pole));

%disp('nasateam/kimura: takes 1:15 (Antarctic) or 2:00 (Arctic)')
%load(strcat('nasateam_l1_pole',num2str(pole)));
%load(strcat('kimura_l1_pole',num2str(pole)));
%outfname=strcat('nasateam_kimuragrid_mean_l1_pole',num2str(pole));
%outfname=strcat('nasateam_kimuragrid_nanmean_l1_pole',num2str(pole));

% ----------------------------------------------
% get ice concentration on velocity spatial grid
% ----------------------------------------------


% get source and target grids in polar stereo
[concx,concy]=ll2ps(conclat,conclon,'TrueLat',slat,'meridian',cmer);
[ivelx,ively]=ll2ps(vellat,vellon,'TrueLat',slat,'meridian',cmer);

% bin ice concentration data onto motion spatial grid at concentration times
concvelgrid=zeros(size(icexvel,1),size(icexvel,2),size(iceconc,3));

tic
% find valid concentration points surrounding each motion point
% (Arakawa A-grid)
concindex=cell(size(icexvel,1),size(icexvel,2));
for i=2:size(icexvel,1)-1
for j=2:size(icexvel,2)-1

  % correct, but would be clearer if i and j were labelled the other way around 
  ivelxmin=mean([ivelx(i,j) ivelx(i,j-1)]);
  ivelxmax=mean([ivelx(i,j) ivelx(i,j+1)]);
  ivelymin=mean([ively(i,j) ively(i-1,j)]);
  ivelymax=mean([ively(i,j) ively(i+1,j)]);
  concindex{i,j}=find((concx>ivelxmin)&(concx<ivelxmax)&...
                      (concy>ivelymin)&(concy<ivelymax));
end
end

% for each velocity point for each concentration day, get mean of valid concs
% (using mean/nanmean maximises/minimises effect of land NaNs in input field)
nt=size(iceconc,3);
nx=size(icexvel,1);
ny=size(icexvel,2);
parfor t=1:nt
for i=2:nx-1
for j=2:ny-1
  temp=iceconc(:,:,t);
  %concvelgrid(i,j,t)=mean(temp(concindex{i,j}));
  concvelgrid(i,j,t)=nanmean(temp(concindex{i,j}));
end
end
end
toc

% ----
% save
% ----

save(outfname,'conctime','vellat','vellon', ...
              'concvelgrid', ...
     '-v7.3');

% ++++++++++++++++++
% END MAIN POLE LOOP
% ++++++++++++++++++
