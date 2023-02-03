%%
%%  Spectral-Structured-Sparse-Bayesian-Learning Main function
%%
% Description
% Inverse solutions are via Spectral Structured Sparse Bayesian Learning (ssSBL),
% the extension of the nominal SBL to the cross-spectrum.
% ssSBL results in two orders of magnitude less distortions than state-of-the-art methods
% for the standard ESI setup, which considers large-scale networks and
% low-density EEG (10-20 system) and is compared to high-density MEG and ECoG.
%
%
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Ariosky Areces Gonzalez
% - Pedro A. Valdes Sosa
%
% Updated: Jan 26, 2023


disp("=====================================================================");
disp("    <<<<< Spectral-Structured-Sparse-Bayesian-Learning >>>>>");
disp("=====================================================================");
disp("-->> Starting process analysis");


addpath('data/');
addpath('functions/');

disp("-->> Getting the data ready. Please wait.");
if(~isfile('data/Svv.mat') || ~isfile('data/Lvj.mat'))
    if ismac
        try
            system('zip -s 0 data/Svv/Svv_files.zip --out Svv.zip');
            system('zip -s 0 data/Lvj/Lvj_files.zip --out Lvj.zip');
        catch
            warning("Please Download The Unarchiver from https://apps.apple.com/us/app/the-unarchiver/id425424353?mt=12");
            return;
        end
    elseif isunix
        
        status = system('7z x data/Svv/Svv_files.zip');
        system('7z x data/Lvj/Lvj_files.zip');
        if(~isequal(status,0))
            warning("Please Download 7-zip from https://www.7-zip.org/");
            return;
        end
    elseif ispc
        try
            system(['"C:\Program Files\7-Zip\7z.exe" e -o"', pwd, '" "data\Svv\Svv_files.zip"']);
            system(['"C:\Program Files\7-Zip\7z.exe" e -o"', pwd, '" "data\Lvj\Lvj_files.zip"']);
        catch
            warning("Please Download 7-zip from https://www.7-zip.org/");
            return;
        end
    else
        disp('Platform not supported')
    end
    unzip('Svv.zip');
    movefile('Svv.mat','data/Svv.mat');
    delete('Svv.zip');
    
    unzip('Lvj.zip');
    movefile('Lvj.mat','data/Lvj.mat');
    delete('Lvj.zip');
end
load('data/mycolormap.mat');
load('data/labels.mat');
load('data/Lvj.mat');
load('data/Svv.mat');
Sc    = load('data/Sc.mat');
%
%% Estimating source directions
% This step estimates for each source the optimal orientation at all 
% frequencies. The orientation is an important anatomical constraint 
% for the sources in order to lock the phase of their components in 
% the x, and y, and z axis. Defining this constraint is very common 
% the Grid orientation but this is not strictily valid.
% Therefore, we propose to ssSBL user to first implementing a joint-MAP at 
% all frequencies to estimate the source orientations.
%
%% Initial ssSBL solution for all frequencies using 2D rotational invariance priors around surface orientation 
Svv_mean       = squeeze(mean(Svv(:,:,:),3)); % averaging cross-spectrum across freqeuncies
VN             = blk_diag(Sc.VertNormals', 1); % 2D rotational invariance prior 
Lvj3D          = Lvj; % 3D Lead Field
Lvj            = Lvj3D*VN; % projecting the 3D Lead Field
[Tjv,s2jj]     = ssSBL(Lvj,Svv_mean); % joint-MAP at all frequencies 
%
%% Estimating source orientations
Tjv3D          = VN*Tjv; % recovering 3D inverse operator
SvvTvj3D       = Svv_mean*Tjv3D'; % recovering 3D source spectrum
for g = 1:size(Lvj3D,2)
    s2jj3D(g) = abs(Tjv3D(g,:)*SvvTvj3D(:,g));
end
SourceNormals  = reshape(s2jj3D,3,length(Lvj3D)/3)'; % recovering source orientations
SourceNormals  = sqrt(SourceNormals);
SourceNormals  = SourceNormals./repmat(sqrt(sum(SourceNormals.^2,2)),1,3);
SourceNormals  = repmat(sign(sum(SourceNormals.*Sc.VertNormals,2)),1,3).*SourceNormals;
%
%% Final ssSBL solution for a freqeuncy (alpha peak in this example) using 2D rotational invariance priors around source orientation
SN         = blk_diag(SourceNormals', 1); % 2D rotational invariance prior 
Lvj        = Lvj3D*SN; % projecting the 3D Lead Field
[Tjv,s2jj] = ssSBL(Lvj,squeeze(Svv(:,:,105))); % joint-MAP at a frequency
%
%% Plotting results
figure_cortical_map = figure('Color','w','Name','cortical spectral topography','NumberTitle','off'); 
plot_sources(s2jj/max(s2jj),Sc,cmap); % ploting spectrum
%
%%
function plot_sources(sources,Sc,cmap)
sources                 = sqrt(sources);
smoothValue             = 0.86;
SurfSmoothIterations    = 30;
Vertices                = tess_smooth(Sc.Vertices, smoothValue, SurfSmoothIterations, Sc.VertConn, 1);
patch('Faces',Sc.Faces,'Vertices',Vertices,'FaceVertexCData',Sc.SulciMap*0.06+...
    log(1+sources),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
set(gca,'xcolor','w','ycolor','w','zcolor','w');
az = 0; 
el = 0;
view(az,el);
rotate3d on;
colormap(gca,cmap);
title('cortical spectral topography','Color','k','FontSize',16);
end