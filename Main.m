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

load('data/Lvj.mat');
load('data/Svv.mat');
load('data/mycolormap.mat');
load('data/labels.mat');
Sc = load('data/Sc.mat');

%% Estimating source directions
% This step estimates for each source the optimal orientation at all 
% frequencies. The orientation is an important anatomical constraint 
% for the sources in order to lock the phase of their components in 
% the x, and y, and z axis. Defining this constraint is very common 
% the Grid orientation but this is not strictily valid.
% Therefore, we propose to ssSBL user to first implementing a joint-MAP at 
% all frequencies to estimate the source orientations.

Svv_mean       = squeeze(mean(Svv(:,:,:),3)); % averaging cross-spectrum
[Tjv3,s2jj3]   = sSSBLpp(Lvj,Svv_mean); % joint-MAP at all freqeuncies
s2jj3          = reshape(s2jj3,3,length(Lvj)/3)'; % Obtaining source directions
plot_sources(sqrt(sum(s2jj3,2)),Sc,cmap)
hold on
SourceNormals  = sqrt(s2jj3);
SourceNormals  = SourceNormals./repmat(sqrt(sum(SourceNormals.^2,2)),1,3);
VertNormals    = Sc.VertNormals; % Projecting SourceOrient to GridOrient
SourceNormals  = sign(SourceNormals.*VertNormals).^SourceNormals;
quiver3(Sc.Vertices(:,1),Sc.Vertices(:,2),Sc.Vertices(:,3),SourceNormals(:,1),SourceNormals(:,2),SourceNormals(:,3),'Color','b','AutoScaleFactor',3);
quiver3(Sc.Vertices(:,1),Sc.Vertices(:,2),Sc.Vertices(:,3),VertNormals(:,1),VertNormals(:,2),VertNormals(:,3),'Color','k','AutoScaleFactor',3);
SourceNormals  = blk_diag(SourceNormals', 1);
Lvj            = Lvj

%% Calling Main fuction
%%
% 
[Tjv,s2jj3]   = sSSBLpp(Lvj,squeeze(Svv(:,:,105)));
s2jj                                    = sum(reshape(abs(s2jj3),3,length(Lvj)/3),1)';
stat                                    = sqrt(2)*s2jj./sqrt(var(s2jj));
indms                                   = find(stat > 1);
sps2jj                                  = zeros(length(stat),1);
sps2jj(indms)                           = s2jj(indms);
%%
%% Plotting results
%%
figure_stat             = figure('Color','w','Name','ssSBL-stat-alpha-band','NumberTitle','off'); 
h                                       = histogram(stat,1000,'Normalization','pdf');
bins                                    = h.BinEdges;
pdf                                     = (1/(sqrt(2*pi)))*(bins.^(-1/2)).*exp(-bins/2);
hold on;
plot(bins,pdf,'LineWidth',2,'Color','r')
legend('Empirical','Chi2');
title('ssSBL-stat','Color','k','FontSize',16);

%%
function plot_sources(sources,Sc,cmap)
sources                 = sources/max(sources(:));
figure_activation       = figure('Color','w','Name','ssSBL-activation-alpha-band','NumberTitle','off'); hold on;
% smoothValue             = 0.66;
% SurfSmoothIterations    = 10;
% Vertices                = tess_smooth(Sc.Vertices, smoothValue, SurfSmoothIterations, Sc.VertConn, 1);
patch('Faces',Sc.Faces,'Vertices',Sc.Vertices,'FaceVertexCData',Sc.SulciMap*0.06+...
    log(1+sources),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
set(gca,'xcolor','w','ycolor','w','zcolor','w');
az = 0; el = 0;
view(az,el);
rotate3d on;
colormap(gca,cmap);
title('ssSBL-activation','Color','k','FontSize',16);
end