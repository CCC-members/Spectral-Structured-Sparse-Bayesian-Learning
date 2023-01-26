

%%
%%  Spectral-Structured-Sparse-Bayesian-Learning Main function
%%
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
        end
    elseif isunix
        system('7z x data/Svv/Svv_files.zip');
        system('7z x data/Lvj/Lvj_files.zip');
    elseif ispc
        try
            system(['"C:\Program Files\7-Zip\7z.exe" e -o"', pwd, '" "data\Svv\Svv_files.zip"']);
            system(['"C:\Program Files\7-Zip\7z.exe" e -o"', pwd, '" "data\Lvj\Lvj_files.zip"']);
        catch
            warning("Please Download 7-zip from https://www.7-zip.org/");
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

%%
%% Calling Main fuction
%%
[s2j,sigma2j,Tjv,Svv,scaleJ,scaleLvj]   = sSSBLpp(Svv(:,:,1),Lvj);
s2j                                     = sum(reshape(abs(s2j),3,length(Ke)/3),1)';
stat                                    = sqrt(2)*s2j./sqrt(var(s2j));
indms                                   = find(stat > sssblpp_th);
h                                       = histogram(stat,1000,'Normalization','pdf');
bins                                    = h.BinEdges;
pdf                                     = (1/(sqrt(2*pi)))*(bins.^(-1/2)).*exp(-bins/2);
J                                       = s2j;
J                                       = J*scaleSvv/scaleKe^2;
Jsp                                     = zeros(length(stat),1);
Jsp(indms)                              = J(indms);

%%
%% Plotting results
%%
figure_stat             = figure('Color','w','Name','ssSBL-stat-alpha-band','NumberTitle','off'); hold on;
plot(bins,pdf,'LineWidth',2,'Color','r')
legend('Empirical','Chi2');
title('ssSBL-stat','Color','k','FontSize',16);

sources_iv              = sqrt(abs(J));
sources_iv              = sources_iv/max(sources_iv(:));
figure_activation       = figure('Color','w','Name','ssSBL-activation-alpha-band','NumberTitle','off'); hold on;
smoothValue             = 0.66;
SurfSmoothIterations    = 10;
Vertices                = tess_smooth(Sc.Vertices, smoothValue, SurfSmoothIterations, Sc.VertConn, 1);
patch('Faces',Sc.Faces,'Vertices',Vertices,'FaceVertexCData',Sc.SulciMap*0.06+...
    log(1+sources_iv),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
set(gca,'xcolor','w','ycolor','w','zcolor','w');
az = 0; el = 0;
view(az,el);
rotate3d on;
colormap(gca,cmap);
title('ssSBL-activation','Color','k','FontSize',16);
    
disp("=====================================================================");
disp("-->> Process finished.");
disp("=====================================================================");