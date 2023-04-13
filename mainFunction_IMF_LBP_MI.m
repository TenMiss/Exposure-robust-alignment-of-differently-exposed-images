% Image alignment + linear optimization
%
% Wu Shiqian. 11 Sep 2010

%% specify the directory that contains your range of differently exposed
clear all; close all; clc
[filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Open images');
if isequal([filename,pathname],[0,0]);
    return
end
thumbs_file = fullfile(pathname,'thumbs.db');
delete(thumbs_file)
[filenames, exposures, numExposures] = readDir(pathname);
% Pre-define parameters
a = imread(filenames{1});
[row,col,h] = size(a);
R_W = zeros(256,1);
for i=1:127
    R_W(i) = i;
end
for i=128:255
    R_W(i)=255-i;
end
I = cell(1,numExposures);
for i=1:numExposures
    a = imread(filenames{i});
    if size(a,3)==3
        a = rgb2gray(a);
    end
    I{i} = a;    
end
Y_NorRef = zeros(row,col);
Y_NorCur = zeros(row,col);
IMFRef2Cur = zeros(256,1);
IMFCur2Ref = zeros(256,1);
shifts = zeros(numExposures,2);
alfa = zeros(1,numExposures);
%% Image alignment
Y_Ref = I{1};  % The brightest image is served as reference
T_Ref = exposures(1);
Ref = imread(filenames{1});
PM = cell(1,2);
MI = zeros(1,numExposures);
tic
for i = 2:numExposures
    Y_Cur = I{i};
    T_Cur = exposures(i);     
    [IMFRef2Cur, IMFCur2Ref]=HistogramBasedIMF(Y_Ref,Y_Cur,row,col,IMFRef2Cur, IMFCur2Ref, T_Ref, T_Cur);
    a_Ref = sum(sum(R_W(Y_Ref+1)));
    a_Cur = sum(sum(R_W(Y_Cur+1)));
    if a_Cur>a_Ref
        PM{1} = Y_Ref;
        PM{2} = UniDirectionalNormalization(Y_Cur, row, col, IMFCur2Ref);
    else
        PM{1} = UniDirectionalNormalization(Y_Ref, row, col, IMFRef2Cur);
        PM{2} = Y_Cur;
    end
    [sft, beta] = imfLBPalignment(PM);
    shifts(i,:) = sft
    alfa(i) = beta
    %%%align two images by using the computed parameters    
    tmp = imread(filenames{i});
    if any(abs(shifts(i,:))>0.25)
        xx = floor(shifts(i,2)*4+0.5)/4;%%%complexity can be reduced by changing 4 as 2
        yy = floor(shifts(i,1)*4+0.5)/4;%%%complexity can be reduced by changing 4 as 2
        tmp = shift(tmp,xx,yy);
    end
    if abs(alfa(i))>0.1
        tmp = imrotate(tmp,-alfa(i),'bicubic','crop');
    end
    MI(i) = mutualInformation(Ref,tmp,50)
    figure, imshow(Ref);
    hold on;
    h2 = imshow(tmp);
    title('overlay two aligned images'); 
    alpha(h2,0.6);
    [pathstr, namestr, extstr] = fileparts(filenames{i}); 
    saveas(h2,['Overlay_' num2str(i) '.tiff'],'tiff');
    hold off
end
Time = toc
fprintf('Finish image alignment\n')




