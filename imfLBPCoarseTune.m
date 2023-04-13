function [rotation, translation,sigbeta,Rshift,Cshift] = imfLBPCoarseTune(T,Unregistered,ANG2,rmaxshift,cmaxshift,THbeta,DH,dbeta)
% Detect the initial rotation and translation parameters
% Input: 
%       'im' is the input images, im{1} is reference and im{2}is the 
%            unregistered image
%       'ANG2' is a 2-element vector which is the angle range to be searched
%       'rmaxshift' is a scalar which is the row shift range to be searched 
%       'cmaxshift' is a scalar which is the column shift range to be searched       
% Output:
%       rotation is a 2-element vector: [row_alfa, col_alfa]
%       translation = [row_shift,column_shift]
%       sigbeta = 1 indicates the determined angle is sure

% Copyright: Wu Shiqian
% Institute for Infocomm Research
% 8 August 2010                                                           
if nargin<6
    THbeta = 2;
    DH =1;
    dbeta = 1;      
end
H_TH = 256;
T1 = T;
U1 = Unregistered;
if ~isa(T1, 'double') 
    T1 = double(T1);
    U1 = double(U1); 
end
[Tr,Tc] = size(T1);
angs = ANG2(1):dbeta:ANG2(2);
rlags = -rmaxshift:rmaxshift;
clags = -cmaxshift:cmaxshift;
rnlags = length(rlags);
cnlags = length(clags);
nAngs = length(angs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FT = LBPCode(T1);%T1;
%FU = LBPCode(U1); %U1;
FT = LBP8Feature(T1,0);
FU = LBP8Feature(U1,0);
FTHC0 = histc(FT,1:H_TH,1);
FTHR0 = histc(FT,1:H_TH,2);
FTHC = FTHC0;
FTHR = FTHR0;
for WH =1:DH
    ctmp = circshift(FTHC0,[0,WH]);    
    rtmp = circshift(FTHR0,[WH,0]);    
    FTHC = FTHC + ctmp;
    FTHR = FTHR + rtmp;
end
FTHC(:,1:DH) = [];
FTHR(1:DH,:)= [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cshift = repmat(900,[nAngs,3]);
Rshift = repmat(900,[nAngs,3]);
for i = 1:nAngs
    beta = angs(i);
    tmp = imrotate(FU,beta,'nearest'); 
    [r1,c1] = size(tmp);
    deltC = floor(0.5*(c1-Tc));
    deltR = floor(0.5*(r1-Tr));
    FU1 = tmp(deltR+1:deltR+Tr,deltC+1:deltC+Tc); 
    FUHC0 = histc(FU1,1:H_TH,1);
    FUHR0 = histc(FU1,1:H_TH,2);
    FUHC = FUHC0;
    FUHR = FUHR0;
    for WH =1:DH
        ctmp = circshift(FUHC0,[0,WH]);
        rtmp = circshift(FUHR0,[WH,0]);        
        FUHC = FUHC + ctmp;
        FUHR = FUHR + rtmp;
    end
    FUHC(:,1:DH) = [];
    FUHR(1:DH,:)= [];
    Csearch = zeros(1,cnlags);
    Rsearch = zeros(1,rnlags);    
    for j = 1:cnlags
        delt = clags(j);
        DHj = abs(circshift(FUHC,[0,delt])-FTHC);
        PDH = DHj(:,cmaxshift+1:end-cmaxshift);
        Csearch(j) = sum(PDH(:));
    end
    CminV = min(Csearch);
    num = find(Csearch == CminV); 
    
    if length(num)==1
        Cshift(i,:) = [beta, clags(num), CminV];    
    else
        mNo = fix(0.5*(length(num)+1));
        MN = num(mNo);
        Cshift(i,:) = [beta,clags(MN), CminV];  
        
    end
    %%%%% find shifts in row (Y) direction
    for j = 1:rnlags
        delt = rlags(j);
        DHj = abs(circshift(FUHR,[delt,0])-FTHR);
        PDH = DHj(cmaxshift+1:end-cmaxshift,:);
        Rsearch(j) = sum(PDH(:));               
    end
    RminV = min(Rsearch);
    num = find(Rsearch == RminV); 
    
    if length(num)==1
        Rshift(i,:) = [beta, rlags(num), RminV];    
    else
        mNo = fix(0.5*(length(num)+1));
        MN = num(mNo);
        Rshift(i,:) = [beta,rlags(MN), RminV];  
        %disp('detect many in Rshift')
    end
end
if size(Cshift,1)==1
    alfa1 = Cshift(1);
    Xshift = Cshift(2);
else
    tmp = Cshift(:,3);
    tmp1 = circshift(tmp,1);
    tmp2 = circshift(tmp,-1);
    tmp3 = tmp+tmp1+tmp2;
    tmp3(1) = inf; tmp3(end) = inf;
    num = find(tmp3 == min(tmp3));
    alfa1 = Cshift(num,1);
    Xshift = Cshift(num,2);
    %figure,plot(Rshift(:,1),tmp3)
end
if size(Rshift,1)==1
    alfa2 = Rshift(1);
    Yshift = Rshift(2);
else
    tmp = Rshift(:,3);
    tmp1 = circshift(tmp,1);
    tmp2 = circshift(tmp,-1);
    tmp3 = tmp+tmp1+tmp2;
    tmp3(1) = inf; tmp3(end) = inf;
    num = find(tmp3 == min(tmp3));
    alfa2 = Rshift(num,1);
    Yshift = Rshift(num,2);
    %figure,plot(Rshift(:,1),tmp3)
end
if length(alfa1)>1 || length(alfa2)>1
    Malfa1 = alfa1 * ones(1,length(alfa2));
    Malfa2 = ones(length(alfa1),1) * alfa2';
    DM = abs(Malfa1 -Malfa2);
    minDM = min(DM(:));
    [n1,n2] = find(DM == minDM);
    alfa11 = alfa1(n1(1));
    Xshift1 = Xshift(n1(1));
    alfa22 = alfa2(n2(1));
    Yshift1 = Yshift(n2(1)); 
    alfa1 = alfa11;
    Xshift = Xshift1;
    alfa2 = alfa22;
    Yshift = Yshift1; 
end
translation = [Yshift,Xshift];
rotation = [alfa2,alfa1];
if abs(alfa1-alfa2)<= THbeta
    sigbeta =1;    
else
    sigbeta =0;        
end