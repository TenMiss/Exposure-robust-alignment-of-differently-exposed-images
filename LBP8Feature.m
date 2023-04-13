function Ic = LBP8Feature(X,TH)
% Extract local feture of image X
% Input: X is an image 
% Output:
%        I is an image with the same size of X
%        Each element in I contains the bit information
% Copyright: Wu Shiqian
% Institute for Infocomm Research
% 3 August 2010
[row,col,h]=size(X);
if h > 1
    X = rgb2gray(X);
end
A =  repmat(0,[row+2 col+2]); 
A1 = A; A2 = A; A3 = A; A4 = A;
A5 = A; A6 = A; A7 = A; A8 = A;
A(2:row+1,2:col+1) = X; 
A1(3:row+2,3:col+2) = X;  A2(3:row+2,2:col+1) = X; 
A3(3:row+2,1:col) = X;    A4(2:row+1,1:col) = X; 
A5(1:row,1:col) = X;      A6(1:row,2:col+1) = X; 
A7(1:row,3:col+2) = X;    A8(2:row+1,3:col+2) = X; 
Ic = ((A1-A)>=TH) + 2*((A2-A)>=TH) + 4*((A3-A)>=TH) + 8*((A4-A)>=TH) ...
    + 16*((A5-A)>=TH) + 32*((A6-A)>=TH) + 64*((A7-A)>=TH) + 128*((A8-A)>=TH);
Ic = Ic(2:row+1,2:col+1);
