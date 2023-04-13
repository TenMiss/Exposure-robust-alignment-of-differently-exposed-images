function mi = mutualInformation(I1,I2,pad)
%
% Calculate mutual information mi between two images I1 and I2 (grayscale)
% which is used to verify the alignment performance
%
% Inputs:
%         I1, I2  - Images containing points that we wish to match,
%                           which can be intensity /color images;
%         pad     - pad the border that is not used for calculation (eg
%                      after translation)
%
% Outputs:
%          mi     - mutual information between two images I1 and I2 
%
% CopyRight by Wu shiqian
% Institute for Infocomm Research
% 21 Nov, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rows cols color]=size(I1);
if color > 1
    I1 = rgb2gray(I1);
    I2 = rgb2gray(I2);
end
I1 = I1(1+pad:end-pad,1+pad:end-pad);
I2 = I2(1+pad:end-pad,1+pad:end-pad);
I1 = double(I1);
I2 = double(I2);
N = 256; 
h = zeros(N,N);
[rows cols]=size(I1);
for r = 1:rows;
    for c = 1:cols;
        h(I1(r,c)+1,I2(r,c)+1) = h(I1(r,c)+1,I2(r,c)+1)+1;
    end
end
b= h./(rows*cols); 
y_marg=sum(b,1); 
x_marg=sum(b,2); 

Hy=0;
for i=1:N;  
    if( y_marg(i)>0 )        
        Hy = Hy -(y_marg(i)*(log2(y_marg(i)))); %marginal entropy for image 1
    end
end
Hx=0;
for i=1:N; 
    if( x_marg(i)>0 )        
        Hx = Hx -(x_marg(i)*(log2(x_marg(i)))); %marginal entropy for image 2
    end
end
Hxy = -sum(sum(b.*(log2(b+(b==0))))); 
nm = min(Hy,Hx);   
mi =(Hx+Hy-Hxy)/nm; 
