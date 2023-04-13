function [delta, phi] = imfLBPalignment(im)
% Copyright: Wu Shiqian
% 11 Sep 2010
sampling = 2;
I = im{1};
[ri,ci,h] = size(I);
if h>1
    im{1} = rgb2gray(im{1});
    im{2} = rgb2gray(im{2});
end
A = [1,2,4,8,16,32,64,128,256,512];
B = ceil(min([ri,ci])/32);
D = B - A;
num = find(D>=0);
scale = 1./A(num);
n = length(num); %%%control the number of scales
im0 = cell(1,n);
im1 = cell(1,n);
A = uint8(im{1});
B = uint8(im{2});
%W = ones(ri,ci);
im0{1} = A;
im1{1} = B;
for i=2:n
    im0{i} = imresize(A,scale(i),'bilinear');
    im1{i} = imresize(B,scale(i),'bilinear');
end
%%%%%%%%%%%%%%%%%%%%%
gama =[];
stot = zeros(1,3);
if(1)
    for k = n:-1:n-2
        f0 = im0{k};
        f1 = im1{k};   
        [rr,cc]= size(f1);
        [ROTS,dxy,sigbeta] = imfLBPCoarseTune(f0,f1,[-20,20],fix(rr/4),fix(cc/4),3,1,1);    
    %%%%[-5,5] is the current setting; [-20, 20] is the original setting; (fix(rr/4),fix(cc/4)) 
    %%%is the original searching range and (fix(rr/16),fix(cc/16)) is the current searching range
        if sigbeta==1   
            disp('Find the initial angle')
            break        
        else
            gama = [gama;ROTS];
        end
    end
    if sigbeta == 1
        nn = k-1;
        stot = [dxy(1),dxy(2),-mean(ROTS)];
        A1 = shift(im1{k-1},2*stot(2),2*stot(1)); 
        A1 = imrotate(A1,-stot(3),'bilinear','crop');
        srow = ceil(2*stot(1));
        scol = ceil(2*stot(2));
        A0 = im0{k-1};
        if srow > 0
            A0(1:srow,:)=[]; 
            A1(1:srow,:)=[];
        else
            A0(end-srow+1:end,:) =[];
            A1(end-srow+1:end,:) =[]; 
        end
        if scol > 0
            A0(:,1:scol)=[];
            A1(:,1:scol)=[];
        else
            A0(:,end-scol+1:end) =[];
            A1(:,end-scol+1:end) =[]; 
        end
        im0{k-1} = A0;
        im1{k-1} = A1;
    else
        nn = n-1;
        ROTS = median(gama(:));
        stot = [0,0,-mean(ROTS)];
        im1{n-1} = imrotate(im1{n-1},-stot(3),'bilinear','crop');
    end
else
    nn = n;
end
for pyrlevel = nn:-1:1 %%%control the number of iteration, original setting is 1
    f0 = im0{pyrlevel};
    f1 = im1{pyrlevel}; 
    f0 = LBP8Feature(f0,0);
    f1 = LBP8Feature(f1,0);    
    [y0,x0]=size(f0);
    xmean=x0/2; ymean = y0/2;
    x=kron((-xmean:xmean-1),ones(y0,1));
    y=kron(ones(1,x0),(-ymean:ymean-1)');
    sigma=1;
    g3 = exp(-(y.^2+x.^2)/(2*sigma^2))/2/pi/sigma^2;
    g1 = -g3.*y; 
    g2 = -g3.*x; 
    a=real(ifft2(fft2(f1).*fft2(g2))); 
    c=real(ifft2(fft2(f1).*fft2(g1))); 
    b=real(ifft2(fft2(f1).*fft2(g3)))-real(ifft2(fft2(f0).*fft2(g3))); 
    R=c.*x-a.*y; 
    [row, col, h] = size(a);
    a11 = 0;
    a12 = 0;
    a13 = 0;
    a22 = 0;
    a23 = 0;
    a33 = 0;
    b1 = 0;
    b2 = 0;
    b3 = 0;
    for i = 1: sampling: row
        for j = 1: sampling : col
            a11 = a11+a(i,j)*a(i,j);
            a12 = a12+a(i,j)*c(i,j);
            a13 = a13+a(i,j)*R(i,j);
            a22 = a22+c(i,j)*c(i,j);
            a23 = a23+c(i,j)*R(i,j);            
            a33 = a33+R(i,j)*R(i,j);
            b1 = b1+a(i,j)*b(i,j);
            b2 = b2+b(i,j)*c(i,j);
            b3 = b3+b(i,j)*R(i,j);
        end
    end
    a21 = a12;
    a31 = a13;
    a32 = a23;
            
%     a11 = sum(sum(a.*a)); a12 = sum(sum(a.*c)); a13 = sum(sum(R.*a));
%     a21 = sum(sum(a.*c)); a22 = sum(sum(c.*c)); a23 = sum(sum(R.*c)); 
%     a31 = sum(sum(R.*a)); a32 = sum(sum(R.*c)); a33 = sum(sum(R.*R));
%     b1 = sum(sum(a.*b)); b2 = sum(sum(c.*b)); b3 = sum(sum(R.*b));
    Ainv = [a11 a12 a13; a21 a22 a23; a31 a32 a33]^(-1);
    s = Ainv*[b1; b2; b3];
    st = s;
    it=1;       
    while ((abs(s(1))+abs(s(2))+abs(s(3))*9/pi>0.01)&&it<=20) %%%original setting is ((abs(s(1))+abs(s(2))+abs(s(3))*180/pi/20>0.01)&&it<=20)
       if abs(st(3)*180/pi)>=0.5
            tmp = shift(im0{pyrlevel},-st(1),-st(2));
            tmp = imrotate(tmp,-st(3)*180/pi,'bilinear','crop');        
            ff0 = LBP8Feature(tmp,0);
        else
            ff0 = shift(f0,-st(1),-st(2));
        end
        b = real(ifft2(fft2(f1).*fft2(g3)))-real(ifft2(fft2(ff0).*fft2(g3)));
        [row, col, h] = size(a);        
        b1 = 0;
        b2 = 0;
        b3 = 0;
        for i = 1: sampling: row
            for j = 1: sampling : col
                b1 = b1+a(i,j)*b(i,j);
                b2 = b2+b(i,j)*c(i,j);
                b3 = b3+b(i,j)*R(i,j);
            end
        end        
%        s = Ainv*[sum(sum(a.*b)); sum(sum(c.*b)); sum(sum(R.*b))];
        s = Ainv*[b1; b2; b3];
        st = st+s;
        it = it+1;        
    end
    st(3)=-st(3)*180/pi;
    st = st';
    st(1:2) = st(2:-1:1);
    stot = [2*stot(1:2)+st(1:2) stot(3)+st(3)]; 
    if pyrlevel>1
        A1 = shift(im1{pyrlevel-1},2*stot(2),2*stot(1)); 
        A1 = imrotate(A1,-stot(3),'bilinear','crop');
        srow = ceil(2*stot(1));
        scol = ceil(2*stot(2));
        A0 = im0{pyrlevel-1};
        if srow > 0
            A0(1:srow,:)=[]; 
            A1(1:srow,:)=[]; 
        else
            A0(end-srow+1:end,:) =[];
            A1(end-srow+1:end,:) =[]; 
        end
        if scol > 0
            A0(:,1:scol)=[];
            A1(:,1:scol)=[];
        else
            A0(:,end-scol+1:end) =[];
            A1(:,end-scol+1:end) =[];
        end
        im0{pyrlevel-1} = A0;
        im1{pyrlevel-1} = A1;
    end
end
phi = stot(3);
delta = stot(1:2);
