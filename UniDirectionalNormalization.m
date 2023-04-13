function [ImgCur]=UniDirectionalNormalization(ImgCur,row,col,IMFCur2Ref)
for ii=1:row
    for jj=1:col
        ImgCur(ii,jj) = IMFCur2Ref(ImgCur(ii,jj)+1)-1;
    end
end


