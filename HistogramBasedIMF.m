function [IMFRef2Cur, IMFCur2Ref] = HistogramBasedIMF(ImgRef,ImgCur,row,col,IMFRef2Cur, IMFCur2Ref, T_Ref, T_Cur)
Ref_His = zeros(256,1);
Cur_His = zeros(256,1);
ARef_His = zeros(256,1);
ACur_His = zeros(256,1);
for ii=1:row
    for jj=1:col
        intR = ImgRef(ii, jj);
        intX = ImgCur(ii, jj);
        Ref_His(intR+1) = Ref_His(intR+1)+1;
        Cur_His(intX+1) = Cur_His(intX+1)+1;
    end
end 

%%%%Stop here.
ARef_His(1) = Ref_His(1);
ACur_His(1) = Cur_His(1);
for ii=2:256
    ARef_His(ii) = ARef_His(ii-1)+Ref_His(ii);
    ACur_His(ii) = ACur_His(ii-1)+Cur_His(ii);
end
IMFRef2Cur  = HistogramOptimization(ARef_His, ACur_His, IMFRef2Cur);
IMFCur2Ref  = HistogramOptimization(ACur_His, ARef_His, IMFCur2Ref);

% %%%an image to be detected is darker than its reference image 
if T_Ref>T_Cur
    bb = 256;
    for hh=256:-1:255
        Pix_Sum = 0;
        Vua_Sum = 0;
        for aa=bb:-1:128
            if IMFCur2Ref(aa)==hh
                Pix_Sum = Pix_Sum+Cur_His(aa);
                Vua_Sum = Vua_Sum+Cur_His(aa)*aa;
            else
                bb = aa;
                break;
            end
        end
        if Pix_Sum>0
            IMFRef2Cur(hh) = floor(Vua_Sum/Pix_Sum+0.5);
        end
    end
    bb = 1;
    for hh=1:5
        Pix_Sum = 0;
        Vua_Sum = 0;
        for aa=bb:128
            if IMFRef2Cur(aa)==hh
                Pix_Sum = Pix_Sum+Ref_His(aa);
                Vua_Sum = Vua_Sum+Ref_His(aa)*aa;
            else
                bb = aa;
                break;
            end
        end
        if Pix_Sum>0
            IMFCur2Ref(hh) = floor(Vua_Sum/Pix_Sum+0.5);
        end
    end
%%%an imag to be detected is brighter than its reference image    
else
    bb = 256;
    for hh=256:-1:255
        Pix_Sum = 0;
        Vua_Sum = 0;
        for aa=bb:-1:128
            if IMFRef2Cur(aa)==hh
                Pix_Sum = Pix_Sum+Ref_His(aa);
                Vua_Sum = Vua_Sum+Ref_His(aa)*aa;
            else
                bb = aa;
                break;
            end
        end
        if Pix_Sum>0
            IMFCur2Ref(hh) = floor(Vua_Sum/Pix_Sum+0.5);
        end
    end
    bb = 1;
    for hh=1:5
        Pix_Sum = 0;
        Vua_Sum = 0;
        for aa=bb:128
            if IMFCur2Ref(aa)==hh
                Pix_Sum = Pix_Sum+Cur_His(aa);
                Vua_Sum = Vua_Sum+Cur_His(aa)*aa;
            else
                bb = aa;
                break;
            end
        end
        if Pix_Sum>0
            IMFRef2Cur(hh) = floor(Vua_Sum/Pix_Sum+0.5);
        end
    end    
end


clear ARef_His;
clear ACur_His;
clear Ref_His;
clear Cur_His;
