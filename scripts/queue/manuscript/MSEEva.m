function [ y ] = MSEEva( testx,Econshort,KernPara,coinv,FSF,FS,fbasis,kernind )
%MSEEva returns the expected MSE of a single SK model

msep1 = KernPara(2)^2; conL = size(Econshort,2);
if size(testx,1) == conL
    fi1 = Kernval(testx',Econshort,KernPara,kernind);
    etav = fbasis(testx) - FS*fi1;
    mse = msep1-sum(fi1'*coinv .*fi1', 2)+sum(etav'*FSF.* etav', 2);
    y = mse';  
elseif size(testx,2) == conL
    fi1 = Kernval(testx,Econshort,KernPara,kernind);
    etav = fbasis(testx) - FS*fi1;
    y = msep1-sum(fi1'*coinv .*fi1', 2)+sum(etav'*FSF.* etav', 2);
end


end

