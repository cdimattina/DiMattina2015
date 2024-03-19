function Y = gDoublePrimeMat(xAug,thetaMat)
    Y = gPrimeMat(xAug,thetaMat)*(1-2*gMat(xAug,thetaMat));
end

