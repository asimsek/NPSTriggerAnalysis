def getCrossSection(histo, delLumi, scale=1, removeOutliers=0.98):
    margin = 0.05
    if removeOutliers>0.5:
        raise Exception("removeOutliers cannot be larger than 0.5. %f"%removeOutliers)
    print("getCrossSection START")
    print(histo, delLumi, scale, removeOutliers)
    npointsMedian = 10000000
    average = 0
    count = 0
    nhisto = histo.Clone(histo.GetName()+delLumi.GetName())
    if histo.Integral() == 0: raise Exception("getCrossSection: %s is empty"%histo.GetName())
    if delLumi.Integral() == 0: raise Exception("getCrossSection: %s is empty"%delLumi.GetName())
    maxAllowedValue = histo.GetMaximum()
    minAllowedValue = histo.GetMinimum()
    if removeOutliers>0: ##compute the quantile of histo/delLumi using only npointsMedian points
        y1 = histo.GetArray()
        y2 = delLumi.GetArray()
        ys = [1E-10] ## to keep the same number of bins with rate>0
        npointsMedian = min(npointsMedian,histo.GetNbinsX()-1)
        jump = max(1., float(histo.GetNbinsX())/ npointsMedian)
        for x in range(1, npointsMedian):
            i = int(x*jump)
            if y1[i]>0 and y2[i]>0:
                ys.append(y1[i]/y2[i])
        maxAllowedIdx = min(int(len(ys) * (1.0-removeOutliers)), len(ys)-1)
        minAllowedIdx = max(int(len(ys) * (removeOutliers)), 0)
#        minAllowedIdx = 0 ## do not apply the lower cut
        if len(ys)>1:
            marginVal = margin * abs(sorted(ys)[maxAllowedIdx])
            maxAllowedValue = sorted(ys)[maxAllowedIdx] + marginVal
            minAllowedValue = sorted(ys)[minAllowedIdx] - marginVal
        else: ##only 1 point
            val = histo.GetMaximum()/delLumi.GetMaximum()
            marginVal = abs(val) * margin
            maxAllowedValue = val + marginVal
            minAllowedValue = val - marginVal
        print(ys)
        print(sorted(ys))
        print("getCrossSection",histo.GetName(),delLumi.GetName(),minAllowedValue, maxAllowedValue, removeOutliers, len(ys), maxAllowedIdx, minAllowedIdx,  histo.Integral(), delLumi.Integral(), jump)
    for i in range(len(histo)-1):
         val = float(histo[i]) 
         lum = float(delLumi[i])
         if lum>0 and val>=0:
             nhisto.SetBinContent(i, val/lum*scale)
             nhisto.SetBinError(i, val**0.5/lum*scale)
         else:
             nhisto.SetBinContent(i, 0)
             nhisto.SetBinError(i, 0)
             if lum<0: print("getCrossSection: lum<0 in %s bin %d"%(delLumi.GetName(), i))
             if val<0: print("getCrossSection: val<0 in %s bin %d"%(histo.GetName(), i))
    nhisto.SetMaximum(maxAllowedValue*scale)
    nhisto.SetMinimum(minAllowedValue*scale)
    print("minAllowedValue=",minAllowedValue, "maxAllowedValue=",maxAllowedValue*scale, "removeOutliers=",removeOutliers)
    print("getCrossSection STOP", nhisto.GetMaximum(), nhisto.GetMinimum())
    return nhisto 
