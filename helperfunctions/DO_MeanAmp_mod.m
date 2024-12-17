function StatMat = vs_MeanAmp(chanlocs,preStim,srate,DatMat,Lat,Mlen,Chan,JK);

mBeg = Lat - Mlen;
mEnd = Lat + Mlen;
dBeg = round((mBeg - preStim)/1000 * srate);
dEnd = round((mEnd - preStim)/1000 * srate);

StatMat = [];

for i = 1 : length(Chan)
   ch(i) = DO_channum(chanlocs,Chan(i))
end;   

for sub = 1 : size(DatMat,1)
   for bed = 1 : size(DatMat,2)
      if length(Chan) > 1 
      for k = 1 : length(Chan) 
          StatMat(sub,bed,k) = mean(DatMat(sub,bed,ch(k),dBeg:dEnd),4);
      end;    
      elseif length(Chan) == 1   
          StatMat(sub,bed) = mean(DatMat(sub,bed,ch(1),dBeg:dEnd),4);
      end;
   end;
end;   

if JK == 1
    for bed = 1:size(StatMat,2)
    StatMat_all(bed) = nanmean(StatMat(:,bed),1);
        for sub = 1 : size(StatMat,1)
            StatMat2(sub,bed) = (size(StatMat,1) * StatMat_all(bed)) - ((size(StatMat,1) - 1) * StatMat(sub,bed));
        end % forsub
    end % for bed

    StatMat = StatMat2;
end % if JK
