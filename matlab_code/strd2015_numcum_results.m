function [medcumsummod] = strd2015_numcum_results(mod, pop2015, xlsxfilename, xlsxsheetname)


%mod(isnan(mod))=0;

%t1 =   find(mod(end,:,1),1,'first') ;


mod2 = zeros(20, 102, length(mod(1,1,:)) );

% fist >0 (to identify when incidence starts, for HIV+ no incidence before
% the epidemic)



for i=1:length(mod(1,1,:))
    
    mod2(3:15,:,i) = mod(:,:,i);
    mod2(16:20,:,i) = repmat(mod(end,102,i), 5, 102);
    
end

%%
modf = mod2;

for i=1:length(mod(1,1,:))
   for c = 2:length(mod2(16,:,1))-1
   
     modf(16,c+1,i) = modf(16-1,c,i) ;  
     
   end
end

for i=1:length(mod(1,1,:))
   for c = 3:length(mod2(16,:,1))-1
   
     modf(17,c+1,i) = modf(17-1,c,i) ;  
     
   end
end

for i=1:length(mod(1,1,:))
   for c = 4:length(mod2(16,:,1))-1
   
     modf(18,c+1,i) = modf(18-1,c,i) ;  
     
   end
end

for i=1:length(mod(1,1,:))
   for c = 5:length(mod2(16,:,1))-1
   
     modf(19,c+1,i) = modf(19-1,c,i) ;  
     
   end
end

for i=1:length(mod(1,1,:))
   for c = 6:length(mod2(16,:,1))-1
   
     modf(20,c+1,i) = modf(20-1,c,i) ;  
     
   end
end

%%

for i=1:length(mod(1,1,:))

    stdmod(:,:,i) = pop2015 .* modf(:,:,i);
    
    summod(i,:) = sum(stdmod(:,:,i),1) ;
    cumsummod(i,:) = cumsum(summod(i,:));
    
end

medcumsummod = array2table(median(cumsummod));

writetable(medcumsummod, xlsxfilename, 'Sheet', xlsxsheetname, 'Range', 'B4:CY4', 'WriteVariableNames',0)


end