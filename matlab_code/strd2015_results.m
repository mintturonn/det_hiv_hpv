function [summod] = strd2015_results(mod, pop2015, lgtvar)


%mod(isnan(mod))=0;

%t1 =   find(mod(end,:,1),1,'first') ;

clngth = 171-lgtvar+1;
mod2 = zeros(20, clngth, length(mod(1,1,:)) );

% fist >0 (to identify when incidence starts, for HIV+ no incidence before
% the epidemic)



for i=1:length(mod(1,1,:))
    
    mod2(3:15,:,i) = mod(:,lgtvar:end,i);
    mod2(16:20,:,i) = repmat(mod(end,lgtvar,i), 5, clngth);
    
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
    
    summod(i,:) = 100000.*sum(stdmod(:,:,i),1) ./ sum(pop2015);
    
end


end