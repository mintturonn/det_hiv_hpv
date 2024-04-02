function [tab] = res_ccctable(res, xlsxfilename)

tab.cccall = zeros(16, 102);
tab.cccneg = zeros(16, 102);
tab.cccpos = zeros(16, 102);
tab.cccart = zeros(16, 102);
tab.cccnoart = zeros(16, 102);

tab.popall = zeros(16, 102);
tab.popneg = zeros(16, 102);
tab.poppos = zeros(16, 102);
tab.popart = zeros(16, 102);
tab.popnoart = zeros(16, 102);

for i = 1:length(res)
    
  ccc(:,:,i) = cumsum(res{i}.ccinc_neg(end-102+1:end,:)' + res{i}.ccinc_art(end-102+1:end,:)' +  res{i}.ccinc_pos(end-102+1:end,:)', 2) ;
  cccneg(:,:,i) = cumsum(res{i}.ccinc_neg(end-102+1:end,:)', 2) ;
  cccpos(:,:,i) = cumsum(res{i}.ccinc_art(end-102+1:end,:)' + res{i}.ccinc_pos(end-102+1:end,:)', 2) ;
  cccart(:,:,i) = cumsum(res{i}.ccinc_art(end-102+1:end,:)', 2) ;
  cccnart(:,:,i) = cumsum(res{i}.ccinc_pos(end-102+1:end,:)', 2) ;
  
end


for i = 1:length(res)
    
  pop(:,:,i) = cumsum(res{i}.psize_age_f(end-102+1:end,:)', 2);
  popneg(:,:,i) = cumsum(res{i}.psize_age_f_neg(end-102+1:end,:)', 2);
  poppos(:,:,i) = cumsum(res{i}.psize_age_f_art(end-102+1:end,:)' + res{i}.psize_age_f_pos(end-102+1:end,:)', 2);
  popart(:,:,i) = cumsum(res{i}.psize_age_f_art(end-102+1:end,:)', 2);
  popnart(:,:,i) = cumsum(res{i}.psize_age_f_pos(end-102+1:end,:)', 2);
  
end

for i=1:13
    for t=1:102

    tab.cccall(2+i,t) = median(ccc(i,t,:));
    tab.cccneg(2+i,t) = median(cccneg(i,t,:));
    tab.cccpos(2+i,t) = median(cccpos(i,t,:));
    tab.cccart(2+i,t) = median(cccart(i,t,:));
    tab.cccnoart(2+i,t) = median(cccnart(i,t,:));
    
    tab.popall(2+i,t) = median(pop(i,t,:));
    tab.popneg(2+i,t) = median(popneg(i,t,:));
    tab.poppos(2+i,t) = median(poppos(i,t,:));
    tab.popart(2+i,t) = median(popart(i,t,:));
    tab.popnoart(2+i,t) = median(popnart(i,t,:));

    end
end



tabnew_cccall = array2table(tab.cccall);
tabnew_cccneg = array2table(tab.cccneg);
tabnew_cccpos = array2table(tab.cccpos);
tabnew_cccart = array2table(tab.cccart);
tabnew_cccnoart = array2table(tab.cccnoart);

tabnew_popall = array2table(tab.popall);
tabnew_popneg = array2table(tab.popneg);
tabnew_poppos = array2table(tab.poppos);
tabnew_popart = array2table(tab.popart);
tabnew_popnoart = array2table(tab.popnoart);

writetable(tabnew_popall,  xlsxfilename, 'Sheet', "All - Pop (N)",   'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_popneg,  xlsxfilename, 'Sheet', "HIV- (N)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_poppos,  xlsxfilename, 'Sheet', "HIV+ (N)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_popart,  xlsxfilename, 'Sheet', "HIV+ ART (N)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_popnoart,xlsxfilename, 'Sheet',"HIV+ no ART (N)",'Range', 'B6:CY21', 'WriteVariableNames',0)

writetable(tabnew_cccall,  xlsxfilename, 'Sheet', "Pop(All) CCC",   'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_cccneg,  xlsxfilename, 'Sheet', "HIV- (CCC)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_cccpos,  xlsxfilename, 'Sheet', "HIV+ (CCC)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_cccart,  xlsxfilename, 'Sheet', "HIV+ ART (CCC)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_cccnoart,xlsxfilename, 'Sheet',"HIV+ no ART (CCC)",'Range', 'B6:CY21', 'WriteVariableNames',0)

end