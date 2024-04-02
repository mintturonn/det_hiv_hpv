function [tab] = res_poptable(res, xlsxfilename, yrlength)

tab.popall = zeros(16, yrlength);
tab.popneg = zeros(16, yrlength);
tab.poppos = zeros(16, yrlength);
tab.popart = zeros(16, yrlength);
tab.popnoart = zeros(16, yrlength);

for i = 1:length(res)
    
  pop(:,:,i) = res{i}.psize_age_f(end-yrlength+1:end,:)';
  popneg(:,:,i) = res{i}.psize_age_f_neg(end-yrlength+1:end,:)';
  poppos(:,:,i) = res{i}.psize_age_f_art(end-yrlength+1:end,:)' + res{i}.psize_age_f_pos(end-yrlength+1:end,:)';
  popart(:,:,i) = res{i}.psize_age_f_art(end-yrlength+1:end,:)';
  popnart(:,:,i) = res{i}.psize_age_f_pos(end-yrlength+1:end,:)';
  
end



for i=1:13
    for t=1:102

    tab.popall(2+i,t) = median(pop(i,t,:));
    tab.popneg(2+i,t) = median(popneg(i,t,:));
    tab.poppos(2+i,t) = median(poppos(i,t,:));
    tab.popart(2+i,t) = median(popart(i,t,:));
    tab.popnoart(2+i,t) = median(popnart(i,t,:));

    end
end

    % last age group assumed the same as the one before

for t=1:102

    tab.popall(16,t) = median(pop(13,t,:));
    tab.popneg(16,t) = median(popneg(13,t,:));
    tab.poppos(16,t) = median(poppos(13,t,:));
    tab.popart(16,t) = median(popart(13,t,:));
    tab.popnoart(16,t) = median(popnart(13,t,:));

end

tabnew_popall = array2table(tab.popall);
tabnew_popneg = array2table(tab.popneg);
tabnew_poppos = array2table(tab.poppos);
tabnew_popart = array2table(tab.popart);
tabnew_popnoart = array2table(tab.popnoart);

writetable(tabnew_popall,  xlsxfilename, 'Sheet', "All - Pop (P-Y)",   'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_popneg,  xlsxfilename, 'Sheet', "HIV- (P-Y)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_poppos,  xlsxfilename, 'Sheet', "HIV+ (P-Y)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_popart,  xlsxfilename, 'Sheet', "HIV+ ART (P-Y)",     'Range', 'B6:CY21', 'WriteVariableNames',0)
writetable(tabnew_popnoart,xlsxfilename, 'Sheet',"HIV+ no ART (P-Y)",'Range', 'B6:CY21', 'WriteVariableNames',0)

end