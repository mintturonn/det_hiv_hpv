function [tab] = res_poptable_prctile(res, percent, xlsxfilename, yrlength)

p=percent;

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
    for t=1:101

    tab.popall(2+i,t) = prctile(pop(i,t,:), p);
    tab.popneg(2+i,t) = prctile(popneg(i,t,:), p);
    tab.poppos(2+i,t) = prctile(poppos(i,t,:), p);
    tab.popart(2+i,t) = prctile(popart(i,t,:), p);
    tab.popnoart(2+i,t) = prctile(popnart(i,t,:), p);

    end
end

    % last age group assumed the same as the one before

for t=1:101

    tab.popall(16,t) = prctile(pop(13,t,:), p);
    tab.popneg(16,t) = prctile(popneg(13,t,:), p);
    tab.poppos(16,t) = prctile(poppos(13,t,:), p);
    tab.popart(16,t) = prctile(popart(13,t,:), p);
    tab.popnoart(16,t) = prctile(popnart(13,t,:), p);

end

tabnew_popall = array2table(tab.popall);
tabnew_popneg = array2table(tab.popneg);
tabnew_poppos = array2table(tab.poppos);
tabnew_popart = array2table(tab.popart);
tabnew_popnoart = array2table(tab.popnoart);

writetable(tabnew_popall,  xlsxfilename, 'Sheet', "All - Pop (P-Y)",   'Range', 'B6:CX21', 'WriteVariableNames',0)
writetable(tabnew_popneg,  xlsxfilename, 'Sheet', "HIV- (P-Y)",     'Range', 'B6:CX21', 'WriteVariableNames',0)
writetable(tabnew_poppos,  xlsxfilename, 'Sheet', "HIV+ (P-Y)",     'Range', 'B6:CX21', 'WriteVariableNames',0)
writetable(tabnew_popart,  xlsxfilename, 'Sheet', "HIV+ ART (P-Y)",     'Range', 'B6:CX21', 'WriteVariableNames',0)
writetable(tabnew_popnoart,xlsxfilename, 'Sheet',"HIV+ no ART (P-Y)",'Range', 'B6:CX21', 'WriteVariableNames',0)

end