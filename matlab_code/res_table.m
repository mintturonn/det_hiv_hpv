function [tab] = res_table(stdres, res, xlsxfilename, xlsxsheetname)

tab = NaN(18, 102);

tab(1,:) = median(stdres(1:end,:));

% 2 first age groups 0
tab(3:4,:) = 0;

for i=1:13
    for t=1:102
tab(4+i,t) = 100000.*median(res(i,t,:));
    end
end

% last age group assumed the same as the one before
for t=1:102
    tab(18,t) = 100000.*median(res(13,t,:));
end

tabnew = array2table(tab);


writetable(tabnew, xlsxfilename, 'Sheet', xlsxsheetname, 'Range', 'B4:CY21', 'WriteVariableNames',0)

end

