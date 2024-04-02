function [tab] = res_table_prctile(stdres, res, percent, xlsxfilename, xlsxsheetname)

tab = NaN(18, 101);

p=percent;

tab(1,:) = prctile(stdres(1:end,:), p);

% 2 first age groups 0
tab(3:4,:) = 0;

for i=1:13
    for t=1:101
tab(4+i,t) = 100000.*prctile(res(i,t,:), p);
    end
end

% last age group assumed the same as the one before
for t=1:101
    tab(18,t) = 100000.*prctile(res(13,t,:), p);
end

tabnew = array2table(tab);


writetable(tabnew, xlsxfilename, 'Sheet', xlsxsheetname, 'Range', 'B4:CX21', 'WriteVariableNames',0)

end

