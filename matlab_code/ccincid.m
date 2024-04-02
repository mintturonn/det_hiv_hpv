function [cc_incid_all, cc_incid_neg, cc_incid_posall, cc_incid_art, cc_incid_pos] = ccincid(res)

for i = 1:length(res)
    
    cc_incid_art(:,:,i) = (res{i}.ccinc_art ./  res{i}.psize_age_f_art)' ;
    cc_incid_pos(:,:,i) = (res{i}.ccinc_pos ./ res{i}.psize_age_f_pos)' ;
    cc_incid_posall(:,:,i) =  ( (res{i}.ccinc_art + res{i}.ccinc_pos )./  (res{i}.psize_age_f_art + res{i}.psize_age_f_pos) )' ;
    cc_incid_neg(:,:,i) = (res{i}.ccinc_neg ./ res{i}.psize_age_f_neg)' ;
                                             
    cc_incid_all(:,:,i) = ( (res{i}.ccinc_neg + res{i}.ccinc_art +  res{i}.ccinc_pos ) ./ ...
                            (res{i}.psize_age_f_neg + res{i}.psize_age_f_art + res{i}.psize_age_f_pos) )';                    
                        
end

end