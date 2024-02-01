function [output] = read_c_output(file)
  function ret = readMatrix(f)
    d1 = fread(f, 1, 'int');
    d2 = fread(f, 1, 'int');
    ret = zeros(d1, d2);
    for i = 1:d1
      for j = 1:d2
        ret(i,j) = fread(f, 1, 'double');
      end
    end
  end

  function ret = read8D(f)
    ret = zeros(5,5,5,4,2,3,13,2);
    for a = 1:5
      for b = 1:5
        for c = 1:5
          for d = 1:4
            for e = 1:2
              for f1 = 1:3
                for g = 1:13
                  for h = 1:2
                    ret(a,b,c,d,e,f1,g,h) = fread(f, 1, 'double');
                  end
                end
              end
            end
          end
        end
      end
    end
  end

  function ret = readMatrix3D(f)
    d1 = fread(f, 1, 'int');
    d2 = fread(f, 1, 'int');
    d3 = fread(f, 1, 'int');
    ret = zeros(d1,d2,d3);
    for i = 1:d1
      for j = 1:d2
        for k = 1:d3
          ret(i,j,k) = fread(f, 1, 'double');
        end
      end
    end
  end

  function ret = readMatrix4D(f)
    d1 = fread(f, 1, 'int');
    d2 = fread(f, 1, 'int');
    d3 = fread(f, 1, 'int');
    d4 = fread(f, 1, 'int');
    ret = zeros(d1,d2,d3,d4);
    for i = 1:d1
      for j = 1:d2
        for k = 1:d3
          for m = 1:d4
            ret(i,j,k,m) = fread(f, 1, 'double');
          end
        end
      end
    end
  end

  function ret = readMatrix5D(f)
    d1 = fread(f, 1, 'int');
    d2 = fread(f, 1, 'int');
    d3 = fread(f, 1, 'int');
    d4 = fread(f, 1, 'int');
    d5 = fread(f, 1, 'int');
    ret = zeros(d1,d2,d3,d4,d5);
    for i = 1:d1
      for j = 1:d2
        for k = 1:d3
          for m = 1:d4
            for p = 1:d5
              ret(i,j,k,m,p) = fread(f, 1, 'double');
            end
          end
        end
      end
    end
  end

  function ret = readArray(f)
    count = fread(f, 1, 'int');
    ret = fread(f, count, 'double');
  end

  function ret = readIntArray(f)
    count = fread(f, 1, 'int');
    ret = fread(f, count, 'int');
  end

  f = fopen(file);
  output.popbyrisk_f = readMatrix(f);
  output.popbyrisk_m = readMatrix(f);
  
  output.hivart_9to14 = readMatrix(f);
  output.hivart_15to49 = readMatrix(f);
  output.hivart_15to74 = readMatrix(f);
  
  output.hiv_9to14 = readMatrix(f);
  output.hiv_15to49 = readMatrix(f);
  output.hiv_15to74 = readMatrix(f);
  output.art_15plus = readMatrix(f);
  
  output.hiv_fage = readMatrix(f);
  output.hiv_mage = readMatrix(f);
  output.art_age = readMatrix(f);
  
  output.pop1980 = read8D(f);
  output.pop2018 = read8D(f);
  output.pop2120 = read8D(f);
  
  output.psize = readMatrix(f);
  output.psize_age_f = readMatrix(f);
  output.psize_age_f_neg = readMatrix(f);
  output.psize_age_f_pos = readMatrix(f);
  output.psize_age_f_art = readMatrix(f);
  output.psize_age_m = readMatrix(f);
  output.psize_age_m_neg = readMatrix(f);
  output.psize_age_m_pos = readMatrix(f);
  output.psize_age_m_art = readMatrix(f);
  output.vacc_age_f = readMatrix(f);
  output.vacc_age_m = readMatrix(f);
  output.vacc_age_f_pos = readMatrix(f);
  output.vacc_age_m_pos = readMatrix(f);
  output.deatht = readMatrix(f);
  output.hivdftot = readMatrix(f);
  output.hivdmtot = readMatrix(f);
  output.ccdftot = readMatrix(f);
  output.ccdftot_neg = readMatrix(f);
  output.ccdftot_pos = readMatrix(f);
  output.ccdftot_art = readMatrix(f);
  output.hpv1cum_pop = readMatrix(f);
  output.hpv2cum_pop = readMatrix(f);
  output.hpv3cum_pop = readMatrix(f);
  output.hpv1hncum_pop = readMatrix(f);
  output.hpv2hncum_pop = readMatrix(f);
  output.hpv3hncum_pop = readMatrix(f);
  
  output.hivsus_pop = readMatrix5D(f);
  
  output.hpv1hncum = readMatrix(f);
  output.hpv2hncum = readMatrix(f);
  output.hpv3hncum = readMatrix(f);
  output.hpv1cum = readMatrix(f);
  output.hpv2cum = readMatrix(f);
  output.hpv3cum = readMatrix(f);
  
  output.ccinc = readMatrix(f);
  output.ccinc_neg = readMatrix(f);
  output.ccinc_pos = readMatrix(f);
  output.ccinc_art = readMatrix(f);
  output.ccinc_nvt = readMatrix(f);
  output.ccinc_nvt_neg = readMatrix(f);
  output.ccinc_nvt_pos = readMatrix(f);
  output.ccinc_nvt_art = readMatrix(f);
  
  output.cin = readMatrix(f);
  output.cin_neg = readMatrix(f);
  output.cin_pos = readMatrix(f);
  output.cin_art = readMatrix(f);
  
  output.hivinctot = readMatrix4D(f);
  
  output.hpv1618_f = readMatrix3D(f);
  output.hpv1618_m = readMatrix3D(f);
  output.nvthpv_f = readMatrix3D(f);
  output.nvthpv_m = readMatrix3D(f);
  output.hpv9vt_f = readMatrix3D(f);
  output.hpv9vt_m = readMatrix3D(f);
  output.hpv_f = readMatrix3D(f);
  output.hpv_m = readMatrix3D(f);
  output.cinpr_f = readMatrix3D(f);
  output.ccpr_f = readMatrix3D(f);
  
  output.screencovf = readMatrix3D(f);
  output.hiv_fsw = readArray(f);
  output.hiv_mcli = readArray(f);
  
  output.errmessage = reshape(readIntArray(f), 1, 3);
  output.testtime = readArray(f);
  output.testtime1yr = readArray(f);
  fclose(f);
end
