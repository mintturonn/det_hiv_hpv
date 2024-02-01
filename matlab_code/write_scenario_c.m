function write_scenario_c(vaccscen, screenscen, output_fn)
  out = fopen(output_fn, 'w');
  
  fwrite(out, vaccscen.setup, "int");
  if (vaccscen.setup >= 1)
    fwrite(out, vaccscen.agef1, "double");
    fwrite(out, vaccscen.agef2hp, "double");
    fwrite(out, vaccscen.agef3hp, "double");
    fwrite(out, vaccscen.agef4hp, "double");
  end
  
  fwrite(out, screenscen.setup, "int");
  if (screenscen.setup > 0)
    fwrite(out, screenscen.yr1, "int");
    fwrite(out, screenscen.yr2, "int");
    fwrite(out, screenscen.yr3, "int");
    fwrite(out, screenscen.cintx(1), "double");
    fwrite(out, screenscen.cintx(2), "double");
    fwrite(out, screenscen.cctx(1), "double");
    fwrite(out, screenscen.cctx(2), "double");
  end
  
  if (screenscen.setup == 1)
    for si = 1:3
      fwrite(out, screenscen.scen_15to24(si), "double");
      fwrite(out, screenscen.scen_25to34(si), "double");
      fwrite(out, screenscen.scen_35to49(si), "double");
    end
    
  elseif (screenscen.setup == 2)
    for si = 1:3
      fwrite(out, screenscen.scen_15to24hn(si), "double");
      fwrite(out, screenscen.scen_25to34hn(si), "double");
      fwrite(out, screenscen.scen_35to49hn(si), "double");
      fwrite(out, screenscen.scen_15to24hp(si), "double");
      fwrite(out, screenscen.scen_25to34hp(si), "double");
      fwrite(out, screenscen.scen_35to49hp(si), "double");
    end   
  end
      
  fclose(out);
    
end