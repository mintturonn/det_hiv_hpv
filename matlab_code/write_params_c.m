function write_params_c(bp0, dp0, hp0, pp0, out_file)

  function process(out, array)
    fwrite(out, length(array), "int");
    for i=1:length(array)
      fwrite(out, array(i), "double");
    end
  end

  out = fopen(out_file, 'w');
  process(out, bp0);
  process(out, dp0);
  process(out, hp0);
  process(out, pp0);
  fclose(out);
end