function write_popfile_c(popn, out_file)
  out = fopen(out_file, "w");
  for d1=1:13
    fwrite(out, popn(d1), 'double');
  end
  fclose(out);
end
