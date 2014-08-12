begin

  ; User-defined Parameter
  ; luse_g16 = False => use g37 grid
  luse_g16 = False
  longdate = systemfunc("date")

  gx1v6_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx1v6_090205.nc"
  gx3v7_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx3v7_120309.nc"

  if (luse_g16) then
    fgrid = addfile(gx1v6_gridfile,"r")
    fin   = addfile("/glade/p/cgd/oce/people/yeager/OBS/phc/PHC2_TEMP_gx1v6.nc","r")
    outfile = "datasets/PHC2_TEMP_mod_gx1v6.nc"
  else
    fgrid = addfile(gx3v7_gridfile,"r")
    fin   = addfile("/glade/p/cgd/oce/people/yeager/OBS/phc/PHC2_TEMP_gx3v7.nc","r")
    outfile = "datasets/PHC2_TEMP_mod_gx3v7.nc"
  end if

  ; Compute Parameters
  dimSizes = new(3, "integer")
  dimSizes(0:1) = fgrid->grid_dims
  dimSizes(2)   = 60

  ; Variable setup
  VarNames = "TEMPERATURE"+ispan(1,12,1)
  VarNames(0:8) = "TEMPERATURE0"+ispan(1,9,1)

  ; Output
  system("/bin/rm -f "+outfile)
  fout = addfile(outfile, "c")
  setfileoption(fout, "DefineMode", True)

  fAtt = True
  fAtt@title = "Monthly Levitus Climatology (from Steve Yeager)"
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"X", "Y", "depth"/)
  dimUnlim = (/False, False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  do i=0,11
    filevardef(fout, VarNames(i),  typeof(fin->TEMP), (/"depth", "Y", "X"/))
    filevarattdef(fout, VarNames(i), fin->TEMP)
  end do
  setfileoption(fout, "DefineMode", False)

  do i=0,11
    fout->$VarNames(i)$ = fin->TEMP(i,:,:,:)
  end do

end

