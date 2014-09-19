begin

  ; User-defined Parameter
  ; luse_g16 = False => use g37 grid
  luse_g16 = False
  lcheck_mask = True
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
    if (lcheck_mask) then
      fmask = addfile("KMT/gx3v7/KMT_140919.nc","r")
    end if
  end if

  ; Compute Parameters
  dimSizes = new(3, "integer")
  dimSizes(0:1) = fgrid->grid_dims
  dimSizes(2)   = 60

  ; Variable setup
  VarNames = "TEMPERATURE"+ispan(1,12,1)
  VarNames(0:8) = "TEMPERATURE0"+ispan(1,9,1)
  FILE_TEMP = fin->TEMP
  delete(fin)

  if (lcheck_mask) then
    do k=0,59
      if (any(ismissing(FILE_TEMP(0,k,:,:)).and.(k.lt.fmask->KMT))) then
        print((/"Problem in level "+(k+1)/))
        do month=0,11
          level_temp = where(ismissing(FILE_TEMP(month,k,:,:)).and.(k.lt.fmask->KMT),0.,FILE_TEMP(month,k,:,:))
          FILE_TEMP(month,k,:,:) = (/level_temp/)
        end do
      end if
    end do
    delete(fmask)
  end if

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
    filevardef(fout, VarNames(i),  typeof(FILE_TEMP), (/"depth", "Y", "X"/))
    filevarattdef(fout, VarNames(i), FILE_TEMP)
  end do
  setfileoption(fout, "DefineMode", False)

  do i=0,11
    fout->$VarNames(i)$ = FILE_TEMP(i,:,:,:)
  end do

end

