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
    fin   = addfile("/glade/p/cgd/oce/people/yeager/OBS/phc/PHC2_SALT_gx1v6.nc","r")
    outfile = "datasets/PHC2_SALT_mod_gx1v6.nc"
  else
    fgrid = addfile(gx3v7_gridfile,"r")
    fin   = addfile("/glade/p/cgd/oce/people/yeager/OBS/phc/PHC2_SALT_gx3v7.nc","r")
    outfile = "datasets/PHC2_SALT_mod_gx3v7.nc"
    if (lcheck_mask) then
      fmask = addfile("KMT/gx3v7/KMT_140919.nc","r")
    end if
  end if

  ; Compute Parameters
  dimSizes = new(3, "integer")
  dimSizes(0:1) = fgrid->grid_dims
  dimSizes(2)   = 60

  ; Variable setup
  VarNames = "SALINITY"+ispan(1,12,1)
  VarNames(0:8) = "SALINITY0"+ispan(1,9,1)
  FILE_SALT = fin->SALT
  delete(fin)

  if (lcheck_mask) then
    do k=0,59
      if (any(ismissing(FILE_SALT(0,k,:,:)).and.(k.lt.fmask->KMT))) then
        print((/"Problem in level "+(k+1)/))
        do month=0,11
          level_salt = where(ismissing(FILE_SALT(month,k,:,:)).and.(k.lt.fmask->KMT),34.7,FILE_SALT(month,k,:,:))
          FILE_SALT(month,k,:,:) = (/level_salt/)
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
    filevardef(fout, VarNames(i),  typeof(FILE_SALT), (/"depth", "Y", "X"/))
    filevarattdef(fout, VarNames(i), FILE_SALT)
  end do
  setfileoption(fout, "DefineMode", False)

  do i=0,11
    fout->$VarNames(i)$ = FILE_SALT(i,:,:,:)
  end do

end

