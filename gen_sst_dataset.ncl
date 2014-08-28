begin

  ; User-defined Parameter
  ; luse_g16 = False => use g37 grid
  luse_g16 = False
  longdate = systemfunc("date")

  gx1v6_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx1v6_090205.nc"
  gx3v7_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx3v7_120309.nc"
  r2d       = 90.d/acos(0.d)

  if (luse_g16) then
    fgrid = addfile(gx1v6_gridfile,"r")
    fin   = addfile("/glade/p/cgd/oce/people/yeager/OBS/phc/PHC2_TEMP_gx1v6.nc","r")
    outfile = "datasets/PHC2_SST_mod_gx1v6.nc"
  else
    fgrid = addfile(gx3v7_gridfile,"r")
    fin   = addfile("/glade/p/cgd/oce/people/yeager/OBS/phc/PHC2_TEMP_gx3v7.nc","r")
    outfile = "datasets/PHC2_SST_mod_gx3v7.nc"
  end if
  TEMP_tmp = fin->TEMP
  TEMP_tmp@long_name = "Sea Surface Temperature"
  TEMP_tmp@coordinates = "longitude latitude"
  TEMP_tmp@_FillValue = 273.15
  lons = fgrid->grid_center_lon
  lats = fgrid->grid_center_lat
  if (lons@units.eq."radians") then
    lons=lons*r2d
    lons@units="degrees_east"
  end if
  if (lats@units.eq."radians") then
    lats=lats*r2d
    lats@units="degrees_north"
  end if

  ; Compute Parameters
  dimSizes = new(3, "integer")
  dimSizes(0:1) = fgrid->grid_dims
  dimSizes(2)   = 60

  ; Variable setup
  VarNames = "SST"+ispan(1,12,1)
  VarNames(0:8) = "SST0"+ispan(1,9,1)

  ; Output
  system("/bin/rm -f "+outfile)
  fout = addfile(outfile, "c")
  setfileoption(fout, "DefineMode", True)

  fAtt = True
  fAtt@title = "Monthly Levitus Climatology (from Steve Yeager)"
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"X", "Y"/)
  dimUnlim = (/False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "latitude",  typeof(lats), (/"Y", "X"/))
  filevardef(fout, "longitude", typeof(lons), (/"Y", "X"/))
  filevarattdef(fout, "latitude",  lats)
  filevarattdef(fout, "longitude", lons)
  do i=0,11
    filevardef(fout, VarNames(i),  typeof(TEMP_tmp), (/"Y", "X"/))
    filevarattdef(fout, VarNames(i), TEMP_tmp)
  end do
  setfileoption(fout, "DefineMode", False)

  fout->latitude  = reshape(lats, (/dimSizes(1), dimSizes(0)/))
  fout->longitude = reshape(lons, (/dimSizes(1), dimSizes(0)/))
  do i=0,11
    fout->$VarNames(i)$ = TEMP_tmp(i,0,:,:)
  end do

end

