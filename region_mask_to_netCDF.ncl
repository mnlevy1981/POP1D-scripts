begin

  luse_g16 = False
  setfileoption("bin", "ReadByteOrder", "BigEndian")
  qm        = str_get_dq()
  r2d       = 90.d/acos(0.d)
  curdate   = systemfunc("date "+qm+"+%y-%m-%d %H:%M:%S"+qm)
  longdate  = systemfunc("date -d "+qm+curdate+qm)
  shortdate = systemfunc("date -d "+qm+curdate+qm+" +%y%m%d")

  if (luse_g16) then
    grid = "gx1v6"
  else
    grid = "gx3v7"
  end if
  outfile = "region_masks/"+grid+"/region_mask_"+shortdate+".nc"
  print((/outfile/))

  ; Read in region mask (binary file)
  if (luse_g16) then
    fgrid = addfile("/glade/p/cesmdata/cseg/mapping/grids/gx1v6_090205.nc", "r")
    mask_file = "/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx1v6/grid/region_mask_20090205.ieeei4"
  else
    fgrid = addfile("/glade/p/cesmdata/cseg/mapping/grids/gx3v7_120309.nc", "r")
    mask_file = "/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx3v7/grid/region_mask_20090831.ieeei4"
  end if
  dimSizes=fgrid->grid_dims
  region_mask=fbindirread(mask_file, 0, (/dimSizes(1), dimSizes(0)/), "integer")
  region_mask@coordinates = "longitude latitude"
  region_mask@_FillValue = default_fillvalue("integer")
  region_mask=where(region_mask.eq.0, default_fillvalue("integer"), region_mask)
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


  ; Output region mask (integer)
  system("/bin/rm -f "+outfile)
  fout = addfile(outfile, "c")
  setfileoption(fout, "DefineMode", True)

  fAtt = True
  fAtt@title = "Region Mask for "+grid
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"X", "Y"/)
  dimUnlim = (/False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "latitude",    typeof(lats), (/"Y", "X"/))
  filevardef(fout, "longitude",   typeof(lons), (/"Y", "X"/))
  filevardef(fout, "region_mask", "integer",                      (/"Y", "X"/))

  filevarattdef(fout, "latitude",  lats)
  filevarattdef(fout, "longitude", lons)
  filevarattdef(fout, "region_mask", region_mask)

  setfileoption(fout, "DefineMode", False)

  fout->latitude  = reshape(lats, (/dimSizes(1), dimSizes(0)/))
  fout->longitude = reshape(lons, (/dimSizes(1), dimSizes(0)/))
  fout->region_mask = region_mask

end
