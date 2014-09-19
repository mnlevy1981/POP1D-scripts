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
  outfile = "KMT/"+grid+"/KMT_"+shortdate+".nc"
  print((/outfile/))

  ; Read in region mask (binary file)
  if (luse_g16) then
    fgrid = addfile("/glade/p/cesmdata/cseg/mapping/grids/gx1v6_090205.nc", "r")
    mask_file = "/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx1v6/grid/region_mask_20090205.ieeei4"
  else
    fgrid = addfile("/glade/p/cesmdata/cseg/mapping/grids/gx3v7_120309.nc", "r")
    mask_file = "/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx3v7/grid/topography_20100105.ieeei4"
  end if
  dimSizes=fgrid->grid_dims
  KMT=fbindirread(mask_file, 0, (/dimSizes(1), dimSizes(0)/), "integer")
  KMT@coordinates = "longitude latitude"
  KMT@_FillValue = default_fillvalue("integer")
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
  fAtt@title = "Number of levels for "+grid
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"X", "Y"/)
  dimUnlim = (/False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "latitude",    typeof(lats), (/"Y", "X"/))
  filevardef(fout, "longitude",   typeof(lons), (/"Y", "X"/))
  filevardef(fout, "KMT",         "integer",    (/"Y", "X"/))

  filevarattdef(fout, "latitude",  lats)
  filevarattdef(fout, "longitude", lons)
  filevarattdef(fout, "KMT"      , KMT)

  setfileoption(fout, "DefineMode", False)

  fout->latitude  = reshape(lats, (/dimSizes(1), dimSizes(0)/))
  fout->longitude = reshape(lons, (/dimSizes(1), dimSizes(0)/))
  fout->KMT = KMT

end
