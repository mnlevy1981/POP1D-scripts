begin

  ; User-defined Parameters
  ; luse_g16 = False => use g37 grid
  luse_g16 = False
  min_lon =   0.d
  max_lon = 360.d
  min_lat =  90.d
  max_lat = -90.d

  ; Constants
  pi_dble   = acos(0.d)*2.d
  qm        = str_get_dq()
  curdate   = systemfunc("date "+qm+"+%y-%m-%d %H:%M:%S"+qm)
  longdate  = systemfunc("date -d "+qm+curdate+qm)
  shortdate = systemfunc("date -d "+qm+curdate+qm+" +%y%m%d-%H%M%S")

  gx1v6_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx1v6_090205.nc"
  gx3v7_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx3v7_120309.nc"

  if (luse_g16) then
    fin = addfile(gx1v6_gridfile,"r")
    outfile = "param_sets/gx1v6_s_restore_"+shortdate+".nc"
  else
    fin = addfile(gx3v7_gridfile,"r")
    outfile = "param_sets/gx3v7_s_restore_"+shortdate+".nc"
  end if

  ; Compute Parameters
  dimSizes = fin->grid_dims
  MAX_LEV = new(dimSizes(0)*dimSizes(1), "integer")
  MAX_LEV@coordinates = "longitude latitude"
  MAX_LEV(ind(fin->grid_imask.eq.1)) = 60

  RTAU = new(dimSizes(0)*dimSizes(1), "double")
  RTAU@coordinates = "longitude latitude"
  RTAU(ind(fin->grid_imask.eq.1)) = 0.1d

  ; Output
  system("/bin/rm -f "+outfile)
  fout = addfile(outfile, "c")
  setfileoption(fout, "DefineMode", True)

  fAtt = True
  fAtt@title = "Salinity Restoring Parameters"
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"X", "Y"/)
  dimSizes = fin->grid_dims
  dimUnlim = (/False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "latitude",  typeof(fin->grid_center_lat), (/"Y", "X"/))
  filevardef(fout, "longitude", typeof(fin->grid_center_lon), (/"Y", "X"/))
  filevardef(fout, "S_RESTORE_MAX_LEVEL", "integer",          (/"Y", "X"/))
  filevardef(fout, "S_RESTORE_RTAU",      "double",           (/"Y", "X"/))

  filevarattdef(fout, "latitude",  fin->grid_center_lat)
  filevarattdef(fout, "longitude", fin->grid_center_lon)
  filevarattdef(fout, "S_RESTORE_MAX_LEVEL", MAX_LEV) 
  filevarattdef(fout, "S_RESTORE_RTAU",      RTAU) 

  setfileoption(fout, "DefineMode", False)

  fout->latitude  = reshape(fin->grid_center_lat, (/dimSizes(1), dimSizes(0)/))
  fout->longitude = reshape(fin->grid_center_lon, (/dimSizes(1), dimSizes(0)/))
  fout->S_RESTORE_MAX_LEVEL = reshape(MAX_LEV, (/dimSizes(1), dimSizes(0)/))
  fout->S_RESTORE_RTAU      = reshape(RTAU,    (/dimSizes(1), dimSizes(0)/))

end

