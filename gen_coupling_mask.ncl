; This script generates a mask (alpha) that ranges from 0 to 1 (inclusive)
; Rather than passing the coupler POP's SST field, we will instead pass
; alpha*POP_SST + (1-alpha)*DOCN_SST
;
; alpha = 1 => full POP cell (active ocean)
; alpha = 0 => full data ocean cell (non-active ocean)

begin

  ; User-defined Parameters
  ; luse_g16 = False => use g37 grid
  luse_g16 = False
  min_lon = 287.d
  max_lon = 355.d
  min_lat = -17.d
  max_lat =  50.d

  ; Constants
  pi_dble   = acos(0.d)*2.d
  r2d       = 180.d/pi_dble
  qm        = str_get_dq()
  longdate  = systemfunc("date")

  gx1v6_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx1v6_090205.nc"
  gx3v7_gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx3v7_120309.nc"

  if (luse_g16) then
    fin = addfile(gx1v6_gridfile,"r")
    outfile = "masks/gx1v6_coupling_mask.nc"
  else
    fin = addfile(gx3v7_gridfile,"r")
    outfile = "masks/gx3v7_coupling_mask.nc"
  end if

  ; Compute Parameters
  dimSizes = fin->grid_dims
  alpha = new(dimSizes(0)*dimSizes(1), "float")
  alpha@long_name = "Weight used when averaging POP_SST and DOCN_SST"
  alpha@coordinates = "longitude latitude"
  alpha@_FillValue = default_fillvalue("float")
  alpha = 0.
  inbox = ind((r2d*fin->grid_center_lon.ge.min_lon).and. \
              (r2d*fin->grid_center_lon.le.max_lon).and. \
              (r2d*fin->grid_center_lat.ge.min_lat).and. \
              (r2d*fin->grid_center_lat.le.max_lat))
  alpha(inbox) = 1.
  alpha(ind(fin->grid_imask.eq.0)) = alpha@_FillValue

  ; Output
  system("/bin/rm -f "+outfile)
  fout = addfile(outfile, "c")
  setfileoption(fout, "DefineMode", True)

  fAtt = True
  fAtt@title = "Potential Temperature Restoring Parameters"
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"X", "Y"/)
  dimSizes = fin->grid_dims
  dimUnlim = (/False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "latitude",  typeof(fin->grid_center_lat), (/"Y", "X"/))
  filevardef(fout, "longitude", typeof(fin->grid_center_lon), (/"Y", "X"/))
  filevardef(fout, "alpha",     "float",                      (/"Y", "X"/))

  filevarattdef(fout, "latitude",  fin->grid_center_lat)
  filevarattdef(fout, "longitude", fin->grid_center_lon)
  filevarattdef(fout, "alpha",     alpha)

  setfileoption(fout, "DefineMode", False)

  fout->latitude  = reshape(fin->grid_center_lat, (/dimSizes(1), dimSizes(0)/))
  fout->longitude = reshape(fin->grid_center_lon, (/dimSizes(1), dimSizes(0)/))
  fout->alpha     = reshape(alpha,                (/dimSizes(1), dimSizes(0)/))

end

