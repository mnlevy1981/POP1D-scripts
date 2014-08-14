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

  ; Box(es) to turn mask on in
  min_lon = (/281.d, 300.d,   0.d/)
  max_lon = (/343.d, 360.d,  20.d/)
  min_lat = (/ 9.1d, -15.d, -15.d/)
  max_lat = (/ 50.d,  9.1d,   9.d/)
  ; Transition region(s)
  north_buf = (/ 5.d, 0.d, 0.d/)
  east_buf  = (/15.d, 0.d, 0.d/)
  south_buf = (/ 0.d, 5.d, 5.d/)
  west_buf  = (/10.d, 0.d, 0.d/)
  ne_buf = north_buf*east_buf.ne.0
  nw_buf = north_buf*west_buf.ne.0
  sw_buf = south_buf*west_buf.ne.0
  se_buf = south_buf*east_buf.ne.0
  ; Nothing in Hudson Bay
  nw_buf(0) = False
  
  num_blocks = dimsizes(min_lon)
  if (any((/dimsizes(max_lon), dimsizes(min_lat), dimsizes(max_lat)/).ne.num_blocks)) then
    print((/"ERROR: min_lon, max_lon, min_lat, and max_lat must all be same dimension"/))
    status_exit(1)
  end if

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
  lons = fin->grid_center_lon
  lats = fin->grid_center_lat
  if (lons@units.eq."radians") then
    lons=lons*r2d
    lons@units="degrees_east"
  end if
  if (lats@units.eq."radians") then
    lats=lats*r2d
    lats@units="degrees_north"
  end if

  ; Compute Parameters
  dimSizes = fin->grid_dims
  gridmask = fin->grid_imask
  delete(fin)
  alpha = new(dimSizes(0)*dimSizes(1), "float")
  alpha@long_name = "Weight used when averaging POP_SST and DOCN_SST"
  alpha@coordinates = "longitude latitude"
  alpha@_FillValue = default_fillvalue("float")
  alpha = 0.
  do i=0,num_blocks-1
    inbox = ind((lons.ge.min_lon(i)).and. \
                (lons.le.max_lon(i)).and. \
                (lats.ge.min_lat(i)).and. \
                (lats.le.max_lat(i)))
    alpha(inbox) = 1.
    delete(inbox)
    ; Transition Regions
    ; North
    if (north_buf(i).gt.0.d) then
      transition = ind((lats.le.max_lat(i)+north_buf(i)).and. \
                       (lats.gt.max_lat(i)).and.              \
                       (lons.ge.min_lon(i)).and.              \
                       (lons.le.max_lon(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (lats(transition) - max_lat(i))/north_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if
    ; East
    if (east_buf(i).gt.0.d) then
      transition = ind((lons.le.max_lon(i)+east_buf(i)).and. \
                       (lons.gt.max_lon(i)).and.             \
                       (lats.ge.min_lat(i)).and.             \
                       (lats.le.max_lat(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (lons(transition) - max_lon(i))/east_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if
    ; South
    if (south_buf(i).gt.0.d) then
      transition = ind((lats.ge.min_lat(i)-south_buf(i)).and. \
                       (lats.lt.min_lat(i)).and.              \
                       (lons.ge.min_lon(i)).and.              \
                       (lons.le.max_lon(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (min_lat(i) - lats(transition))/south_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if
    ; West
    if (west_buf(i).gt.0.d) then
      transition = ind((lons.ge.min_lon(i)-west_buf(i)).and. \
                       (lons.lt.min_lon(i)).and.             \
                       (lats.ge.min_lat(i)).and.             \
                       (lats.le.max_lat(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (min_lon(i) - lons(transition))/west_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if

    ; NW
    if (nw_buf(i)) then
      transition = ind((lons.ge.min_lon(i)-west_buf(i)).and. \
                       (lons.lt.min_lon(i)).and.             \
                       (lats.le.max_lat(i)+north_buf(i)).and. \
                       (lats.gt.max_lat(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (min_lon(i) - lons(transition))/west_buf(i)) * \
              dble2flt(1.d - (lats(transition) - max_lat(i))/north_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if
    ; NE
    if (ne_buf(i)) then
      transition = ind((lons.le.max_lon(i)+east_buf(i)).and. \
                       (lons.gt.max_lon(i)).and.             \
                       (lats.le.max_lat(i)+north_buf(i)).and. \
                       (lats.gt.max_lat(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (lons(transition) - max_lon(i))/east_buf(i)) * \
              dble2flt(1.d - (lats(transition) - max_lat(i))/north_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if
    ; SE
    if (se_buf(i)) then
      transition = ind((lons.le.max_lon(i)+east_buf(i)).and. \
                       (lons.gt.max_lon(i)).and.             \
                       (lats.ge.min_lat(i)-south_buf(i)).and. \
                       (lats.lt.min_lat(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (lons(transition) - max_lon(i))/east_buf(i)) * \
              dble2flt(1.d - (min_lat(i) - lats(transition))/south_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if
    ; SW
    if (sw_buf(i)) then
      transition = ind((lons.ge.min_lon(i)-west_buf(i)).and. \
                       (lons.lt.min_lon(i)).and.             \
                       (lats.ge.min_lat(i)-south_buf(i)).and. \
                       (lats.lt.min_lat(i)))
      if (.not.any(ismissing(transition))) then
        tmp = dble2flt(1.d - (min_lon(i) - lons(transition))/west_buf(i)) * \
              dble2flt(1.d - (min_lat(i) - lats(transition))/south_buf(i))
        alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
        delete(tmp)
      end if
      delete(transition)
    end if
  end do
  alpha(ind(gridmask.eq.0)) = alpha@_FillValue

  ; Output
  system("/bin/rm -f "+outfile)
  fout = addfile(outfile, "c")
  setfileoption(fout, "DefineMode", True)

  fAtt = True
  fAtt@title = "Potential Temperature Restoring Parameters"
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"X", "Y"/)
  dimUnlim = (/False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "latitude",  typeof(lats), (/"Y", "X"/))
  filevardef(fout, "longitude", typeof(lons), (/"Y", "X"/))
  filevardef(fout, "alpha",     "float",      (/"Y", "X"/))

  filevarattdef(fout, "latitude",  lats)
  filevarattdef(fout, "longitude", lons)
  filevarattdef(fout, "alpha",     alpha)

  setfileoption(fout, "DefineMode", False)

  fout->latitude  = reshape(lats,  (/dimSizes(1), dimSizes(0)/))
  fout->longitude = reshape(lons,  (/dimSizes(1), dimSizes(0)/))
  fout->alpha     = reshape(alpha, (/dimSizes(1), dimSizes(0)/))

end

