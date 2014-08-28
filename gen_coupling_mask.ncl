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
  ; luse_region_mask = False => use bounding box[es]
  luse_region_mask = True

  ; Constants
  qm        = str_get_dq()
  longdate  = systemfunc("date")

  if (luse_g16) then
    outfile = "masks/gx1v6_coupling_mask.nc"
  else
    outfile = "masks/gx3v7_coupling_mask.nc"
  end if
  if (luse_region_mask) then
    docn_region = 6 ; Atlantic Domain!
    if (luse_g16) then
      mask_file = "region_masks/gx1v6/region_mask_140818.nc"
    else
      mask_file = "region_masks/gx3v7/region_mask_140818.nc"
    end if
    fin = addfile(mask_file, "r")
    lons = fin->longitude
    lats = fin->latitude
    alpha = where(fin->region_mask.eq.docn_region, 0.d, 1.d)
    alpha@_FillValue = -1.d
    dimSizes = dimsizes(alpha)
    ; Transition region: average neighboring alpha values
    alpha_tmp = new(dimSizes, "double")
    do k=1,3
      do i=0,dimSizes(0)-1
        do j=0,dimSizes(1)-1
          idims = mod(dimSizes(0)+ispan(i-1,i+1,1),dimSizes(0))
          jdims = mod(dimSizes(1)+ispan(j-1,j+1,1),dimSizes(1))
          if ((.not.ismissing(alpha(i,j))).and.                           \
              (.not.all(ismissing(alpha(idims, jdims))))) then
            alpha_tmp(i,j) = avg(alpha(idims, jdims))
          end if
        end do
      end do
      alpha = (/alpha_tmp/)
      ; Fix Panama (no overflow to Pacific!)
      if ((docn_region.eq.6).and.(.not.luse_g16)) then
        alpha(61, 86) = 1.d
        alpha(62, 87) = 0.d
      end if
    end do
  else
    ; Box(es) to turn mask on in
    west_edge = (/281.d, 300.d,   0.d/)
    east_edge = (/345.d, 360.d,  20.d/)
    south_edge = (/ 9.1d, -15.d, -15.d/)
    north_edge = (/ 50.d,  9.1d,   9.d/)
    ; Transition region(s)
    north_buf = (/ 5.d, 0.d, 0.d/)
    east_buf  = (/10.d, 0.d, 0.d/)
    south_buf = (/ 0.d, 5.d, 5.d/)
    west_buf  = (/10.d, 0.d, 0.d/)
    ne_buf = north_buf*east_buf.ne.0
    nw_buf = north_buf*west_buf.ne.0
    sw_buf = south_buf*west_buf.ne.0
    se_buf = south_buf*east_buf.ne.0
    ; Nothing in Hudson Bay
    nw_buf(0) = False
  
    num_blocks = dimsizes(west_edge)
    if (any((/dimsizes(east_edge), dimsizes(south_edge), dimsizes(north_edge)/).ne.num_blocks)) then
      print((/"ERROR: west_edge, east_edge, south_edge, and north_edge must all be same dimension"/))
      status_exit(1)
    end if

    ; Constants
    pi_dble   = acos(0.d)*2.d
    r2d       = 180.d/pi_dble

    if (luse_g16) then
      gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx1v6_090205.nc"
    else
      gridfile = "/glade/p/cesmdata/cseg/mapping/grids/gx3v7_120309.nc"
    end if
    fin = addfile(gridfile, "r")

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
    dimSizes = fin->grid_dims((/1,0/))
    gridmask = fin->grid_imask
    delete(fin)

    alpha = new(dimSizes(0)*dimSizes(1), "double")
    alpha@long_name = "Weight used when averaging POP_SST and DOCN_SST"
    alpha@coordinates = "longitude latitude"
    alpha@_FillValue = -1.d
    alpha = 0.d
    do i=0,num_blocks-1
      inbox = ind((lons.ge.west_edge(i)).and. \
                  (lons.le.east_edge(i)).and. \
                  (lats.ge.south_edge(i)).and. \
                  (lats.le.north_edge(i)))
      alpha(inbox) = 1.d
      delete(inbox)
      ; Transition Regions
      ; North
      if (north_buf(i).gt.0.d) then
        transition = ind((lats.le.north_edge(i)+north_buf(i)).and. \
                         (lats.gt.north_edge(i)).and.              \
                         (lons.ge.west_edge(i)).and.              \
                         (lons.le.east_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = 1.d - (lats(transition) - north_edge(i))/north_buf(i)
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if
      ; East
      if (east_buf(i).gt.0.d) then
        transition = ind((lons.le.east_edge(i)+east_buf(i)).and. \
                         (lons.gt.east_edge(i)).and.             \
                         (lats.ge.south_edge(i)).and.             \
                         (lats.le.north_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = 1.d - (lons(transition) - east_edge(i))/east_buf(i)
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if
      ; South
      if (south_buf(i).gt.0.d) then
        transition = ind((lats.ge.south_edge(i)-south_buf(i)).and. \
                         (lats.lt.south_edge(i)).and.              \
                         (lons.ge.west_edge(i)).and.              \
                         (lons.le.east_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = 1.d - (south_edge(i) - lats(transition))/south_buf(i)
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if
      ; West
      if (west_buf(i).gt.0.d) then
        transition = ind((lons.ge.west_edge(i)-west_buf(i)).and. \
                         (lons.lt.west_edge(i)).and.             \
                         (lats.ge.south_edge(i)).and.             \
                         (lats.le.north_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = 1.d - (west_edge(i) - lons(transition))/west_buf(i)
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if

      ; NW
      if (nw_buf(i)) then
        transition = ind((lons.ge.west_edge(i)-west_buf(i)).and. \
                         (lons.lt.west_edge(i)).and.             \
                         (lats.le.north_edge(i)+north_buf(i)).and. \
                         (lats.gt.north_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = (1.d - (west_edge(i) - lons(transition))/west_buf(i)) * \
                (1.d - (lats(transition) - north_edge(i))/north_buf(i))
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if
      ; NE
      if (ne_buf(i)) then
        transition = ind((lons.le.east_edge(i)+east_buf(i)).and. \
                         (lons.gt.east_edge(i)).and.             \
                         (lats.le.north_edge(i)+north_buf(i)).and. \
                         (lats.gt.north_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = (1.d - (lons(transition) - east_edge(i))/east_buf(i)) * \
                (1.d - (lats(transition) - north_edge(i))/north_buf(i))
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if
      ; SE
      if (se_buf(i)) then
        transition = ind((lons.le.east_edge(i)+east_buf(i)).and. \
                         (lons.gt.east_edge(i)).and.             \
                         (lats.ge.south_edge(i)-south_buf(i)).and. \
                         (lats.lt.south_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = (1.d - (lons(transition) - east_edge(i))/east_buf(i)) * \
                (1.d - (south_edge(i) - lats(transition))/south_buf(i))
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if
      ; SW
      if (sw_buf(i)) then
        transition = ind((lons.ge.west_edge(i)-west_buf(i)).and. \
                         (lons.lt.west_edge(i)).and.             \
                         (lats.ge.south_edge(i)-south_buf(i)).and. \
                         (lats.lt.south_edge(i)))
        if (.not.any(ismissing(transition))) then
          tmp = (1.d - (west_edge(i) - lons(transition))/west_buf(i)) * \
                (1.d - (south_edge(i) - lats(transition))/south_buf(i))
          alpha(transition) = dim_max_n((/alpha(transition),tmp/), 0)
          delete(tmp)
        end if
        delete(transition)
      end if
    end do
    alpha(ind(gridmask.eq.0)) = alpha@_FillValue
    ; Note: I had the mask backwards! Want 0 in Atlantic, 1 elsewhere!
    alpha = 1.d - alpha
  end if

  ; Output
  system("/bin/rm -f "+outfile)
  fout = addfile(outfile, "c")
  setfileoption(fout, "DefineMode", True)

  fAtt = True
  fAtt@title = "Potential Temperature Restoring Parameters"
  fAtt@creation_date = longdate
  fileattdef(fout, fAtt)

  dimNames = (/"Y", "X"/)
  dimUnlim = (/False, False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "latitude",  typeof(lats), (/"Y", "X"/))
  filevardef(fout, "longitude", typeof(lons), (/"Y", "X"/))
  filevardef(fout, "alpha",     "double",     (/"Y", "X"/))

  filevarattdef(fout, "latitude",  lats)
  filevarattdef(fout, "longitude", lons)
  filevarattdef(fout, "alpha",     alpha)

  setfileoption(fout, "DefineMode", False)

  fout->latitude  = reshape(lats,  dimSizes)
  fout->longitude = reshape(lons,  dimSizes)
  fout->alpha     = reshape(alpha, dimSizes)

end

