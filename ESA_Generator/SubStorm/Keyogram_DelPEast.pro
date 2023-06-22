
;;;------------- SIMPLE Code to access and plot
;;;------------- 3-D GITM/GDC Data Cubes
;;; Put your custom directory here for the DataCube:
;restore, '../ProcessedDataCubes/GITM_DataCube.save'
restore, '/Users/jmbell4/GDC_CodeTesting/SmallerGITMSims/SMC/ProcessedDataCubes/GITM_ExtendedDataCube.save'
;;; KEY VARIABLES:
;;; nAlts:  GITM Alt Points
;;; nLats:  GITM Lat Points
;;; nLons:  GITM Lon Points
;;; nGITMFiles:  Number of Simulation Files
;;;
;;; GITMDataCubeAlts  = fltarr(nGITMFiles,nAlts)
;;; GITMDataCubeLats  = fltarr(nGITMFiles,nLats)
;;; GITMDataCubeLons  = fltarr(nGITMFiles,nLons)
;;;-----
;;; --- Neutral Temps
;;; GITMDataCubeTemps = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; --- Neutral Densities
;;; GITMDataCubeRho   = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeN2    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeO2    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeO     = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeHe    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeMBar  = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeNDenTotal  = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; --- Neutral Winds
;;; GITMDataCubeUn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeVn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeWn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; --- Ion Winds
;;; GITMDataCubeUi    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeVi    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeWi    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeMBarI = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeIDenTotal  = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; ---Time Variables
;;; GITMDataCubeUT          = strarr(nGITMFiles)
;;; GITMDataCubeDate        = strarr(nGITMFiles)
;;; GITMDataCubeLocalTime   = strarr(nGITMFiles,nLons)
;;; GITMDataCubeJulianDate  = strarr(nGITMFiles)
;;; GITMDataCubeJulian2000  = strarr(nGITMFiles)
;;; -----
;;; FORCING TERMS
;;; Ion Drag (m/s^2) on the neutrals
;;; GITMDataCubeIonDragUn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeIonDragVn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeIonDragWn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; Coriolis (m/s^2) on the neutrals
;;; GITMDataCubeCorUn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeCorVn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeCorWn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; Viscosity--Vertical (m/s^2) on the neutrals
;;; GITMDataCubeViscUn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeViscVn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeViscWn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; Pressure Gradient -- (m/s^2) on the neutrals
;;; GITMDataCubeDelPUn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeDelPVn    = fltarr(nGITMFiles,nLons,nLats,nAlts)
;;; GITMDataCubeDelPWn    = fltarr(nGITMFiles,nLons,nLats,nAlts)

SelectLats  = [50.0, 60.0, 70.0, 80.0]
SelectLTs   = [6.0, 12.0, 18.0, 24.0]
;;; identify indicies closest to our Lons
nSelLTs = n_elements(SelectLTs)
nSelLats = n_elements(SelectLats)
iLTIndex = intarr(nTimes, nSelLTs)
iLatIndex = intarr(nTimes, nSelLats)

for iGITM = 0, nTimes - 1 do begin
   for i = 0, nSelLTs - 1 do begin
      DeltaLTs = abs( reform(GITMDataCubeLocalTimes[iGITM,*]) - SelectLTs[i])
       MinDelta = min(DeltaLTs,MinIndex)
       iLTIndex[iGITM,i] = MinIndex
   endfor ;i = 0, n_elements(SelectLons) do begin

   for i = 0, nSelLats - 1 do begin
      DeltaLats = abs( reform(GITMDataCubeLats[iGITM,*]) - SelectLats[i])
       MinDelta = min(DeltaLats,MinIndex)
       iLatIndex[iGITM,i] = MinIndex
   endfor ;i = 0, n_elements(SelectLons) do begin

endfor 

;;; Next, Sub-select the data based upon our chosen "orbits"

GITMOrbitTemp       = fltarr(nTimes,nSelLTs,nSelLats,nAlts)
GITMOrbitRho        = fltarr(nTimes,nSelLTs,nSelLats,nAlts)
GITMOrbitUn         = fltarr(nTimes,nSelLTs,nSelLats,nAlts)
GITMOrbitVn         = fltarr(nTimes,nSelLTs,nSelLats,nAlts)
GITMOrbitWn         = fltarr(nTimes,nSelLTs,nSelLats,nAlts)

;;; GITMDataCubeDelPUn    = fltarr(nTimes,nLons,nLats,nAlts)
GITMOrbitDelPUn         = fltarr(nTimes,nSelLTs,nSelLats,nAlts)
GITMOrbitDelPVn         = fltarr(nTimes,nSelLTs,nSelLats,nAlts)
GITMOrbitDelPWn         = fltarr(nTimes,nSelLTs,nSelLats,nAlts)

GITMOrbitTime       = strarr(nTimes)
GITMOrbitDate       = strarr(nTimes)
GITMOrbitLT       = fltarr(nTimes,nSelLTs,nSelLats)

;;; nLats:  GITM Lat Points
;;; nTimes:  Number of Simulation Files
for iGITM = 0, nTimes - 1 do begin
    GITMOrbitTime[iGITM] = GITMDataCubeUT[iGITM]
    GITMOrbitDate[iGITM] = GITMDataCubeDate[iGITM]
    for iOrbitLT = 0, nSelLTs - 1 do begin
      for iOrbitLat = 0, nSelLats - 1 do begin

          GITMOrbitTemp[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeTemp[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

          GITMOrbitRho[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeRho[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

          GITMOrbitUn[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeUn[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

          GITMOrbitVn[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeVn[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

          GITMOrbitWn[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeWn[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

          GITMOrbitLT[iGITM,iOrbitLT,iOrbitLat] = $
            GITMDataCubeLocalTimes[iGITM,iLTIndex[iGITM,iOrbitLT]]

          GITMOrbitDelPUn[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeDelPUn[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

          GITMOrbitDelPVn[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeDelPVn[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

          GITMOrbitDelPWn[iGITM,iOrbitLT,iOrbitLat,0:nAlts-1] = $
            GITMDataCubeDelPWn[iGITM,iLTIndex[iGITM,iOrbitLT], $
                  iLatIndex[iGITM,iOrbitLat], 0:nAlts-1]

      endfor ;iOrbitLat = 0, nSelLats - 1 do begin
    endfor ;iOrbitLT = 0, nSelLTs - 1 do begin
endfor ;iGITM = 0, nTimes - 1 do begin

LocalTimeString = ['Dawn', 'Noon', 'Dusk', 'Midnight']
LocalTimeColorString = ['orange', 'firebrick', 'blue', 'black']
LatString = strarr(nSelLats)
for iLat = 0, nSelLats - 1 do begin
    if(SelectLats[iLat] ge 0.0) then begin
      LatString[iLat] = strmid(string(SelectLats[iLat]),6,4) + ' N'
    endif else begin
      LatString[iLat] = strmid(string(SelectLats[iLat]),6,4) + ' S'
    endelse
endfor 

;;; Create 1-page per time-stamp
nPages = nSelLats
for iPage = 0, nPages - 1 do begin
  
    iSelLat = iPage
    iSelLT  = 0
    SubAlts = reform(GITMDataCubeAlts[0,*])
    SubAltIndex = where(SubAlts ge 200.0 and SubAlts le 505.0)

    GITMAlts = SubAlts
    GITMTimes = GITMOrbitTime[0:nTimes-1]

    ; Calibrate our Colorbar
    ;sample_var = reform(GITMOrbitTemp[0:nTimes-1,*,*,SubAltIndex])
    ;levels = 20
    ;min_data = min(sample_var)
    ;max_data = max(sample_var)
    ;nCLevels = floor( (max_data - min_data)/levels) + 1
    ;ctable = COLORTABLE(73, ncolors = levels, /REVERSE)

    ; Set overall panel sizes
    panel_height = 0.18
    panel_width  = 0.70

    ; set global positions
    plot_x0 = 0.125
    plot_y0 = 0.150
    plot_y1 = 0.90

    deltay = 0.025

    ; cycle through and plot the panels
    panel_x0 = plot_x0
    panel_x1 = panel_x0 + panel_width
    ct = COLORTABLE(72, /reverse)
    npanels = nSelLTs

    ;;; Create time axes using GITM Julian Dates    
    ;;; This is an idl-specific method
    date_time = TIMEGEN(nTimes,STEP_SIZE = 5, UNITS='Minutes', $
                          START=JULDAY(2,23,1997,0,0,0), $
                          FINAL=JULDAY(2,23,1997,6,0,0))
    date_label = LABEL_DATE(DATE_FORMAT = ['%H:%I'])
    

    event_start = JULDAY(2, 23, 1997, 1,45,0)
    event_end   = JULDAY(2, 23, 1997, 4,45,0)


    for ipanel = 0, npanels - 1 do begin
       panel_y0 =  plot_y0 + (npanels-1 - ipanel)*(panel_height) + (npanels-2-ipanel)*deltay
       panel_y1 = panel_y0 + panel_height

       ;contour_var = reform(GITMOrbitDelPUn[0:nTimes-1,ipanel,iSelLat,SubAltIndex])
       contour_var = reform(GITMOrbitDelPUn[0:nTimes-1,ipanel,iSelLat,SubAltIndex])
       ;contour_var2 = reform(GITMOrbitDelPVn[0:nTimes-1,ipanel,iSelLat,SubAltIndex])
       ;contour_var3 = reform(GITMOrbitDelPWn[0:nTimes-1,ipanel,iSelLat,SubAltIndex])

    stop
       y_var = GITMAlts[SubAltIndex]
       ;x_var = indgen(nTimes)
       ;x_var = GITMDataCubeUT
       x_var = date_time

       min_x = min(x_var)
       max_x = max(x_var)

       min_y = min(y_var)
       max_y = max(y_var)

       if(ipanel eq 0) then begin
          panel0 = CONTOUR( contour_var, x_var, y_var, $
                           xstyle = 1, xrange = [min_x, max_x], $
                           ystyle = 1, yrange = [min_y, max_y], $
                           font_size = 13, /FILL, $ 
                           xtickformat = 'LABEL_DATE', $
                           xtickunits = 'TIME', $
                           RGB_TABLE = ct, $
                           xshowtext = 0, $
                           position = [panel_x0, panel_y0, panel_x1, panel_y1], $
                           title = ' Latitude: ' + LatString[iPage])

          pline = polyline([event_start, event_start], [min_y, max_y], $
                   target = panel0, color = !COLOR.Black, linestyle = 2, thick = 2.0, /data)

          pline = polyline([event_end, event_end], [min_y, max_y], $
                   target = panel0, color = !COLOR.Black, linestyle = 3, thick = 2.0, /data)
       endif else begin
          if(ipanel eq npanels - 1) then begin
             panel0 = CONTOUR( contour_var, x_var, y_var, $
                              xstyle = 1, xrange = [min_x, max_x], $
                              ystyle = 1, yrange = [min_y, max_y], $
                              font_size = 13, /FILL, $ 
                              xtickformat = 'LABEL_DATE', $
                              xtickunits = 'TIME', $
                              RGB_TABLE = ct, /CURRENT, $
                              position = [panel_x0, panel_y0, panel_x1, panel_y1])

          pline = polyline([event_start, event_start], [min_y, max_y], $
                   target = panel0, color = !COLOR.Black, linestyle = 2, thick = 2.0, /data)

          pline = polyline([event_end, event_end], [min_y, max_y], $
                   target = panel0, color = !COLOR.Black, linestyle = 3, thick = 2.0, /data)
          endif else begin
             panel0 = CONTOUR( contour_var, x_var, y_var, $
                              xstyle = 1, xrange = [min_x, max_x], $
                              ystyle = 1, yrange = [min_y, max_y], $
                              font_size = 13, /FILL, $ 
                              xtickformat = 'LABEL_DATE', $
                              xtickunits = 'TIME', $
                              xtitle = 'UT (hh:mm)', $
                              RGB_TABLE = ct, /CURRENT, $
                              xshowtext = 0, $
                              position = [panel_x0, panel_y0, panel_x1, panel_y1])

          pline = polyline([event_start, event_start], [min_y, max_y], $
                   target = panel0, color = !COLOR.Black, linestyle = 2, thick = 2.0, /data)
          pline = polyline([event_end, event_end], [min_y, max_y], $
                   target = panel0, color = !COLOR.Black, linestyle = 3, thick = 2.0, /data)
          endelse
       endelse
       text = TEXT(panel_x0, panel_y1 + 0.00, LocalTimeString[ipanel],$
                   font_size = 14, $
                   color = LocalTimeColorString[ipanel])


    endfor ;ipanel = 0, npanels - 1 do begin
    cb_x0 = panel_x1 + 0.110
    cb_x1 = panel_x1 + 0.16
    cb_y0 = 0.25
    cb_y1 = 0.75
    cb = COLORBAR( title = '$(\nabla P/\rho)_{East}$ (m/s$^{2}$)', orientation = 1, $
          font_size = 12, $
          position = [cb_x0, cb_y0, cb_x1, cb_y1])

    text = TEXT(panel_x0 - 0.125, 0.50 - 0.05, 'Altitude (km)', $
                font_size = 14, orientation = 90)
  filename = 'GDC_DelPEast_KeyOGram.pdf'; + '_Page_' + strmid(string(iPage),7,1)+'.pdf' 

  if(iPage eq 0) then begin
      panel0.SAVE, filename, border = 5, $
            resolution = 300, /TRANSPARENT, /append
  endif else begin
      panel0.SAVE, filename, border = 5, $
             resolution = 300, /TRANSPARENT, /append
  endelse
  
endfor ;iPanels = 0, nTimes - 1 do begin
panel0.SAVE, filename, /close


end
