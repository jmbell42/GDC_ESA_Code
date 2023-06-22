;;;------------- SIMPLE Code to access and plot
;;;------------- 3-D GITM/GDC Data Cubes
;;; Some Constants
AMU = 1.66054e-27    ;; kg
kb  = 1.38065e-23    ;; J/K
MassN2 =  28.0*AMU
MassO2 =  32.0*AMU
MassO  =  16.0*AMU
MassHe =   4.0*AMU
restore, '../../GitmSimulations/SubStorm/ProcessedDataCubes/GITM_ExtendedDataCube.save'
;;; KEY VARIABLES:
;;; nAlts:  GITM Alt Points
;;; nLats:  GITM Lat Points
;;; nLons:  GITM Lon Points
;;; nTimes:  Number of Timestamps during the Simulation
;;;
;;; GITMDataCubeAlts  = fltarr(nTimes,nAlts)
;;; GITMDataCubeLats  = fltarr(nTimes,nLats)
;;; GITMDataCubeLons  = fltarr(nTimes,nLons)
;;;-----
;;; GITMDataCubeTemps = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeRho   = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeN2    = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeO2    = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeO     = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeHe    = fltarr(nTimes,nLons,nLats,nAlts)
;;;
;;; GITMDataCubeUn    = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeVn    = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeWn    = fltarr(nTimes,nLons,nLats,nAlts)
;;;
;;; GITMDataCubeUi    = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeVi    = fltarr(nTimes,nLons,nLats,nAlts)
;;; GITMDataCubeWi    = fltarr(nTimes,nLons,nLats,nAlts)
;;;
;;; GITMDataCubeUT    = strarr(nTimes)
;;; GITMDataCubeDate  = strarr(nTimes)

print, 'GITM, lons, lats, alts = ', $
        nLons, nLats, nAlts
GITMTemp = dblarr(nTimes,nLons,nLats,nAlts)
GITMRho  = dblarr(nTimes,nLons,nLats,nAlts)
;;;
GITMN2  = dblarr(nTimes,nLons,nLats,nAlts)
GITMO2  = dblarr(nTimes,nLons,nLats,nAlts)
GITMO   = dblarr(nTimes,nLons,nLats,nAlts)
GITMHe  = dblarr(nTimes,nLons,nLats,nAlts)

GITMNe  = dblarr(nTimes,nLons,nLats,nAlts)
;;;
GITMUn  = dblarr(nTimes,nLons,nLats,nAlts)
GITMVn  = dblarr(nTimes,nLons,nLats,nAlts)
GITMWn  = dblarr(nTimes,nLons,nLats,nAlts)

GITMUi  = dblarr(nTimes,nLons,nLats,nAlts)
GITMVi  = dblarr(nTimes,nLons,nLats,nAlts)
GITMWi  = dblarr(nTimes,nLons,nLats,nAlts)

GITMLats = dblarr(nLats)
GITMLons = dblarr(nLons)
GITMAlts = dblarr(nAlts)
;;;;
;;;; GITMDataCubeUn[i,0:nLons-1,0:nLats-1,0:nAlts-1] 
;;;; SelectLons and SelectAlts are arrays of Longitudes and Altitudes
;;;; that we will fly through
;;;; This code imagines hypothetical satellites flying at each longitude
;;;; at each of the constant altitudes

SelectAlts = [400.0]

;;; identify indicies closest to our Lons
nSelLons = nLons
nSelAlts = n_elements(SelectAlts)

iAltIndex = intarr(nTimes, nSelAlts)

;;;; Find the Indices in GITM that are the closest to our simple emphemeris
;;;; points.  For a real orbit, you will need the actual physical location
;;;; as a function of time.
for iGITM = 0, nTimes - 1 do begin
   for i = 0, nSelAlts - 1 do begin
   DeltaAlts = abs( reform(GITMDataCubeAlts[*]) - SelectAlts[i])
    MinDelta = min(DeltaAlts,MinIndex)
    iAltIndex[iGITM,i] = MinIndex
    print, ' Select Alts = ', GITMDataCubealts[iAltIndex[iGITM,i]]
    AltRef = GITMDataCubeAlts[iAltIndex[iGITM,i]]
   endfor ;i = 0, n_elements(SelectAlts) do begin
endfor 


GITMLats = GITMDataCubeLats
GITMLons = GITMDataCubeLons
;;; Next, Sub-select the data based upon our chosen "orbits"
;;; These variables are our "along track" fields for each satellite 
;;; We can then plot these trajectories for all latitudes and times
GITMOrbitTemp = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitRho  = fltarr(nTimes,nLons,nLats,nSelAlts)

GITMOrbitN2  = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitO2  = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitO   = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitHe  = fltarr(nTimes,nLons,nLats,nSelAlts)

GITMOrbitNe  = fltarr(nTimes,nLons,nLats,nSelAlts)

GITMOrbitUn  = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitVn  = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitWn  = fltarr(nTimes,nLons,nLats,nSelAlts)

GITMOrbitUi  = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitVi  = fltarr(nTimes,nLons,nLats,nSelAlts)
GITMOrbitWi  = fltarr(nTimes,nLons,nLats,nSelAlts)

for iGITM = 0, nTimes - 1 do begin
    for iOrbitAlt = 0, nSelAlts - 1 do begin
         GITMOrbitTemp[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeTemp[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitRho[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeRho[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitN2[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeN2[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitO2[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeO2[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitO[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeO[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitHe[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeHe[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitNe[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeIDenTotal[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         ;;; NEUTRAL WINDS
         GITMOrbitUn[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeUn[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitVn[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeVn[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitWn[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeWn[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         ;;; ION WINDS
         GITMOrbitUi[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeUi[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitVi[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeVi[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]

         GITMOrbitWi[iGITM,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            GITMDataCubeWi[iGITM,0:nLons-1, $
                  0:nLats-1, iAltIndex[iGITM,iOrbitAlt]]
    endfor ;iOrbitAlt = 0, nSelAlts - 1 do begin
endfor ;iGITM = 0, nTimes - 1 do begin

GITMDataCubeNDenTotal = $
    GITMDataCubeN2 + GITMDataCubeO2 + GITMDataCubeO + GITMDataCubeHe
GITMDataCubePressure = $
    GITMDataCubeNDenTotal*kb*GITMDataCubeTemp


nGITMTimes = nTimes 
nGITMLons = nLons
nGITMLats = nLats
nGITMAlts = nAlts


restore, '../../MsisSimulationFiles/SubStorm/ProcessedDataCubes/MSIS_ExtendedDataCube.save'
;;; KEY VARIABLES:
;;; nAlts:  MSIS Alt Points
;;; nLats:  MSIS Lat Points
;;; nLons:  GITM Lon Points
;;; nTimes:  Number of Timestamps during the Simulation
;;;
;;; MSISDataCubeAlts  = fltarr(nTimes,nAlts)
;;; MSISDataCubeLats  = fltarr(nTimes,nLats)
;;; MSISDataCubeLons  = fltarr(nTimes,nLons)
;;;-----
;;; MSISDataCubeTemps = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeRho   = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeN2    = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeO2    = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeO     = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeHe    = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeUn    = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeVn    = fltarr(nTimes,nLons,nLats,nAlts)
;;; MSISDataCubeWn    = fltarr(nTimes,nLons,nLats,nAlts)
;;;
;;; MSISDataCubeUT    = strarr(nTimes)
;;; MSISDataCubeDate  = strarr(nTimes)
print, 'MSIS, Lons, Lats, Alts = ', $
        nLons, nLats, nAlts

MSISTemp = fltarr(nTimes,nLons,nLats,nAlts)
MSISRho  = fltarr(nTimes,nLons,nLats,nAlts)
;;;
MSISN2  = fltarr(nTimes,nLons,nLats,nAlts)
MSISO2  = fltarr(nTimes,nLons,nLats,nAlts)
MSISO   = fltarr(nTimes,nLons,nLats,nAlts)
MSISHe  = fltarr(nTimes,nLons,nLats,nAlts)

MSISNe  = fltarr(nTimes,nLons,nLats,nAlts)
;;;
MSISUn  = fltarr(nTimes,nLons,nLats,nAlts)
MSISVn  = fltarr(nTimes,nLons,nLats,nAlts)
MSISWn  = fltarr(nTimes,nLons,nLats,nAlts)
;;;;
;;;; SelectLons and SelectAlts are arrays of Longitudes and Altitudes
;;;; that we will fly through
;;;; This code imagines hypothetical satellites flying at each longitude
;;;; at each of the constant altitudes

;;; identify indicies closest to our Lons
nSelLons = nLons
nSelAlts = n_elements(SelectAlts)

iAltIndex = intarr(nTimes, nSelAlts)

;;;; Find the Indices in MSIS that are the closest to our simple emphemeris
;;;; points.  For a real orbit, you will need the actual physical location
;;;; as a function of time.
for iMSIS = 0, nTimes - 1 do begin
   for i = 0, nSelAlts - 1 do begin
   DeltaAlts = abs( reform(MSISDataCubeAlts[*]) - SelectAlts[i])
    MinDelta = min(DeltaAlts,MinIndex)
    iAltIndex[iMSIS,i] = MinIndex
   endfor ;i = 0, n_elements(SelectAlts) do begin
endfor 

for iMSIS = 0, nTimes - 1 do begin
   for i = 0, nSelAlts - 1 do begin
    print, ' MSIS Select Alts = ', MSISDataCubealts[iAltIndex[iMSIS,i]]
   endfor ;i = 0, n_elements(SelectAlts) do begin
endfor 

;;; Next, Sub-select the data based upon our chosen "orbits"
;;; These variables are our "along track" fields for each satellite 
;;; We can then plot these trajectories for all latitudes and times
MSISOrbitTemp = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitRho  = fltarr(nTimes,nLons,nLats,nSelAlts)

MSISOrbitN2  = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitO2  = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitO   = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitHe  = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitNe  = fltarr(nTimes,nLons,nLats,nSelAlts)

MSISOrbitUn  = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitVn  = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitWn  = fltarr(nTimes,nLons,nLats,nSelAlts)

MSISOrbitUi  = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitVi  = fltarr(nTimes,nLons,nLats,nSelAlts)
MSISOrbitWi  = fltarr(nTimes,nLons,nLats,nSelAlts)


for iMSIS = 0, nTimes - 1 do begin
    for iOrbitAlt = 0, nSelAlts - 1 do begin
         MSISOrbitTemp[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeTemp[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitRho[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeRho[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitN2[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeN2[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitO2[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeO2[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitO[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeO[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitHe[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeHe[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitNe[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeIDenTotal[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         ;;; NEUTRAL WINDS
         MSISOrbitUn[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeUn[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitVn[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeVn[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitWn[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeWn[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitUi[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeUn[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitVi[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeVn[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

         MSISOrbitWi[iMSIS,0:nLons-1,0:nLats-1,iOrbitAlt] = $
            MSISDataCubeWn[iMSIS,0:nLons-1, $
                  0:nLats-1, iAltIndex[iMSIS,iOrbitAlt]]

    endfor ;iOrbitAlt = 0, nSelAlts - 1 do begin
endfor ;iMSIS = 0, nTimes - 1 do begin

MSISDataCubeNDenTotal = $
    MSISDataCubeN2 + MSISDataCubeO2 + MSISDataCubeO + MSISDataCubeHe
MSISDataCubePressure = $
    MSISDataCubeNDenTotal*kb*MSISDataCubeTemp
;; For the Ion Densities, we will perform a fit versus pressure
;;; where we have Ln[Ni] = F[Ln(Pn)] with a polynomial fit
;;; We can then plot these trajectories for all latitudes and times
ForwardModelTempShapeFactor = fltarr(nTimes,nLons,nLats,nAlts)

for iAlt = 0, nAlts-1 do begin
    ForwardModelTempShapeFactor[0:nTimes-1,0:nLons-1,0:nLats-1,iAlt] = $
       MSISDataCubeTemp[0:nTimes-1,0:nLons-1,0:nLats-1,iAlt]/$  
          MSISOrbitTemp[0:nTimes-1,0:nLons-1,0:nLats-1,0] 
endfor ;iAlt = 0, nAlts-1 do begin


;----------------------------------------------------------------------------------

;;; BEGIN FORWARD MODELING
;;; Step 1:  We extend the data upward using hydrostatic assumption
;;; and the shape factor from MSIS
nTimes = nGITMTimes 
nLons  = nGITMLons 
nLats  = nGITMLats 
nAlts  = nGITMAlts 

AltMax = 565.0
DeltaAlt = 1.0
nUpwardAlts = floor((AltMax-AltRef)/DeltaAlt) + 1
UpwardAlts  = dblarr(nUpwardAlts)
UpwardTemp  = dblarr(nTimes,nLons,nLats,nUpwardAlts)

UpwardN2   = dblarr(nTimes,nLons,nLats,nUpwardAlts)
UpwardO2   = dblarr(nTimes,nLons,nLats,nUpwardAlts)
UpwardO    = dblarr(nTimes,nLons,nLats,nUpwardAlts)
UpwardHe   = dblarr(nTimes,nLons,nLats,nUpwardAlts)

UpwardNe   = dblarr(nTimes,nLons,nLats,nUpwardAlts)

UpwardUn   = dblarr(nTimes,nLons,nLats,nUpwardAlts)
UpwardVn   = dblarr(nTimes,nLons,nLats,nUpwardAlts)
UpwardWn   = dblarr(nTimes,nLons,nLats,nUpwardAlts)

UpwardUi   = dblarr(nTimes,nLons,nLats,nUpwardAlts)
UpwardVi   = dblarr(nTimes,nLons,nLats,nUpwardAlts)
UpwardWi   = dblarr(nTimes,nLons,nLats,nUpwardAlts)

;;; Step 1: Set up the Forward Model Grid Upward
;;;         Calculate the Gravitational Force
;;;         Establish Constants (Mass, kb, etc)

;;; Calculate Gravity (with centrifugal force)
;;; Establish Uniform Grid
UpwardGravity = dblarr(nTimes,nLons,nLats,nUpwardAlts)
for iAlt = 0, nUpwardAlts - 1 do begin
     UpwardAlts[iAlt] = AltRef + DeltaAlt*(iAlt)
endfor
 for iLat = 0, nLats -1 do begin
    for iLon = 0, nLons -1 do begin
      for iTime = 0, nTimes -1 do begin
         for iAlt = 0, nUpwardAlts - 1 do begin
          ; for now ignore centrifugal component
          RE = 6371.0 
          R = RE + UpwardAlts[iAlt] ; km
          Omega = 2.0*!pi/86400.0  ; rads/s
          UpwardGravity[iTime,iLon,iLat,iAlt] = $
             -9.8*(RE/R)^2.0 + $
             R*(Omega*Omega)*cos(GITMLats[iLat]*!pi/180.0)
         endfor ;iTime = 0, nTimes -1 do begin
       endfor ;iLon = 0, nLons -1 do begin
    endfor ;iLat = 0, nLats -1 do begin
endfor 

;;; Step2:  Use MSIS model shape function to map our
;;; "measured" temperatures onto our UpwardGrid
;;; -(a) Interpolate MSIS ShapeFactor onto grid
;;;  (b) Apply shape function to fill in the temperatures
UpwardShapeFunction   = dblarr(nTimes,nLons,nLats,nUpwardAlts)

 for iLat = 0, nLats -1 do begin
    for iLon = 0, nLons -1 do begin
      for iTime = 0, nTimes -1 do begin
          ; result = SPLINE(X,Y,T,/DOUBLE)
          ; X = MSISDataCubeAlts  
          ; Y = ForwardModelTempShapeFactor
          ; T = UpwardAlts
          UpwardShapeFunction[iTime,iLon,iLat,*] = $
             SPLINE(MSISDataCubeAlts[0:nAlts-1],$
          ForwardModelTempShapeFactor[iTime,iLon,iLat,0:nAlts-1], $
                  UpwardAlts[0:nUpwardAlts-1], /DOUBLE)
      endfor ;iTime = 0, nTimes -1 do begin
    endfor ;iLon = 0, nLons -1 do begin
 endfor ;iLat = 0, nLats -1 do begin

 
; Fill the Baseline data with the "measured data"  going 400 upward to 565
UpwardTemp[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitTemp[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardN2[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitN2[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardO2[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitO2[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardO[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitO[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardHe[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitHe[0:nTimes-1,0:nLons-1,0:nLats-1]

UpwardUn[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitUn[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardVn[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitVn[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardWn[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
   GITMOrbitWn[0:nTimes-1,0:nLons-1,0:nLats-1]

UpwardUi[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitUi[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardVi[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
      GITMOrbitVi[0:nTimes-1,0:nLons-1,0:nLats-1]
UpwardWi[0:nTimes-1,0:nLons-1,0:nLats-1,0] = $
   GITMOrbitWi[0:nTimes-1,0:nLons-1,0:nLats-1]


;;; Next, Fill the Temperature Field with the Shape Function (why only temp?)
;;; Calculated above
for iTime = 0, nTimes -1 do begin
  for iLon = 0, nLons -1 do begin
    for iLat = 0, nLats -1 do begin
                  UpwardTemp[iTime,iLon,iLat,1:nUpwardAlts-1] = $
                  UpwardTemp[iTime,iLon,iLat,0]*$
          UpwardShapeFunction[iTime,iLon,iLat,1:nUpwardAlts-1]

                  UpwardUn[iTime,iLon,iLat,1:nUpwardAlts-1] = $
                  UpwardUn[iTime,iLon,iLat,0]
;
                  UpwardVn[iTime,iLon,iLat,1:nUpwardAlts-1] = $
                  UpwardVn[iTime,iLon,iLat,0]
;
                  UpwardWn[iTime,iLon,iLat,1:nUpwardAlts-1] = $
                  UpwardWn[iTime,iLon,iLat,0]
 
                  UpwardUi[iTime,iLon,iLat,1:nUpwardAlts-1] = $
                  UpwardUi[iTime,iLon,iLat,0]
;
                  UpwardVi[iTime,iLon,iLat,1:nUpwardAlts-1] = $
                  UpwardVi[iTime,iLon,iLat,0]
;
                  UpwardWi[iTime,iLon,iLat,1:nUpwardAlts-1] = $
                  UpwardWi[iTime,iLon,iLat,0]

    endfor ;iLat = 0, nLats -1 do begin
  endfor ;iLon = 0, nLons -1 do begin
endfor ;iTime = 0, nTimes -1 do begin
  
; Next, Integrate Upward to Calculate the Mass Densities
for iTime = 0, nTimes -1 do begin
  for iLon = 0, nLons -1 do begin
    for iLat = 0, nLats -1 do begin
       for iAlt = 1, nUpwardAlts-1 do begin
           ; Note: we use the mean values of the scale height
           dAlt = 1000.0*( UpwardAlts[iAlt] - UpwardAlts[iAlt-1])
           MeanTemp = 0.5*(UpwardTemp[iTime,iLon,iLat,iAlt] + $
                           UpwardTemp[iTime,iLon,iLat,iAlt-1])
           MeanGravity = -0.5*(UpwardGravity[iTime,iLon,iLat,iAlt] + $
                               UpwardGravity[iTime,iLon,iLat,iAlt-1])

           ScaleHeight = kb*MeanTemp/MeanGravity/MassN2
           UpwardN2[iTime,iLon,iLat,iAlt] = $
           UpwardN2[iTime,iLon,iLat,iAlt-1]*$
             (UpwardTemp[iTime,iLon,iLat,iAlt-1]/$
              UpwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp( -1.0*dAlt/ScaleHeight)

           ScaleHeight = kb*MeanTemp/MeanGravity/MassO2
           UpwardO2[iTime,iLon,iLat,iAlt] = $
           UpwardO2[iTime,iLon,iLat,iAlt-1]*$
             (UpwardTemp[iTime,iLon,iLat,iAlt-1]/$
              UpwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp( -1.0*dAlt/ScaleHeight)

           ScaleHeight = kb*MeanTemp/MeanGravity/MassO
           UpwardO[iTime,iLon,iLat,iAlt] = $
           UpwardO[iTime,iLon,iLat,iAlt-1]*$
             (UpwardTemp[iTime,iLon,iLat,iAlt-1]/$
              UpwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp( -1.0*dAlt/ScaleHeight)

           ScaleHeight = kb*MeanTemp/MeanGravity/MassHe
           UpwardHe[iTime,iLon,iLat,iAlt] = $
           UpwardHe[iTime,iLon,iLat,iAlt-1]*$
             (UpwardTemp[iTime,iLon,iLat,iAlt-1]/$
              UpwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp( -1.0*dAlt/ScaleHeight)

       endfor 

    endfor ;iLat = 0, nLats -1 do begin
  endfor ;iLon = 0, nLons -1 do begin
endfor ;iTime = 0, nTimes -1 do begin

UpwardNDenTotal = UpwardN2 + UpwardO2 + UpwardO + UpwardHe
UpwardPressure = UpwardNDenTotal*kb*UpwardTemp

LogUpwardPressure = alog(UpwardPressure)

; Get the empirical pressure, log(P), and Log(Ne)
MSISPressure = MSISDataCubeTemp*kb*$
              ( MSISDataCubeN2 + MSISDataCubeO2 + $
                MSISDataCubeO  + MSISDataCubeHe )
LogMSISPressure = alog(MSISPressure         )
LogMSISNe       = alog(MSISDataCubeIDenTotal)

;--------------
; Need to fit the upper part of the ionosphere only
; Step1: Find the peak of the ionosphere
MaxIndex = intarr(nTimes, nLons, nLats)

UpwardLogNe   = dblarr(nTimes,nLons,nLats,nUpwardAlts)
   UpwardNe   = dblarr(nTimes,nLons,nLats,nUpwardAlts)

for iTime = 0, nTimes - 1 do begin
  for iLon = 0, nLons - 1 do begin
    for iLat = 0, nLats - 1 do begin
        MaxValue = MAX(LogMSISNe[iTime,iLon,iLat,0:nAlts-1], iMaxIndex_)
        MaxIndex[iTime,iLon,iLat] = iMaxIndex_
        ; Create our Interpolation Arrays at this location and time
        LocalLogMSISP  = REVERSE(reform(LogMSISPressure[iTime,iLon,iLat,iMaxIndex_:nAlts-1]))
        LocalLogMSISNe = REVERSE(reform(      LogMSISNe[iTime,iLon,iLat,iMaxIndex_:nAlts-1]))
        LocalAlts      = reform(GITMDataCubeAlts[iMaxIndex_:nAlts-1])

        LocalUpwardLogP = REVERSE(reform(LogUpwardPressure[iTime,iLon,iLat,0:nUpwardAlts-1]))
        ; Simple Test for now
        LocalUpwardLogNe = $
            REVERSE(SPLINE(LocalLogMSISP, LocalLogMSISNe, LocalUpwardLogP, /DOUBLE))
        LocalUpwardNe = exp(LocalUpwardLogNe)
        ;for iAlt = 0,nUpwardAlts-1 do begin
        ;    Value = LocalUpwardNe[iAlt]
        ;    LocalUpwardNe[iAlt] = max(1.0, Value)
        ;endfor ;iAlt = 0,nUpwardAlts-1 do begin
       
        ;; Need to scale these values relative to the 400 km value
        ;; Truth Value is GITMOrbitNe, Approximate value is the UpwardNe
        ScaleFactor = GITMOrbitNe[iTime,iLon,iLat,0]/LocalUpwardNe[0]
        LocalUpwardNe = ScaleFactor*LocalUpwardNe

        UpwardNe[iTime,iLon,iLat,0:nUpwardAlts-1] = LocalUpwardNe[0:nUpwardAlts-1]

    endfor ;iLat = 0, nLats - 1 do begin
  endfor ;iLon = 0, nLons - 1 do begin
endfor ;iTime


;;; Next, we fit a polynomial of the Empirical Ni and Pn
;;; use that polynomial to then forward model the Ni using our
;;; Forward model Pressures

;;; Ions Step1:  Fit Ln(Pn) and Ln(Ni)
;;; We use a 4th order polynomial fit
AltMin = 200.0
DeltaAlt = 1.0
nDownwardAlts = floor((AltRef-AltMin)/DeltaAlt) + 2
DownwardAlts  = dblarr(nDownwardAlts)
DownwardTemp  = dblarr(nTimes,nLons,nLats,nDownwardAlts)

DownwardN2   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardO2   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardO    = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardHe   = dblarr(nTimes,nLons,nLats,nDownwardAlts)

DownwardNe  = dblarr(nTimes,nLons,nLats,nDownwardAlts)

DownwardUn   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardVn   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardWn   = dblarr(nTimes,nLons,nLats,nDownwardAlts)

DownwardUi   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardVi   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardWi   = dblarr(nTimes,nLons,nLats,nDownwardAlts)

;;; Step 1: Set up the Forward Model Grid Downward
;;;         Calculate the Gravitational Force
;;;         Establish Constants (Mass, kb, etc)

AMU = 1.66054e-27    ;; kg
kb  = 1.38065e-23    ;; J/K
MassN2 =  28.0*AMU
MassO2 =  32.0*AMU
MassO  =  16.0*AMU
MassHe =   4.0*AMU

;;; Calculate Gravity (with centrifugal force)
;;; Establish Uniform Grid
DownwardGravity = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardAlts[nDownwardAlts-1] = AltRef
for iAlt = nDownwardAlts - 2, 0, -1 do begin
     DownwardAlts[iAlt] = DownwardAlts[iAlt+1] - DeltaAlt
endfor

;print, DownwardAlts
 for iLat = 0, nLats -1 do begin
    for iLon = 0, nLons -1 do begin
      for iTime = 0, nTimes -1 do begin

         for iAlt = 0, nDownwardAlts - 1 do begin
          ; for now ignore centrifugal component
          RE = 6371.0 
          R = RE + DownwardAlts[iAlt] ; km
          Omega = 2.0*!pi/86400.0  ; rads/s
          DownwardGravity[iTime,iLon,iLat,iAlt] = $
             -9.8*(RE/R)^2.0 + $
             R*(Omega*Omega)*cos(GITMLats[iLat]*!pi/180.0)
         endfor ;iTime = 0, nTimes -1 do begin

       endfor ;iLon = 0, nLons -1 do begin
    endfor ;iLat = 0, nLats -1 do begin
endfor 

;;; Step2:  Use MSIS model shape function to map our
;;; "measured" temperatures onto our DownwardGrid
;;; -(a) Interpolate MSIS ShapeFactor onto grid
;;;  (b) Apply shape function to fill in the temperatures
DownwardShapeFunction   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
DownwardNeShapeFunction   = dblarr(nTimes,nLons,nLats,nDownwardAlts)

 for iLat = 0, nLats -1 do begin
    for iLon = 0, nLons -1 do begin
      for iTime = 0, nTimes -1 do begin
          ; result = SPLINE(X,Y,T,/DOUBLE)
          ; X = MSISDataCubeAlts  
          ; Y = ForwardModelTempShapeFactor
          ; T = DownwardAlts
          DownwardShapeFunction[iTime,iLon,iLat,*] = $
             SPLINE(MSISDataCubeAlts[0:nAlts-1],$
          ForwardModelTempShapeFactor[iTime,iLon,iLat,0:nAlts-1], $
                  DownwardAlts[0:nDownwardAlts-1], /DOUBLE)

      endfor ;iTime = 0, nTimes -1 do begin
    endfor ;iLon = 0, nLons -1 do begin
 endfor ;iLat = 0, nLats -1 do begin

;for iAlt = 0, nAlts - 1 do begin
;print, MSISDataCubeAlts[iAlt], ForwardmodelTempShapeFactor[0,1,1,iAlt]
;endfor ;iAlt = 0, nAlts - 1 do begin

;for iAlt = 0, nDownwardAlts - 1 do begin
;   print, DownwardAlts[iAlt], DownwardShapeFunction[0,1,1,iAlt]
;endfor ;iAlt = 0, nDownwardAlts - 1 do begin

; Fill the Baseline data with the "measured data"
iRef = nDownwardAlts-1
DownwardTemp[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitTemp[0:nTimes-1,0:nLons-1,0:nLats-1]

DownwardN2[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitN2[0:nTimes-1,0:nLons-1,0:nLats-1]
DownwardO2[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitO2[0:nTimes-1,0:nLons-1,0:nLats-1]
DownwardO[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitO[0:nTimes-1,0:nLons-1,0:nLats-1]
DownwardHe[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitHe[0:nTimes-1,0:nLons-1,0:nLats-1]

DownwardNe[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitNe[0:nTimes-1,0:nLons-1,0:nLats-1]

DownwardUn[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitUn[0:nTimes-1,0:nLons-1,0:nLats-1]
DownwardVn[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitVn[0:nTimes-1,0:nLons-1,0:nLats-1]
DownwardWn[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitWn[0:nTimes-1,0:nLons-1,0:nLats-1]

DownwardUi[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitUi[0:nTimes-1,0:nLons-1,0:nLats-1]
DownwardVi[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitVi[0:nTimes-1,0:nLons-1,0:nLats-1]
DownwardWi[0:nTimes-1,0:nLons-1,0:nLats-1,iRef] = $
      GITMOrbitWi[0:nTimes-1,0:nLons-1,0:nLats-1]


;;; Next, Fill the Temperature Field with the Shape Function
;;; Calculated above
for iTime = 0, nTimes -1 do begin
  for iLon = 0, nLons -1 do begin
    for iLat = 0, nLats -1 do begin
                  DownwardTemp[iTime,iLon,iLat,1:nDownwardAlts-1] = $
                  DownwardTemp[iTime,iLon,iLat,iRef]*$
          DownwardShapeFunction[iTime,iLon,iLat,1:nDownwardAlts-1]
 
                  DownwardUn[iTime,iLon,iLat,1:nDownwardAlts-1] = $
                  DownwardUn[iTime,iLon,iLat,iRef]

                  DownwardVn[iTime,iLon,iLat,1:nDownwardAlts-1] = $
                  DownwardVn[iTime,iLon,iLat,iRef]

                  DownwardWn[iTime,iLon,iLat,1:nDownwardAlts-1] = $
                  DownwardWn[iTime,iLon,iLat,iRef]

                  DownwardUi[iTime,iLon,iLat,1:nDownwardAlts-1] = $
                  DownwardUi[iTime,iLon,iLat,iRef]

                  DownwardVi[iTime,iLon,iLat,1:nDownwardAlts-1] = $
                  DownwardVi[iTime,iLon,iLat,iRef]

                  DownwardWi[iTime,iLon,iLat,1:nDownwardAlts-1] = $
                  DownwardWi[iTime,iLon,iLat,iRef]
    endfor ;iLat = 0, nLats -1 do begin
  endfor ;iLon = 0, nLons -1 do begin
endfor ;iTime = 0, nTimes -1 do begin
  
; Next, Integrate Downward to Calculate the Mass Densities
for iTime = 0, nTimes -1 do begin
  for iLon = 0, nLons -1 do begin
    for iLat = 0, nLats -1 do begin
       ;for iAlt = 1, nDownwardAlts-1 do begin
       ; Here, we integrate downward
       for iAlt = iRef-1, 1, -1 do begin
           ; Note: we use the mean values of the scale height
           dAlt = 1000.0*( DownwardAlts[iAlt+1] - DownwardAlts[iAlt  ])
           MeanTemp = 0.5*(DownwardTemp[iTime,iLon,iLat,iAlt+1] + $
                           DownwardTemp[iTime,iLon,iLat,iAlt  ])
           MeanGravity = -0.5*(DownwardGravity[iTime,iLon,iLat,iAlt+1] + $
                               DownwardGravity[iTime,iLon,iLat,iAlt  ])

           ScaleHeight = kb*MeanTemp/MeanGravity/MassN2
           DownwardN2[iTime,iLon,iLat,iAlt] = $
           DownwardN2[iTime,iLon,iLat,iAlt+1]*$
             (DownwardTemp[iTime,iLon,iLat,iAlt+1]/$
              DownwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp(  1.0*dAlt/ScaleHeight)

           ScaleHeight = kb*MeanTemp/MeanGravity/MassO2
           DownwardO2[iTime,iLon,iLat,iAlt] = $
           DownwardO2[iTime,iLon,iLat,iAlt+1]*$
             (DownwardTemp[iTime,iLon,iLat,iAlt+1]/$
              DownwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp(  1.0*dAlt/ScaleHeight)

           ScaleHeight = kb*MeanTemp/MeanGravity/MassO
           DownwardO[iTime,iLon,iLat,iAlt] = $
           DownwardO[iTime,iLon,iLat,iAlt+1]*$
             (DownwardTemp[iTime,iLon,iLat,iAlt+1]/$
              DownwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp(  1.0*dAlt/ScaleHeight)

           ScaleHeight = kb*MeanTemp/MeanGravity/MassHe
           DownwardHe[iTime,iLon,iLat,iAlt] = $
           DownwardHe[iTime,iLon,iLat,iAlt+1]*$
             (DownwardTemp[iTime,iLon,iLat,iAlt+1]/$
              DownwardTemp[iTime,iLon,iLat,iAlt  ])*$
           exp(  1.0*dAlt/ScaleHeight)

       endfor 

    endfor ;iLat = 0, nLats -1 do begin
  endfor ;iLon = 0, nLons -1 do begin
endfor ;iTime = 0, nTimes -1 do begin

DownwardNDenTotal = DownwardN2 + DownwardO2 + DownwardO + DownwardHe
GITMDataCubeNDenTotal = $
    GITMDataCubeN2 + GITMDataCubeO2 + GITMDataCubeO + GITMDataCubeHe
DownwardPressure = DownwardNDenTotal*kb*DownwardTemp

LogDownwardPressure = alog(DownwardPressure)
MaxIndex = intarr(nTimes, nLons, nLats)

DownwardLogNe   = dblarr(nTimes,nLons,nLats,nDownwardAlts)
   DownwardNe   = dblarr(nTimes,nLons,nLats,nDownwardAlts)

for iTime = 0, nTimes - 1 do begin
  for iLon = 0, nLons - 1 do begin
    for iLat = 0, nLats - 1 do begin
        MaxValue = MAX(LogMSISNe[iTime,iLon,iLat,0:nAlts-1], iMaxIndex_)
        MaxIndex[iTime,iLon,iLat] = iMaxIndex_
        ; Create our Interpolation Arrays at this location and time
        LocalLogMSISP  = REVERSE(reform(LogMSISPressure[iTime,iLon,iLat,iMaxIndex_:nAlts-1]))
        LocalLogMSISNe = REVERSE(reform(      LogMSISNe[iTime,iLon,iLat,iMaxIndex_:nAlts-1]))
        LocalAlts      = reform(GITMDataCubeAlts[iMaxIndex_:nAlts-1])

        LocalDownwardLogP = REVERSE(reform(LogDownwardPressure[iTime,iLon,iLat,0:nDownwardAlts-1]))
        ; Simple Test for now
        LocalDownwardLogNe = $
            REVERSE(SPLINE(LocalLogMSISP, LocalLogMSISNe, LocalDownwardLogP, /DOUBLE))
        LocalDownwardNe = exp(LocalDownwardLogNe)
        ;for iAlt = 0,nDownwardAlts-1 do begin
        ;    Value = LocalDownwardNe[iAlt]
        ;    LocalDownwardNe[iAlt] = max(1.0, Value)
        ;endfor ;iAlt = 0,nDownwardAlts-1 do begin
        ;; Need to scale these values relative to the 400 km value
        ;; Truth Value is GITMOrbitNe, Approximate value is the DownwardNe
        ScaleFactor = GITMOrbitNe[iTime,iLon,iLat]/$
                           LocalDownwardNe[nDownwardAlts-1]
        LocalDownwardNe = ScaleFactor*LocalDownwardNe

        DownwardNe[iTime,iLon,iLat,0:nDownwardAlts-1] = LocalDownwardNe[0:nDownwardAlts-1]

    endfor ;iLat = 0, nLats - 1 do begin
  endfor ;iLon = 0, nLons - 1 do begin
endfor ;iTime

nMergedAlts = nUpwardAlts + nDownwardAlts - 1 

MergedAlts = fltarr(nMergedAlts)

MergedTotalDensity     = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedTemperature      = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedN2Density        = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedO2Density        = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedODensity         = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedHeDensity        = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedTotalIDensity    = fltarr(nTimes, nLons, nLats, nMergedAlts)

MergedUn     = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedVn     = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedWn     = fltarr(nTimes, nLons, nLats, nMergedAlts)

MergedUi     = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedVi     = fltarr(nTimes, nLons, nLats, nMergedAlts)
MergedWi     = fltarr(nTimes, nLons, nLats, nMergedAlts)

MergedAlts[0:nDownwardAlts-2] = DownwardAlts[0:nDownwardAlts-2]
MergedAlts[nDownwardAlts-1:nMergedAlts-1] = UpwardAlts[0:nUpwardAlts-1]


MergedTotalDensity[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
 DownwardNDenTotal[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedTotalDensity[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
   UpwardNDenTotal[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedTemperature[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
     DownwardTemp[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedTemperature[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
       UpwardTemp[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedN2Density[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
     DownwardN2[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedN2Density[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
       UpwardN2[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedO2Density[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
     DownwardO2[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedO2Density[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
       UpwardO2[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedODensity[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
     DownwardO[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedODensity[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
       UpwardO[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedHeDensity[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
     DownwardHe[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedHeDensity[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
       UpwardHe[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedTotalIDensity[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
         DownwardNe[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedTotalIDensity[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
           UpwardNe[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedUn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
         DownwardUn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedUn[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
           UpwardUn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedVn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
         DownwardVn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedVn[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
           UpwardVn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedWn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
         DownwardWn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedWn[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
           UpwardWn[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedUi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
         DownwardUi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedUi[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
           UpwardUi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedVi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
         DownwardVi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedVi[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
           UpwardVi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]

MergedWi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
         DownwardWi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
MergedWi[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
           UpwardWi[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]



p1 = PLOT(GITMDataCubeNDenTotal[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
           /XLOG, xrange = [1.0e+11, 1.0e+16], xstyle = 1, $
           xtitle = 'Density (m$^{-3}$)', ytitle = 'Altitude (km)', $
           title = 'Total Neutral Density', FONT_SIZE = 12, $ 
           thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardNDenTotal[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardNDenTotal[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedTotalDensity[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalNeutralDensity.png'


p1 = PLOT(GITMDataCubeTemp[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
           xrange = [600.0, 1000.0], xstyle = 1, $
           xtitle = 'Temperature (K)', ytitle = 'Altitude (km)', $
           title = ' Neutral Temperature', FONT_SIZE = 12, $ 
           thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardTemp[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardTemp[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedTemperature[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.85,0.80])
p1.SAVE, 'Temperature.png'


p1 = PLOT(GITMDataCubeIDenTotal[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
           /XLOG, xrange = [1.0e+09, 1.0e+12], xstyle = 1, $
           xtitle = 'Density (m$^{-3}$)', ytitle = 'Altitude (km)', $
           title = 'Total Ion Density', FONT_SIZE = 12, $ 
           thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardNe[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardNe[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Downward' )

p4 = PLOT(MergedTotalIDensity[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )
leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalIonDensity.png'


p1 = PLOT(GITMDataCubeUn[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
          xrange = [-100.0, 100.0], xstyle = 1, $
          xtitle = 'EastWinds (m/s)', ytitle = 'Altitude (km)', $
          title = ' Neutral Winds ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardUn[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardUn[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedUn[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalNeutralEastWinds.png'


p1 = PLOT(GITMDataCubeVn[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
          xrange = [-200.0, 100.0], xstyle = 1, $
          xtitle = 'EastWinds (m/s)', ytitle = 'Altitude (km)', $
          title = ' Neutral Winds ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardVn[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardVn[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedVn[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalNeutral_NorthWinds.png'



p1 = PLOT(GITMDataCubeWn[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
          xrange = [-100.0, 100.0], xstyle = 1, $
          xtitle = 'Upward Winds (m/s)', ytitle = 'Altitude (km)', $
          title = ' Neutral Winds ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardWn[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardWn[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedWn[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalNeutralUpwardWinds.png'


p1 = PLOT(GITMDataCubeUi[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
          xrange = [-100.0, 100.0], xstyle = 1, $
          xtitle = 'EastWinds (m/s)', ytitle = 'Altitude (km)', $
          title = ' Ion Winds ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardUi[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardUi[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedUi[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalIonEastWinds.png'


p1 = PLOT(GITMDataCubeVi[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
          xrange = [-100.0, 100.0], xstyle = 1, $
          xtitle = 'EastWinds (m/s)', ytitle = 'Altitude (km)', $
          title = ' Ion Winds ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardVi[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardVi[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedVi[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalIon_NorthWinds.png'



p1 = PLOT(GITMDataCubeWi[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [250.0, 520.0],ystyle = 1, $
          xrange = [-100.0, 100.0], xstyle = 1, $
          xtitle = 'Upward Winds (m/s)', ytitle = 'Altitude (km)', $
          title = ' Ion Winds ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p3 = PLOT(UpwardWi[0,1,1,*], UpwardAlts[*], /OVERPLOT, $
          color = 'red', thick = 2.0, linestyle = 2, name = 'Forward Model Up' )

p2 = PLOT(DownwardWi[0,1,1,*], DownwardAlts[*], /OVERPLOT, $
          color = 'blue', thick = 2.0, linestyle = 2, name = 'Forward Model Down' )

p4 = PLOT(MergedWi[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p3,p2,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'TotalIon_UpwardWinds.png'
;;;; CALCULATE ION-NEUTRAL DRAG

   ;--  Forcing Constants and Derived Quanities
   nSpecies = 6
   inO_  = 1-1
   inO2_ = 2-1
   inN2_ = 3-1
   inN_  = 4-1
   inNO_ = 5-1
   inHe_ = 6-1

   nIons = 10
   inOP_  = 1-1
   inO2P_ = 2-1
   inN2P_ = 3-1
   inNP_  = 4-1
   inNOP_ = 5-1
   inO_2DP_ = 6-1 
   inO_2PP_ = 7 -1
   inHP_ = 8 -1
   inHeP_ = 9-1
   ine_ = 10-1

   nDims = 3
   inEast_  = 0
   inNorth_ = 1
   inUp_    = 2

   AMU = 1.66054e-27     ; kg
   kb  = 1.38064852e-23  ; J/K
   Mass = fltarr(nSpecies)
   mu   = fltarr(nSpecies)

   ; Mass is the species mass in kg
   ; mu is the species mass in amu
   Mass(inO_)  = 16.0 * AMU
     mu(inO_)  = 16.0 
   Mass(inO2_) =2.0* 16.0 * AMU
     mu(inO2_) =2.0* 16.0 
   Mass(inN2_) =2.0* 14.0 * AMU
     mu(inN2_) =2.0* 14.0 
   Mass(inN_)  =1.0* 14.0 * AMU
     mu(inN_)  =1.0* 14.0 
   Mass(inNO_)  =Mass(inO_) + Mass(inN_)
     mu(inNO_)  =mu(inO_) + mu(inN_)
   Mass(inHe_)  =4.0* AMU
     mu(inHe_)  =4.0  

   ; Neutral Polarizability:  10^-24 cm^3
   Polarize = fltarr(nSpecies) 
   Polarize(inO_  ) = 0.802
   Polarize(inN_  ) = 0.802*(7.5/5.4)  ; Rees
   Polarize(inN2_ ) = 1.74
   Polarize(inNO_ ) = 0.5*(Polarize(inN_) + Polarize(inO_))  ; Approx
   Polarize(inO2_ ) = 1.74*(10.7/11.8) ;; Rees
   Polarize(inHe_ ) = 0.802

   MassI = fltarr(nIons)
     muI = fltarr(nIons)
   MassI(inOP_ ) = Mass(inO_ )
     muI(inOP_ ) = mu(inO_ )
   MassI(inO2P_) = Mass(inO2_)
     muI(inO2P_) = mu(inO2_)
   MassI(inN2P_) = Mass(inN2_)
     muI(inN2P_) = mu(inN2_)
   MassI(inNP_ ) = Mass(inN_ )
     muI(inNP_ ) = mu(inN_ )
   MassI(inNOP_) = Mass(inNO_)
     muI(inNOP_) = mu(inNO_)
   MassI(inO_2DP_ ) = Mass(inO_ )
     muI(inO_2DP_ ) = mu(inO_ )
   MassI(inO_2PP_ ) = Mass(inO_ )
     muI(inO_2PP_ ) = mu(inO_ )
   MassI(inHP_ ) = AMU
     muI(inHP_ ) = 1.0
   MassI(inHeP_ ) = Mass(inHe_)
     muI(inHeP_ ) = mu(inHe_)
   MassI(ine_   ) = 9.109e-31
     muI(ine_   ) = MassI(ine_)/AMU

   MergedNu_nTotal   = fltarr(nTimes, nLons, nLats, nMergedAlts) ; MeanMass Neutrals
   MergedNu_iTotal   = fltarr(nTimes, nLons, nLats, nMergedAlts) ;; MeanMass Ions
   MergedIonDragUn  = fltarr(nTimes, nLons, nLats, nMergedAlts ) ; 
   MergedIonDragVn  = fltarr(nTimes, nLons, nLats, nMergedAlts ) ;
   MergedIonDragWn  = fltarr(nTimes, nLons, nLats, nMergedAlts ) ; 

;MergedTotalIDensity[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2] = $
;         DownwardNe[0:nTimes-1,0:nLons-1,0:nLats-1,0:nDownwardAlts-2]
;MergedTotalIDensity[0:nTimes-1,0:nLons-1,0:nLats-1,nDownwardAlts-1:nMergedAlts-1] = $
;           UpwardNe[0:nTimes-1,0:nLons-1,0:nLats-1,0:nUpwardAlts-1]
;

   ;; Assume mostly O+ and O
   ;Polarize(inO_  ) = 0.802
   MergedNu_nTotal = MergedTotalIDensity*(2.0e-15)*sqrt(0.5*Polarize(inO_))

   MergedIonDragUn = 4.0*MergedNu_nTotal*(MergedUi - MergedUn)
   MergedIonDragVn = 4.0*MergedNu_nTotal*(MergedVi - MergedVn)


p1 = PLOT(GITMDataCubeIonDragUn[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [200.0, 520.0],ystyle = 1, $
          xrange = [-0.05, 0.15], xstyle = 1, $
          xtitle = 'IonDrag Accel (m/s$^{2}$)', ytitle = 'Altitude (km)', $
          title = 'Ion-Neutral Drag East ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p4 = PLOT(MergedIonDragUn[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'ForwardModel_EastIonDragAccel.png'


p1 = PLOT(GITMDataCubeIonDragVn[0,1,1,*], GITMDataCubeAlts, $
          color = 'black', yrange = [150.0, 520.0],ystyle = 1, $
          xrange = [-0.05, 0.15], xstyle = 1, $
          xtitle = 'IonDrag Accel (m/s$^{2}$)', ytitle = 'Altitude (km)', $
          title = 'Ion-Neutral Drag North ', FONT_SIZE = 12, $ 
          thick = 2.0 , linestyle = 0, name = 'Truth ')

p4 = PLOT(MergedIonDragVn[0,1,1,*], MergedAlts[*], /OVERPLOT, $
          color = 'grey', thick = 2.0, linestyle = 1, name = 'Merged Model' )

leg = legend(target = [p1,p4], /auto_text_color, pos = [0.75,0.80])
p1.SAVE, 'ForwardModel_NorthIonDragAccel.png'




end
