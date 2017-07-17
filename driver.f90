PROGRAM driver

  USE model_global
  USE model_Parameters  !ONLY: IND_*
  USE model_Rates,       ONLY: Update_SUN, Update_RCONST
  USE model_integrator,  ONLY: integrate!, IERR_NAMES
  USE model_monitor,     ONLY: spc_names
  USE model_Init
  USE model_Util
  USE constants

  IMPLICIT NONE
  REAL(dp) :: ENDSTATE(NVAR), total, RATIO, TNOX, TNOX_OLD
  REAL(dp) :: STARTSTATE(NVAR), TIMESCALE
  REAL(dp) :: RH
  REAL(dp) :: RSTATE(20)
  REAL(dp) :: DIURNAL_OLD(NVAR,3000), DIURNAL_NEW(NVAR,3000)
  REAL(dp) :: DIURNAL_RATES(NREACT, 3000)
  INTEGER  :: ERROR, IJ
  LOGICAL :: SCREWED

! Photolysis calculation variables
  REAL(dp) :: Alta

  Integer  :: Daycounter, CONSTNOXSPEC, JK,counter

  REAL(dp) :: NOXRATIO, NEW_TIME
  REAL(dp) :: Fracdiff, SpeedRatio, oldfracdiff, FRACCOUNT
 
! find the correct half hour of the model integration
  INTEGER :: HALFHOUR
  REAL    :: OBS_TEMP(48)=(/&
       -4.90988,-5.998,-6.76708,-8.16098,-8.60077,-8.97174,-9.09161,-9.34328,-9.2665,-9.42021,&
       -9.32238,-9.36907,-9.49675,-9.27415,-9.07073,-9.18774,-10.3613,-10.3786,-10.3822,-10.5816,&
       -10.5055,-10.5738,-10.9845,-11.2199,-11.5564,-11.5178,-11.4663,-11.4704,-11.4965,-11.5318,&
       -11.6473,-11.7589,-11.2993,-10.784,-10.0277,-9.14637,-7.99273,-7.53081,-6.87639,-6.05381,&
       -5.63831,-5.54391,-4.88845,-4.97976,-4.86069,-4.50494,-4.76706,-4.63537/)

  REAL    :: EMISS(48)=(/&
       0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,&
       0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,&
       0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,&
       0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,&
       0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9/)

! 2012 HONO
  REAL     :: HONO(48)=(/&
       28.88,27.85,24.45,24.59,29.24,28.68,30.13,31.81,32.57,&
       41.21,43.53,44.03,47.46,47.28,49.,54.26,59.4,58.79,58.81,&
       59.1, 61.15,60.21,62.51,63.9, 62.79, 62.91, 67.86,68.44, 69.7,&
        74.19, 76.06 ,92.93, 92.57, 86.27, 68.17, 45.34, 37.58, 30.85, 29.,&
        26.8, 30.18, 38.34, 33.73, 33.6, 29.07, 26.01, 29.28,32.42/)

  REAL     :: CLNO2(48)=(/&
       21.0994,26.3102,45.4727,78.7776,153.944,191.638,233.375,272.866,275.795,267.478,&
       272.051,282.215,259.541,230.965,235.998,259.276,301.291,297.618,300.514,340.369,&
       301.005,279.349,261.297,279.124,296.919,291.245,287.455,276.348,279.633,292.379,&
       244.924,192.209,160.123,119.415,93.7459,59.6054,39.4277,26.4976,19.0286,15.7525,&
       13.7339,13.5311,14.6732,15.9641,15.8982,15.7743,15.1901,16.8845/)

! Have called C10 aromatics 5-ethyl-m-xylene (DIME35EB) as this is the only C10 in MCMv3.2
  REAL      ::  DIME35EB(48)=(/&
       69.452,70.3056,67.4846,74.1145,74.8584,73.826,82.6164,88.9278,76.7496,83.3863,&
       75.3736,67.7071,67.4484,64.812,64.2898,67.9836,70.5326,65.4262,69.1544,76.5433,&
       85.9401,75.6455,79.2941,75.1509,82.3132,80.2801,75.915,87.1978,81.9256,81.6976,&
       104.801,104.91,100.605,108.964,101.727,98.7474,102.016,103.402,94.3561,84.8273,&
       74.6906,80.2924,81.3013,77.1724,76.3041,74.0445,64.804,66.351/)

! Have grouped C11 and C12 aromatics as 3,5-diethyl toluene (DIET35TOL) as this is the only C10+ in MCMv3.2
  REAL      ::  DIET35TOL(48)=(/&
       36.8734,37.0257,34.21,35.1549,37.7419,34.8994,37.9065,40.8461,36.5408,37.1806,&
       34.3297,32.7743,32.0106,31.163,29.9589,32.8109,32.9365,32.3456,34.3427,36.1953,&
       38.1519,35.3669,34.9246,34.7808,35.198,35.4073,34.135,38.581,36.941,35.414,&
       44.537,44.8959,45.9063,47.5709,46.9998,46.6037,48.8989,48.9082,46.6074,44.5145,&
       41.5593,44.9039,46.0723,43.1969,40.4023,39.1257,34.5272,35.63/)

!       
! end mje 

  STEPMIN = 0.0_dp
  STEPMAX = 0.0_dp
  RTOL(:) = 1.0e-7_dp
  ATOL(:) = 1.0_dp
  counter=0
  LAST_POINT=.FALSE.
 

!  IF YOU WANT TO CONSTRAIN THE NOX THEN
   CONSTRAIN_NOX=.FALSE.

!  If we are running a constrained run we want one file with the final points calculated
   IF (CONSTRAIN_RUN .EQ. .TRUE. .AND. OUTPUT_LAST .EQ. .FALSE.) THEN 
      call newinitsavedata(1)
   ENDIF


CALL INITVAL(0)
!This is the loop of different points in the Init_cons.dat file

!$OMP PARALLEL
!$OMP DO

do counter=1,LINECOUNT-3
!100 counter=counter+1

! If we are running a non-constrained run then we want one file per input in the Init_cons.dat file
   IF (CONSTRAIN_RUN .EQ. .FALSE. .OR. OUTPUT_LAST .EQ. .TRUE.) THEN 
      CALL NEWINITSAVEDATA(COUNTER)
   ENDIF
   
! Read in the next initial conditions
     WRITE (6,*) 'Reading in point', counter
!$OMP CRITICAL
     CALL InitVal(counter)
!$OMP END CRITICAL 
! Set up the output files file

     M   = CFACTOR
     O2 = 0.21 * CFACTOR
     N2 = 0.78 * CFACTOR

     write (6,*) 'Starting Jday:',jday

! tstart is the starting time, variations due to day of year are dealt with somewhere else 
     tstart = (mod(jday,1.))*24.*60.*60.              ! time start

! convert tstart to local time
     tstart = tstart+LON/360.*24*60*60.

! tend is the end time. IntTime is determined from the Init_cons.dat file
     tend = tstart + IntTime    

!dt is the output timestep and the timestep between times rate constants and notably photolysis rates are calcualted
     dt = 600.                       

     write (6,*) 'Starting time:',tstart
     write (6,*) 'Ending time:', tend
     write (6,*) 'Time step:', dt
     
!    Set up the photolysis rates
!    First calculate pressure altitude from altitude
     WRITE (6,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
     WRITE (6,*) 'Using TUV to calculate photolysis rates as a function of SZA'
        alta=(1-(press/1013.25)**0.190263)*288.15/0.00198122*0.304800/1000.
        write (6,*) 'Aerosol surface area', SAREA
        write (6,*) 'Aerosol particle radius 1', RP1
        write (6,*) 'Altitude =', alta
        write (6,*) 'Pressure =', Press
        write (6,*) 'Temperature =', Temp
        write (6,*) 'Latitude =', Lat
        write (6,*) 'Lon =', Lon
        write (6,*) 'Local Time =', Tstart/(60.*60.)
        write (6,*) 'SZA =',ZENANG(int(jday),Tstart/(60.*60.),lat)*180./(4*ATAN(1.))
        if (o3col .eq. 0) then 
           o3col=323.
           write (6,*) 'Ozone column not specified using 323 Dobsons'
        else
           write (6,*) 'Ozone column =', o3col
        endif
        
        if (albedo .eq. 0) then 
	   albedo=0.1
	   write (6,*) 'Albedo not specified using 0.1'
	else
	   write (6,*) 'Albedo =', albedo
	endif       
!    Calculate the photolysis rates for the run
!$OMP CRITICAL 
	IF (JREPEAT .EQ. 0 .OR. COUNTER .EQ. 1) THEN 
        call set_up_photol(O3col,Albedo, alta, temp, bs,cs,ds,szas,svj_tj)
	ELSE
	WRITE (6,*) 'Using previously calculated photolysis params'
	ENDIF
!$OMP END CRITICAL
     WRITE (6,*) 'Photolysis rates calculated'
     WRITE (6,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
     time = tstart

     OLDFRACDIFF=0.
    
! If NOx is being constrained calculate the total NOx in the model 
     IF (CONSTRAIN_NOX) THEN 
      TNOX_OLD=0.
      DO JK=1,NVAR
         TNOX_OLD=TNOX_OLD+C(JK)*NOX(JK)
      ENDDO
     ENDIF

! Define the initial state of the model 
     DO I=1,NVAR
        STARTSTATE(I)=C(I)
        DO IJ=1,3000
           DIURNAL_OLD(I,IJ)=0.
        ENDDO
     ENDDO

! Calculate clear sky photolysis rates
     JFACTNO2=1.
     JFACTO1D=1.

! Update the rate constants
     CALL Update_RCONST()

     WRITE (6,*) 'JO1D Calc=', J(1)
     WRITE (6,*) 'JO1D Measre =', JO1D
! Calcualte correction factors for the model photolysis rates
     IF (JO1D .NE. 0. .AND. J(1) .GT. 0.) THEN
        JFACTO1D=JO1D/J(1)
     ENDIF

     IF (JNO2 .NE. 0. .AND. J(4) .GT. 0.) THEN
        JFACTNO2=JNO2/J(4)
     ENDIF

     IF (JNO2 .EQ. 0. .AND. JO1D .NE. 0.) THEN 
        JFACTNO2=JFACTO1D
     ENDIF

     IF (JO1D .EQ. 0. .AND. JNO2 .NE. 0.) THEN 
        JFACTO1D=JFACTNO2
     ENDIF
! Fix JFACTO1D and JFACTNO2 to observed
      JFACTO1D=0.5167336
      JFACTNO2=0.6643193

     WRITE (6,*) 'Correction JO1D and JNO2 by', JFACTO1D,JFACTNO2

! If we are running a free running model output the initial condition so T=0 of the output file gives the initial condition
     IF (CONSTRAIN_RUN .EQ. .FALSE.) THEN 
           CALL NEWSAVEDATA()
     ENDIF

! Set up a counter to count the number time that the model has been run for 
     Daycounter=0

! This is the main loop for integrations
     time_loop: DO WHILE (time < TEND)

!mje decide on the deposition rate
        WRITE (6,*) TIME, TIME/(60*60), MOD(TIME/(60*60),24.)

         DEPOSITION=(1e-6+1.8e-7)*1.5
        IF (MOD(TIME/(60*60),24.) .GT. 10. .AND. &
             MOD(TIME/(60*60),24.) .LT. 19.) THEN 
           DEPOSITION=(1e-6+1.8e-7)*10.*1.5
        ENDIF

        O3DEPOSITION=(1.e-6+1.8e-7-3.e-7)*1.5
        IF (MOD(TIME/(60*60),24.) .GT. 10. .AND. &
             MOD(TIME/(60*60),24.) .LT. 19.) THEN 
           O3DEPOSITION=(1.e-6+1.8e-7-3.e-7)*10.*1.5
        ENDIF

        HNO3DEPOSITION=(1e-6+1.8e-7)*1.5
        IF (MOD(TIME/(60*60),24.) .GT. 10. .AND. &
             MOD(TIME/(60*60),24.) .LT. 19.5) THEN 
                   HNO3DEPOSITION=(1e-6+1.8e-7)*10.*1.5
        ENDIF

        PNADEPOS = DEPOSITION    
        N2O5DEPOSITION=0. 

        N2O5HYD=(1.2e-3)
        IF (MOD(TIME/(60*60),24.) .GT. 1. .AND. &
             MOD(TIME/(60*60),24.) .LT. 15.5) THEN 
           N2O5HYD=(9.e-3)
        ENDIF

        KNOSCAV=1.0e-15
        KNO2SCAV=1.0e-16
        NO2SCAVDEP=0.
        NOSCAVDEP=0.

! Update the rate constants
        CALL Update_RCONST()
                
! Integrate the model forwards 1 timestep
        CALL INTEGRATE( TIN = time, TOUT = time+DT, RSTATUS_U = RSTATE, &
             ICNTRL_U = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /),&
             IERR_U=ERROR)

        IF (ERROR .NE. 1) THEN 
           WRITE (6,*) 'Integration error.'
           WRITE (6,*) 'Skipping point'
           DO I=1,NVAR
              C(I)=0.
           ENDDO
           GOTO 1000
        ENDIF
! Traps for NaN
        SCREWED=.FALSE.
        DO I=1,NVAR
           IF (ISNAN(C(I))) SCREWED=.TRUE.
        ENDDO

        IF (SCREWED) THEN 
           SCREWED=.FALSE.
           DO i=1,NVAR
              C(I)=0.
           ENDDO
           GOTO 1000
        ENDIF
! Update the time to reflect the integration has taken place and 
        time = RSTATE(1)
 Daycounter=Daycounter+1

! mje for Pete update
! tstart is local time, but Pete's numbers are GMT
    HALFHOUR=INT(MOD((TIME-LON/360.*24*60*60)*2./(60.*60.),48.))+1
    TEMP=OBS_TEMP(HALFHOUR)+273.15
    C(ind_EMISS)=EMISS(HALFHOUR)*1e-8*CFACTOR*1.2
    C(ind_HONO)=HONO(HALFHOUR)*1e-12*CFACTOR
    C(ind_CLNO2)=CLNO2(HALFHOUR)*1e-12*CFACTOR*0.6
    C(ind_DIME35EB)=DIME35EB(HALFHOUR)*1e-12*CFACTOR
    C(ind_DIET35TOL)=DIET35TOL(HALFHOUR)*1e-12*CFACTOR

!    WRITE (6,*) HALFHOUR, TIME,LON/360*24.*60.*60.
!

! If we are constraining NOx then:
        IF (CONSTRAIN_NOX) THEN 
           WRITE (6,*) 'Constraining NOx'

! Calcualte the total NOx in the box
           TNOX=0
           DO I=1,NVAR
              IF (NOX(I) .NE. 0) THEN 
                 TNOX=TNOX+C(I)*NOX(I)
              ENDIF
           ENDDO
           
! Update all NOx variables so that the total NOx in the box is the same as it was
           DO I=1,NVAR
              IF (NOX(I) .NE. 0) THEN 
                 C(I)=C(I)*TNOX_OLD/TNOX
              ENDIF
           ENDDO
        ENDIF

! If constrain species concentALrations if necessary
        DO I=1,NVAR
           IF (CONSTRAIN(I) .GT. 0) THEN             
              C(I)=CONSTRAIN(I)
           ENDIF
        ENDDO
        
! If we are not doing a constrained run then output the concentrations
        IF (CONSTRAIN_RUN .EQ. .FALSE.) THEN 
           CALL NEWSAVEDATA()
        ENDIF

! If we are doing a constrained run we need to store the diurnal profile of all the species
        IF (CONSTRAIN_RUN .EQ. .TRUE.) THEN
           DO I=1,NVAR
              DIURNAL_NEW(I,DAYCOUNTER)=C(I)
           ENDDO

           DO I=1,NREACT
              DIURNAL_RATES(I,DAYCOUNTER)=RCONST(I)
           ENDDO
           
! Are we at the end of a day?
! If so we need to 
!   1) fiddle with the NOX to ensure it has the right concentrations see if we have reached a steady state
           IF (DAYCOUNTER*DT .GE. 24.*60.*60.) THEN

! Sort out the NOx. Need to increase the NOx concentration so that the constrained species is right
! What is  the constrained NOx species? Put result into CONSTNOXSPEC
           DO I=1,NVAR
              IF (NOX(I) .NE. 0) THEN 
                 IF (CONSTRAIN(I) .LT. 0) THEN 
                    CONSTNOXSPEC=I
                 ENDIF
              ENDIF
           ENDDO

! Calculate the ratio between the value we the constrained NOx species and what we have
! Remember the constrained NOx species is given by the negative constrained value
    NOXRATIO=-CONSTRAIN(CONSTNOXSPEC)/C(CONSTNOXSPEC)
 
    ! Multiply all the NOx species by the ratio so 
    DO I=1,NVAR
       IF (NOX(I) .NE. 0) THEN 
          C(I)=C(I)*NOXRATIO
       ENDIF
    ENDDO
    
! Update the total amount of NOx in box
    TNOX_OLD=TNOX_OLD*NOXRATIO
 


! Lets see how much the diurnal ratios have changed since the last itteration

! Frac diff is our metric for how much it has changed 
              FRACDIFF=0.
              FRACCOUNT=0.
! Add up for all species and for each time point in the day
              DO I=1,NVAR
                 DO JK=1,DAYCOUNTER

!If there is a concentration calculated
                    IF (DIURNAL_NEW(I,JK) .GT. 1.e2 .AND. &
                         TRIM(SPC_NAMES(I)) .NE. 'DUMMY') THEN 

!Calculate the absolute value of the fractional difference and add it on
! Increment the counter to calculate the average
                       FRACDIFF=FRACDIFF+&
                            ABS(DIURNAL_OLD(I,JK)-DIURNAL_NEW(I,JK))/&
                            DIURNAL_NEW(I,JK)
                       FRACCOUNT=FRACCOUNT+1
                    ENDIF
                     
                 ENDDO
              ENDDO

! Calculate the average fractional difference
              FRACDIFF=FRACDIFF/FRACCOUNT

! Output the diagnostic
              WRITE (6,*) 'Fraction difference in the diurnal profile:', FRACDIFF


! Store the new diurnal profile as the old one so we can compare with the next day

              DO I=1,NVAR
                 DO JK=1,DAYCOUNTER
                    DIURNAL_OLD(I,JK)=DIURNAL_NEW(I,JK)
                 ENDDO
              ENDDO

! reset the day counter to 0
             

! if the system has converged then end the simulation for this point

              IF (FRACDIFF .LE. 1e-3) THEN
                 GOTO 1000
              ENDIF
              DAYCOUNTER=0
              OLDFRACDIFF=FRACDIFF
           ENDIF
        ENDIF
     ENDDO time_loop


1000  IF (CONSTRAIN_RUN .EQ. .TRUE. .AND. OUTPUT_LAST .EQ. .FALSE.) THEN 
   call newsavedata()
ENDIF

IF (OUTPUT_LAST .EQ. .TRUE.) THEN 
   DO I=1,DAYCOUNTER
      NEW_TIME=I*DT
      WRITE (10,999) NEW_TIME,LAT, LON, PRESS, TEMP,H2O, CFACTOR, &
           JFACTNO2, JFACTO1D,RO2, &
           (DIURNAL_NEW(JK,I),JK=1,NVAR)
      WRITE (12,999) NEW_TIME,LAT, LON, PRESS, TEMP, CFACTOR,& 
           (DIURNAL_RATES(JK,I),JK=1,NREACT)
   ENDDO
999 FORMAT(E24.16,100000(1X,E24.16))
ENDIF

IF (CONSTRAIN_RUN .EQ. .FALSE.) THEN
   ! close output file
   call newclosedata()
ENDIF
WRITE (6,*) 'Outputed point', counter

enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

!if (.NOT. LAST_POINT) goto 100

IF (CONSTRAIN_RUN .EQ. .TRUE.) THEN 
   call Newclosedata()
ENDIF

END PROGRAM driver


