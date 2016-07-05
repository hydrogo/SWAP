C *******************************************************************
C          Physically based model for permafrost
c	Subroutine model_minus (npar,ncell_bas,ncl_riv,river)          ! ***
      Subroutine model_minus (npar,ncell_bas,river)          ! ***
c made from gswpper6.for 
c calculation from 1 June 
c---------------------------------------------------------------------------
      IMPLICIT NONE
      
      REAL e2,secd,tabs0
      parameter (tabs0=273.16)
      PARAMETER (E2=0.01) 
      PARAMETER (secd=86400.)

      character*1 ichar1
      character*2 ichar2  
      character*3 ichar3,riv  
      character*4 ichar4  
      character*5 ichar5,river
     
c monthly fields        
      REAL LAI_M(12),EFTR_M(12),albedo_m(12),d00m(12)
      real grnfr_m(12),Z0VEGM(12) 
      REAL DBMES(12),DBMES_N(12),DBMES_S(12)
      
c daily fields
      REAL  LAID(366),eftrd(366),grnfrd(366),z0vegd(366)
      real d00d(366),albedod(366)
      REAL SDBMES(366)
c 1.CONSTANTS FOR EVERYTHING
C   a) prescribed parameters and other values                        
      REAL CP,RV,KAPA,SIG,L,G,LPL,KINV,ZCB,ZRB,CIC,RL
      REAL ALBSNOW,RS0,SER,ZCRT,ev,emkl,DEXP,exp2,HPOD,LPOD
      REAL LEVEL0,zwind0,ZEROT1,TSTEP,sncover
      real FUNRO              
      integer MON(12)
      integer iyrfir,numyr,firday,TAYL,ncell
c   b) calculated at the beginning of the program     
      real DELTAT 
      integer NSTEP,mmm
c   c) at the beginning of the program after reading parameters

c 2.FIXED FIELDS needed at each time_step (should be massives (15238))
c   a) given by organizers of GSWP   
      REAL h0,nb,wzav,por,bpar,K0,fi0,hroot,lat,elevat 
c   b) derived by ourselves previously     
      REAL albzim,prles,emkh,sai,hveg,leaf,sryt,amt,xi,yi,ksiot00  
c   c) calculated at the beginning of the program after reading the above parameters
      REAL vdhgr,trmn,zrt,wh,uss,umg,lc,sigk0,hk,POR2,level,zwind 
c 3.FORCINGS (3-hour values) should be read at every time_step
      REAL SW,RDOWN,rainf,snowf,T45,U45,PRES,Q45 
	Real SW_d,RDOWN_d,q45_d,U45_d
C 4.VARIABLES NEEDED FOR THE NEXT time_STEP for each point
C  a) current values (different for the points, so should be massives (15238),initialization before the annual cycle)      
      REAL KSIOT,INTERPR,nhour,TQ,PRHOUR,W0,W1M,w2,u1
      real ice1,sumqbz     
      REAL PSN,HSNOW,WLSN,SUMSN,WLES,SNLES,RSPR,rslespr,SUMTP
      REAL HYD1,HYD2,HYD11,HYD22,hyd1pr,sumhyd,hsr,H2SOIL,H2PSOIL
      REAL t0_d,t2_d,t2daypr,T0DAYPR,L1 
      REAL hhyd0(24)             
      integer NLETO,NZIMA,K80,KODZIM  
c  b) daily values for GSWP2 output (different for the points, so should be massives, initialization before  srok-cycle)   
      REAL swnet_d,lwnet_d,qle_d,qh_d,qg_d,qf_d,qv_d,hprecp_d,precip_d
      REAL ec_d,stok_d,drain_d,snmlt_d,snowfr_d,vd_d,dvsoil_d,dswe_d
      REAL dint_d,snowt_d,vegt_d,avgt_d,radt_d,epot_d,inter_d,et_d
      REAL canint_d,evapsn_d,radtmin,radtmax,esnow_d,subsl_d,ep_d 
C  c) the same for all points (not massives)
      INTEGER IMDAY_1,NYEAR,JDY      
c 5. INTERNAL VARIABLES (need not for the next time_step)
c  a) daily outputs at the end of a day (not massives)
      REAL albedo_d,swe_d,sweveg_d,WK_d,w2_d,soilw_d
      REAL snalb_d,snow_d,hgr_d,snfr_d,ksi_d,ksiot_d
c  b) the rest variables         
      REAL pr,totalpr,L3,L2,XX,YY,ZZ,LAI,LSAI,EF,Z0VEG
      REAL INTERVT,INTER,INTER2,KOND,HGR,ETMMDEC,EPMMDEC,ADED,X1,X2            
      REAL WK,AMB,d,W2PR,PS,RAD,STOK,uzimpr,icepr,tqpr,EP00,ZC3,WH5,B      
      REAL ZC2,ZA3,delsoilm,TC45,A3,A4,T0,EFIZL,P,U2,T2,Q2,EP0,ALBEDO                   
      REAL zimw,SNOWMELT,ET00,DRAIN2,VD,ESNOWMM,w2zon,DELWLSN                   
      REAL GMLT,MLES,SOILWET     
      REAL DRAIN1,STOKCAL,DELSUM,ET,EP,ESNOW,EC,ECMMDEC,XDV
      REAL DPDED,RO2,D0DED,DT0DED,DD,BPRINT,BKEY3,ksipr,rspr0              
      REAL ep00pr,TOTHYDR,DRHYDR,wlsnpr,tles,RADLES,exples   
      REAL ELES,trles,sumsnpr,subles,wlespr,rsles0,SNLESPR,DELWL,DELSN
      REAL elespr,BALANCE,XBAL,pr_int,prh_int,sumqbzpr    
      REAL ecanop,prsolid,xbal_all,ml,bkey2,qf,qv,del_e,eles_wt
      REAL trles_wt,sumtppr,prliq,HL,BZV,DELSWE,DELINTER
      REAL snowt,avgsurft,radt,snfr,snalb,sliqfr,qle,swnet,lwnet,epot
      REAL SUBSOIL,ESOIL,psum,del_les,canopint,dv2,dv3,bbb
      REAL D00,ALBLET,prhour0,snowfrz,zerot2,zerot3,delp,hg,hg1           
      integer ISTEP,kodsubl,KODSUBp,KEYITER,KODHYD,kodhyd2,keyrecal
      INTEGER ncl,land,ico,jco,year_print,jdy_pr  
      REAL lon,surrun_transf,drain_transf,k_lw,k_sw 
c 6. TIME AND CELL etc.CYCLES: I_CYCLE (DAILY), II_CYCLE (SROK), ICELL_CYCLE (POINTS)     
      INTEGER icell,iyr_1,I,II,KODEND,KODPRINT,IMDAY,IYR,IM
c	integer npar,ncell_bas,ncl_riv(ncell_bas),ibas,icel_cal
      integer npar,ncell_bas,ibas
      real opt_par(100),k0_opt,hroot_opt,krain,ksnow,h0_opt, kalblet

c *********    NEW  ***********   NEW   ***********   
c 7. Variables connected with permafrost_program
c  a) needed at the next time step  
      REAL ice2,UU2,KSI_UP,t_frost,ksiotpr,hsig0,b0,tga,nb2,wzav2,k02
      integer marfir,kksiot,k_snowt,kod_exp, kod_prec,kod_calibr,kod_cal 
c  b) not needed at the next time step  
      real t_frost0,TZV,za2,y1,y2,ksiot0,tayot

                
      COMMON /METEO/ U45,T45,Q45,RAD,PRES
      COMMON /PARAM1/KAPA,SIG,SER,L,LPL,CP,L1,L2,L3,WH,KINV,RV,ZA2
      COMMON /EVAP/ ETMMDEC,EPMMDEC,LAI,LSAI
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy  
      common /year1/ NB,k0,HROOT,bpar,wzav,fi0,POR,por2,umg,a4,u2,lc,
     !              h2soil,W0,WK,w2,PS,PSN,STOK,
     !              EP00,PR,T2,Q2,kodhyd,kodhyd2,hpod,emkl,h0,
     !              pr_int,ksiot,w1m,lpod,t0,hgr,uu2,
     !              ice2,RL,u1,ice1,kodzim,drain2,tstep,drain1,TZV,
     !              ZEROT1,interpr,k02,nb2,wzav2,epot,ef,ev,w2pr,ksiot0

      common /Cleto/ nhour,prhour,tq,hk,sigk0,trmn,leaf,
     !              ET00,inter,inter2,kond
      common /Czima1/ ZA3,ZC2,ZRB,USS,RSpr,IM,
     !              t0daypr,SUMTP,HSNOW,
     !              WLSN,DELWLSN,gmlt,SUMSN,DELSUM,B,VD,
     !              SNOWMELT,ESNOWMM,KEYITER,
     !              MLES,ELES,WLES,SNLES,rslespr,WLSNPR,trles,
     !              sncover,tga,b0,marfir,sumqbz,emkh,prh_int,prsolid,
     !              KSI_UP,T_FROST,hsig0,ML,del_les,snlespr,
     !              ZC3,snowfrz,del_e,zimw,ksiot00

      common /hydro/ xi,yi,lat
	common /opt_par/ opt_par,kod_exp, kod_prec,albsnow,iyrfir,numyr,
     !	firday,kod_calibr,kod_cal

C input of the months' length      
      DATA MON/31,28,31,30,31,30,31,31,30,31,30,31/
C input of the data for the calculation of G by Budyko
      DATA DBMES_N/-0.23,-0.15,0.08,0.15,0.23,0.23,0.19,0.12,-0.08,
     ! -0.12,-0.19,-0.23/        ! for the northern hemisphere
      DATA DBMES_S/0.19,0.12,-0.08,-0.12,-0.19,-0.23,-0.23,-0.15,0.08,
     ! 0.15,0.23,0.23/    ! for the southern hemisphere
     
C input of the physical constants
C     [cp]=[J/kg/K], [RV]=[J/kg/K;m**2/s**2/k], [sig]=[W/m**2/grad**4],
C     [L]=[J/kg],[G]=[m/s**2]
C          lpl  - (kal/g) specific heat of ice fusion (plavleniya)
C          kinv - (m**2/s) coefficient of kinematic vyazkosti
C          zcb  - (kal/g/grad) cpecific heat capacity of water
C          zrb  - density of water (g/cm^3) 
C          cic  - specific heat capacity of ice (kal/g/grad)         

C          RL - ice density (g/cm^3);
      DATA CP,RV,KAPA,SIG,L,G,LPL,KINV,ZCB,ZRB,cic,RL/1005.,
     *287.05,0.43,0.000000057,2501000.,9.81,80.,0.000013,1.,1.,
     *0.5,0.9/ 
        FUNRO(XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      
       tstep=3.0    
       nstep=24./tstep
       deltat=tstep*3600.
C  PILPS' parameters  which are the same for all grid cells
      LEVEL0=2.          ! the reference height for atmospheric forcings (air temperature and air humidity) (m)
      zwind0=10.         ! the reference height of wind speed measurements (m) 
c      albsnow=0.75       ! albedo of fresh snow cover
      RS0=0.15           ! snowpack density for fresh snow (g/cm^3)
      SER=1.             ! coefficient of "serosti"
      ZCRT=0.22          ! (kal/g/grad) specific heat capacity of dry soil
      zerot1=90.         ! transition from positive to negative air temperatures at autumn (number of day) 
      zerot2=111.        ! many-year transition of monthly air temperature across mean annual temperature in spring (number of day) 
      zerot3=283.        ! many-year transition of monthly air temperature across mean annual temperature in autumn (number of day) 
      ev=1               ! parameter for the calculation of transpiration
      HPOD=0.            ! thickness of podstilki
      LPOD=0.00015       ! TEPLOPROVODNOST' PODSTILKI
      emkl=0.1           ! coefficient for the calculation of interception capacity for liquid precipitations
      TAYL=3             ! for the calculation of water table depth  (in days)
      DEXP=0.2           ! parameter for heat exchange in the forest
      exp2=1.            ! parameter for the forest     
ccc      iyrfir=1983        ! the first full year of the calculations (for the calculations started in summer)
ccc      numyr=14          ! number of years for the calculation
ccc      firday=182         ! the first day for the calculations (if they are started in summer)
c        numyr=10        ! for Kolyma
c         iyrfir=1970    ! for Kolyma
      mmm=TAYL*nstep
	riv=river
c
c        write (*,*) 'Input the number of grid cell'
c      read(*,*) icel_cal
c       ncell=15238                ! number of grid cells  
c        ncell=2
c open input files with parameters which are different for each grid cell
        open (7, file='c:\d\rivers\fixed_fields\
     !fixed_param_'//riv//'.txt')
        open (10,file='c:\d\rivers\monthly_fields\monthly_param_'//
     !	  riv//'.txt') 
c *************************** the beginning of the cycle by cells*******************************
c      icel_cal=0
	ncell=ncell_bas
      do icell=1,ncell 
c input model parameters which are different for the cells
c  
      read (7,*) ncl,land,ico,jco,lon,lat,nb,wzav,por,bpar,k0,fi0,
     !      amt,sryt,albzim,prles,sai,leaf,emkh,xi,yi,hroot,h0,hveg,
     !      elevat,ksiot00
      read (10,*) land,(albedo_m(i),i=1,12),(D00M(I),I=1,12),
     !            (grnfr_m(i),i=1,12),(Z0VEGM(i),i=1,12),
     !            (LAI_M(I),I=1,12),(EFTR_M(i),i=1,12) 

c       do ibas=1,ncell_bas  
c	    if(ncl_riv(ibas).eq.icell) icel_cal=icell
c	 enddo 
	
	  
c	if (icell.LT.icel_cal.or.icell.GT.icel_cal.or.sryt.ge.-1.) then
      if (sryt.ge.-1.) then
	    GOTO  1879  
	  else   
	  
c********************************************************************************
       if (por-nb.le..01) nb=por-.01
	 if(h0.lt.hroot) h0=hroot+0.1
	 krain=1
	 ksnow=1
	 k_sw=1
	 k_lw=1	
C		calibration of albzim, albsnow, albedo								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if (kod_cal.ne.1) 	then
		
c		albsnow=opt_par(5)
		albzim=opt_par(4)*albzim
	    kalblet=opt_par(6)
	    albedo_m=kalblet*albedo_m
	    nb=opt_par(3)*nb
    				if (por-nb.le..01) nb=por-.01
          k0_opt=opt_par(1)
	    k0=k0*exp(k0_opt)
	    hroot_opt=opt_par(2)
          hroot=hroot_opt*hroot
   	    h0=hroot*opt_par(9)
                  if(h0.lt.hroot) h0=hroot+0.1
          krain=opt_par(14)
	    ksnow=opt_par(15)
	    k_sw=opt_par(7)
      	k_lw=opt_par(8)
c       
c          write(*,*) ncl,k0,k0_opt,h0,h0_opt
        endif
       endif	     

c*****************************************************************************************     
c open input file with forcings  for i-cell AND OUTPUT FILES  
c open output files for  i-cell
    
	if (land.gt.9999) then
	  write(ICHAR5,'(I5)') land
        OPEN (55, FILE='c:\d\rivers\res\point '//ichar5//'.txt')
      else
	  write(ICHAR4,'(I4)') land
        OPEN (55, FILE='c:\d\rivers\res\point 0'//ichar4//'.txt')
      endif
      Write(55,*) 'Year	Day	Evap	Rainf	Snowf	Qs_transf	Qsb_tr
     !ansf   SW   RDOWN    T2    U2    Q2'
	
c	if (icell.eq.icel_cal) then     
       if (ncl.lt.10) then
          WRITE(ICHAR1,'(I1)') ncl
          open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !		ichar1//'.dat')
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar1//'.dat')
			if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !ichar1//'.dat')
c            OPEN (4, FILE='c:\d\rivers\res\evap.'//ichar1)  
c            OPEN (33, FILE='c:\d\rivers\res\Eb.'//ichar1)
c            OPEN (444, FILE='c:\d\rivers\res\WB.'//ichar1)
c            OPEN (12, FILE='c:\d\rivers\res\d_o1o378.'//ichar1)
c            OPEN (15, FILE='c:\d\rivers\res\d_o2o4.'//ichar1)
c            OPEN (19, FILE='c:\d\rivers\res\d_o5o6.'//ichar1) 
c            OPEN (333, FILE='c:\d\rivers\res\balan_w.'//ichar1) 
c            OPEN (6, FILE='c:\d\rivers\res\balan_t.'//ichar1) 
c            OPEN (55, FILE='c:\d\rivers\res\point.'//ichar1)
       else
         if (ncl.ge.10000) then
           WRITE(ICHAR5,'(I5)') ncl
           open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !		 ichar5//'.dat')
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar5//'.dat')
			if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !   ichar5//'.dat')
c             OPEN (4, FILE='c:\d\rivers\res\evap.'//ichar5)  
c             OPEN (33, FILE='c:\d\rivers\res\Eb.'//ichar5)
c             OPEN (444, FILE='c:\d\rivers\res\WB.'//ichar5)
c             OPEN (12, FILE='c:\d\rivers\res\d_o1o378.'//ichar5)
c             OPEN (15, FILE='c:\d\rivers\res\d_o2o4.'//ichar5)
c             OPEN (19, FILE='c:\d\rivers\res\d_o5o6.'//ichar5) 
c             OPEN (333, FILE='c:\d\rivers\res\balan_w.'//ichar5) 
c             OPEN (6, FILE='c:\d\rivers\res\balan_t.'//ichar5) 
c		   OPEN (55, FILE='c:\d\rivers\res\point.'//ichar5)   
         else 
           if (ncl.ge.10.and.ncl.lt.100) then
              WRITE(ICHAR2,'(I2)') ncl
              open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !			ichar2//'.dat') 
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar2//'.dat')
     		if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !   ichar2//'.dat')  
c                OPEN (4, FILE='c:\d\rivers\res\evap.'//ichar2)  
c                OPEN (33, FILE='c:\d\rivers\res\Eb.'//ichar2)
c                OPEN (444, FILE='c:\d\rivers\res\WB.'//ichar2)
c                OPEN (12, FILE='c:\d\rivers\res\d_o1o378.'//ichar2)
c                OPEN (15, FILE='c:\d\rivers\res\d_o2o4.'//ichar2)
c                OPEN (19, FILE='c:\d\rivers\res\d_o5o6.'//ichar2) 
c                OPEN (333, FILE='c:\d\rivers\res\balan_w.'//ichar2) 
c                OPEN (6, FILE='c:\d\rivers\res\balan_t.'//ichar2) 
c                OPEN (55, FILE='c:\d\rivers\res\point.'//ichar2) 
           else
              if (ncl.ge.100.and.ncl.lt.1000) then
                WRITE(ICHAR3,'(I3)') ncl
              open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !			ichar3//'.dat') 
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar3//'.dat')
     		if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !  ichar3//'.dat')			        
c                  OPEN (4, FILE='c:\d\rivers\res\evap.'//ichar3)  
c                  OPEN (33, FILE='c:\d\rivers\res\Eb.'//ichar3)
c                  OPEN (444, FILE='c:\d\rivers\res\WB.'//ichar3)
c                  OPEN (12, FILE='c:\d\rivers\res\d_o1o378.'//ichar3)
c                  OPEN (15, FILE='c:\d\rivers\res\d_o2o4.'//ichar3)
c                  OPEN (19, FILE='c:\d\rivers\res\d_o5o6.'//ichar3) 
c                  OPEN (333, FILE='c:\d\rivers\res\balan_w.'//ichar3) 
c                  OPEN (6, FILE='c:\d\rivers\res\balan_t.'//ichar3) 
c               	OPEN (55, FILE='c:\d\rivers\res\point.'//ichar3)
              else
                 WRITE(ICHAR4,'(I4)') ncl
              open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !			ichar4//'.dat') 
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar4//'.dat')
     		if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !  ichar4//'.dat')			     
c                  OPEN (4, FILE='c:\d\rivers\res\evap.'//ichar4)  
c                  OPEN (33, FILE='c:\d\rivers\res\Eb.'//ichar4)
c                  OPEN (444, FILE='c:\d\rivers\res\WB.'//ichar4)
c                  OPEN (12, FILE='c:\d\rivers\res\d_o1o378.'//ichar4)
c                  OPEN (15, FILE='c:\d\rivers\res\d_o2o4.'//ichar4)
c                  OPEN (19, FILE='c:\d\rivers\res\d_o5o6.'//ichar4) 
c                  OPEN (333, FILE='c:\d\rivers\res\balan_w.'//ichar4) 
c                  OPEN (6, FILE='c:\d\rivers\res\balan_t.'//ichar4) 
c	            OPEN (55, FILE='c:\d\rivers\res\point.'//ichar4)
              endif
           endif
         endif
       endif      
c 9999 continue     
c making headers for output files
c       write(4,4455) '      PotEvap       Ecanop        Tveg          Es
c     !oil       RootMoist      SubSnow       subsoil       ec_mm       
c     !  subles   SWet  SMoist2    WaterT     SubSurf        eles'
c      write(33,4132) 'nyear  jdy  i   Swnet     Lwnet   Qle      Qh     
c     !     Qg        Qf      Qv     SnowT     VegT     AvgSrT   RadT   n 
c     !z nl kd t2daypr t0daypr  lai   lsai  ef    z0veg  d00  albl   bbb'    
c      write(444,1111) 'jdy    i      Evap           Qs           Qsb    
c     !       Qsm           Qfz           Qst         DSM    DelSWE   Del
c     !In    balance        Snowf         Rainf         totalpr       can
c     !opInt      xdv           dv2            dv3        etmmdec   epmmd
c     !ec   inter    esnowmm      eles     trles      SWE     SWEVeg  Alb
c     !edo  SnFr    ksi_up ksiot  SnDep   SAlbedo'
 1111 Format(A354)    
      
 4132 format(A181)
 4455 format(A175)   
c      endif
	            
c h0 - the depth of unpermeable layer (m)
c NB - field capacity
c WZAV - wilting point;   
C POR - porosity (m^3/m^3)
C BPAR - Clapp-Hornberger b parameter;
C K0 - hydraulic conductivity at saturation (m/s);
c FI0 - saturated matric potential (m)
c HROOT - root-zone depth (m);
c Z0VEG and Z0SNOW - roughness coefficients of vegetation and snow (m); 
c LAI_M = total_LAI  
c grnfr_m - greenness fraction (green_LAI/total_LAI 
C d00 - zero displacement height  (m)
c albedo_m - snow_free albedo
c hveg - the height of vegetation (m)
C AMT - amlitude of air temperature (in C),
c SRYT - soil temperature at the depth where seasonal variations are absent (in C);
C LEAF -size of grass's leaves in (cm)
c WH - min. amount of unfrozen water 
C LPOD - TEPLOPROVODNOST' PODSTILKI,
C HPOD - thickness of podstilki
C ZEROT1 - transition from positive to negative air temperatures at autumn;
C LEVEL0 - the reference height for atmospheric forcings (m);
c ZWIND0 - the height if wind speed measurements (m)
C ZCRT (kal/g/grad) specific heat capacity of dry soil;
c HGR - depth to water table (m)
c K02 - hydraulic conductivity at saturation in the 2nd zone (m/s);
c 
c transformation of input parameters 
       tayot=-25.2+2.42*lat+0.0141*elevat      ! calendar day of the beginning of soil thawing
       if (tayot.lt.0.) tayot=0.
       if (tayot.gt.firday) tayot=firday
       zerot1=75.+0.16*tayot                   ! calendar day of the beginning of soil heating
       if (zerot1.ge.(tayot-5)) zerot1=tayot-5.
c       write(*,*) tayot,zerot1
c       stop
      
        do i=1,12
          if (z0vegm(i).eq.0.) z0vegm(i)=0.0024
        enddo  
       
       if (prles.eq.0.) then
          emkh=0.
          sai=0.
          hveg=0.2
      	d00m=0.
      	z0vegm=hveg*2./3.
       endif 
          
       por2=por
       k02=k0   
       wzav2=wzav
       IF ((POR-NB).LT.0.1) NB=POR-0.1  
       nb2=nb

       tzv=sryt
       if (sryt.gt.0.) sryt=0.   ! for the first approximation for the territories without permafrost (really there should be another program)
       hroot=hroot*1000.      ! translation from m to mm
       h0=h0*1000.            ! translation from m to mm
c CALCULATION OF the rest MODEL PARAMETERS:
       level=level0+hveg    ! m
       zwind=zwind0+hveg    ! m  
       vdhgr=por-nb       
C calculation of TRMN=ET00/EP00  
        IF (LEAF.EQ.0.) THEN
          TRMN=0.
        ELSE
          TRMN=0.825*0.5/SQRT(LEAF)+0.784
        ENDIF  
C calculation of soil parameters
         ZRT=(1.-POR)*2.65
C          2.65 - density of hard fraction of the soil (g/cm^3) 
C          ZRT  - density of dry soil (g/cm^3)  
c translation of FI0 into positive value
         fi0=fi0/(-1.)
C  calculation by Kalyuzhnyi and Pavlova
c         WH=1.04*NB-0.06*ZRT 
c         wh=0.18  
          wh=(nb+wzav)/2.
         USS=WZAV/1.35
         UMG=USS
         lc=(0.102*exp(4.7*umg)+0.45*zrt-0.35)*0.0024
C            UMG - maximum hygroscopicity
C            USS - static water  
         k02=k02*1000.*deltat     ! mm/time_step
C calculation of standard deviation of saturated hydraulic conductivity (SIGK0) 
C  translation of m/s into mm/min
         K0=K0*60.*1000.                        
         SIGK0=0.9*K0**0.6                 
         IF ((SQRT(3.)*SIGK0).GT.K0) SIGK0=K0/SQRT(3.)
C  translation of mm/min into mm/hour
         K0 = K0 * 60.
         SIGK0 = SIGK0 * 60.                   
C calculation of the effective capilar. potential HK (mm)
         X1 = 2. + 3. / BPAR
         X2 = FI0 ** X1
         X1 = 1. - X1
         HK = -X2 * FI0 ** X1 / X1
         HK = HK * 1000.  
C ACCOUNTING FOR  SOUTHERN OR NORTHERN hemisphere (for the calculation of G by Budyko)
         DO I = 1,12
            IF (lat.lt.0.) THEN
               dbmes(i)=dbmes_S(i)
            ELSE
               dbmes(i)=DBMES_N(I)
            ENDIF   
         enddo
ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii   
C initial conditions    
         W0 = nb
         W1M=nb  
         w2=w0
          u1=w0
          uu2=w0
          ice2=0.
          ice1=0.
         psn=0.
c         HSNOW=100000.  ! cm
         hsnow=0.
         KSI_UP=0.
c         ksiot=90.           ! cm  initialization for the 1st July   
         ksiot=ksiot00         ! cm  initialization for the 1st  July
         ksiotpr=ksiot*10.     ! mm daily ksiot at the previous day
         sumsn=hsnow*10*rs0   ! mm
         rspr=rs0
         WLSN=0.   
         snles=0.
         wles=0.
         interpr=0.
         t2daypr=tabs0+5.
         t0daypr=tabs0+5.  
         sumhyd=0.
         hsr=sumhyd/mmm 
         h2soil=h0
         h2psoil=h2soil
         T_frost=0.
         hsig0=0.  
         if (ksiot00.lt.5.) hsig0=50.
         do i=1,mmm
            hhyd0(i)=0.
         enddo          
         kodhyd2=0     
            k80=0
             HYD1=1.e-5
             HYD2=1.e-5
             HYD11=1.e-5
             HYD22=1.e-5  
            hyd1pr=1.e-5 

           
           TQ = -(DELTAT/3600./2.)
           NHOUR = 0.
           PRHOUR =NHOUR
             SUMTP = 0.
           
             sumqbz=0. 
             NLETO=0
             NZIMA=0  
              
        marfir=mon(1)+mon(2)+mon(3)+15
        sncover=1.
        kksiot=0
c formula by Bracht         
         L1 = 0.0049 * RS0 ** 2 
C*****************************************************************
      NYEAR=1
      kodend=1                    
      kodprint=0
	jdy_pr=0
      year_print=iyrfir-1
c      kodprint=1
    
  100   CONTINUE
      if (kodend.eq.4) kodprint=1     
cc        if (nyear.ge.4) kodprint=1
c         WRITE(*,*) '************'
c         WRITE(*,*) NYEAR, kodend,kodprint
c         WRITE(*,*) '____________' 
         IMDAY=365
         MON(2)=28 
         jdy=firday       
C IYR - year of run
         IYR=iyrfir-1+NYEAR  
         iyr_1=iyr-1
c         write(*,*)  iyr   
C LEAP YEARS  
         if (abs(iyr/4.-int(iyr/4.)).eq.0.) then
            IMDAY=366
            MON(2)=29
         END IF 
         imday_1=365
         if (abs(iyr_1/4.-int(iyr_1/4.)).eq.0.) then
            IMDAY_1=366
            jdy=firday+1
         endif                                    
C calculation of daily values of LAI, LSAI,  B:
         call messyt(mon,dbmes,SDBMES,IMDAY)
         CALL MESSYT(MON,LAI_M,LAID,IMDAY) 
         CALL MESSYT(MON,grnfr_m,grnfrd,IMDAY)
         CALL MESSYT(MON,EFTR_M,eftrd,IMDAY) 
         call messyt(mon,z0vegm,z0vegd,imday) 
         call messyt(mon,D00M,D00D,imday)  
         call messyt(mon,albedo_m,albedod,imday)  
C number of time steps in a year         
         IM=IMDAY*nstep
c ____________________________________________________________________
C CALCULATION DAY BY DAY
      DO 30 II=1,IMDAY    
      
      
       if (jdy.gt.IMDAY_1) jdy=jdy-IMDAY_1 
         IF (NYEAR.EQ.NUMYR.AND.JDY.EQ.1) GOTO 1333  
c        IF (NYEAR.EQ.NUMYR.AND.ii.EQ.1) stop             
c  =0 - for the calculation of daily totals:
      k_snowt=0
      swnet_d=0.
	sw_d=0.
	rdown_d=0.
	u45_d=0.
	q45_d=0.
      lwnet_d=0.
      qle_d=0.
      qh_d=0.
      qg_d=0.
      qf_d=0.
      qv_d=0.
      snowt_d=0.
      vegt_d=0.
      avgt_d=0.
      radt_d=0.
      t2_d=0.
      t0_d=0.
      precip_d=0.
      hprecp_d=0. 
      ec_d=0.   
      et_d=0.
      ep_d=0.
      inter_d=0.
      stok_d=0.
      drain_d=0.
	surrun_transf=0.
	drain_transf=0.
      snow_d=0.
      WK_d=0.
      esnow_d=0.
      albedo_d=0.
      swe_d=0.
      snmlt_d=0.
      snowfr_d=0.
      ksi_d=0.
      ksiot_d=0.   
      w2_d=0. 
      vd_d=0.
      hgr_d=0.
      sweveg_d=0.
      dvsoil_d=0.
      dswe_d=0.
      dint_d=0.
      soilw_d=0.
      epot_d=0.
      subsl_d=0.
      snfr_d=0.
      snalb_d=0.
      canint_d=0.
      radtmax=0.
      radtmin=1000.
C calculation of the ground heat flux by Budyko
C  relation between annual amplitude of air temperature and 
          AMB = (.184 * AMT * 1000. / 30.4)/2.
C  calculation of the ground heat flux B (cal/cm**2/day)
          b = sdbmes(jdy)* AMB
C  transformation of (cal/cm**2/day) into (W/m**2)
          b = b / 1440. * 697.37 
C time-step value of LAI, LSAI, EF, Z0VEG, D00, ALBLET:    
          lai=LAID(jdy)*grnfrd(jdy)   ! calculation of green_LAI
          LSAI=LAID(jdy)+sai          ! total_LAI + SAI
          ef=eftrd(jdy)  
          z0veg=z0vegd(jdy)
          D00=D00D(JDY)
          alblet=albedod(jdy)         
ccccc     zaglushka
c       b=1.2                     
c       lai=0.8
c       lsai=1.2
c       ef=1.
c       z0veg=1.
c       d00=4.8
c       alblet=0.12      
cccccccccccccccccccccccccccccccc 
      bbb=b
c
          IF (T2daypr.GT.tabs0) THEN 
             NLETO=NLETO+1
          ELSE
             NLETO=0
          ENDIF      
          IF (T2daypr.Le.tabs0) THEN
             NZIMA=NZIMA+1
          ELSE
             NZIMA=0
          ENDIF                  
          if ((ksiotpr).gt.hroot) then
             kksiot=kksiot+1
          else
             if (kksiot.lt.5) kksiot=0
          endif   

C-------------------------------------------------------------------------        
C CALCULATION STEP BY STEP 
        DO 200 I=1,nstep  
c	  if (ncl.gt.13500) write(*,*) nyear,ii,i 
c        write(*,*) nyear,ii,i
c for the calculation of snow coverage for 1 March      
       if (jdy.eq.marfir) then
c           y1=(sumsn+wlsn)*2.
c           y2=(sumsn+wlsn)*0.
            y1=(sumsn+wlsn)*1.
            y2=(sumsn+wlsn)*1.
           b0=y1-y2
           tga=y1-y2  
       endif    
        istep=(jdy-1)*nstep+i    ! number of current time step in a year  
        if (((ksiot.eq.0.).and.(jdy.gt.240)).or.(ksiot00.lt.5.)) then    
           t_frost=t_frost+deltat 
        else      
           t_frost=0. 
           hsig0=0.
        endif 
        
        if ((sumsn+wlsn).eq.0.) rspr=rs0
        if ((snles+wles).eq.0.) rslespr=rs0
        INTER2=0.  
        inter=0.
        KOND=0.
        b=bbb  
C calculation of soil parameters which depend on soil moisture
C CALCULATION BY (3.11), P.111 FROM Kalyuzhnyi,Palova,Lavrov (W/m/grad)           
       L3=0.102*EXP(4.7*W1M)+0.45*ZRT-0.35 
       if (l3.le.0.) write(*,*) 'l3<0',l3 
         L3=L3*0.0024    ! translation into kal/cm/c/grad                                                        
         ZC3 = ZCRT * ZRT + ZCB * ZRB * W1M
C   WH5 - the amount of unfrozen water under the soil temperature =-5C (Kalyuzhnyi and Pavlova)  
         WH5=0.94*WZAV+0.017*ZRT
         ZC2=ZCB*U1*ZRB+CIC*ICE1*RL+ZCRT*ZRT+(WH-WH5)*LPL*ZRB/5.
         L2 = 1.15 * L3                                  
         ZA3 = L3 / ZC3 
         za2=l2/zc2   
C  L3 - teploprovodnost' of unfrozen soil    
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C reading values of forcing factors
C precipitatiom - kg/m^2/s^1=mm/s; radiation - W/m**2;temperature -K;
C wind speed - m/s at 10 m; pressure - Pa; air specific humidity - kg/kg
        READ (1,1030) SW,RDOWN,rainf,snowf,T45,U45,PRES,Q45   
	  if(kod_exp.eq.2) read(2,*) rainf,snowf
	   
 1030   format (2F10.3,2e14.4,3f10.3,e14.4)
ccc      rainf=rainf*krain
ccc	  snowf=snowf*ksnow
		rdown=rdown*k_lw
	    sw=sw*k_sw
	  

	if(kod_prec.eq.2) then 
              read (60,*) totalpr
      else
            totalpr=rainf+snowf
      endif 

              if (T45-tabs0.lt.0.) then
		        totalpr=totalpr*ksnow
		    else
		        totalpr=totalpr*krain
	        end if
	
		   
c	write (*,*) ncell, ii, tabs0, t45, totalpr, kod_prec
         
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
c transformation of forcings
c the wind speed at 10 m                             
        IF (U45.LT.0.5) U45=0.5  
c translation of pressure from Pa into gPa=mb               
        pres=pres/100.   
c translation of precipitation from mm/s into mm/time_step
        totalpr=totalpr*deltat    
        PR=totalpr 
c determination of KODZIM. IF KODZIM=1 then ZIMA-subprogram is used,
c                          if KODZIM=0 then LETO-subprogram is used
      if (jdy.le.210) then                           
         IF (((HSNOW.EQ.0.).AND.(kksiot.gt.5).and.(NZIMA.LE.7)).and.
     !      (NLETO.ge.1)) THEN
             kodzim=0
         else
             kodzim=1
         endif
      else
         IF ((((T2daypr.lt.TABS0).or.(HSNOW.gt.0.).or.(KSI_UP.gt.0.).
     !    or.(NZIMA.GE.7)).and.(NLETO.LE.7)).or.(hsig0.gt.0.))  THEN 
             kodzim=1
         else
             kodzim=0
         endif
      endif           
C calculation of albedo, sumsn - is  snow pack (mm)         
c        albsnow=0.834-21.988*(rspr-0.1)**3 
       albsnow=opt_par(5)-21.988*(rspr-0.1)**3
       if (albsnow.lt.0.1) albsnow=0.1 

c         albsnow=0.8
       IF ((sumsn+wlsn). LT. 10.0) THEN 
        ALBEDO = alblet + (ALBSNOW - alblet) * SQRT((sumsn+wlsn)*0.1)
       ELSE
         ALBEDO=ALBSNOW
        ENDIF   
c        if (albedo.lt.0.) then
c          write(*,*) alblet,albsnow,rspr                                   
c          stop  
c        endif  
C**************************************************        
        IF (KODZIM.EQ.1) THEN
           exples=EXP(-PRLES*LSAI)
           radles=SW*(1.- ALBzim)+RDOWN 
           rad=(SW*(1.-albedo)+rdown)*exples 
        else
           exples=1
           albedo=alblet
           RAD = SW * (1. - ALBedo) + RDOWN
           radles=0.        
        ENDIF       
C**************************************************     
        TC45 = T45 - TABS0  
        A3 = (0.623 * 6.1 / PRES) * EXP(17.1 * TC45 / (235. + TC45))
        A4 = A3 * 17.1 * 235. / (235. + TC45) ** 2 
C iterations' block 
C//////////////////////////////////////////////////////////////////////////      
       if (kodzim.eq.1) zimw=1/(1+8.*ice1)**2*(u1-wzav)/(2./3.*nb) 
      CALL ITBLOCK1(KEYITER,HSNOW,KSI_up,KSIOT,A3,A4,
     *D00,LEVEL,TABS0,B,EFIZL,P,U2,T2,Q2,T0,G,EP0,
     *zwind,e2,WLSN,kodzim,radles,tles,LSAI,expLES,d,
     !MLES,ELES,WLES,SNLES,DEXP,exp2,z0veg,ef,ev,trles,LPOD,HPOD,
     ! zimw,SRYT,T_FROST,HSIG0,BZV,HL,kodsubp,kodsubl,ksiot00,
     ! zerot2,zerot3) 
c//////////////////////////////////////////////////////////////////////////
C translation of W/m^2 into mm/day     
        EP00=(EP0/697.37)*24. 
        eles=(eles/697.37)*24. 
        trles=(trles/697.37)*24.
C translating mm/day  into mm/time_step
        EP00=EP00/nstep       
        eles=eles/nstep 
        trles=trles/nstep 
        kodhyd=0
c the depth of the second soil zone        
c        hgr=h0-hsr*1000./VDHGR           
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
      IF (KODZIM.EQ.1) THEN 
         STOK=0.
         SNOWMELT=0.  
         ET00=0.
         nhour=0.
         TQ = -(DELTAT/3600./2.)
c--------------------------------
        w2pr=w2
        uzimpr=UU2      !!!!!!!!!!!
        icepr=ice2      !!!!!!!!!!!!
        ksipr=ksi_UP 
        ksiot0=ksiot
        sumsnpr=sumsn
        sumtppr=sumtp
       
        sumqbzpr=sumqbz
        wlsnpr=wlsn 
        wlespr=wles        
        SNLESPR=SNLES
        rspr0=rspr 
        elespr=eles
        rsles0=rslespr       
        t_frost0=t_frost
  998   CONTINUE       
        CALL ZIMA1(EPMMDEC,PRES,T45,LSAI,DELWL,DELSN,tabs0) 
          ETMMDEC=0.  
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ELSE
c §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
         wlespr=wles
         del_les=0.    
         del_e=0.                                         
         kodsubp=0
         kodsubl=0 
         snles=0.
         wles=0.
         snlespr=0.
         sumqbz=0.
         DELSUM=0.
         DELWLSN=0.
         wlsn=0.
         hsnow=0.
         sumsn=0.    
         snowmelt=0. 
         snowfrz=0.
         vd=0.
         esnowmm=0.  
         gmlt=0.  
         prsolid=0. 
         ksi_up=0.
C---------------------
         w2pr=w2
         ep00pr=ep00 
         prhour0=prhour
         tqpr=tq         
         ksiot0=ksiot
        
  999    continue  
        CALL LETO1(PRES,zrb,tabs0,wlespr)  
         hsig0=0.
       
         DELWL=0.
         DELSN=0.
c DRAIN2 - is drainage from the 2nd layer 
c DRAIN1 - is drainage drom the 1st layer    
c §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§           
      END IF
        stokcal=stok
c hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
c the calculation of hydrograph:
       
  888  call hydrogr(stok,DRAIN2,tothydr,kodhyd,hyd1,hyd2) 
         if (kodhyd.eq.1) then  
          IF (KODZIM.EQ.1) THEN 
              pr=totalpr
              w2=w2pr   
              UU2=uzimpr
              ice2=icepr         
              ksi_UP=ksipr
              ksiot=ksiot0
              sumsn=sumsnpr
              sumtp=sumtppr
             
              sumqbz=sumqbzpr 
              wlsn=wlsnpr   
              wles=wlespr  
              eles=elespr
              SNLES=SNLESPR   
              rspr=rspr0 
              rslespr=rsles0
              t_frost=t_frost0
             goto 998
          else   
             w2=w2pr 
             inter2=0.
             ep00=ep00pr  
             ksiot=ksiot0
             prhour=prhour0
             nhour=nhour-DELTAT/3600. 
           
             pr=totalpr
             tq=tqpr
             goto 999  
          endif   
         endif              
       call hydrogr(0.,DRAIN2,drhydr,kodhyd,hyd11,hyd22) 
      stok=tothydr-drhydr
      if (stok.lt.0.) then
          stok=0.
          drhydr=tothydr
      endif                              
C RECALCULATION OF WATER BALANCE COMPONENTS      
c     ********************************
          k80=k80+1
          sumhyd=sumhyd+hyd1-hhyd0(k80)
          hhyd0(k80)=hyd1
          hsr=sumhyd/mmm 
          if (k80.eq.mmm) k80=0
c hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
C translation of mm/time-step into W/m^2 (for ET,EP,Esnow and INTERVT)         
        IF (KOND.LT.0.) EPMMDEC=KOND 
        ET=ETMMDEC*nstep/24.*697.37
        EP=EPMMDEC*nstep/24.*697.37
        esnow=esnowmm*nstep/24.*697.37     
        INTERVT=INTER*nstep/24.*697.37   
        eles_wt=eles*nstep/24.*697.37
        trles_wt=trles*nstep/24.*697.37  
        del_e=del_e*nstep/24.*697.37                  
        EC=ET+EP+INTERVT+esnow
        ECMMDEC=ETMMDEC+EPMMDEC+INTER+esnowmm
        XDV=(WK-W0)*HROOT                
c w1000 - soil moisture in 1 m
        if (hroot.lt.1000.) then
           w1M=(wk*hroot+w2*(1000.-hroot))/1000.
        else
           w1m=wk
        endif      
c calculation of water table        
        delp=por2-ice2
        if (nb2.ge.delp) then 
          if (uu2.lt.delp) then
              hg=0.
          else
              hg=h2soil-hroot
          endif
        else
          if (uu2.le.nb2) then
              hg=0.
          else
              hg=(uu2-nb2)/(delp-nb2)*(h2soil-hroot)
          endif                      
        endif
        
      if (hg.ge.(h2soil-hroot)) then
           delp=por-ice1 
        if (nb.ge.delp) then 
          if (u1.lt.delp) then
              hg1=0.
          else
              hg1=hroot
          endif
        else
          if (u1.le.nb) then
              hg1=0.
          else
              hg1=(u1-nb)/(delp-nb)*hroot
          endif                      
        endif                           
      else
         hg1=0.  
      endif   
      hgr=h2soil-hg-hg1
      
C RECALCULATION of T0, B, P and EFIZL (U2 in m/s; D0DED and DD in cm/s)            
        DPDED=0.76*U2/(0.82*U2**0.5+1.)  
C  air density in kg/m^3 
        RO2=FUNRO(T2,Q2,PRES) 
C  coefficient for calculation of P in W/m^2                               
        aded=ro2*cp/100.
        keyrecal=0 
        if (kodzim.eq.1.and.ep00.eq.ecmmdec) goto 111
           d0ded=5.1*u2**(2./3.)+0.000001
           IF (LEAF.EQ.0.) THEN
              DT0DED=0.
           ELSE
              DT0DED=1.09*U2**(2./3.)/SQRT(LEAF)
           ENDIF   
           DD=D0DED*LSAI*DT0DED/(LSAI*DT0DED+D0DED)
           DD=DD*(1.-EXP(-0.45*LSAI))+DPDED*EXP(-0.45*LSAI)
        IF (KODZIM.EQ.0) THEN  
c summer
           T0=(RAD-B-EC+3.*T45**4*SIG*SER+ADED*DD*T2)/(4.*T45**3*SER*SIG
     *      +ADED*DD) 
           P=ADED*DD*(T0-T2) 
        else
c winter            
            T0=(RAD-B-EC+3.*T45**4*SIG*SER+ADED*DD*T2*exples+
     !         sig*ser*(4*t45**3*tles-3*t45**4)*exp2*(1-exples)+
     !         cp*ro2*d*(1-exples)*tles)
            t0=t0/(4.*T45**3*SER*SIG+ADED*DD*exples+cp*ro2*d*(1-exples))  
           P=ADED*DD*(T0-T2) 
           p=p*exples           
        endif   
      EFIZL=SIG*SER*(4.*T45**3*T0-3.*T45**4)    
       keyrecal=1
  111 continue   
C the end of the recalculation
c         
C check for heat balance OF THE SURFACE      
        if (kodzim.eq.1) then                      
          BALANCE=rad-efizl+sig*ser*(4*t45**3*tles-3*t45**4)*exp2
     !             *(1-exples)
     !        -B-EC-P+cp*ro2*d*(tles-t0)*(1-exples)
        else
           BALANCE=RAD-EFIZL-B-EC-P
        endif   
c  check for water balance in the whole system  
        dv2=(w2*(h2soil-hroot)-w2pr*(h2psoil-hroot))  
        dv3=0.
        IF (KODZIM.EQ.0) THEN 
           XBAL_all=totalpr+INTERPR+wlespr-ECMMDEC-XDV-INTER2-tothydr-
     !     dv2-dv3-(hyd1-hyd1pr)*1000. 
        ELSE   
           XBAL_all=totalpr+snlespr+wlespr+interpr-(wles+snles)-DELSUM-
     !      DELWLSN-XDV-tothydr-ECMMDEC-eles-trles-dv2-dv3-
     !       (hyd1-hyd1pr)*1000.
        ENDIF   
c***********************************************************************          
c printing descripances in heat and water balances                  
c      if (abs(balance).gt.3.) then
c        WRITE(6,1040) nYEAR,II,I,EC,P,BALANCE,T0,T2,Tles,ET,EP,
c     *           B,EFIZL,RAD,U45,keyiter,WK,WH,ps,HSNOW,KSI_up,KSIOT
c      endif
 1040  FORMAT(I5,I4,I2,3F8.3,3F7.2,5F8.3,F5.2,I3,3F6.3,3f7.1)            
c water balance of the whole system
c       if (abs(xbal_all).gt.0.0001.and.kod_calibr.ne.1) then
c        WRITE(333,1051)  nYEAR,II,I,totalpr,DELWL,DELSN,DELSUM,DELWLSN,
c     !   XDV,tothydr,ECMMDEC,trles,dv2,dv3,(hyd1-hyd1pr)*1000.,XBAL_all
c       endif        
 1051  format(I5,I4,I2,13F12.5)    
C*********************************************************************** 
       ec=ec+eles_wt+trles_wt
       ecmmdec=ecmmdec+(eles+trles)
C calculation of daily values           
        T0_d=T0_d+T0/nstep 
        t2_d=t2_d+T45/nstep 
	u45_d=u45_d+u45/nstep
	q45_d=q45_d+q45/nstep
c gswp2 outputs at 3-hour  time step
      IF (KODPRINT.EQ.1) THEN 
c printing PILPS' output field:      
c ******************    ENERGY BALANCE COMPONENTS  *****************     
       if (kodzim.eq.1) then
          bprint=b                              
          bkey3=0.  
          bkey2=0.
c the situations with snowmelt (3) or freezing water in snow cover (2):
          if (keyiter.eq.3) then
             bkey3=gmlt  
             bprint=b-gmlt 
          else    
             if (keyiter.eq.2) then
                bkey2=gmlt
                bprint=b-gmlt         
             endif
          endif    
c melting and freezing  of snow on trees crowns (ml) and on the ground (bkey3+bkey2):       
         qf=ml+bkey3+bkey2            ! energy of fusion 
         qv=0.
         subsoil=0.
         IF (KODSUBp.EQ.1) THEN
            SUBSOIL=EP
         ENDIF   
         if (kodsubl.eq.1) then
           qv=eles_wt+esnow         ! energy of sublimation
         else
           qv=esnow
         endif     
         qv=qv+subsoil
         SWNET=SW*((1.-albedo)*exples+(1.-albzim)*(1.-exples)) 
         LWNET=RDOWN-EFIZL*exples-sig*ser*(4*t45**3*tles-3*t45**4)*
     !     (1.-exples)
         Psum=P+HL*(1.-EXPLES)
         B=BPRINT+BZV*(1.-EXPLES)+del_les+del_e 
         BALANCE=SWNET+LWNET-B-qf-EC-Psum              
         if (keyrecal.eq.1) then
            b=b+balance
            BALANCE=SWNET+LWNET-B-qf-EC-Psum  
         endif   
       else     
          qf=0.
          qv=0.
          swnet=SW*(1.-albedo)
          Lwnet=RDOWN-EFIZL 
          psum=p
          BALANCE=SWNET+LWNET-B-EC-Psum                                    
       endif
       qle=ec-qv  
c daily averaged values
        sw_d=sw_d+sw/nstep
	  rdown_d=rdown_d+rdown/nstep
	  swnet_d=swnet_d+swnet/nstep
        lwnet_d=lwnet_d+lwnet/nstep
        qle_d=qle_d+qle/nstep
        qh_d=qh_d+psum/nstep
        qg_d=qg_d+b/nstep
        qf_d=qf_d+qf/nstep
        qv_d=qv_d+qv/nstep       
C ******************  SURFACE STATE VARIABLES   ************************  
C  SWE  on the ground (sumsn+wlsn)   
c  SWEVeg on trees crowns (SNLES+WLES)
         if (sumsn.gt.0.) then   
            if (t0.gt.tabs0) t0=tabs0
            snowt=t0
            snfr=sncover
            snalb=albedo   
            sliqfr=wlsn/(sumsn+wlsn)   
         else
            snowt=-9999.
            snfr=0.              
            snalb=-9999.
            sliqfr=0.
         endif       
         
         vegt_d= vegt_d+tles/nstep 
         
         IF (KODZIM.EQ.0) then
             tles=-9999.  
             avgsurft=t0
             radt=t0
         else    
             avgsurft=tles*(1-exples)+t0*exples
             radT=tles**4*(1-exples)+t0**4*exples  
             radt=sqrt(sqrt(radt))
         endif 
        
         if (radt.lt.radtmin) radtmin=radt
         if (radt.gt.radtmax) radtmax=radt
c eb.* output files                   
cc      WRITE(33,444) nyear,jdy,i,SWNET,LWNET,qle,Psum,B,qf,qv,snowt,TLES,
cc     ! avgsurft,radt,nzima,nleto,kodzim,t2daypr,t0daypr
cc     !,lai,lsai,ef,z0veg,sumsn,wlsn  
  444 format (2i5,i3,7F9.3,f11.3,3F9.3,3i3,2f8.3,6f6.2)                           
c DAILY AVERAGED VALUES
        if (snowt.gt.-99.) then
            k_snowt=k_snowt+1
            snowt_d=snowt_d+snowt
        endif
        
         avgt_d=avgt_d+avgsurft/nstep
         radt_d=radt_d+radt/nstep
C VALUES AT THE END OF DAILY TIME STEP:   
         sweveg_d=snles+wles
         if (snles.eq.0.) sweveg_d=0.
         albedo_d=(ALBZIM*(1-exples)+ALBEDO*exples) 
         swe_d=sumsn+wlsn                        
C ******************   WATER BALANCE COMPONENTS  *********************
        if (prsolid.gt.0) then 
          PRSOLID=totalpr 
          prliq=0.
          hprecp_d=hprecp_d+totalpr    ! daily snowfall
        else
          PRLIQ=totalpr
          prsolid=0.
          precip_d=precip_d+totalpr    ! daily rainfall
        endif  
        DELSWE=DELSUM+DELWLSN
        DELINTER=snles-snlespr+wles-wlespr+INTER2-INTERPR
        Canopint=inter2+snles+wles
        delsoilm=XDV+DV2+dv3
        BALANCE=PRSOLID+PRLIQ-ECMMDEC-stokcal-drain2-delsoilm
     !   -DELINTER-DELSWE 
C daily  VALUES:
        ec_d=ec_d+ecmmdec 
        stok_d=stok_d+stokcal
        drain_d=drain_d+drain2
        surrun_transf=surrun_transf+stok
        drain_transf=drain_transf+drhydr
        snmlt_d=snmlt_d+SNOWMELT         
        snowfr_d=snowfr_d+snowfrz
        vd_d=vd_d+vd
        dvsoil_d=dvsoil_d+delsoilm
        dswe_d=dswe_d+delswe
        dint_d=dint_d+delinter
        if (hgr.lt.0.) then
          write(*,*) 'hgr',hgr
c          stop
          hgr=0.
        endif  
        hgr_d=hgr/1000.          ! water table depth (m)    
c        hgr_d=h2soil/1000.          ! water table depth (m)   
c ******************** SUBSURFACE STATE VARIABLES  *******************
       w2zon=w2*(h2soil-hroot)+por2*(h0-h2soil)       ! average layer soil moisture in the 2nd zone (mm)           
        SOILWET=(wk*hroot+w2zon-wzav*h0)/
     !    (por*hroot+por2*(h0-hroot)-wzav*h0)    
C VALUES AT THE END OF DAILY TIME STEP:        
         WK_d=wk*hroot     ! average layer soil moisture in the root zone (mm)  
         w2_d=w2zon
         soilw_d=soilwet
C ********************  EVAPORATION COMPONENTS (mm/time_step or mm) ***********************
        epot_d=epot_d+epot  
        canint_d=canint_d+canopint
        subsoil=0.
        subles=0.
        evapsn_d=0. 
      IF (KODZIM.EQ.1) THEN
        IF (KODSUBp.EQ.1) THEN
            SUBSOIL=EPMMDEC        ! sublimation from soil
            ESOIL=0.               ! soil evaporation
        ELSE
            ESOIL=EPMMDEC
        ENDIF    
        
        if (KODSUBl.eq.1) then
           ecanop=0.                  ! Ecanop (evaporation of canopy interception)
           subles=eles                ! sublimation of canopy interception)
        else   
           ecanop=eles 
        endif    
C daily VALUES (in mm/3hour):          
          et_d=et_d+trles       
          ep_d=ep_d+esoil      
          esnow_d=esnow_d+esnowmm 
          subsl_d=subsl_d+subsoil+subles
          inter_d=inter_d+ecanop  
C VALUES AT THE END OF TIME STEP:  
c evap.* output files  
cc        WRITE(4,479) epot/deltat,Ecanop/DELTAT,TRLES/DELTAT,
cc     !   ESOIL/DELTAT,wk*HROOT,ESNOWMM/DELTAT,SUBSOIL/DELTAT,
cc     !   ECMMDEC/DELTAT,subles/deltat,SOILWET,w2zon,hgr/1000.,
cc     !   (subsoil+subles)/deltat,eles/deltat,ice2,ice1,uu2,u1,nb

cc        WRITE(4,479) wk,w2,hgr/1000.,ECMMDEC/DELTAT,totalpr/deltat,
cc     ! (stokcal+drain2)/deltat,(stok+drhydr)/deltat,radt,ksi_up/100.,
cc     ! ksiot/100.,(sumsn+wlsn)   
     
      ELSE 
           inter_d=inter_d+INTER
           et_d=et_d+ETMMDEC 
           ep_d=ep_d+EPMMDEC
           esnow_d=0. 
           subsl_d=0.  
cc        WRITE(4,479) epot/deltat,INTER/DELTAT,ETMMDEC/DELTAT,
cc     !   EPMMDEC/DELTAT,wk*HROOT,ESNOWMM/DELTAT,SUBSOIL/DELTAT,
cc     !   ECMMDEC/DELTAT,subles/deltat,SOILWET,w2zon,hgr/1000.,
cc     ! (subsoil+subles)/deltat,eles/deltat,ice2,ice1,uu2,u1,nb
cc      WRITE(4,479) wk,w2,hgr/1000.,ECMMDEC/DELTAT,totalpr/deltat,
cc     ! (stokcal+drain2)/deltat,(stok+drhydr)/deltat,radt,ksi_up/100.,
cc     ! ksiot/100.,(sumsn+wlsn)   

      ENDIF  
	       
  479 format (3f6.3,4e14.4,f7.2,2f6.3,f8.2) 
c  479 format (9e14.3,f6.3,F9.3,f10.5,2e14.4,5e12.4)    
C *******************  COLD SEASON PROCESSES *****************************
        snfr_d=snfr
        ksiot_d=ksiot_d+ksiot/100./nstep 
        snow_d=HSNOW/100.     ! snow depth in m
        snalb_d=snalb    
        ksi_d=0.
        if(ksi_up.gt.0.) ksi_d=ksi_up/100.
c wb.* output files    
cc      WRITE(444,478) jdy,i,ECMMDEC/deltat,
cc     !  stokcal/deltat,drain2/deltat,SNOWMELT/deltat,snowfrz/deltat,
cc     ! vd/deltat,delsoilm,DELSWE,DELINTER,BALANCE,
cc     ! prsolid/deltat,prliq/deltat,totalpr/deltat,canopint,xdv,dv2,dv3,
cc     ! etmmdec,epmmdec,inter,esnowmm,eles,trles,(sumsn+wlsn),
cc     ! (snles+wles),(ALBZIM*(1-exples)+ALBEDO*exples),
cc     ! snfr,ksi_up/100.,ksiot/100.,hsnow/100.,snalb,snles,wles,emkl*lsai
cc     ! ,emkh*lsai,hsig0,t_frost,kodzim,w2,kksiot,ksiotpr      
  478   FORMAT(2i4,6e14.3,3f8.2,5e14.3,3e14.3,6f10.3,2f10.3,5f7.3,e14.4
     !  ,2f10.4,3f7.3,e14.4,i3,f10.5,i3,f7.2)          
c------------------------------------------------------------------------ 
c **************************************************************
      endif  
 1098  FORMAT(3I5,i3,F8.3,10F10.3,3f6.2)       
        INTERPR=INTER2
        W0=WK                        
        PSN=PS   
        hyd1pr=hyd1 
         
  200 CONTINUE        
C      ------------------(i)-----the end of time-step calculations-------------------- 
      t0daypr=t0_d 
      t2daypr=t2_d
      ksiotpr=ksiot_d*1000.       ! mm
      if (k_snowt.gt.0) then
         snowt_d=snowt_d/k_snowt
      else 
         snowt_d=-9999.
      endif
      if (kodzim.eq.0) vegt_d=-9999.      
c      
      IF (KODPRINT.EQ.1) THEN 
	if (jdy.lt.jdy_pr) year_print=iyr_1+1                             
c d_010378.* output daily files      
c      write(12,584) year_print,jdy,swnet_d,lwnet_d,qle_d,qh_d,
c     ! qg_d,qf_d,qv_d,snowt_d,vegt_d,avgt_d,
c     ! radt_d,albedo_d,swe_d,sweveg_d,radtmax,radtmin
c     ! ,snfr_d,snalb_d,snow_d,ksi_d,ksiot_d    
  584 format(i5,i4,7f8.2,4f10.2,f5.2,2f8.2,2f7.2,f5.2,4f10.3)  
c d_0204.* output daily files
c      write(15, 585) year_print,jdy,hprecp_d/secd,precip_d/secd,
c     ! ec_d/secd,stok_d/secd,drain_d/secd,snmlt_d/secd,
c     ! snowfr_d/secd,vd_d/secd,dvsoil_d,dswe_d,dint_d,
c     ! WK_d,w2_d,soilw_d
      write(55,587) year_print,jdy,ec_d/secd,precip_d/secd,hprecp_d/secd
ccc     !	 ,surrun_transf/secd,drain_transf/secd,sw_d,rdown_d,t2_d,u45_d
ccc     !      ,q45_d
  
     !,surrun_transf/secd,drain_transf/secd,sw_d,rdown_d,t2_d,u45_d,
     !q45_d,swe_d,ksi_d,ksiot_d,snow_d

  587 format(i5,i4,14e14.3)



  585 format(i5,i4,8e14.3,2f10.3,e14.4,2f8.2,f5.2)
c d0506.* output daily files
c      write(19,586) year_print,jdy,epot_d/secd,inter_d/secd,et_d/secd
c     ! ,ep_d/secd,WK_d,canint_d,evapsn_d,esnow_d/secd,
c     ! subsl_d/secd,hgr_d
  586 format(i5,i4,4e14.3,f10.2,4e14.3,f10.3)
      jdy_pr=jdy
	ENDIF   
         jdy=jdy+1  
      
   30 CONTINUE 
C      ------------------(ii)-----the end of daily calculations------------------------
C      goto 1254
C  EQUILIBRIUM (the 4th year will be printed - the last year)      
      if (kodend.le.3) then
         kodend=kodend+1    
         close(1)  
	   close(2)
	   close(60)                          
c      if (icell.eq.8194.or.icell.eq.8195) then     ! for VALDAI     
c       if (icell.eq.icel_cal) then     ! for VALDAI
        if (ncl.lt.10) then
          WRITE(ICHAR1,'(I1)') ncl
          open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !		ichar1//'.dat')
      if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar1//'.dat')
				if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !  ichar1//'.dat')
        else
          if (ncl.ge.10000) then
            WRITE(ICHAR5,'(I5)') ncl
            open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !		  ichar5//'.dat') 
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar5//'.dat')
    			if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !  ichar5//'.dat')	 	   
          else 
            if (ncl.ge.10.and.ncl.lt.100) then
              WRITE(ICHAR2,'(I2)') ncl
              open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !			ichar2//'.dat')  
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar2//'.dat')
     			if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !  ichar2//'.dat') 		 
            else
               if (ncl.ge.100.and.ncl.lt.1000) then
                WRITE(ICHAR3,'(I3)') ncl
              open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !			ichar3//'.dat') 
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar3//'.dat')
     			if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !  ichar3//'.dat')			        
               else
                 WRITE(ICHAR4,'(I4)') ncl
              open (1,file='c:\d\rivers\forcing_data\'//river//'\F'//
     !			ichar4//'.dat') 
	if(kod_exp.eq.2) open (2,file='c:\d\rivers\forcing_data\p3_'//
     !	ichar4//'.dat')
     			if(kod_prec.eq.2)
     !    open (60,file='c:\d\rivers\forcing_data\'//river//'\P'//
     !  ichar4//'.dat')		     
               endif
            endif
          endif
        endif     
c       endif
       GOTO 100
      endif 
 1254 continue
C__________________________________
        NYEAR=NYEAR+1
C   *********************************************      
      IF (NYEAR.LE.NUMYR) GOTO 100  
1333    CLOSE (1)  
        close(2)
	  close(60)
        CLOSE (33)
        CLOSE (6)
c        close (4)
        CLOSE (12)
        CLOSE (15)
c        close (333)
        close (444)
        close (19)
	close (55)
      if (kod_calibr.eq.2) WRITE(*,*) 'end minus',ncl 
 1879 continue      
      enddo
	close(1)
	close(2)
	close(60)
	close(7)
	close(10)
c      STOP
      return
 1880 END                                   
 
C \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  
C//////////////////////////////////////////////////////////////////////////      
C iterations' block                                         
      SUBROUTINE itblock1(KEYITER,HSNOW,KSI,KSIOT,A3,A4,
     *D00,LEVEL,TABS0,B,EFIZL,P,U2,T2,Q2,T0,G,EP0,
     *zwind,e2,WLSN,kodzim,radles,tles,LSAI,expLES,d,
     !MLES,ELES,WLES,SNLES,DEXP,exp2,z00,ef,ev,trles,LPOD,HPOD
     !,zimw,sryt,T_FROST,HSIG0,BZV,HL,kodsubp,kodsubl,ksiot00,
     ! zerot2,zerot3)
      
      IMPLICIT NONE
                                          
      REAL KSI,KSIOT,LEVEL,LZV,KINV,L,LPL,KAPA,L1,L2,L3,deltat,radles
      REAL B,Q2,WLSN,T0,EP0,T2,G,U2,ZWIND,FF,HSNOW,P,CP,RAD
      REAL TABS0,WH,RO,SIG,D00,RV,SER,FF10,A3,A4,A5,E2,PRES,Q45,EFIZL
      REAL Z00,T45,U45,FUNRO,XX,YY,ZZ,SUMFF,SUMFFU2,SUMFF10,UZV
      REAL X30,X20,X40,X60,X70,X50,X71,X72,TZV,QZV,TETA0,TETA45,Q0
      REAL DEL2,T0ST,ZX,FF0,FFU2,TC2,Q2SAT,DQ2,DE2
      real tles,bles,T00PR,LSAI,RADSUR,EFIZLES,exples,d,tp,bzv,exp2
      REAL BALLES,BALSUR,HL,PL0,snles,wles,mles,eles,DEXP,tlespr,cles
      real z0veg,trk,trles,ef,ev,LPOD,HPOD,zimw,hsig,TZVMN
      REAL HSIG0,T_FROST,a2,par,sryt,lcan,le,zerot2,zerot3,ksiot00
      
      INTEGER ISTEP,I,KODZIM,II,KEYITER,NYEAR,ITER,KEYKOD,K57,NSTEP,LES
      integer kodsubp,kodsubl,jdy
      
      COMMON /METEO/ U45,T45,Q45,RADsur,PRES
c      COMMON /PARAM/ KAPA,SIG,SER,L,LPL,CP,L1,L2,L3,WH,KINV,RV  
      COMMON /PARAM1/KAPA,SIG,SER,Le,LPL,CP,L1,L2,L3,WH,KINV,RV,A2
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy
      FUNRO (XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))
                             
        TZVMN=sryt+tabs0
        A5 = (1000. / PRES) ** .288                      
C Setting up initial (at the first iteration) values of FF 
c cles - teploemkost' drevesiny 2700 J/kg/grad; biomass=2x10^5 kg/ga=20 kg/m^2
        cles=3000.*40.
        cles=cles/deltat
        lcan=le   
        l=le
        kodsubp=0
        kodsubl=0
        if (snles.gt.0.) then
           lcan=Le+LPL*4.19*1000.  
           kodsubl=1
        endif   
         IF((HSNOW.GT.0.).OR.(KSI.GT.0.).or.(KSIOT.EQ.0.)) then   
           L=Le+LPL*4.19*1000. 
           kodsubp=1
        endif
        RO = FUNRO(T45, Q45, PRES)
        T0=T45
        T00PR=T0
        tles=t45 
        tlespr=tles
        tp=t45
        rad=radsur
        eles=0.
        trles=0.
        mles=0. 
        bles=0. 
c FF for T and q at 2 m, FFu2 and ff10 for wind speed at 2 and 10 m       
        FF10 = LOG((zwind - D00) / Z00)
        FF = LOG((LEVEL - D00) / Z00)
        FFu2 = LOG((level - D00) / Z00)   
        if (kodzim.eq.1.and.exples.ne.1) THEN
           les=1                 
           TLES=T0
           KEYKOD=1  
        else
           les=0.
           keykod=1000   ! there is no forest   
           hl=0.
           bzv=0.
        ENDIF   
          
        KEYITER = 0  
 3747   if (kodzim.eq.1) then
          if (les.eq.1) then
            tp=t0
            rad=radles 
            bles=-sig*ser*(4*T45**3*Tp-3*T45**4)     
          else
            TLES=T0
            rad=radsur
            EFIZLES=SIG*SER*exp2*(1-EXPLES)*(4*T45**3*TLES-3*T45**4)
            RAD=RAD+EFIZLES              
          endif 
        endif   
       
        SUMFF10=FF10
        SUMFF=FF
        sumffu2=ffu2
        ITER=1
C calculation of the dynamical speed
 4300   UZV=U45*KAPA/FF10 
        U2=UZV*FFU2/KAPA
        d=0.0076*u2*exp(-DEXP*lsai)/(1+0.82*sqrt(u2))   
        if (les.eq.1) goto 4004
          X30=SIG*SER*T45**4-RAD-KAPA*RO*L*UZV*(Q45-A3)/FF
     !        -cp*ro*d*(1-exples)*(tles-t45)
          X20=KAPA*RO*L*UZV*A4/A5+KAPA*RO*CP*UZV*exples+
     !        4.*SIG*SER*FF*T45**3./A5+cp*ro*d*(1-exples)*ff/a5
          X40=697.37*60.
          X60=T45-FF*X30/A5/X20
          X70=FF/A5/X20     
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      if (ksiot00.lt.5.) then                                  
          if (jdy.gt.zerot2.and.jdy.lt.zerot3) then
            hsig=sqrt(12.*a2*(jdy-zerot2)*deltat+hsig0*hsig0)
          else
             if(jdy.le.zerot2) then 
               hsig=sqrt(12.*a2*(365-zerot3+jdy)*deltat+hsig0*hsig0)
             else
               hsig=sqrt(12.*a2*(jdy-zerot3)*deltat+hsig0*hsig0)  
             endif
          endif     
          keyiter=9 
          X50=HSNOW/L1+hsig/2./L2+HPOD/LPOD  
          X71=X40*X70/X50
          t0=(x60+x71*TZVMN)/(1+x71)
          IF (T0.GE.TABS0) THEN
             T0=TABS0
             KEYITER=3
          ENDIF     
      else         
       if ((kodzim.eq.1).and.(ksiot.eq.0.).and.(t0.lt.tzvmn).and.
     !   hsig0.gt.0.) then  
           hsig=sqrt(12.*a2*t_frost+hsig0*hsig0)
           keyiter=8  
           X50=HSNOW/L1+hsig/2./L2+HPOD/LPOD  
           X71=X40*X70/X50
           t0=(x60+x71*TZVMN)/(1+x71)
           IF (T0.GE.TABS0) THEN
              T0=TABS0
              KEYITER=3
           ENDIF               
       else   
C winter 
        k57=0
        IF ((HSNOW.GT.0.).AND.(KEYITER.EQ.3)) THEN
          T0=TABS0
        ELSE 
          IF (HSNOW.GT.0.) THEN
            IF (WLSN.GT.0.) THEN   
               X50=HSNOW/2./L1 
               KEYITER=2           
            ELSE                               
                X50=(HSNOW+1.)/L1+HPOD/LPOD 
               KEYITER=1          
            endif                                
            X71=X40*X70/X50
            T0=(X60+X71*TABS0)/(1.+X71)        
            IF (T0.GE.TABS0) THEN
              T0=TABS0
              KEYITER=3
            ENDIF                 
          ELSE              
  554       IF ((KSI.GT.1.0).OR.(HPOD.GT.0.)) THEN             
              X50=HSNOW/L1+ksi/L2+HPOD/LPOD 
              X71=X40*X70/X50
              T0=(X60+X71*tabs0)/(1.+X71)
              KEYITER=4            
              IF (T0.GT.TABS0) then
                 if ((ksiot.gt.0.).OR.(HPOD.GT.0.)) then
                     X50=KSIOT/L3+HPOD/LPOD
                     X72=X70*X40/X50
                     T0=(X60+X72*TABS0)/(1.+X72)
                     keyiter=6  
                 ELSE 
                     goto 678              
                 ENDIF    
               ENDIF  
            ELSE                   
c by Budyko            
  678         TZV=SIG*SER*T45**4+B-RAD-KAPA*RO*L*UZV*(Q45-A3)/FF
     !            -cp*ro*d*(1-exples)*(tles-t45)
               TZV=TZV/(KAPA*RO*L*UZV*A4/A5+KAPA*RO*CP*UZV*exples+
     *           4.*SIG*SER*FF*T45**3./A5+cp*ro*d*(1-exples)*ff/a5)
               T0=T45-TZV*FF/A5
               QZV=(Q45-A3)/FF+A4*TZV/A5
               IF (t0.ge.tabs0) then
                 keyiter=7
               ELSE
                 keyiter=5
               ENDIF
               GOTO 4375
            ENDIF
          ENDIF       
        ENDIF                    
      endif
      endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
        TETA0=T0*A5
        TETA45=T45*A5
        TZV=(TETA45-TETA0)/FF
        QZV=(Q45-A3-A4*T0+A4*T45)/FF
 4375     Q0=A3+A4*(T0-T45)
C     calculation of the air density
          RO=FUNRO(T0,Q0,PRES) 
c*********************************************forest***************
 4004 if (les.eq.1) then
       IF (SNles.eq.0.) THEN          
          if (wles.eq.0.) then
              trk=zimw*ef*ev
          else
             trk=1.  
          endif                   
          
          TZV=(1+exp2)*SIG*SER*T45**4+Bles-RAD+cp*ro*d*(t45-tp)+
     !      cles*(t45-tlespr)-trk*kapa*ro*Lcan*UZV*(Q45-A3)/FF
          TZV=TZV/(KAPA*RO*CP*UZV+4.*(1+exp2)*SIG*SER*FF*T45**3./A5+
     !         (cp*ro*d+cles)*ff/a5+trk*KAPA*RO*Lcan*UZV*A4/A5)
          T0=T45-TZV*FF/A5 
          QZV=(Q45-A3)/FF+A4*TZV/A5   
          Q0=A3+A4*(T0-T45)
          RO=FUNRO(T0,Q0,PRES) 
          trles=-trk*kapa*ro*uzv*qzv  
          if ((wles.gt.0.).or.(trles.lt.0.)) then    
             eles=trles
             trles=0.
          endif 
       else
          trk=1.     
         IF (Wles.eq.0.) THEN 
            TZV=(1+exp2)*SIG*SER*T45**4+Bles-RAD+cp*ro*d*(t45-tp)
     !          -trk*kapa*ro*Lcan*UZV*(Q45-A3)/FF+cles*(t45-tlespr)
            TZV=TZV/(KAPA*RO*CP*UZV+4.*(1+exp2)*SIG*SER*FF*T45**3./A5
     !          +(cp*ro*d+cles)*ff/a5+trk*KAPA*RO*Lcan*UZV*A4/A5)
            T0=T45-TZV*FF/A5  
            if (t0.ge.tabs0) goto 222 
         else
  222       T0=TABS0                                         
            tzv=a5*(t45-t0)/ff
         ENDIF 
         QZV=(Q45-A3)/FF+A4*TZV/A5   
         Q0=A3+A4*(T0-T45)  
         RO=FUNRO(T0,Q0,PRES)   
         ELES=-trk*kapa*ro*uzv*qzv            
         IF(T0.EQ.TABS0) THEN  
           Mles=rad-(1+exp2)*sig*ser*t0**4-bles+kapa*ro*uzv*(cp*tzv+
     !          Lcan*qzv*trk)-cp*ro*d*(t0-tp)-cles*(t0-tlespr)
         ELSE
           mles=0.
         ENDIF   
       endif
       HL=-KAPA*CP*RO*UZV*TZV        
      ENDIF                 
c*************************************************************************
C     calculation of the scale LZV
          LZV=(TZV/T45+0.61*QZV)*G*KAPA**2
          LZV=UZV**2/LZV  
      IF (ITER.EQ.1) GOTO 4750
                   
      del2=t0-t0st 
  
      IF(((ABS(DEL2).LT.E2).and.(iter.ge.3)).OR.(ITER.GT.5)) GOTO 4771 
C     calculation of the universal functions 
 4750      zx=zwind
        CALL UNFUN(ZX,D00,Z00,LZV,FF,UZV,KINV,ff0,par)    
           ff10=ff-ff0
          ZX=LEVEL       
        CALL UNFUN(ZX,D00,Z00,LZV,FF,UZV,KINV,ff0,par)
           ffu2=ff-ff0

            t0st=t0
            ITER=ITER+1
            SUMFF=SUMFF+FF 
            sumff10=sumff10+ff10
            sumffu2=sumffu2+ffu2
              FF=SUMFF/ITER 
              ff10=sumff10/iter 
              ffu2=sumffu2/iter           
              GOTO 4300
 4771 if (kodzim.eq.1) then
       IF (LES.EQ.1) THEN
            LES=0  
            GOTO 3747
       ELSE  
           DEL2=Tles-T00PR
       IF (((ABS(DEL2).gT.e2).or.(KEYKOD.LT.3)).and.(keykod.lt.15)) THEN
              t00pr=tles
              LES=1 
              KEYKOD=KEYKOD+1
              GOTO 3747                           
       ENDIF
       ENDIF               
      endif                                         
C calculation of fluxes 
              P=-KAPA*CP*RO*UZV*TZV  
              p=p*exples
              EP0=-KAPA*RO*UZV*QZV  
              EFIZL=SIG*SER*(4.*T45**3*T0-3.*T45**4)
C calculation of air temperature, wind speed and deficit of water vapour pressure at 2m
               T2=T45
               Q2=Q45
               U2=UZV*FFU2/KAPA
               TC2=T2-TABS0  
               Q2SAT=(0.623*6.1/PRES)*EXP(17.1*TC2/(235.+TC2))
C  deficit of specific humidity (kg/kg)
               DQ2=Q2SAT-Q2                 
C  deficit of water vapor pressure (mb)    
              DE2=DQ2*PRES/0.623 
       t0=anint(t0*100.)/100.                                               
              if (keyiter.eq.8.or.keyiter.eq.9) B=(T0-TZVMN)*X40/X50
              IF ((KEYITER.EQ.1).OR.(KEYITER.EQ.4).or.(keyiter.eq.2).
     !           OR.(KEYITER.EQ.6)) B=(T0-TABS0)*X40/X50
              IF (KEYITER.EQ.3) THEN
                 B=RAD-EFIZL-P-L*ep0+cp*ro*d*(1-exples)*(tles-t0)
                 IF (B.LT.0.) then
                    IF (WLSN.GT.0.) then
                      KEYITER=2
                    else  
                      KEYITER=1
                    endif  
                 endif
              END IF  
                 
       if ((kodzim.eq.1).and.(keykod.lt.1000)) then  
         bles=sig*ser*(4*T45**3*T0-3*T45**4) 
         PL0=cp*ro*d*(tles-t0)
         bzv=cles*(tles-tlespr)  
         efizles=sig*ser*(4*T45**3*Tles-3*T45**4) 
         eles=eles*lcan
         trles=trles*lcan
         BALLES=RADLES-(1+exp2)*efiZLES+BLES-HL-PL0-bzv-eles-trles
     !         -mles
         BALSUR=RADSUR-BLES+EFIZLES*exp2*(1-EXPLES)-L*EP0-P-B+
     !          PL0*(1-exples)                                            
           mles=mles*(1.-EXPLES)
           eles=eles*(1.-EXPLES)
           trles=trles*(1.-exples)    
      else 
          balsur=RAD-EFIZL-P-L*ep0-b 
       endif
       ep0=ep0*L
       tlespr=tles 
       
      RETURN
      END                                                            
C //////////////////////////////////////////////////////////////////////////
     
C //////////////////////////////////////////////////////////////////////////                                        
      SUBROUTINE LETO1(PRES,rb,tabs0,wlespr)   
  
      IMPLICIT NONE
      
      REAL M1,NB,LAI,k0,k00,lc,l,lleto,m,LSAI,kw2,por2
      REAL ET,EP,XX,YY,ZZ,RV,FUNRO,T2,U2,CPLETO,WZAV,HROOT,TSTEP
      REAL EFF,w2,PR,WK,PS,BPAR,WN,Q2UP,UMG,H2SOIL,PSN,POR,Q22,A4,EP00
      REAL ET00,PRES,Q2DOWN,FI0,EV,ETDAY,FLAI,FLSAI,WKP,F0,ET0
      REAL EP0,FI00,RO,CP,DP,DWV,QPS,DIF,APS,SBR1,SBR2,SBR,ww,h0
      REAL Q3UP,Q2SUM,Q3DOWN,DRAIN,k02,deltat,pr_int
      real aa,bb,epc0,prmmhour,k0zv,emk,leaf,inter,stok,nhour,prhour,tq
      real y5,hk,dw,x5,x6,sigk0,x7,q,trmn,inter2,kond,hpod,emkl
      real kapa,sig,ser,lpl,rb,l1,l2,l3,wh,kinv,ksiot,w1m,lpod,t0
      real tabs0,tc,ll,hgr,ICE_2,ice2,uu2,unn,rl,ksiotmm,iice
      real u1,ice1,ice_1,w3,wlespr,w2pr                
      real podzstok,ksimm,drain1,del1,del2,nb2,wzav2 
      REAL ZA2,TZV,ZEROT1,QT,S,interpr,epot,tcot,ksiot0
      
      INTEGER I,II,ISTEP,NYEAR,NSTEP,kodhyd,kodhyd2,kodzim,jdy               
   
      COMMON /PARAM1/KAPA,SIG,SER,L,LPL,CPleto,L1,L2,L3,WH,KINV,RV,ZA2  
      
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy
      COMMON /EVAP/ ET,EP,LAI,LSAI      
      common /Cleto/ nhour,prhour,tq,hk,sigk0,trmn,leaf,
     !              ET00,inter,inter2,kond   
       
      common /year1/ NB,k0,HROOT,bpar,wzav,fi0,POR,por2,umg,a4,u2,lc,
     !              h2soil,WN,WK,w2,PS,PSN,STOK,
     !              EP00,PR,T2,Q22,kodhyd,kodhyd2,hpod,emkl,h0,
     !              pr_int,ksiot,w1m,lpod,t0,hgr,uu2,
     !              ice2,rl,u1,ice1,kodzim,podzstok,tstep,drain1,TZV,
     !              ZEROT1,interpr,k02,nb2,wzav2,epot,eff,ev,w2pr,ksiot0 
                              
        FUNRO (XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))
c  the units: mm/time_step, kal/g/grad, g/cm^3 etc.
C  all the components of water balance are in mm/time_step (here, mm/3hour)
c  WK, WN and NB are full (not available) values of soil moisture
C  PS - the thickness of drying layer, cm
c  PSN - the value of PS from previous time-step, cm
C      constants:
         sbr=0.
         DW = POR - NB
         M1=0.45                               
         prmmhour=pr
         kond=0.                            
         tc=t0-tabs0
C interception capacity in mm     
            EMK=emkl*LsAI      
            IF (LEAF.EQ.0.) EMK=0.      
c functions of LAI and LSAI
         IF (LAI.GT.10.) LAI=10.
         IF (LSAI.GT.10.) LSAI=10.
         FLAI=EXP(-M1*LAI)
         FLSAI=EXP(-M1*LSAI)
C calculation of intercepted precipitation
               INTER=INTERpr+PRMMHOUR+wlespr
               IF (INTER.GE.EMK) THEN
                  PRMMHOUR=INTER-EMK
                  INTER=EMK
               ELSE
                  PRMMHOUR=0.
               END IF
               pr_int=pr-prmmhour
C calculation of runoff in mm/hour
               STOK = 0.
               NHOUR=NHOUR+DELTAT/3600.          
C translation from mm/time_step to mm/hour
               PRMMHOUR=PRMMHOUR/tstep
               IF (PRMMHOUR.GT.0.) THEN
                  IF ((NHOUR-PRHOUR).GT.24.) TQ=-DELTAT/3600./2.
                  TQ=TQ+DELTAT/3600. 
                  Y5=SQRT(TQ)  
                  K0ZV=SQRT(HK*DW/TQ/8.+PRMMHOUR)-SQRT(2.*HK*DW)/4./Y5
                  K0ZV=K0ZV*K0ZV
                  X5=SQRT(3.)
                  X6=K0-X5*SIGK0
                  IF (X6.LT.0.) X6=0.
                  IF (K0ZV.LE.X6) THEN
                     STOK=0.
                  ELSE
                     X7=K0+X5*SIGK0
                     BB=1./2./X5/SIGK0
                     AA=SQRT(2.*HK*DW)*BB/2./Y5
                     IF ((K0ZV.GT.X6).AND.(K0ZV.LE.X7)) THEN
                        Q=2.*AA/3.*(K0ZV**(3./2.)-X6**(3./2.))+BB/2.*(k0
     *                    ZV**2-X6**2)+PRMMHOUR*BB*(X7-K0ZV)
                     ELSE
                        Q=2.*AA/3.*(X7**(3./2.)-X6**(3./2.))+BB/2.*(X7**
     *                    2- X6**2)
                     END IF
                     STOK=PRMMHOUR-Q  
                     if(stok.lt.0.) then
C                        write(*,*) 'STOK<0',ii,i,stok
C                        stop
                         STOK=0.
                     endif   
                  END IF                 
                  PRHOUR=NHOUR
               END IF
C the end of stok calculation           
c  translation of STOK and PRMMHOUR from mm/hour into mm/time-step
           STOK=STOK*tstep
           PRMMHOUR=PRMMHOUR*tstep  
c calculation of precipitation at the land surface            
           PR=PRMMHOUR-STOK  
C calculation of potential evaporation from the canopy
           epc0=ep00*(1.-FLSAI)  
C calculation of potential transpiration
           IF (EP00.EQ.0.) THEN
              ET00=0.
           ELSE
              ET00=EP00*TRMN
           END IF
C calculation of evaporation (condensation) of intercepted precipitation
      IF (EMK.EQ.0.) THEN 
         INTER=0.
         INTER2=0.
         IF (EP00.LT.0.) THEN
           KOND=EP00
           PR=PR-EP00
           EP00=0.
         ENDIF 
      else 
        IF (EPc0.GT.INTER) THEN
           EPc0=EPc0-INTER
           ET00=EPc0/(1.-FLSAI)*TRMN
        ELSE
           INTER2=INTER-EPc0
              INTER=EPc0
              ET00=0.
              if (ep00.lt.0.) then
                 kond=ep00*FLSAI
                 pr=pr-kond
                 ep00=0.
              endif
        END IF          
      endif           
c potential transpiration should be in mm/day for the calculation of WKP
         etday=et00*nstep
C the calculation of ET0 and EP0
         WKP=(6.+0.42*etday)/100.
         F0=EV*EFF*(WN-wzav)/WKP
         IF (F0.LT.0.) F0=0.
         IF (F0.GT.1.) F0=1.
         ET0=ET00*(1.-FLAI)
         EP0=EP00*FLSAI        
         Epot=inter+ET0+EP0  
c CALCULATION OF PARAMETERS FOR DRYING LAYER         
c  k0 -saturated hydraulic conductivity (mm/hour)
C      translation of K0 from mm/hour into mm/time_step
         k00=k0*tstep                 
c  fi0 - matric potential at saturation (m)        
c      translation of FI0 from (m) into (mm) and negative value
         fi00=fi0*(-1000.)
C  air density         
         RO=FUNRO(T2,Q22,PRES)
         RO=RO/1000.
c  L   - [J/kg]        
c      translation of L from (J/kg) into (kal/g)        
         Lleto=L*0.2387/1000.
c  CP  - [J/kg/K]                               
c      translation of CP from (J/kg/K) into (kal/g/grad)       
         CP=CPleto*0.2387/1000.
c  [DP] in (cm/c) ; 
         DP=0.76*U2/(0.82*U2**0.5+1.) 
c    coefficient of diffusivity of water vapour [DWV] in (cm^2/c) 
         dwv=0.66*(por-umg)*(0.24+0.42*u2)   
c water influx to the lower boundary of drying layer         
         QPS=2.*fi00*bpar*k00*(WN/por)**(bpar+3.)/hroot/(bpar+3.)
c calculation of heat flux from unfrozen soil to permafrost &&&&&&&&&&&&&&&&&&&&  new &&&&&&&&&&&&&&     
         
        IF (ISTEP.GT.ZEROT1*NSTEP) THEN
          QT=9./4.*KSIOT0**2+12.*ZA2*(DELTAT*(iSTEP-ZEROT1*nstep)) 
c       if (qt.lt.0.) write(*,*) '777',ii,jdy,i,qt,ksiot,za2,istep,zerot1 
          QT=2.*L2*TZV/(SQRT(QT)-3./2.*KSIOT0)     
          LL=LPL*RB*(W1m-WH)
          IF (ksiot0.ge.100.) LL=lpl*rb*((w2+W1M)/2.-WH)  
          IF (LL.LE.0.) LL=lpl*rb*0.1   
c          if (ll.eq.0.) ll=0.000001
          S=QT*DELTAT/LL 
       ELSE
          S=0.
          LL=Lpl*RB*(W1m-WH)    
          IF (ksiot0.ge.100.) LL=lpl*rb*((w2+w1m)/2.-WH)          ! ********************   
          IF (LL.LE.0.) LL=lpl*rb*0.1     
c           if (ll.eq.0.) ll=0.000001
       ENDIF                          
C       WRITE(222,444) JDY,I,ISTEP,ZEROT1,QT,LL,KSIOT,S,ZA2
C  444  FORMAT(3I6,F6.1,5E14.4)    
C calculation of soil thawing depth  KSIOT   
c            SUMTC=SUMTC+TC   
            tcot=tc
            if (tcot.lt.0.) tcot=0.
c            IF (SUMTC.LE.0.) SUMTC=0.  
c        if(ii.eq.74.and.i.eq.2) write(*,*) 'l3',l3,sumtc,ll,w1m,wh    
              KSIOT=S-L3/LPOD*HPOD+SQRT((ksiot0+L3/LPOD*HPOD)**2+
     !                2.*L3*tcot*DELTAT/LL+s*s)  
c        write(222,222) jdy,i,s,qt,ll,istep,zerot1
c  222  format(2i5,3e14.4,i5,f6.1)              
        if (ksiot.gt.300.) ksiot=300.
c        if(ii.eq.74.and.i.eq.2) write(*,*) 'ksiot',ksiot   
c calculation of ICE
         unn=wh
         ICE_2=(W2-UNN)*RB/RL  
         ice_1=(wn-unn)*rb/rl    
         if (ice_2.lt.0.) ice_2=0.
         if (ice_1.lt.0.) ice_1=0.
c calculation of ICE2
c ksiotmm [mm], ksiot [cm]         
         ksiotmm=ksiot*10.
         if (ksiotmm.lt.h2soil) then    
            if (ksiotmm.lt.hroot) then
              ice1=ice_1*(hroot-ksiotmm)/hroot
              ice2=ice_2
            else
              ice1=0.
              ice2=ice_2*(h2soil-ksiotmm)/(h2soil-hroot)
            endif 
         else
            ice1=0.
            ice2=0.
         endif
        if (ice1.lt.0.) ice1=0.
         if (ice2.lt.0.) ice2=0.
        uu2=w2-ice2*RL/RB   
        u1=wn-ice1*rl/rb       
c        if(ii.eq.74.and.i.eq.2) write(*,*) 'wn',wn,uu2,hroot,h2soil   
C           the calculation of Q2UP                            
             ww=(wn*hroot+UU2*(h2soil-hroot))/h2soil
             iice=ice2*(h2soil-hroot)/h2soil
        if (ww.lt.0) write(*,*) ii,i,wn,w2,h2soil,hroot,ww           
              kw2=k00*(ww/por)**(2.*bpar+3.)/(1.+8.*iice)**2
              dif=kw2*fi00/por**(-bpar)*(-bpar)*ww**(-bpar-1)/(1.+
     !            8.*iice)
              q2up=dif*(uu2-wn)/(h2soil/2.)
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   end new &&&&&&&&&&&&
         IF (Q2UP.LT.0.) Q2UP=0.
         QPS=QPS-q2up
         APS=(LLETO*A4*RO*CP/(LC+RO*CP*DWV)+CP/DWV)*DP/(CP+LLETO*A4)
c QPS in mm/time_step, PSN in cm, APS in 1/cm
        M=1./(1.+APS*PSN)         
C the calculation of ET and EP (mm/time_step) 
         ET=F0*ET0 
         EP=M*EP0 
         IF (WN.LT.UMG) EP=0.
c calculation of PS (cm) 
         if (wn.eq.umg) then
            PS=PSN+0.1*(EP+QPS)/(WN-UMG+0.0000001)  
         else   
            PS=PSN+0.1*(EP+QPS)/(WN-UMG)  
         endif
           
            if (ps.lt.psn) ps=psn
            
         if (wn.eq.umg) then
            PS=PS+0.1*(-pr-q2up)/(WN-UMG+0.0000001)
         else   
            PS=PS+0.1*(-pr-q2up)/(WN-UMG) 
         endif    
         
         IF (PS.LT.0.) PS=0 
         if (ps.gt.(hroot/10.)) then
            ps=hroot/10.
            if (EP.GT.Q2up) EP=Q2up
         endif 
c calculation of WK (in the 1st zone)          
         WK=WN+(PR-ET-EP+Q2up)/HROOT 
         if (wk.lt.0.) then
            et=0.
            WK=WN+(PR-ET-EP+Q2up)/HROOT 
            if (wk.lt.0.) then
              ep=0.
              WK=WN+(PR-ET-EP+Q2up)/HROOT 
            endif    
         endif
c calculation of SBR from the 1st zone      
       IF ((Por-ICE1).GT.NB) THEN
           if ((Wk-ICE1*RL/RB).GT.NB) then
              if ((Wk-ICE1*RL/RB).GT.(Por-ICE1)) then 
                 SBR1=Wk-ICE1*RL/RB-Por+ICE1
                 W3=Por-ICE1+ICE1*RL/RB
                 Wk=ICE1*RL/RB+NB+(Por-ICE1-NB)*EXP(-4.*K0*((nb-wzav)/
     !              (por-wzav))**3/hroot/(por-wzav)/(1.+8.*ICE1)**2)
                 SBR2=W3-Wk
                 SBR=SBR1+SBR2       
              ELSE                  
                 W3=Wk
               Wk=ICE1*RL/RB+NB+(W3-ICE1*RL/RB-NB)*EXP(-4.*K0*((nb-wzav)
     !              /(por-wzav))**3/hroot/(por-wzav)/(1.+8.*ICE1)**2)
                 SBR=W3-Wk
              ENDIF   
           ELSE     
                 SBR=0.
           ENDIF      
        ELSE 
           IF ((ICE1+(Wk-ICE1*RL/RB)).GT.Por) THEN
                 W3=ICE1*RL/RB+(Por-ICE1)
                 SBR=Wk-W3
                 Wk=W3
           ELSE
                 SBR=0.
           ENDIF         
        END IF      
          
       q2down=sbr*hroot
       q2sum=q2down-q2up
       DRAIN1=Q2SUM 
       q3up=0.
       
         w2=w2pr+(q2sum+q3up)/(h2soil-hroot) 
         if (w2.lt.0.) then
            q2sum=0.
            w2=w2pr+(q2sum+q3up)/(h2soil-hroot)   
         endif           
         
C calculation of SBR from the 2nd zone
c calculation of drainage out of the 2nd zone (with taking into account
c change in the water table  
         ksimm=0.
         DEL1=0.
         DEL2=0. 
         sbr=0.
      if (ksimm.gt.h2soil) then
        sbr=0.
        IF ((w2-ICE2*RL/RB).gt.(por2-ice2)) then
           DEL2=(w2-(por2-ice2+ice2*RL/RB))*(H2SOIL-HROOT)
           sbr=del2                                                   
           w2=por2-ice2+ICE2*RL/RB 
        else                                                                   
           if ((w2-ICE2*RL/RB).GT.NB2) then                                    
              if ((w2-ICE2*RL/RB).GT.(Por2-ICE2)) then                        
                 SBR1=w2-ICE2*RL/RB-Por2+ICE2                                  
                 W3=Por2-ICE2+ICE2*RL/RB                                      
            w2=ICE2*RL/RB+NB2+(Por2-ICE2-NB2)*EXP(-4.*K02*((nb2-wzav2)/       
     !              (por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)
                 SBR2=W3-w2                                                     
                 SBR=SBR1+SBR2                                                  
              ELSE                                                              
                 W3=w2                                                          
                w2=ICE2*RL/RB+NB2+(W3-ICE2*RL/RB-NB2)*EXP(-4.*K02*((nb2-          
     !       wzav2)/(por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)        
                 SBR=W3-w2                                                      
              ENDIF                                                             
           ELSE                                                                 
             SBR=0.                                                          
           ENDIF                                                                 
        ENDIF                                                                    
      else   
        IF ((Por2-ICE2).GT.NB2) THEN
           if ((w2-ICE2*RL/RB).GT.NB2) then
              if ((w2-ICE2*RL/RB).GT.(Por2-ICE2)) then 
                 SBR1=w2-ICE2*RL/RB-Por2+ICE2
                 W3=Por2-ICE2+ICE2*RL/RB
             w2=ICE2*RL/RB+NB2+(Por2-ICE2-NB2)*EXP(-4.*K02*((nb2-wzav2)/
     !              (por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)
                 SBR2=W3-w2
                 SBR=SBR1+SBR2       
              ELSE                  
                 W3=w2
                w2=ICE2*RL/RB+NB2+(W3-ICE2*RL/RB-NB2)*EXP(-4.*K02*((nb2-
     !       wzav2)/(por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)
                 SBR=W3-w2
              ENDIF
           ELSE     
                 SBR=0.
           ENDIF      
        ELSE 
           IF ((ICE2+(w2-ICE2*RL/RB)).GT.Por2) THEN
                 W3=ICE2*RL/RB+(Por2-ICE2)
                 SBR=w2-W3
                 w2=W3
           ELSE
                 SBR=0.
           ENDIF         
        END IF 
      endif  
c       if(ii.eq.74.and.i.eq.2) write(*,*) 'sbr',sbr   
       q3down=sbr*(h2soil-hroot)
       drain=q3down-q3up
       if (ice1.lt.0.) ice1=0.
       if (ice2.lt.0.) ice2=0.
        u1=wk-ice1*rl/rb
        uu2=w2-ice2*RL/RB             

c        
        podzstok=drain
        if (wk.lt.0.) then
          write(*,*) 'wk<0',wk
          stop
        endif  
  771 FORMAT(3I4,f7.2,i3,11F6.3,3F7.1)  
       
      RETURN
      END
C////////////////////////////////////////////////////////////////////////////
      SUBROUTINE ZIMA1(EP,PRES,TMLL,lsai,DELWL,DELSN,tabs0)     
    
      IMPLICIT NONE
                                                            
      REAL L1,L2,L3,L,LZV,LL,NB,KSI0,KSIOT,KFILTR,K0,ICE_1,LLET
      REAL ICE_2,tstep
      REAL LSUM,lc,m,ksimm,ksiotmm,ice2,ice1,KW2,ww,iice,w3,por2
      REAL STOK,Q2,CPZIM,WLSN,T0DAYPR,ZEROT1,U2,RSPR,W0,WZAV
      REAL H,DELSUM,ESNOW,HROOT,P,QVPIT,RB,W1M,U,EP,DELWLSN
      REAL SUMTP,TABS0,VD,BZ,TCCC,RL,WH,PR,WK,PS,BPAR,QT,RV,US,UMG,h0
      REAL H2SOIL,PSN,TBYD,GMLT,A3,A4,EP00,C2,SUM,SNOWMELT,eles0
      REAL TMLL,PRES,H2,Q2DOWN,DELTAT,FI0,FUNRO,XX,YY,ZZ,k02
      REAL T0DAY,TC,TML,TZV,TBYDC,T0,FI00,RO,CP,DP,DWV,APS,BZIM,GWLSN
      REAL WLSNPR,HT,HSUM0,ESOIL,QBZIM,TCC,RS,HH,DV1,DV,ESN,WLSN0,TP1
      REAL S,S1,UNN,XXUU,U1,DIF,Q2UP,QPS,SBR1,SBR2
      REAL SBR,Q3UP,Q2SUM,W2,UU2,Q3DOWN,DRAIN,DLW,DRAIN1,DEL2,DEL1
      REAL KAPA,SIG,SER,KINV,MLES,ELES,WLES,SNLES,ML,lsai,ESNW,emkh
      real rsles,rslespr,delwl,delwl2,delwl3,delsn,delsn2,EMK,R77,HPOD
      REAL PRECIP,trles,LPOD,sncover,kondliq,tga,b0,sumqbz,sumq,emkl                                
      real pr_int,prh_int,prsolid,KSI_UP,T_FROST,hsig0
      real hgr,podzstok,del_les
      real snlespr,nb2,wzav2,interpr,eleskond,snowfrz,tcot
      real del_e,zimw,ef,ev,ZA2,ZC3,epot,w2pr,ksiot00,ksiot0

      INTEGER I,II,ISTEP,NYEAR,NSTEP,IM,KEYITER,kodhyd,kodhyd2
      integer marfir,kodzim,jdy                                       
               
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy 
      COMMON /PARAM1/KAPA,SIG,SER,LLET,L,cpzim,L1,L2,L3,WH,KINV,RV,ZA2
   
      common /year1/ NB,KFILTR,HROOT,bpar,wzav,fi0,P,por2,umg,a4,u2,lc,
     !              h2soil,W0,WK,w2,PS,PSN,STOK,
     !              EP00,PRECIP,TBYD,Q2,kodhyd,kodhyd2,hpod,emkl,h0
     !             ,pr_int,ksiot,w1m,lpod,tccc,hgr,uu2,
     !              ice2,rl,u1,ice1,kodzim,podzstok,tstep,drain1,TZV,
     !               ZEROT1,interpr,k02,nb2,wzav2,epot,ef,ev,w2pr,ksiot0
      
      common /Czima1/ A3,C2,RB,US,RS,IM,
     !              t0daypr,SUMTP,H2,
     !              WLSN,DELWLSN,gmlt,SUM,DELSUM,BZ,VD,
     !              SNOWMELT,ESNOW,KEYITER,
     !              MLES,ELES,WLES,SNLES,rsles,WLSNPR,trles,
     !              sncover,tga,b0,marfir,sumqbz,emkh,prh_int,prsolid,  
     !              KSI_UP,T_FROST,hsig0,ML,del_les,snlespr,
     !              ZC3,snowfrz,del_e,zimw,ksiot00

      FUNRO (XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))
C I - step; II- day 
c WLSN  - volume of liquid water in snow cover (MM)
C W0    - total initial moisture in rooting zone
c w2   -total soil moisture in the 2nd zone
C translation of the surface temperature TC, temperature at the GCM's
C lowest level TML, temperature at 2m TBYDC and constant temperature
C of deep soil TZV from K to C
                  
        epot=eles+trles+ep00  
        eleskond=0.  
        qvpit=0.              
        snowfrz=0.
        PR=PRECIP
        rspr=rs
        rslespr=rsles    
        TC=TCCC-TABS0   
c        TCC=TC
c        tcot=tc 
        IF (TC.GT.0.) then
           TCC=0.
           tcot=tc
        else 
           tcc=tc
           tcot=0.
        endif      
       rsles=rslespr*(1.+0.1/10./nstep*(snles+wles)*EXP(.08*TCc-
     !        21.*RSlespr)) 
        ESN=0.11/RSles-0.11   
        ESNW=ESN/(1-ESN)
       emk=emkh*lsai  
        t0day=t0daypr-tabs0   
        
        TML=TMLL-TABS0
        TBYDC=TBYD-TABS0   
        if (jdy.gt.200) sncover=1.
c PSN - dry layer depth at the beginning of the time step (i.e. from previous time step),cm        
C constants and parameters
C K0 - coef. of filtration in mm/time_step;       
        K0=KFILTR*24./nstep 
c  fi0 - matric potential at saturation (m)        
c      translation of FI0 from (m) into (mm) and negative value
        fi00=fi0*(-1000.)        
C LSUM - in kal/g
        LSUM=L+LLET/4.19/1000.
c air density in kg/m^3
        RO=FUNRO(Tbyd,Q2,PRES)
        ro=ro/1000.
c [CP] - [J/kg/K] ; translation of CP from (J/kg/K) into (kal/g/grad)       
        CP=CPzim*0.2387/1000.
c [DP] - [cm/c] ; 
        DP=0.76*U2/(0.82*U2**0.5+1.) 
c coefficient of diffusivity of water vapour [DWV] in (cm^2/c) 
        dwv=0.66*(p-umg)*(0.24+0.42*u2)   
c water influx to the lower boundary of drying layer         
        APS=(Lsum*A4*RO*CP/(LC+RO*CP*DWV)+CP/DWV)*DP/(CP+Lsum*A4)
c PSN in cm, APS in 1/cm
        M=1./(1.+APS*PSN)             
c ****************************************************************
C calculation of the components of water balance of snow (mm/time_step)
C  H - snow pack, VD - 'vodootdacha', DELSUM - change in snow pack
C  ESOIL - soil evaporation (after snow cover loss)  
C  QBZIM - for soil thawing out or heating 
C translation of BZIM from W/m^2 to kal/cm^2/s
        BZIM=BZ/697.37/60.        
        ML=MLES/697.37/60. 
C translation of kal/cm^2/s to mm/time_step
        gwlsn=0.
        BZIM=10.*BZIM/RB/L*DELTAT  
        ML=10.*ML/RB/L*DELTAT 
        qbzim=0.
        IF (KEYITER.EQ.2) then
           qbzim=bzim*(1.-sncover)
           GWLSN=BZIM*sncover               
        endif   
        IF (KEYITER.NE.3) then
           BZIM=0. 
        else
           qbzim=bzim*(1.-sncover)
           bzim=bzim*sncover  
        endif                           
c WLSNPR - volume of liquid water in snow cover at the previous time step 
C recalculation of precipitation PR into hard precipitation HT
          IF (TML.LT.0.) THEN
             HT=PR 
             pr=0. 
          ELSE
             HT=0.
          END IF
          prsolid=ht
c calculation of snow interception             
          delwl=0.
          delwl2=0.
          delwl3=0.
          delsn=0.
          delsn2=0.
          del_e=0. 
        if (emk.eq.0.) then   
           pr=pr+interpr
           goto 876           
        endif
        if((eles.lt.0.).and.(snles.gt.0.)) then      ! condition of sublimation
             eleskond=eles
             eles=0.
             snles=snles-eleskond
          endif                         
c interception of solid and liquid precipitation by a forest
         snles=snles+ht
         if (eles.lt.0.) then
             wles=wles+pr+interpr-eles 
             eles0=0.
             eleskond=eles
             eles=0.
          else
             wles=wles+pr+interpr  
             eles0=eles
          endif   
          pr_int=0.  
          prh_int=0.
c comparison with interception capacity and removing excess of solid and liquid fractions (delsn and delwl) 
       if(ht.gt.0.) then   
          if(snles.gt.emk) then 
             delsn=snles-emk
             snles=emk
             delwl=delsn*wles/(emk+delsn)
             wles=wles-delwl
             prh_int=ht-delsn
          else
             prh_int=ht
          endif  
       else                       
          R77=esnw*SNLES+emkl*lsai
          if (wles.gt.R77) then
             delwl=wles-R77
             wles=R77
             pr_int=pr-delwl
             if (pr_int.lt.0.) pr_int=0.
          else
             pr_int=pr
          endif                             
       endif
c evaporation, melting and freezing of intercepted precipitation              
       if (snles.gt.0.) then
          r77=snles/(snles+wles)
          if (snles.lt.(eles*r77)) then
             eles=snles+wles
             snles=0.
             wles=0.              
          else 
             snles=snles-eles*r77
             wles=wles-eles*(1-r77)
             if(ml.gt.0.) then
               if(ml.gt.snles) ml=snles
             else
               if(abs(mL).gt.wles) ml=-wles  
             endif  
             snles=snles-ml
             wles=wles+ml
          endif      
c comparison with interception capacity         
          R77=ESNW*SNLES+emkl*lsai
          if (wles.gt.R77) then
             delwl2=wles-R77
             wles=R77
          endif      
          if (snles.gt.emk) then
             delsn2=snles-emk
             snles=emk
             delwl3=delsn2*wles/(emk+delsn2)
             wles=wles-delwl3
          endif
          delsn=delsn+delsn2
          delwl=delwl+delwl2+delwl3      
       else
          if (wles.gt.0.) then
             if (wles.ge.eles) then
                wles=wles-eles
             else   
                eles=wles
                wles=0.
             endif    
          endif   
       endif                             
c       
       if (wles.lt.0.) wles=0. 
       if (snles.lt.0.) snles=0.
       
       if ((eles.gt.0.).and.(eles.lt.eles0)) then
          trles=(eles0-eles)*zimw*ef*ev
          del_e=eles0-eles-trles
       endif   

          ht=delsn
          pr=delwl
  876  continue    
c calculation of snow on the ground           
          HSUM0=SUM  
          ESNOW=EP00
          esoil=0.
          vd=0.
          kondliq=0. 
          
          IF (EP00.LT.0.) then
            sum=sum+ht
            if (sum.eq.0.) THEN
               esoil=ep00
               esnow=0.
            else 
               if (tc.ne.0.) then
                  SUM=SUM-EP00
               else
                  kondliq=ep00
               endif   
                  if (sum.gt.bzim) then
                     sum=sum-bzim
                     vd=bzim
                  else
                     vd=sum
                     sum=0.
                  endif
            endif   
          else  
            IF ((SUM+HT).LT.(BZIM+EP00)) THEN
               IF ((SUM+HT).GE.EP00) THEN
                  VD=SUM+HT-EP00
                  SUM=0.
               ELSE
                  VD=0.
                  ESOIL=EP00-SUM-HT
                  SUM=0.    
               ENDIF
            ELSE
               SUM=SUM+HT-BZIM-EP00
               VD=BZIM
            ENDIF 
            ESNOW=ESNOW-ESOIL            
          endif

         QBZIM=qbzim+(BZIM-VD) 
         sumqbz=sumqbz+qbzim
         bzim=vd

c now, BZIM is real energy consumed to snowmelt (mm/time_step)
         H=SUM
         DELSUM=SUM-HSUM0
         SNOWMELT=VD   
         VD=VD+PR-kondliq       
         if (esoil.lt.0.) then 
            vd=vd-esoil
         endif                          
C the end of the calculation 
C
C translation of mm/time_step into kal/cm^2/time_step (soil thawing)          
         sumQ=0.1*sumQBZ*L 
         
         KSI0=KSI_UP 
         if (ksi0.le.0) ksi0=0.0000001
C calculation of snow density
C RSPR -snow density at previous time-step
        RS=RSpr*(1.+0.1/10./nstep*(H+wlsn)*EXP(.2*TCc-21.*0.75*RSpr))    
c        write(*,*) rs
c        stop
        if (rs.gt.0.7) rs=0.7
c       RS=RSpr*(1.+0.1/10./nstep*(H+wlsn)*EXP(.08*TCc-21.*RSpr)) 
C calculation of snow heat conductivity  
c formula by Bracht (Kuchment, 1983, p. 135)          
          L1=0.0049*RS**2
C calculation of snow height (cm) using snow pack (mm)
         H2=(H+wlsnpr)/10./RS 
         if (h.eq.0.) h2=0.
C calculation of 'privedennoi' snow height
         HH=H2*L2/L1+hpod*L2/LPOD 
C calculation of change in volume of liquid water in snow cover
       DV1=0.
       IF ((TC.LT.0.).AND.(WLSNpr.GT.0.)) THEN
c freezing of water in 'wet' snow cover
c here, snowfreezing is pozitive
          DV=GWLSN         
          WLSN=WLSNPR+DV  
          IF (WLSN.LT.0.) THEN 
             DV1=WLSN
             WLSN=0. 
             SUM=SUM+WLSNPR
             snowfrz=wlsnpr
          ELSE
             SUM=SUM-DV
              snowfrz=-dv
          ENDIF
          gwlsn=gwlsn-dv1
       ELSE   
          ESN=0.11/RS-0.11  
          WLSN0=RS*H2*ESN*10./RB                
          WLSN=WLSNpr+VD
          IF (WLSN.LE.WLSN0) THEN
             VD=0.       
          ELSE
             VD=WLSN-WLSN0
             WLSN=WLSN0
          ENDIF  
          GWLSN=0.
       ENDIF       
C calculation of soil surface temperature   
         IF (KSI0.EQ.0.) KSI0=0.00000001
         TP1=KSI0*t0day/(KSI0+HH) 
        
         IF (TP1.GT.0.) TP1=0.
         IF (TP1.LE.-0.1) THEN 
            SUMTP=0.
         ELSE
            SUMTP=SUMTP+DELTAT
         END IF
         
      if (ksiot00.lt.5) then
          hsig0=50.
          ksiot=0.
          ksi_up=0.
      else                   
C heat flux from  UNfrOZEN ZONE  TO PERMAFROST
        S=0.
        IF (ISTEP.GT.ZEROT1*NSTEP) THEN
           QT=9./4.*KSIOT0**2+12.*ZA2*(DELTAT*(iSTEP-ZEROT1*nstep)) 
           QT=2.*L2*TZV/(SQRT(QT)-3./2.*KSIOT0) 
           LZV=L*RB*(W1m-WH)
           IF (ksiot0.ge.100) LZV=l*rb*((w2+w1m)/2.-WH) 
           IF (LZV.LE.0.) LZV=L*RB*0.1  
           S=QT*DELTAT/LZV  
        ENDIF  
     
        LL=L*RB*(W1m-WH)
        IF (ksiot0.ge.100) LL=l*rb*((w2+w1m)/2.-WH) 
         IF (LL.LE.0.) LL=L*RB*0.1 
C calculation of soil thawing depth  KSIOT
c         IF (H2.GT.0.001) THEN    
c             IF (TC.GT.0.) SUMTC=SUMTC-TC
c             if (sumtc.GT.0.) then            
c              KSIOT=s-L3/LPOD*HPOD+SQRT((ksiot0*ksiot0+L3/LPOD*HPOD)**2+
c     !                 2.*L3*TCot*DELTAT/LL+s*s)  
c             else
c                ksiot=0.
c             endif      
c          ELSE
c             IF (SUMTC.LE.0.) SUMTC=0.
c             if (s.lt.0.and.sumtc.eq.0.) s=0. 
            KSIOT=S-L3/LPOD*HPOD+SQRT((ksiot0+HPOD*L3/LPOD)**2+
     !              2.*L3*TCot*DELTAT/LL+s*s)   
c          END IF    
         if (ksiot.gt.300.) ksiot=300. 
                               
        if ((t_frost.gt.0.).and.(t_frost.lt.(deltat*24/tstep*60.))) 
     !      ksiot=0.        
c calculation of soil freezing depth KSI_up
        QT=0.
          if (ksiot.gt.0.) then
            LZV=L*RB*(W1m-WH)+C2*ABS(TP1)/2.
            IF (LZV.EQ.0.) LZV=1E-11
            S=QT*DELTAT/LZV
            S1=(KSI0+HH)**2-(2.*L2*TCC*DELTAT/LZV)*sncover+S**2
            
            if ((tc.gt.0.).and.(hh.eq.0.)) 
     !         S1=(KSI0+HH)**2-(2.*L2*TC*DELTAT/LZV)*sncover+S**2
            IF (S1.GT.0.) GOTO 60
               KSI_UP=0.
               GOTO 70
   60       KSI_UP=-HH-S+SQRT(S1)
   70       CONTINUE 
            KSI_UP=KSI_UP-DV1/(W1m-WH)*0.1
            IF (KSI_UP.LE.0.) KSI_UP=0.
          endif  
       IF ((KSIOT.GT.0.).AND.(KSI_UP.ge.ksiot)) THEN
            HSIG0=KSI_UP 
c            T_FROST=T_FROST+DELTAT  
            KSIOT=0.
            KSI_UP=0. 
c            sumtc=0.
        ENDIF      
       if (ksi_up.gt.ksiot) ksi_up=0.
       if (ksiot.lt.0.) ksiot=0. 
      endif 
C calculation of unfrozen water UNN and ICE in frozen zone
           unn=wh
           ICE_1=(W0-UNN)*RB/RL 
           ICE_2=(W2-UNN)*RB/RL 
           if (ice_1.lt.0.) ice_1=0.
           if (ice_2.lt.0.) ice_2=0. 
c calculation of liquid water U1 in root zone and UU2 in the 2nd zone
        ksimm=ksi_up*10.
        ksiotmm=ksiot*10. 
c***********************************
        if (ksimm.eq.0.) then
           if(ksiotmm.eq.0.) then
             ice1=ice_1
             ice2=ice_2            
           ELSE
             IF (KSIOTMM.GT.H2SOIL) THEN
                ICE1=0.
                ICE2=0.
             ELSE     
                if (ksiotmm.gt.hroot) then
                   ice1=0.
                   ice2=ice_2*(h2soil-ksiotmm)/(h2soil-hroot)
                else
                   ice1=ice_1*(hroot-ksiotmm)/hroot          
                   ice2=ice_2
                endif  
             endif
           ENDIF  
        else    
           IF (KSIMM.GT.H2SOIL) THEN
              ICE1=ICE_1
              ICE2=ICE_2
           ELSE   
             if (ksimm.gt.hroot) then  
               ice1=ice_1   
               if (ksiotmm.gt.h2soil) then
                 ice2=ice_2*(ksimm-hroot)/(h2soil-hroot)
               else
                 ice2=ice_2*(h2soil-ksiotmm+ksimm-hroot)/(h2soil-hroot)
               endif  
             else
               if (ksiotmm.lt.hroot) then
                 ice1=ice_1*(hroot-(ksiotmm-ksimm))/hroot
                 ice2=ice_2
               else
                 ice1=ice_1*ksimm/hroot
                 ice2=ice_2*(h2soil-ksiotmm)/(h2soil-hroot)
               endif
             endif          
           endif
        endif 
        if (ice1.lt.0.) ice1=0.
        if (ice2.lt.0.) ice2=0. 
        if (ice1.gt.(p-wh)) ice1=p-wh 
        if (ice2.gt.(por2-wh)) ice2=por2-wh 
        u1=w0-ice1*RL/RB 
         uu2=w2-ice2*RL/RB 
C calculation of infiltration QVPIT and runoff STOK 
C the calculation of liquid water into wet zone of frozen soil(by Ye.M.Gusev)
           IF(VD.LT.0.) then
              WRITE(*,*) II,I,ICE1,VD,DELSN,DELWL,DELWL2,DELWL3
              write(*,*) 'VD<0'
              stop
           endif
             XXUU=VD*(1.+8.*ICE1)**2/K0  
             XXUU=SQRT(XXUU)
             XXUU=SQRT(XXUU)
             U=US+(P-US)*XXUU
             IF (U.LT.(P-ICE1)) THEN
                QVPIT=VD
                STOK=0.
             ELSE
                QVPIT=K0*((P-ICE1-US)/(P-US))**4/(1.+8.*ICE1)**2
                STOK=VD-QVPIT
                U=P-ICE1
             ENDIF

        ww=(u1*hroot+uu2*(h2soil-hroot))/h2soil
        iice=(ice1*hroot+ice2*(h2soil-hroot))/h2soil
        kw2=k0*(ww/p)**(2.*bpar+3.)/(1.+8.*iice)**2
        dif=kw2*fi00/p**(-bpar)*(-bpar)*ww**(-bpar-1)/(1.+
     !      8.*iice)
        q2up=dif*(uu2-u1)/(h2soil/2.)           
        if (q2up.lt.0.) q2up=0.
C the calculation of  EP (mm/time_step) 
ccc by Budagovsky:
         EP=M*ESOIL
         if (ep.lt.0.) ep=0.     
         IF (W0.LT.UMG) EP=0.
c calculation of PS (cm)
            PS=PSN+0.1*EP/(W0-UMG)  
            if (ps.lt.psn) ps=psn
            PS=PS+0.1*(-vd-q2up)/(W0-UMG) 
         IF (PS.LT.0.) ps=0. 
         if (ps.gt.(hroot/10.)) then
            ps=hroot/10.
            if (EP.GT.Q2up) EP=Q2up
         endif
         qps=0.      
c calculation of WK in the 1st (root) zone         
         Wk=W0+(QVPIT-EP-trles+q2up)/HROOT 
         if (wk.lt.0.) then
            trles=0.
            Wk=W0+(QVPIT-EP-trles+q2up)/HROOT 
            if (wk.lt.0.) then
              ep=0.
              Wk=W0+(QVPIT-EP-trles+q2up)/HROOT 
            endif    
         endif
C calculation of  SBR1 from the root zone
        IF ((P-ICE1).GT.NB) THEN
           if ((Wk-ICE1*RL/RB).GT.NB) then
              if ((Wk-ICE1*RL/RB).GT.(P-ICE1)) then 
                 SBR1=Wk-ICE1*RL/RB-P+ICE1
                 W3=P-ICE1+ICE1*RL/RB
                 Wk=ICE1*RL/RB+NB+(P-ICE1-NB)*EXP(-4.*K0*((nb-wzav)/
     !              (p-wzav))**3/hroot/(p-wzav)/(1.+8.*ICE1)**2)
                 SBR2=W3-Wk
                 SBR=SBR1+SBR2       
              ELSE                  
                 W3=Wk
               Wk=ICE1*RL/RB+NB+(W3-ICE1*RL/RB-NB)*EXP(-4.*K0*((nb-wzav)
     !              /(p-wzav))**3/hroot/(p-wzav)/(1.+8.*ICE1)**2)
                 SBR=W3-Wk
              ENDIF   
           ELSE     
                 SBR=0.
           ENDIF      
        ELSE 
           IF ((ICE1+(Wk-ICE1*RL/RB)).GT.P) THEN
                 W3=ICE1*RL/RB+(P-ICE1) 
c                 write(*,*) '555',rl,rb,ice1,p,wh
c                 stop
                 SBR=Wk-W3
                 Wk=W3
           ELSE
                 SBR=0.
           ENDIF         
        END IF
C CALCULATION OF RETANTION OF WATER IN THE FIRST ZONE BEFORE ITS UNFREEZING
C the calculation of Q3UP, Q2DOWN, Q2SUM 
        q2down=sbr*hroot
        q2sum=q2down-q2up
        DRAIN1=Q2SUM
        q3up=0.
c the calculation of UU2 and ICE2 in the 2nd layer 
         w2=w2pr+(q2sum+q3up)/(h2soil-hroot)  
         if (w2.lt.0.) then
            q2sum=0.
            w2=w2pr+(q2sum+q3up)/(h2soil-hroot)   
         endif   
C calculation of SBR from the 2nd zone  
      DEL1=0.
      DEL2=0.
      sbr=0.
      if (ksimm.gt.h2soil) then
        sbr=0.
        IF ((w2-ICE2*RL/RB).gt.(por2-ice2)) then
           DEL2=(w2-(por2-ice2+ice2*RL/RB))*(H2SOIL-HROOT)
              sbr=del2
           w2=por2-ice2+ICE2*RL/RB 
        else                                                                 
           if ((w2-ICE2*RL/RB).GT.NB2) then                                
              if ((w2-ICE2*RL/RB).GT.(Por2-ICE2)) then                         
                 SBR1=w2-ICE2*RL/RB-Por2+ICE2                                 
                 W3=Por2-ICE2+ICE2*RL/RB                                       
             w2=ICE2*RL/RB+NB2+(Por2-ICE2-NB2)*EXP(-4.*K02*((nb2-wzav2)/       
     !              (por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)
                 SBR2=W3-w2                                                    
                 SBR=SBR1+SBR2                                                 
              ELSE                                                             
                 W3=w2                                                          
               w2=ICE2*RL/RB+NB2+(W3-ICE2*RL/RB-NB2)*EXP(-4.*K02*((nb2-          
     !       wzav2)/(por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)      
                 SBR=W3-w2                                                     
              ENDIF                                                            
           ELSE                                                                
                 SBR=0.                                                        
           ENDIF                                                                
        ENDIF
      else   
        IF ((Por2-ICE2).GT.NB2) THEN
           if ((w2-ICE2*RL/RB).GT.NB2) then
              if ((w2-ICE2*RL/RB).GT.(Por2-ICE2)) then 
                 SBR1=w2-ICE2*RL/RB-Por2+ICE2
                 W3=Por2-ICE2+ICE2*RL/RB
             w2=ICE2*RL/RB+NB2+(Por2-ICE2-NB2)*EXP(-4.*K02*((nb2-wzav2)/
     !              (por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)
                 SBR2=W3-w2
                 SBR=SBR1+SBR2       
              ELSE                  
                 W3=w2
                w2=ICE2*RL/RB+NB2+(W3-ICE2*RL/RB-NB2)*EXP(-4.*K02*((nb2-
     !       wzav2)/(por2-wzav2))**3/hroot/(por2-wzav2)/(1.+8.*ICE2)**2)
                 SBR=W3-w2
              ENDIF   
           ELSE     
                 SBR=0.
           ENDIF      
        ELSE 
           IF ((ICE2+(w2-ICE2*RL/RB)).GT.Por2) THEN
                 W3=ICE2*RL/RB+(Por2-ICE2)
                 SBR=w2-W3
                 w2=W3
           ELSE
                 SBR=0.
           ENDIF         
        END IF 
      endif  
c calculation of drainage out of the 2nd zone (with taking into account
c change in the water table        
       q3down=sbr*(h2soil-hroot)
       drain=q3down-q3up 
C  ADDITION OF LATERAL SBROS
      drain=drain
c      ***********************************************	         
       if (ice1.lt.0.) ice1=0.
       if (ice2.lt.0.) ice2=0.
       u1=wk-ice1*RL/RB 
       uu2=w2-ice2*RL/RB            
c       
        podzstok=drain
       DELSUM=SUM-HSUM0
       DELWLSN=WLSN-WLSNPR
       H=SUM                                 
C***************************       
c calculation of snow coverage after 1 March    
        if ((jdy.ge.marfir).and.(jdy.lt.250)) then  
           if (h.ge.b0/2.) then 
             sncover=1.  
           else  
            sncover=sqrt(2*h/tga)  
          endif
        endif   
C****************************
c GMLT in mm/time_step
      gmlt=bzim
      if (keyiter.eq.2) gmlt=gwlsn  
c translation of mm/time_step into W/m**2
      gmlt=gmlt/10./deltat*697.37*60.*rb*L
      ml=ml/10./deltat*697.37*60.*rb*L 
      del_les=mles-ml
      
      IF (Wk.Le.0.) THEN
        DLW=(0.-Wk)*hroot
        Wk=0.  
         IF (EP.GT.DLW) THEN
            EP=EP-DLW
         ELSE       
            EP=0.       
         ENDIF   
         DLW=0.
      ELSE
         DLW=0.
      ENDIF                
      IF (ESOIL.LT.0.) EP=ESOIL  
       eles=eles+eleskond
      RETURN

      END       
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    