C *******************************************************************
C          Physically based model for seasonally frozen soils
c	Subroutine model_plus (npar,ncell_bas,ncl_riv,river)          ! ***
	Subroutine model_plus (i0,npar,ncell_bas,river)          ! ***
c The program is made from gswpnew6.for 
c changes: 
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
	character*6 ichar6 
c monthly fields                 
      REAL LAI_M(12),grnfr_m(12),z0vegm(12),D00M(12),EFTR_M(12)
      REAL albedo_m(12)
      REAL DBMES(12),DBMES_N(12),DBMES_S(12)
c daily fields
      REAL LAID(366),grnfrd(366),z0vegd(366),D00D(366),eftrd(366)
      REAL SDBMES(366),albedod(366)
     
c 1.CONSTANTS FOR EVERYTHING
C   a) prescribed parameters and other values
      REAL CP,RV,KAPA,SIG,L,G,LPL,KINV,ZCB,ZRB,cic,RL
      REAL albsnow,RS0,SER,ZCRT,ev,emkl,DEXP,exp2,HPOD,LPOD
      REAL LEVEL0,zwind0,zerot1,TSTEP,sncover
      REAL FUNRO   
      integer MON(12),ierr1,ierr2,ierr3
      integer iyrfir,numyr,firday,TAYL,ncell
c   b) calculated at the beginning of the program     
      REAL DELTAT
      INTEGER NSTEP,mmm
c   c) at the beginning of the program after reading parameters
     
c 2.FIXED FIELDS needed at each time_step (should be massives (15238))
c   a) given by organizers of GSWP
      REAL h0,nb,wzav,por,bpar,K0,fi0,hroot,lat,elevat  
c   b) derived by ourselves previously     
      REAL albzim,prles,emkh,sai,hveg,leaf,sryt,amt,xi,yi,ksiot0
c   c) calculated at the beginning of the program after reading the above parameters
      REAL vdhgr,trmn,zrt,wh,uss,umg,lc,sigk0,hk,POR2,level,zwind  
c 3.FORCINGS (3-hour values) should be read at every time_step
      REAL SW,RDOWN,rainf,snowf,T45,U45,PRES,Q45,sw_d,rdown_d
	REAL u45_d,q45_d 
C 4.VARIABLES NEEDED FOR THE NEXT time_STEP for each point
C  a) current values (different for the points, so should be massives (15238),initialization before the annual cycle)
      REAL KSI,KSIOT,INTERPR,nhour,TQ,PRHOUR,W0,W1M,w2,UZIM,ICE,u1
      REAL ice1,sumqbz
      REAL PSN,HSNOW,WLSN,SUMSN,WLES,SNLES,RSPR,rslespr,SUMTP,SUMTC
      REAL HYD1,HYD2,HYD11,HYD22,hyd1pr,sumhyd,hsr,H2SOIL,H2PSOIL
      REAL t0_d,t2_d,t2daypr,T0DAYPR,L1
      REAL hhyd0(24)    
      INTEGER NLETO,NZIMA,K80,KODZERO,KODZIM,k_snowt
c  b) daily values for GSWP2 output (different for the points, so should be massives, initialization before  srok-cycle)   
      REAL swnet_d,lwnet_d,qle_d,qh_d,qg_d,qf_d,qv_d,hprecp_d,precip_d
      REAL ec_d,stok_d,drain_d,snmlt_d,snowfr_d,vd_d,dvsoil_d,dswe_d
      REAL dint_d,snowt_d,vegt_d,avgt_d,radt_d,epot_d,inter_d,et_d
      REAL canint_d,evapsn_d,radtmin,radtmax,esnow_d,subsl_d,ep_d      
C  c) the same for all points (not massives)
      INTEGER IMDAY_1,NYEAR,JDY,i0
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
	REAL GMLT,A2Z,MLES,SOILWET
      REAL DRAIN1,STOKCAL,DELSUM,ET,EP,ESNOW,EC,ECMMDEC,XDV
      REAL DPDED,RO2,D0DED,DT0DED,DD,X40,A1Z,BPRINT,BKEY3,ksipr,rspr0
      REAL ep00pr,TOTHYDR,DRHYDR,sumtcpr,wlsnpr,tles,RADLES,exples
      REAL ELES,trles,sumsnpr,subles,wlespr,rsles0,SNLESPR,DELWL,DELSN
      REAL X10,X30,elespr,BALANCE,XBAL,pr_int,prh_int,sumqbzpr
      REAL ecanop,prsolid,q3hgr,xbal_all,ml,bkey2,qf,qv,del_e,eles_wt
      REAL trles_wt,sumtppr,prliq,snowfrz,HL,BZV,DELSWE,DELINTER
      REAL snowt,avgsurft,radt,snfr,snalb,sliqfr,qle,swnet,lwnet,epot
      REAL SUBSOIL,ESOIL,psum,del_les,canopint,dv2,dv3,bbb
      REAL D00,ALBLET,mps,surrun_transf,drain_transf, kalblet 
      integer ISTEP,kodsubl,KODSUBp,KEYITER,KODHYD,kodhyd2,keyrecal
      integer prhour0,year_print,jdy_pr
      
c      INTEGER ncell_bas,ncl,ncl_riv(ncell_bas),land,ico,jco,ibas,icel_cal
      INTEGER ncell_bas,ncl,land,ibas
      REAL lon,opt_par(100),k0_opt,hroot_opt,krain,ksnow,h0_opt
	real k_lw,k_sw,d0
c 6. TIME AND CELL etc.CYCLES: I_CYCLE (DAILY), II_CYCLE (SROK), ICELL_CYCLE (POINTS)     
      INTEGER icell,iyr_1,I,II,KODEND,KODPRINT,IMDAY,IYR,IM,npar,kod_exp
	INTEGER kod_prec,kod_calibr,kod_cal,nmon,ixm,ia3,ia2,ia33,i0_cal
c      REAL, ALLOCATABLE :: run_day(:,:)
c	integer, ALLOCATABLE :: n1(:),n2(:),n3(:)
c	integer, ALLOCATABLE :: lat_cell(:),lon_cell(:)

	real run_day(1700,12000),pr_day(1700,12000),ec_day(1700,12000)
c      Real, ALLOCATABLE ::  run_day(:,:),ec_day(:,:),pr_day(:,:)
	real clay,sand

	integer n1(1700),n2(1700),n3(1700)
c	integer lat_cell(5000),lon_cell(5000)
	integer ind1(1700),ico(1700),jco(1700)
	real area(1700),sum_cl_s,k_ratio


       equivalence (POR,POR2)
C      
      COMMON /METEO/ U45,T45,Q45,RAD,PRES
      COMMON /PARAM/KAPA,SIG,SER,L,LPL,CP,L1,L2,L3,WH,KINV,RV
      COMMON /EVAP/ ETMMDEC,EPMMDEC,LAI,LSAI
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy  
      common /year/ NB,k0,HROOT,bpar,wzav,fi0,POR,umg,a4,u2,lc,
     !              h2soil,h2psoil,W0,WK,w2,PS,PSN,STOK,drain2,
     !              EP00,PR,T2,Q2,kodhyd,kodhyd2,hpod,emkl,h0,
     !              pr_int,q3hgr,drain1,ef,ev,epot,interpr,mps,
     !              w2pr
      common /Cleto/ nhour,prhour,tq,hk,sigk0,trmn,leaf,
     !              ET00,inter,inter2,kond
      common /Czima/ ZA3,ZC2,ZRB,RL,USS,RSpr,SRYT,ZEROT1,IM,
     !              T0,t0daypr,SUMTC,SUMTP,W1M,KSI,KSIOT,HSNOW,
     !              WLSN,DELWLSN,gmlt,SUMSN,DELSUM,B,VD,
     !              SNOWMELT,uzim,ice,ESNOWMM,KEYITER,kodzero,
     !              MLES,ELES,WLES,SNLES,rslespr,WLSNPR,trles,LPOD,
     !             sncover,sumqbz,emkh,prh_int,prsolid,
     !             ml,del_les,del_e,zimw,kodprint,snowfrz
    
   
      common /hydro/ xi,yi,lat
	common /opt_par/ opt_par,kod_exp, kod_prec,albsnow,iyrfir,numyr,
     !	firday,kod_calibr,kod_cal
c	common /rout/ lat_cell,lon_cell,n1,n2,n3,run_day 
      common /rout/ ind1,area,n1,n2,n3,run_day,ico,jco,ec_day,pr_day  

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
      DATA CP,RV,KAPA,SIG,L,G,LPL,KINV,ZCB,ZRB,cic,RL/1005.,287.05,0.43,
     ! 0.000000057,2501000.,9.81,80.,0.000013,1.,1.,0.5,0.9/ 
        FUNRO(XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
c      ALLOCATE (run_day(ncell,m0),ec_day(ncell,m0),pr_day(ncell,m0),
c     !	  stat=ierr2)

      tstep=24.0           ! time_step (3 hour)
      nstep=24/tstep      ! number of time steps in a day
      deltat=tstep*3600.  ! time_step in sec. 
C  PILPS' parameters  which are the same for all grid cells
      LEVEL0=2.          ! the reference height for atmospheric forcings (air temperature and air humidity) (m)
      zwind0=10.         ! the reference height of wind speed measurements (m) 
c      albsnow=0.75       ! albedo of fresh snow cover
      RS0=0.15           ! snowpack density for fresh snow (g/cm^3)
      SER=1.             ! coefficient of "serosti"
      ZCRT=0.22          ! (kal/g/grad) specific heat capacity of dry soil
      zerot1=290.        ! transition from positive to negative air temperatures at autumn (number of day)
      ev=1               ! parameter for the calculation of transpiration
      HPOD=0.            ! thickness of podstilki
      LPOD=0.00015       ! TEPLOPROVODNOST' PODSTILKI
      emkl=0.1           ! coefficient for the calculation of interception capacity for liquid precipitations
      TAYL=3             ! for the calculation of water table depth  (in days)
      DEXP=0.2           ! parameter for heat exchange in the forest
      exp2=1.            ! parameter for the forest     
ccc      iyrfir=1983        ! the first full year of the calculations (for the calculations started in summer)
ccc      numyr=14           ! number of years for the calculation
c      numyr=2
ccc      firday=182         ! the first day for the calculations (if they are started in summer)
      mmm=TAYL*nstep
c	mmm=TAYL*8
      sncover=1. 
	riv=river 

c	  lat_cell=0
c	  lon_cell=0
       ico=0
	jco=0
	ind1=0
	area=0.
c
c      write(*,*) 'numyr',numyr
c	                                 open (111,file='d:\ISI-MIP2\res\ttt')

c open input files with parameters which are different for each grid cell
      OPEN (7, file='d:\ISI-MIP2\fixed_fields\fixed_param_cl_'//riv//
     !	'.txt')    
      open (10,file='d:\ISI-MIP2\monthly_fields\monthly_param_'//
     !	riv//'.txt') 
c *************************** the beginning of the cycle by cells*******************************
      ncell=ncell_bas   ! number of grid cells  within the river basin

      do icell=1,ncell 
       
c input model parameters which are different for the cells
      read (7,*) ncl,ico(icell),jco(icell),lon,lat,ind1(icell),
     !	area(icell),sand,clay,
     !   amt,sryt,albzim,prles,sai,leaf,emkh,xi,yi,hroot,h0,hveg,elevat
     !   ,ksiot0,d0
      
cc	if (kod_cal.ne.1) then
cc	   sum_cl_s=(clay+sand)*opt_par(12)
cc	    k_ratio=sand/clay*opt_par(3)
cc	    clay=sum_cl_s/(1+k_ratio)
cccc	    clay=sum_cl_s/(1+opt_par(3))
cc	    sand=sum_cl_s-clay
cc	endif
      bpar=0.157*clay*100.-0.003*sand*100.+3.1                         ! Cosby
	fi0=-0.0095*sand*100.+0.0063*(100.-sand*100.-clay*100.)+1.54     ! Cosby
      fi0=-(10**fi0/100.)
	por=(-0.037*clay*100.-0.142*sand*100.+50.5)/100.                    ! Cosby
	k0=-0.0064*clay*100.+0.0126*sand*100.-0.6   ! log(k) Cosby
      k0=10**k0/60./60.                             ! äþéì/ñ
      k0=k0*2.54*0.01                             ! (m/s)
	
	if (kod_cal.ne.1) then
	    k0_opt=opt_par(3)
	    k0=k0*exp(k0_opt)
      endif

	nb=por*(0.1/1000./24./60./60./K0)**(1/(2*bpar+3))      ! fc_Cosby(under k=0.1mm/day)
! 0.1 mm/day=0.1*10^-3/24/60/60 m/s  
 	wzav=por*(-150/fi0)**(-1/bpar)   ! wp_Cosby(fi=-150m)
 
c      write(*,*) bpar,fi0,por,k0,nb,wzav
c	write(*,*) 'par',(opt_par(i),i=1,12)
c	pause

	read (10,*) land,(albedo_m(i),i=1,12),
     !            (grnfr_m(i),i=1,12),(Z0VEGM(i),i=1,12),
     !            (LAI_M(I),I=1,12),(EFTR_M(i),i=1,12) 

	do i=1,12
	   D00M(I)=d0
	   if (kod_cal.ne.1) then
	     albedo_m(i)=albedo_m(i)*opt_par(6)
c	     LAI_M(I)=LAI_M(I)*opt_par(7)
         endif
      enddo



      if(sryt.lt.-1.) then
         GOTO  1879     
      else  

c********************************************************************************     
c       if (por-nb.le..01) nb=por-.01
c	 if(h0.lt.hroot) h0=hroot+0.1
C		calibration of albzim, albsnow, albedo
       krain=1
	 ksnow=1
	 k_sw=1
	 k_lw=1								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if (kod_cal.ne.1) 	then
		
c		albsnow=opt_par(5)
c		albzim=opt_par(4)*albzim
c	    kalblet=opt_par(6)
c	    albedo_m=kalblet*albedo_m
c	    nb=opt_par(3)*nb
c    				if (por-nb.le..01) nb=por-.01
c          k0_opt=opt_par(1)
c	    k0=k0*exp(k0_opt)
c	    hroot_opt=opt_par(2)
          hroot=hroot*opt_par(1)
   	    h0=hroot*opt_par(2)
c                  if(h0.lt.hroot) h0=hroot+0.1
          krain=opt_par(8)
	    ksnow=opt_par(9)
	    k_sw=opt_par(10)
      	k_lw=opt_par(11)
c       
c         write(*,*) (opt_par(i),i=1,15)
c		Write(*,*) albsnow,albzim,albedo_m,nb,k0,hroot,h0
c	stop
        endif
      endif	

c*****************************************************************************************     
      		
c open input file with forcings  for i-cell AND OUTPUT FILES  
c open output files for  i-cell

c	if (land.gt.9999) then
c	  write(ICHAR5,'(I5)') land
c        OPEN (55, FILE='d:\ISI-MIP2\res\point '//ichar5//'.txt')
c      else
c	  write(ICHAR4,'(I4)') land
c        OPEN (55, FILE='d:\ISI-MIP2\res\point 0'//ichar4//'.txt')
c      endif 
c      Write(55,*) 'Year	Day	Evap	Rainf	Snowf	Qs_transf	Qsb_tr
c     !ansf   SW    RDOWN    T2     U2    Q2'
cc       lat_cell(icell)=lat
cc	 lon_cell(icell)=lon
      
	 

       if (ncl.lt.10) then
          WRITE(ICHAR1,'(I1)') ncl
          open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !		ichar1//'.dat')
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar1//'.dat')
			if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !   ichar1//'.dat')
c            OPEN (4, FILE='d:\ISI-MIP2\res\evap.'//ichar1)  
c            OPEN (33, FILE='d:\ISI-MIP2\res\Eb.'//ichar1)
c            OPEN (444, FILE='d:\ISI-MIP2\res\WB.'//ichar1)
c            OPEN (12, FILE='d:\ISI-MIP2\res\d_o1o378.'//ichar1)
c            OPEN (15, FILE='d:\ISI-MIP2\res\d_o2o4.'//ichar1)
c            OPEN (19, FILE='d:\ISI-MIP2\res\d_o5o6.'//ichar1) 
c            OPEN (333, FILE='d:\ISI-MIP2\res\balan_w.'//ichar1) 
c            OPEN (6, FILE='d:\ISI-MIP2\res\balan_t.'//ichar1) 
c            OPEN (55, FILE='d:\ISI-MIP2\res\point.'//ichar1)
         else
         if (ncl.ge.100000) then
           WRITE(ICHAR6,'(I6)') ncl
           open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !		 ichar6//'.dat')
	    else
         if (ncl.ge.10000) then
           WRITE(ICHAR5,'(I5)') ncl
           open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !		 ichar5//'.dat')
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar5//'.dat')
			if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !  ichar5//'.dat')
c             OPEN (4, FILE='d:\ISI-MIP2\res\evap.'//ichar5)  
c             OPEN (33, FILE='d:\ISI-MIP2\res\Eb.'//ichar5)
c             OPEN (444, FILE='d:\ISI-MIP2\res\WB.'//ichar5)
c             OPEN (12, FILE='d:\ISI-MIP2\res\d_o1o378.'//ichar5)
c             OPEN (15, FILE='d:\ISI-MIP2\res\d_o2o4.'//ichar5)
c             OPEN (19, FILE='d:\ISI-MIP2\res\d_o5o6.'//ichar5) 
c             OPEN (333, FILE='d:\ISI-MIP2\res\balan_w.'//ichar5) 
c             OPEN (6, FILE='d:\ISI-MIP2\res\balan_t.'//ichar5) 
c		   OPEN (55, FILE='d:\ISI-MIP2\res\point.'//ichar5)  
         else 
           if (ncl.ge.10.and.ncl.lt.100) then
              WRITE(ICHAR2,'(I2)') ncl
              open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !			ichar2//'.dat')   
      if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar2//'.dat')
			if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !  ichar2//'.dat')
c                OPEN (4, FILE='d:\ISI-MIP2\res\evap.'//ichar2)  
c                OPEN (33, FILE='d:\ISI-MIP2\res\Eb.'//ichar2)
c                OPEN (444, FILE='d:\ISI-MIP2\res\WB.'//ichar2)
c                OPEN (12, FILE='d:\ISI-MIP2\res\d_o1o378.'//ichar2)
c                OPEN (15, FILE='d:\ISI-MIP2\res\d_o2o4.'//ichar2)
c                OPEN (19, FILE='d:\ISI-MIP2\res\d_o5o6.'//ichar2) 
c                OPEN (333, FILE='d:\ISI-MIP2\res\balan_w.'//ichar2) 
c                OPEN (6, FILE='d:\ISI-MIP2\res\balan_t.'//ichar2)
c			  OPEN (55, FILE='d:\ISI-MIP2\res\point.'//ichar2) 
           else
              if (ncl.ge.100.and.ncl.lt.1000) then
                WRITE(ICHAR3,'(I3)') ncl
              open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !			ichar3//'.dat')  
      if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar3//'.dat')
     		if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !  ichar3//'.dat')			       
c                  OPEN (4, FILE='d:\ISI-MIP2\res\evap.'//ichar3)  
c                  OPEN (33, FILE='d:\ISI-MIP2\res\Eb.'//ichar3)
c                  OPEN (444, FILE='d:\ISI-MIP2\res\WB.'//ichar3)
c                  OPEN (12, FILE='d:\ISI-MIP2\res\d_o1o378.'//ichar3)
c                  OPEN (15, FILE='d:\ISI-MIP2\res\d_o2o4.'//ichar3)
c                  OPEN (19, FILE='d:\ISI-MIP2\res\d_o5o6.'//ichar3) 
c                  OPEN (333, FILE='d:\ISI-MIP2\res\balan_w.'//ichar3) 
c                  OPEN (6, FILE='d:\ISI-MIP2\res\balan_t.'//ichar3) 
c               	OPEN (55, FILE='d:\ISI-MIP2\res\point.'//ichar3)
              else
                 WRITE(ICHAR4,'(I4)') ncl
              open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !			ichar4//'.dat') 
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar4//'.dat')	
     		if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !  ichar4//'.dat')	     
c                  OPEN (4, FILE='d:\ISI-MIP2\res\evap.'//ichar4)  
c                  OPEN (33, FILE='d:\ISI-MIP2\res\Eb.'//ichar4)
c                  OPEN (444, FILE='d:\ISI-MIP2\res\WB.'//ichar4)
c                  OPEN (12, FILE='d:\ISI-MIP2\res\d_o1o378.'//ichar4)
c                  OPEN (15, FILE='d:\ISI-MIP2\res\d_o2o4.'//ichar4)
c                  OPEN (19, FILE='d:\ISI-MIP2\res\d_o5o6.'//ichar4) 
c                  OPEN (333, FILE='d:\ISI-MIP2\res\balan_w.'//ichar4) 
c                  OPEN (6, FILE='d:\ISI-MIP2\res\balan_t.'//ichar4) 
c              	OPEN (55, FILE='d:\ISI-MIP2\res\point.'//ichar4)
              endif
           endif
         endif
       endif 
	endif
    
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
c     !edo  SnFr    ksi   ksiot  SnDep   SAlbedo'
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
c Z0VEGm - roughness coefficient of vegetation (m); 
c LAI_M = total_LAI  
c grnfr_m - greenness fraction (green_LAI/total_LAI 
C d00 - zero displacement height  (m)
c albedo_m - snow_free albedo
c hveg - the height of vegetation (m)
C AMT - amlitude of air temperature (in C),
c SRYT - soil temperature at the depth where seasonal variations are absent (in C);
C LEAF -size of grass's leaves in (cm)
c WH - min. amount of unfrozen water 
c 
c calculation of zerot1 (in days from 1 January)
        if (lat.ge.0.) then
           zerot1=373.+1.914*lat-0.05682*lat*lat-0.03234*elevat
           if (zerot1.gt.365.) zerot1=365.
           if (zerot1.lt.0.) zerot1=0.
        else
           zerot1=662.+14.06*lat+0.102*lat*lat-0.0604*elevat+
     !            0.0000015*elevat*elevat   
           if (zerot1.lt.0.) zerot1=0.
           if (zerot1.gt.firday) zerot1=firday
        endif
            
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

c transformation of input parameters             
        IF ((POR-NB).LT.0.1) NB=POR-0.1
        if (sryt.lt.0) sryt=0. ! for the first approximation for the territories with permafrost (really there should be another program)
        hroot=hroot*1000.      ! translation from m to mm
        h0=h0*1000.            ! translation from m to mm
c CALCULATION OF the rest MODEL PARAMETERS:
        level=level0+hveg      ! m
        zwind=zwind0+hveg      ! m
        vdhgr=por-nb           ! water yield of groundwater
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
        fi0=fi0*(-1.)
C  calculation by Kalyuzhnyi and Pavlova
c        WH=1.04*NB-0.06*ZRT    
        wh=(nb+wzav)/2.
        USS=WZAV/1.35
        UMG=USS 
        lc=(0.102*exp(4.7*umg)+0.45*zrt-0.35)*0.0024
C            UMG - maximum hygroscopicity
C            USS - static water  
C calculation of standard deviation of saturated hydraulic conductivity (SIGK0)  for the calculation of runoff
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
         W0=nb
         W1M=nb
         w2=w0
         UZIM=W0
         ice=0.  
         psn=0.
         HSNOW=0.
         KSI=0.                        
         ksiot=0.
         ice1=0.
         u1=w0
         sumsn=0.
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
         do i=1,mmm
            hhyd0(i)=0.0000000001
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
          SUMTC = 0.
          sumqbz=0. 
         NLETO=0
         NZIMA=0
        kodzero=0   ! if =1, then zerot1 is not changed; =0, then zerot1 should be recalculated    
        L1 = 0.0049 * RS0 ** 2   ! formula by Bracht   
C*****************************************************************
     
c      totrunof=0.
	
	NYEAR=1
      kodend=1                    
      kodprint=0
	year_print=iyrfir-1
	jdy_pr=0
	i0_cal=0
c      kodprint=1
      
  100   CONTINUE
       if (kodend.eq.4) kodprint=1
c         WRITE(*,*) '************'
c         WRITE(*,*) NYEAR, kodend,kodprint
c         WRITE(*,*) '____________' 
         IMDAY=365
         MON(2)=28  
         jdy=firday   
C IYR_1 - year of run
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
C _____________________________________________________________________
C CALCULATION DAY BY DAY 
      
      DO 30 II=1,IMDAY    
        if (kodprint.eq.1) i0_cal=i0_cal+1                                         
c      if (ii.eq.10) stop
          if (jdy.gt.IMDAY_1) jdy=jdy-IMDAY_1   
           IF (NYEAR.EQ.NUMYR.AND.JDY.EQ.1) GOTO 1333  
c        IF (NYEAR.EQ.NUMYR.AND.ii.EQ.1) stop                                     
c  =0 - for the calculation of daily totals:
      k_snowt=0
      swnet_d=0.
c	sw_d=0.
c	rdown_d=0.
c	u45_d=0.
c	q45_d=0.
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
          b = sdbmes(jdy) * AMB
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
c       lai=3.
c       lsai=4.
c       ef=1.
c       z0veg=0.2
c       d00=0.3
c       alblet=0.2      
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
C-------------------------------------------------------------------------        
C CALCULATION STEP BY STEP 
      DO 200 I=1,nstep    
        istep=(ii-1)*nstep+i        ! number of current time step in a year      
        if ((sumsn+wlsn).eq.0.) rspr=rs0
        if ((snles+wles).eq.0.) rslespr=rs0
        INTER2=0.        ! liquid water on the canopy at the end of time step (after evaporation of intercepted precipitation) (in LETO)
        inter=0.         ! evaporation of intercepted precipitation during the warm season (in LETO)
        KOND=0.
        b=bbb            
C calculation of soil parameters which depend on soil moisture
C CALCULATION BY (3.11), P.111 FROM Kalyuzhnyi,Palova,Lavrov (W/m/grad)           
        L3=0.102*EXP(4.7*W1M)+0.45*ZRT-0.35    ! teploprovodnost' of unfrozen soil 
        L3=L3*0.0024    !  translation into kal/cm/c/grad                                                        
        ZC3 = ZCRT * ZRT + ZCB * ZRB * W1M
C   WH5 - the amount of unfrozen water under the soil temperature =-5C (Kalyuzhnyi and Pavlova)  
        WH5=0.94*WZAV+0.017*ZRT
        ZC2=ZCB*UZIM*ZRB+CIC*ICE*RL+ZCRT*ZRT+(WH-WH5)*LPL*ZRB/5.
        L2 = 1.15 * L3                                  
        ZA3 = L3 / ZC3  
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C reading values of forcing factors
C precipitatiom - kg/m^2/s^1=mm/s; radiation - W/m**2;temperature -K;
C wind speed - m/s at 10 m; pressure - Pa; air specific humidity - kg/kg
        READ (1,1030) SW,RDOWN,rainf,snowf,T45,U45,PRES,Q45  
      if(kod_exp.eq.2) read(2,*) rainf,snowf

c	write(111,*) nYEAR, jdy,i, rainf,snowf
 1030  format (2F10.3,2e14.4,3f10.3,e14.4)  
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
         IF (((HSNOW.EQ.0.).AND.(KSI.EQ.0.).and.
     *       (NZIMA.LE.7)).and.(NLETO.ge.1)) THEN
             kodzim=0
         else
             kodzim=1
         endif
      else
        IF ((((T2daypr.lt.TABS0).or.(HSNOW.gt.0.).or.(KSI.gt.0.).or.
     *  (NZIMA.GE.7)).and.(NLETO.LE.7)).or.(hsnow.gt.0.).or.(ksi.gt.0.))
     ! THEN
             kodzim=1  
         else
             kodzim=0  
         endif
      endif  
C calculation of albedo, sumsn - is  snow pack (mm)         
c        albsnow=0.834-21.988*(rspr-0.1)**3 
c       albsnow=0.8-21.988*(rspr-0.1)**3 
      albsnow=opt_par(7)-21.988*(rspr-0.1)**3 
    
	 if (albsnow.lt.0.1) albsnow=0.1 
	  
	  IF ((sumsn+wlsn). LT. 1.0) THEN 
          ALBEDO = ALBlet + (ALBSNOW - ALBlet) * SQRT((sumsn+wlsn)*0.1)
        ELSE
          ALBEDO=ALBSNOW
        ENDIF                                      
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
c
      CALL ITBLOCK(KEYITER,HSNOW,KSI,KSIOT,A3,A4,
     *D00,LEVEL,TABS0,B,EFIZL,P,U2,T2,Q2,T0,G,EP0,
     *zwind,e2,WLSN,kodzim,radles,tles,LSAI,expLES,d,
     !MLES,ELES,WLES,SNLES,DEXP,exp2,z0veg,ef,ev,trles,LPOD,HPOD,
     ! zimw,HL,BZV,kodsubp,kodsubl)    
c//////////////////////////////////////////////////////////////////////////
C translation of W/m^2 into mm/day   
        EP00=EP0/697.37*24.                         
        eles=eles/697.37*24. 
        trles=trles/697.37*24.
C translating mm/day  into mm/time_step
        EP00=EP00/nstep       
        eles=eles/nstep 
        trles=trles/nstep
        kodhyd=0
c the depth of the second soil zone        
        hgr=h0-hsr*1000./VDHGR
        if ((kodhyd2.eq.1).and.(hgr.lt.h2soil)) then
           hgr=h2soil
        else
           h2soil=hgr  
        endif
        if (h2soil.le.hroot) then
           h2soil=hroot+0.1
           hgr=h2soil
        endif
        kodhyd2=0         
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
      IF (KODZIM.EQ.1) THEN 
         STOK=0.
         SNOWMELT=0.
         ET00=0.
         nhour=0.
         TQ = -(DELTAT/3600./2.)
c--------------------------------
        w2pr=w2
        uzimpr=uzim
        icepr=ice  
        ksipr=ksi
        sumsnpr=sumsn
        sumtppr=sumtp
        sumtcpr=sumtc 
        sumqbzpr=sumqbz
        wlsnpr=wlsn 
        wlespr=wles        
        SNLESPR=SNLES
        rspr0=rspr 
        elespr=eles
        rsles0=rslespr 
  998     CONTINUE 
        CALL ZIMA(EPMMDEC,PRES,T45,LSAI,DELWL,DELSN,u1,ice1,tabs0)
          ETMMDEC=0.
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ELSE
c §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
         wlespr=wles
         del_les=0.    
         del_e=0.                                         
         kodsubp=0
         kodsubl=0   
         snowmelt=0. 
         snowfrz=0.
         snles=0.
         wles=0.
         eles=0.
         trles=0.
         snlespr=0.
         SUMTC=0.  
         sumqbz=0.
         DELSUM=0.
         DELWLSN=0.
         wlsn=0.
         hsnow=0.
         sumsn=0.
         ksiot=0.
         ksi=0. 
         vd=0.
         esnowmm=0.  
         gmlt=0.  
         prsolid=0.
C---------------------
         w2pr=w2
         ep00pr=ep00 
         prhour0=prhour
         tqpr=tq
  999    continue                       
        CALL LETO(tstep,CP,L,RV,PRES,wlespr)  
         kodzero=0 
         DELWL=0.
         DELSN=0.
c DRAIN2 - is drainage from the 2nd layer 
c §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§           
      END IF
        stokcal=stok
c hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
c the calculation of hydrograph:  
       
      call hydrogr(stok,DRAIN2,tothydr,kodhyd,hyd1,hyd2) 
        if (kodhyd.eq.1) then 
          IF (KODZIM.EQ.1) THEN 
              pr=totalpr
              w2=w2pr   
              uzim=uzimpr
              ice=icepr         
              ksi=ksipr
              sumsn=sumsnpr
              sumtp=sumtppr
              sumtc=sumtcpr 
              sumqbz=sumqbzpr 
              wlsn=wlsnpr   
              wles=wlespr  
              eles=elespr
              SNLES=SNLESPR   
              rspr=rspr0 
              rslespr=rsles0
             goto 998
          else   
             w2=w2pr  
             inter2=0.
             ep00=ep00pr 
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
          k80=k80+1
          sumhyd=sumhyd+hyd1-hhyd0(k80)
          hhyd0(k80)=hyd1
          hsr=sumhyd/mmm 
          if (k80.eq.mmm) k80=0
          hgr=h2soil      
 1234 format (3i5,3f14.4,i5,2f14.4,i5,5f14.4,i5)   
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
c W1M - soil moisture in 1-m soil layer
        if (hroot.lt.1000.) then
           w1M=(wk*hroot+w2*(1000.-hroot))/1000.
        else
           w1m=wk
        endif
C RECALCULATION of T0, B, P and EFIZL (U2 in m/s; D0DED and DD in cm/s)            
        DPDED=0.76*U2/(0.82*U2**0.5+1.)  
C  air density in kg/m^3 
        RO2=FUNRO(T2,Q2,PRES) 
C  coefficient for calculation of P in W/m^2                               
        aded=ro2*cp/100.
        keyrecal=0
      IF (KODZIM.EQ.0) THEN 
C SUMMER
           d0ded=5.1*u2**(2./3.)+0.000001
           IF (LEAF.EQ.0.) THEN
              DT0DED=0.
           ELSE
              DT0DED=1.09*U2**(2./3.)/SQRT(LEAF)
           ENDIF   
           DD=D0DED*LSAI*DT0DED/(LSAI*DT0DED+D0DED)
           DD=DD*(1.-EXP(-0.45*LSAI))+DPDED*EXP(-0.45*LSAI)
           T0=(RAD-B-EC+3.*T45**4*SIG*SER+ADED*DD*T2)/(4.*T45**3*SER*SIG
     *      +ADED*DD)  
       P=ADED*DD*(T0-T2)          
      ELSE   
C WINTER
       if (ep00.eq.ecmmdec) goto 111                   
          X40=697.37*60.
          DD=DPDED   
          
          A1Z=RAD-EC+3.*T45**4*SIG*SER+ADED*DD*T2*exples+
     !        sig*ser*(4*t45**3*tles-3*t45**4)*exp2*(1-exples)+
     !        cp*ro2* d*(1-exples)*tles     
          A2Z=4.*T45**3*SER*SIG+ADED*DD*exples+cp*ro2*d*(1-exples)
             IF ((KSI.GT.0.).OR.(HPOD.GT.0.)) THEN 
               IF (KSIOT.EQ.0.) THEN
                 KEYITER=4
                 X10=X40/(KSI/L2+HPOD/LPOD)
                 T0=(A1Z+X10*TABS0)/(A2Z+X10)
                 B=(T0-TABS0)*X10   
               ELSE
                 KEYITER=6
                 X30=X40/(KSIOT/L3+HPOD/LPOD)
                 T0=(A1Z+X30*TABS0)/(A2Z+X30) 
                 B=(T0-TABS0)*X30 
               ENDIF    
             ELSE   
C by Budyko
               dd=dpded
               T0=(RAD-B-EC+3.*T45**4*SIG*SER+ADED*DD*T2*exples+
     !             sig*ser*(4*t45**3*tles-3*t45**4)*exp2*(1-exples)+
     !              cp*ro2*d*(1-exples)*tles)
     
           t0=t0/(4.*T45**3*SER*SIG+ADED*DD*exples+cp*ro2*d*(1-exples))
               IF (t0.ge.tabs0) then
                 keyiter=7
               ELSE
                 keyiter=5
               ENDIF    
             ENDIF  
        P=ADED*DD*(T0-T2) 
        p=p*exples                                                
        keyrecal=1
  111 continue          
      ENDIF   
        EFIZL=SIG*SER*(4.*T45**3*T0-3.*T45**4) 
C the end of the recalculation
c         
C check for heat balance under the trees      
        if (kodzim.eq.1) then
           BALANCE=RAD-EFIZL+sig*ser*(4*t45**3*tles-3*t45**4)*exp2
     !             *(1-exples)
     !        -B-EC-P+cp*ro2*d*(tles-t0)*(1-exples)
        else
           BALANCE=RAD-EFIZL-B-EC-P
        endif   
c check for water balance in the whole system  
        dv2=(w2*(h2soil-hroot)-w2pr*(h2psoil-hroot))  
        dv3=por2*(h2psoil-h2soil)                            
        IF (KODZIM.EQ.0) THEN 
           XBAL_all=totalpr+INTERPR+wlespr-ECMMDEC-XDV-INTER2-tothydr-
     !     dv2-dv3-(hyd1-hyd1pr)*1000.
        ELSE 
           XBAL_all=totalpr+snlespr+wlespr+interpr-(wles+snles)-delsum-
     !     DELWLSN-XDV-tothydr-ECMMDEC-eles-trles-dv2-dv3-
     !     (hyd1-hyd1pr)*1000.
        ENDIF           
c**********************************************************************
c printing descripances in heat and water balances                  
c      if (abs(balance).gt.3.) 
c     !   WRITE(6,1040) nYEAR,jdy,I,EC,P,BALANCE,T0,T2,Tles,ET,EP,
c     *           B,EFIZL,RAD,U45,keyiter,WK,WH,ps
c      if (abs(xbal_all).gt.0.001.and.kod_calibr.ne.1)
c     !   WRITE(333,1059) nYEAR,jdy,I,totalpr,snlespr,wlespr,interpr,wles
c     !   ,snles,delsum,DELWLSN,XDV,tothydr,ECMMDEC,eles,
c     !        trles,dv2,dv3,(hyd1-hyd1pr)*1000.,xbal_all,kodzim
 1059  format(3i5,17f12.4,i3)
 1040  FORMAT(I5,I4,I2,3F8.3,3F7.2,5F8.3,F5.2,I3,3F6.3)     
c***********************************************************************   
        ec=ec+eles_wt+trles_wt
        ecmmdec=ecmmdec+(eles+trles) 
C calculation of daily values           
        T0_d=T0_d+T0/nstep 
        t2_d=t2_d+T45/nstep
c	  u45_d=u45_d+u45/nstep
c	  q45_d=q45_d+q45/nstep 
c gswp2 outputs at 3-hour  time step
      IF (KODPRINT.EQ.1) THEN    
c printing PILPS' output fild:      
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
c        sw_d=sw_d+sw/nstep
c	  rdown_d=rdown_d+rdown/nstep
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
            snfr=1.
            snalb=albedo   
            sliqfr=wlsn/(sumsn+wlsn)   
         else
            snowt=-9999
            snfr=0.     
            snalb=-9999
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
          precip_d=precip_d+totalpr      ! daily rainfall
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
c        totrunof(kstok)=totrunof(kstok)+stok+drhydr
        surrun_transf=surrun_transf+stok
        drain_transf=drain_transf+drhydr
	  snmlt_d=snmlt_d+SNOWMELT         
        snowfr_d=snowfr_d+snowfrz
        vd_d=vd_d+vd
        dvsoil_d=dvsoil_d+delsoilm
        dswe_d=dswe_d+delswe
        dint_d=dint_d+delinter
        hgr_d=hgr/1000.          ! water table depth (m) 
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
cc     !   (subsoil+subles)/deltat,eles/deltat
cc      WRITE(4,479) wk,w2,hgr/1000.,ECMMDEC/DELTAT,totalpr/deltat,
cc     ! (stokcal+drain2)/deltat,(stok+drhydr)/deltat,radt,ksi/100.,
cc     ! ksiot/100.,(sumsn+wlsn)  

      ELSE 
           inter_d=inter_d+INTER
           et_d=et_d+ETMMDEC 
           ep_d=ep_d+EPMMDEC
           esnow_d=0. 
           subsl_d=0.  
cc       WRITE(4,479) epot/deltat,INTER/DELTAT,ETMMDEC/DELTAT,
cc     !   EPMMDEC/DELTAT,wk*HROOT,ESNOWMM/DELTAT,SUBSOIL/DELTAT,
cc     !   ECMMDEC/DELTAT,subles/deltat,SOILWET,w2zon,hgr/1000.,
cc     ! (subsoil+subles)/deltat,eles/deltat
cc      WRITE(4,479) wk,w2,hgr/1000.,ECMMDEC/DELTAT,totalpr/deltat,
cc     ! (stokcal+drain2)/deltat,(stok+drhydr)/deltat,radt,ksi/100.,
cc     ! ksiot/100.,(sumsn+wlsn)  
      ENDIF 
c  479 format (3f6.3,4e14.4,f7.2,2f6.3,f8.2)        
  479 format (9e14.3,f6.3,F9.3,f10.5,2e14.4)    
C *******************  COLD SEASON PROCESSES *****************************
         snfr_d=snfr
         ksi_d=ksi/100.
         ksiot_d=ksiot/100. 
        snow_d=HSNOW/100.     ! snow depth in m
        snalb_d=snalb
c wb.* output files    
c      WRITE(444,478) jdy,i,ECMMDEC/deltat,
c     !  stokcal/deltat,drain2/deltat,SNOWMELT/deltat,snowfrz/deltat,
c     ! vd/deltat,delsoilm,DELSWE,DELINTER,BALANCE,
c     ! prsolid/deltat,prliq/deltat,totalpr/deltat,canopint,xdv,dv2,dv3,
c     ! etmmdec,epmmdec,inter,esnowmm,eles,trles,(sumsn+wlsn),
c     ! (snles+wles),(ALBZIM*(1-exples)+ALBEDO*exples),
c     ! snfr,ksi/100.,ksiot/100.,hsnow/100.,snalb,snles,wles,emkl*lsai
c     ! ,emkh*lsai  
  478   FORMAT(2i4,6e14.3,3f8.2,5e14.3,3e14.3,6f10.3,2f10.3,5f7.3,e14.4
     !  ,2f10.4,2f7.3)  
c **************************************************************
      ENDIF   
        INTERPR=INTER2
        W0=WK                        
        PSN=PS 
        h2psoil=h2soil 
        hyd1pr=hyd1         
  200 CONTINUE                                                       
C      ------------------(i)-----the end of time-step calculations--------------------     
      t0daypr=t0_d 
      t2daypr=t2_d 
      if (k_snowt.gt.0) then
         snowt_d=snowt_d/k_snowt
      else 
         snowt_d=-9999.
      endif  
      if (kodzim.eq.0) vegt_d=-9999.  
       
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
c	write(55,587) year_print,jdy,ec_d/secd,precip_d/secd,hprecp_d/secd
cccc     !	 ,surrun_transf/secd,drain_transf/secd,sw_d,rdown_d,t2_d,u45_d
cccc     !      ,q45_d

c     !,surrun_transf/secd,drain_transf/secd,sw_d,rdown_d,t2_d,u45_d,
c     !q45_d,swe_d,ksi_d,ksiot_d,snow_d

      nmon=1 
	ia3=jdy 
	  if (abs(year_print/4.-int(year_print/4.)).eq.0.) then
            MON(2)=29
         else
	      MON(2)=28
         END IF 
	do ixm=1,12
          ia3=ia3-mon(ixm)
c	Write (*,*) '***********',year_print,jdy,ia2,ia3,ixm,mon(ixm)
		if (ia3.gt.0) then
	       nmon=nmon+1
          else
	       ia2=nmon
	       ia33=ia3+mon(ixm) 
		   goto 9997
	    endif
      enddo
 9997 continue

      
c 	Write (*,*) year_print,jdy,ia2,ia33
c	pause
c     !	(surrun_transf+drain_transf)
	 n1(i0_cal)=year_print
       n2(i0_cal)=ia2
       n3(i0_cal)=ia33
       run_day(icell,i0_cal)=surrun_transf+drain_transf
	 ec_day(icell,i0_cal)=ec_d
	 pr_day(icell,i0_cal)=precip_d+hprecp_d 
    	
ccc	write(55,587) year_print,jdy,
ccc     !	 surrun_transf/secd,drain_transf/secd
  587 format(i5,2i4,e14.3)
ccc  587 format(i5,i4,5e14.3)
c  587 format(4i5,14e14.3)

  585 format(i5,i4,8e14.3,2f10.3,e14.4,2f8.2,f5.2)
c d0506.* output daily files
c      write(19,586) nyear,jdy,epot_d/secd,inter_d/secd,et_d/secd
c     ! ,ep_d/secd,WK_d,canint_d,evapsn_d,esnow_d/secd,
c     ! subsl_d/secd,hgr_d
  586 format(i5,i4,4e14.3,3f10.2,2e14.3,f10.3)
      jdy_pr=jdy
	ENDIF   
         jdy=jdy+1  
   30 CONTINUE 
C      ------------------(ii)-----the end of daily calculations------------------------
c       goto 1254
C  EQUILIBRIUM (the 4th year will be printed - the last year)      
      if (kodend.le.3) then
         kodend=kodend+1    
         close(1) 
	   close(2)
	   close(60)	                              
c      if (icell.eq.8194.or.icell.eq.8195) then     ! for VALDAI     
c      if (icell.eq.icel_cal) then     ! for VALDAI
       if (ncl.lt.10) then
          WRITE(ICHAR1,'(I1)') ncl
          open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !		ichar1//'.dat')
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar1//'.dat')
				if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !ichar1//'.dat')
       else

	if (ncl.ge.100000) then
           WRITE(ICHAR6,'(I6)') ncl
           open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !		 ichar6//'.dat')
	    else
         if (ncl.ge.10000) then
           WRITE(ICHAR5,'(I5)') ncl
           open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !		 ichar5//'.dat')
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar5//'.dat')
     			if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !  ichar5//'.dat')	   
         else 
           if (ncl.ge.10.and.ncl.lt.100) then
              WRITE(ICHAR2,'(I2)') ncl
              open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !			ichar2//'.dat') 
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar2//'.dat')
     			if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !  ichar2//'.dat')  		
			  
           else
              if (ncl.ge.100.and.ncl.lt.1000) then
                WRITE(ICHAR3,'(I3)') ncl
              open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !			ichar3//'.dat')  
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar3//'.dat')
     			if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     !   ichar3//'.dat')			
			       
              else
                 WRITE(ICHAR4,'(I4)') ncl
              open (1,file='d:\ISI-MIP2\forcing_data\'//river//'\F'//
     !			ichar4//'.dat') 
	if(kod_exp.eq.2) open (2,file='d:\ISI-MIP2\forcing_data\p3_'//
     !	ichar4//'.dat')
     			if(kod_prec.eq.2)
     !    open (60,file='d:\ISI-MIP2\forcing_data\'//river//'\P'//
     ! ichar4//'.dat')			     
              endif
           endif
         endif
       endif  
	 endif   
c      endif
         goto 100
      endif
 1254 continue
C__________________________________
        NYEAR=NYEAR+1
C   *********************************************      
      IF (NYEAR.LE.NUMYR) GOTO 100  
 1333   CLOSE (1) 
        close(60)
	  close(2) 
        CLOSE (33)
        CLOSE (6)
        close (4)
        CLOSE (12)
        CLOSE (15)
ccc        close (333)
        close (444)
        close (19)
	close (55)
	
      if (kod_calibr.eq.2) WRITE(*,*) 'end plus +++',ncl 
c      WRITE(*,*) 'end',ncl  
 1879 continue

      enddo
c      STOP
      close(1)
	close(2)
	close(60)
	close(7)
	close(10)
      return
 1880 END                                   
 
C \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      subroutine messyt(nday,x,xcyt,im)
      implicit none
C transformation of monthly values into daily
      INTEGER DMES,i,j,im,nday(12),IMIN1,IPL1
      REAL x(12),Xcyt(im)
      DO 10 j=1,IM
      dmes=j
C determination of current month (i) and day of month (dmes)
      DO 20 i=1,12
        dmes=dmes-nday(i)
        if (dmes.LE.0) then
           dmes=dmes+nday(i)
           GOTO 30
        endif
   20  CONTINUE
   30  CONTINUE
      if (dmes.EQ.15) then
       xcyt(J)=x(i)
      else
       IMIN1=I-1
       IPL1=I+1
       IF (IMIN1.EQ.0) IMIN1=12
       IF (IPL1.EQ.13) IPL1=1
       if (dmes.LT.15) then
          xcyt(J)=x(i)-(x(i)-x(IMIN1))/(nday(IMIN1)*1.0)*(15-dmes)
       else
          xcyt(J)=x(i)+(x(IPL1)-x(i))/(nday(i)*1.0)*(dmes-15)
       end if
      end if
  100 FORMAT(I4,F6.2,I4,F6.2,I4,F6.2)   
   10 CONTINUE   
      RETURN
      END  
	 
C//////////////////////////////////////////////////////////////////////////      
C iterations' block                                         
      SUBROUTINE itblock(KEYITER,HSNOW,KSI,KSIOT,A3,A4,
     *D00,LEVEL,TABS0,B,EFIZL,P,U2,T2,Q2,T0,G,EP0,
     *zwind,e2,WLSN,kodzim,radles,tles,LSAI,expLES,d,
     !MLES,ELES,WLES,SNLES,DEXP,exp2,z00,ef,ev,trles,LPOD,HPOD
     !,zimw,HL,BZV,kodsubp,kodsubl)
      
      IMPLICIT NONE
                                          
      REAL KSI,KSIOT,LEVEL,LZV,KINV,L,LPL,KAPA,L1,L2,L3,deltat,radles
      REAL B,Q2,WLSN,T0,EP0,T2,G,U2,ZWIND,FF,HSNOW,P,CP,RAD
      REAL TABS0,WH,RO,SIG,D00,RV,SER,FF10,A3,A4,A5,E2,PRES,Q45,EFIZL
      REAL Z00,T45,U45,FUNRO,XX,YY,ZZ,SUMFF,SUMFFU2,SUMFF10,UZV
      REAL X30,X20,X40,X60,X70,X50,X71,X72,TZV,QZV,TETA0,TETA45,Q0
      REAL DEL2,T0ST,ZX,FF0,PAR,FFU2,TC2,Q2SAT,DQ2,DE2
      REAL tles,bles,T00PR,LSAI,RADSUR,EFIZLES,exples,d,tp,bzv,exp2
      REAL BALLES,BALSUR,HL,PL0,snles,wles,mles,eles,DEXP,tlespr,cles
      REAL trk,trles,ef,ev,LPOD,HPOD,zimw,le,lcan
      
      INTEGER ISTEP,I,KODZIM,II,KEYITER,NYEAR,ITER,KEYKOD,K57,NSTEP,LES
      integer kodsubp,kodsubl,jdy
      
      COMMON /METEO/ U45,T45,Q45,RADsur,PRES
      COMMON /PARAM/ KAPA,SIG,SER,Le,LPL,CP,L1,L2,L3,WH,KINV,RV
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy
      FUNRO (XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))

        A5 = (1000. / PRES) ** .288                      
C Setting up initial (at the first iteration) values  
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
        IF((HSNOW.GT.0.).OR.((KSI.GT.0.).AND.(KSIOT.EQ.0.))) then
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
        if (kodzim.eq.1.and.exples.ne.1) then
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

 3747 if (kodzim.eq.1) then
        if (les.eq.1) then
          tp=t0
          rad=radles 
          bles=-sig*ser*(4*T45**3*Tp-3*T45**4)     
        else
          HL=-KAPA*CP*RO*UZV*TZV        
          TLES=T0
          rad=radsur
          EFIZLES=SIG*SER*exp2*(1-EXPLES)*(4*T45**3*TLES-3*T45**4)
          RAD=RAD+EFIZLES              
        endif 
      endif   
c        
        SUMFF10=FF10
        SUMFF=FF
        sumffu2=ffu2
        ITER=1
C calculation of the dynamical speed
 4300   UZV=U45*KAPA/FF10 
        U2=UZV*FFU2/KAPA  
        d=0.0076*u2*exp(-DEXP*lsai)/(1+0.82*sqrt(u2))   
        if (les.eq.1) goto 4004
        IF (kodzim.eq.1) then
C winter                 
        k57=0
        IF ((HSNOW.GT.0.).AND.(KEYITER.EQ.3)) THEN
          T0=TABS0
        ELSE 
          X30=SIG*SER*T45**4-RAD-KAPA*RO*L*UZV*(Q45-A3)/FF
     !        -cp*ro*d*(1-exples)*(tles-t45)
          X20=KAPA*RO*L*UZV*A4/A5+KAPA*RO*CP*UZV*exples+
     !        4.*SIG*SER*FF*T45**3./A5+cp*ro*d*(1-exples)*ff/a5
          X40=697.37*60.
          X60=T45-FF*X30/A5/X20
          X70=FF/A5/X20                    
          IF (HSNOW.GT.0.) THEN
            IF (WLSN.GT.0.) THEN   
               X50=HSNOW/2./L1 
               KEYITER=2                     
            ELSE                               
               X50=HSNOW/L1+KSI/L2+HPOD/LPOD 
               KEYITER=1 
            endif                                
            X71=X40*X70/X50
            T0=(X60+X71*TABS0)/(1.+X71)        
            IF (T0.GE.TABS0) THEN
              T0=TABS0
              KEYITER=3
            ENDIF  
          ELSE              
  554       IF ((KSI.GT.0.).OR.(HPOD.GT.0.)) THEN             
              X50=HSNOW/L1+ksi/L2+HPOD/LPOD 
              X71=X40*X70/X50
              T0=(X60+X71*tabs0)/(1.+X71)
              KEYITER=4                     
              if ((t0.lt.tabs0).and.(ksiot.gt.0.)) then
                 X50=KSIOT/2./L2+HPOD/LPOD 
                 X71=X40*X70/X50 
                 T0=(X60+X71*tabs0)/(1.+X71)
                 keyiter=8
              endif    
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
        TETA0=T0*A5
        TETA45=T45*A5
        TZV=(TETA45-TETA0)/FF
        QZV=(Q45-A3-A4*T0+A4*T45)/FF
      ELSE             
c summer
C calculation of the scale of temperature
             TZV=SIG*SER*T45**4+B-RAD-KAPA*RO*L*UZV*(Q45-A3)/FF
             TZV=TZV/(KAPA*RO*L*UZV*A4/A5+KAPA*RO*CP*UZV+
     *        4.*SIG*SER*FF*T45**3./A5)
C calculation of the surface temperature
             T0=T45-TZV*FF/A5
C calculation of the humidity scale
             QZV=(Q45-A3)/FF+A4*TZV/A5
      END IF                   
C calculation of the saturation humidity at the surface
 4375     Q0=A3+A4*(T0-T45)
C     calculation of the air density
          RO=FUNRO(T0,Q0,PRES) 
c calculation of the forest crowns:
 4004 if (les.eq.1) then  
       IF (SNles.eq.0.) THEN          
          if (wles.eq.0.) then  
             trk=zimw*ef*ev         ! transpiration  (TRLES)
          else 
             trk=1.                 ! evaporation of intercepted precipitation (ELES)  
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
         trk=1.         ! evaporation of intercepted precipitation (ELES) 
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
     !          trk*Lcan*qzv)-cp*ro*d*(t0-tp)-cles*(t0-tlespr)
         ELSE
           mles=0.
         ENDIF   
       endif
       HL=-KAPA*CP*RO*UZV*TZV     
      ENDIF                 
C     calculation of the scale LZV
          LZV=(TZV/T45+0.61*QZV)*G*KAPA**2
          LZV=UZV**2/LZV 
      IF (ITER.EQ.1) GOTO 4750
    
                      
      del2=t0-t0st 
      IF ((ABS(DEL2).LT.E2).and.(iter.ge.3)) GOTO 4771 
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
              les=1
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
              IF ((KEYITER.EQ.1).OR.(KEYITER.EQ.4).or.(keyiter.eq.2).
     !         OR.(KEYITER.EQ.6).or.(keyiter.eq.8)) B=(T0-TABS0)*X40/X50
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
       if (kodzim.eq.1.and.keykod.lt.1000) then  
            bles=sig*ser*(4*T45**3*T0-3*T45**4) 
            PL0=cp*ro*d*(tles-t0)
            bzv=cles*(tles-tlespr)
            efizles=sig*ser*(4*T45**3*Tles-3*T45**4) 
            eles=eles*lcan
            trles=trles*lcan
         BALLES=RADLES-(1+exp2)*efiZLES+BLES-HL-PL0-bzv-eles-
     !    trles-mles
         BALSUR=RADSUR-BLES+EFIZLES*exp2*(1-EXPLES)-EP0*l-P-B+
     !          PL0*(1-exples)         
                                             
           mles=mles*(1.-EXPLES)
           eles=eles*(1.-EXPLES)
           trles=trles*(1.-exples)   
c       write(19,440) ii,i,balsur,balles,radles,(1+exp2)*efiZLES,BLES,HL,
c     ! PL0,bzv,eles,trles,mles,t0-tabs0,tles-tabs0,t45-tabs0,
c     !snles,wles
c  440  format(2i5,16F11.4)       
       endif 
       ep0=ep0*L
       tlespr=tles            
      RETURN
      END
C //////////////////////////////////////////////////////////////////////////
      SUBROUTINE UNFUN(ZX,D00,Z00,LZV,FF,UZV,KINV,ff0,par)
      IMPLICIT NONE
      
      REAL PAR,ZX,D00,PAR0,UZV,FF,Z00,FF0,X1,X2,FUNK1,FUNK0
      REAL LZV,KINV                                        
      
         PAR=(ZX-D00)/LZV
         PAR0=Z00/LZV 
         IF (PAR.GT.0.) then 
            if (par.ge.0.5) then 
                par=0.5 
                lzv=(zx-d00)/par
                par0=z00/lzv
             endif   
             FF=LOG((ZX-D00)/Z00)+10.*((ZX-D00)-Z00)/LZV  
         endif
         IF ((PAR.GE.-0.07).AND.(PAR.LE.0.)) FF=LOG((ZX-D00)/Z00)
         IF (PAR.LT.-0.07) THEN
            X1=ABS(PAR)
            X2=ABS(PAR0)
            FUNK1=0.25-1.2*1./(X1**(1./3.))
            FUNK0=0.25-1.2*1./(X2**(1./3.))
            FF=FUNK1-FUNK0
         END IF          
C taking into account the layer between the surface and z0     
         FF0=0.13*(Z00*UZV/KINV)**0.45
         FF=FF+FF0
      RETURN
      END
C //////////////////////////////////////////////////////////////////////////                                        
      SUBROUTINE LETO(tstep,CPLETO,L,RV,PRES,wlespr)   
  
      IMPLICIT NONE
      
      REAL M1,NB,LAI,k0,k00,lc,l,lleto,m,LSAI,kw2,KW3,por2
      REAL ET,EP,XX,YY,ZZ,RV,FUNRO,H2PSOIL,T2,U2,CPLETO,WZAV,HROOT
	REAL TSTEP,ET0,epot,tq
      REAL EFF,w2,PR,WK,PS,BPAR,WN,Q2UP,UMG,H2SOIL,PSN,POR,Q22,A4,EP00
      REAL ET00,PRES,Q2DOWN,PODZSTOK,FI0,EV,ETDAY,FLAI,FLSAI,WKP,F0
      REAL EP0,FI00,RO,CP,DP,DWV,QPS,DIF,APS,SBR1,SBR2,W1,SBR,ww,h0
      REAL DIF3,Q3UP,Q2SUM,Q3DOWN,DRAIN,fik,k02,q3hgr,deltat,pr_int
      REAL aa,bb,epc0,prmmhour,k0zv,emk,leaf,inter,stok,nhour,prhour
      REAL y5,hk,dw,x5,x6,sigk0,x7,q,trmn,inter2,kond,hpod,emkl,drain1 
      REAL interpr,wlespr,w2pr
      EQUIVALENCE (POR,POR2)
      
      INTEGER I,II,ISTEP,NYEAR,NSTEP,kodhyd,kodhyd2,jdy
      
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy
      COMMON /EVAP/ ET,EP,LAI,LSAI      
      common /Cleto/ nhour,prhour,tq,hk,sigk0,trmn,leaf,
     !              ET00,inter,inter2,kond
       
      common /year/ NB,k0,HROOT,bpar,wzav,fi0,POR,umg,a4,u2,lc,
     !              h2soil,h2psoil,WN,WK,w2,PS,PSN,STOK,podzstok,
     !              EP00,PR,T2,Q22,kodhyd,kodhyd2,hpod,emkl,h0,
     !              pr_int,q3hgr,drain1,eff,ev,epot,interpr,m,
     !              w2pr

      FUNRO (XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))
c  the units: mm/time_step, kal/g/grad, g/cm^3 etc.
C  all the components of water balance are in mm/time_step (here, mm/3hour)
c  WK, WN and NB are full (not available) values of soil moisture
C  PS - the thickness of drying layer, cm
c  PSN - the value of PS from previous time-step, cm
C      constants:
         DW=POR-NB
         if ((kodhyd2.eq.1).and.(q3hgr.ge.0.))  h2soil=h2psoil      
         M1=0.45                               
         prmmhour=pr
         kond=0.
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
                     if (stok.lt.0) stok=0.
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
         k02=k00*(POR2/POR)**(2*bpar+2.)
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
C           the calculation of Q2UP                            
                   ww=(wn*hroot+w2*(h2soil-hroot))/h2soil
        if (ww.lt.0) write(*,*) ii,i,wn,w2,h2soil,hroot,ww           
                   kw2=k00*(ww/por)**(2.*bpar+3.)
                   dif=kw2*fi00/por**(-bpar)*(-bpar)*ww**(-bpar-1)
                   q2up=dif*(w2-wn)/(h2soil/2.)
         IF (Q2UP.LT.0.) Q2UP=0.
         QPS=QPS-q2up
         APS=(LLETO*A4*RO*CP/(LC+RO*CP*DWV)+CP/DWV)*DP/(CP+LLETO*A4)
c QPS in mm/time_step, PSN in cm, APS in 1/cm
         M=1./(1.+APS*(PSN+HPOD))  
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
            PS=PSN+0.1*(EP+QPS)/(WN-UMG+0.0000001)
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
C calculation of SBR from the 1st zone
       sbr1=0.
       sbr2=0.
       if (WK.GT.POR) THEN
          SBR1=WK-POR
          WK=POR
       ENDIF                         
       IF (WK.GT.NB) THEN
            W1=WK
            WK=NB+(W1-NB)*EXP(-4.*K00*((nb-wzav)/(por-wzav))**3/
     *         hroot/(por-wzav))
            SBR2=W1-WK
       ENDIF  
       sbr=sbr1+sbr2
c calculation of Q3UP, Q2DOWN, Q2SUM and soil moisture (w2) in the 2nd zone             
          ww=(w2*(h2soil-hroot)+por2*(h0-h2soil))/(h0-hroot) 
      if (ww.lt.0.) write(*,*) II,I,w2,h2soil,hroot,por2,h0,ww,k02,bpar
          kw3=k02*(ww/por2)**(2.*bpar+3.)
          dif3=kw3*fi00/por2**(-bpar)*(-bpar)*ww**(-bpar-1)
          fik=fi00*(ww/por2)**(-bpar)
       q2down=sbr*hroot
       q2sum=q2down-q2up
       drain1=q2sum
c        
        if (kodhyd2.eq.1) then
           q3up=0.
        else 
           if (kodhyd.eq.1) then
              if (q3up.gt.0.) then
                 q3up=0.
              else
                 h2soil=h2psoil 
                 kodhyd2=1 
              endif
           else 
              q3up=dif3*(por2-w2)/((h2soil-hroot)/2.)
              if (q3up.le.0.) q3up=0. 
           endif     
        endif   
        q3hgr=(por2-w2)*(h2psoil-h2soil)    
       
        w2=w2pr+(q2sum+q3up)/(h2soil-hroot)   
       
        if (w2.lt.0.) then
            q2sum=0.
            w2=w2pr+(q2sum+q3up)/(h2soil-hroot)   
         endif        
C calculation of SBR from the 2nd zone
       sbr1=0.
       sbr2=0.
       if (W2.GT.POR2) THEN
          SBR1=W2-POR2
          W2=POR2
       ENDIF                         
       IF (W2.GT.NB) THEN
            W1=W2
            W2=NB+(W1-NB)*EXP(-4.*K02*((nb-wzav)/(por2-wzav))**3/
     *         (h2soil-hroot)/(por2-wzav))
            SBR2=W1-W2
       ENDIF      
       sbr=sbr1+sbr2
c calculation of drainage out of the 2nd zone (with taking into account
c change in the water table
       q3down=sbr*(h2soil-hroot)
       drain=q3down-q3up
       podzstok=drain-q3hgr   
c the depth of the 2nd soil zone on the previous time step    
      RETURN
      END
C////////////////////////////////////////////////////////////////////////////
      SUBROUTINE ZIMA(EP,PRES,TMLL,lsai,DELWL,DELSN,u1,ice1,tabs0)     
    
      IMPLICIT NONE
                                                            
      REAL L1,L2,L3,L,LZV,LL,NB,KSI,KSI0,KSIOT,KFILTR,K0,ICE,LLET
      REAL LSUM,lc,m,ksimm,ksiotmm,ice2,ice1,KW2,KW3,ww,iice,w3,por2
      REAL H2PSOIL,STOK,Q2,CPZIM,WLSN,T0DAYPR,ZEROT1,U2,RSPR,W0,WZAV
      REAL H,DELSUM,ESNOW,HROOT,SUMTC,P,QVPIT,RB,W1M,U,EP,DELWLSN
	REAL q3hgr,DEL1
      REAL SUMTP,TABS0,VD,BZ,TCCC,RL,WH,PR,WK,PS,BPAR,QT,RV,US,UMG,h0
      REAL H2SOIL,PSN,TBYD,GMLT,A3,A4,EP00,C2,SUM,SNOWMELT,eles0
      REAL TMLL,PRES,H2,Q2DOWN,DELTAT,PODZSTOK,FI0,FUNRO,XX,YY,ZZ,k02
      REAL T0DAY,TC,TML,TZV,TBYDC,T0,FI00,RO,CP,DP,DWV,APS,BZIM,GWLSN
      REAL WLSNPR,HT,HSUM0,ESOIL,QBZIM,TCC,RS,HH,DV1,DV,ESN,WLSN0,TP1
      REAL DQBZIM,S,S1,QT1,TSR,UNN,XXUU,U1,DIF,Q2UP,QPS,SBR1,SBR2,fik
      REAL SBR,DIF3,Q3UP,Q2SUM,W2,UU2,Q3DOWN,DRAIN,DLW,DRAIN1,DEL2
      REAL KAPA,SIG,SER,KINV,MLES,ELES,WLES,SNLES,ML,lsai,ESNW,emkh
      REAL rsles,rslespr,delwl,delwl2,delwl3,delsn,delsn2,EMK,R77,HPOD
      REAL PRECIP,trles,LPOD,sncover,kondliq,sumqbz,sumq,emkl                                
      REAL pr_int,prh_int,prsolid,del_les,ef,ev,del_e,zimw,epot    
      REAL snowfrz,interpr,eleskond,w2pr
      
      INTEGER I,II,ISTEP,NYEAR,NSTEP,IM,KEYITER,KODZERO,kodhyd,kodhyd2
      integer jdy,kodprint   
      EQUIVALENCE(P,POR2)
      
      COMMON /TIME/ I,II,ISTEP,NYEAR,deltat,NSTEP,jdy 
      COMMON /PARAM/KAPA,SIG,SER,LLET,L,cpzim,L1,L2,L3,WH,KINV,RV
      common /year/ NB,KFILTR,HROOT,bpar,wzav,fi0,P,umg,a4,u2,lc,
     !              h2soil,h2psoil,W0,WK,w2,PS,PSN,STOK,podzstok,
     !              EP00,PRECIP,TBYD,Q2,kodhyd,kodhyd2,hpod,emkl,h0
     !              ,pr_int,q3hgr,drain1,ef,ev,epot,interpr,m,w2pr
      
      common /Czima/ A3,C2,RB,RL,US,RS,TZV,ZEROT1,IM,
     !              TCCC,t0daypr,SUMTC,SUMTP,W1M,KSI,KSIOT,H2,
     !              WLSN,DELWLSN,gmlt,SUM,DELSUM,BZ,VD,
     !              SNOWMELT,uu2,ice2,ESNOW,KEYITER,kodzero,
     !              MLES,ELES,WLES,SNLES,rsles,WLSNPR,trles,LPOD,
     !              sncover,sumqbz,emkh,prh_int,prsolid,
     !              ml,del_les,del_e,zimw,kodprint,snowfrz
   
     
      FUNRO (XX,YY,ZZ)=ZZ*100./(RV*XX*(1.+0.61*YY))
C I - step; II- day 
c WLSN  - volume of liquid water in snow cover (MM)
C W0    - total initial moisture in rooting zone
c w2   -total soil moisture in the 2nd zone
c KSI   - soil freezing depth (cm)
C translation of the surface temperature TC, temperature at the GCM's
C lowest level TML, temperature at 2m TBYDC from K to C
C TZV - constant temperature of deep soil 
        epot=eles+trles+ep00  
        eleskond=0.  
        QVPIT=0. 
        snowfrz=0.
        if ((kodhyd2.eq.1).and.(q3hgr.ge.0.)) h2soil=h2psoil
        if (ksi.eq.0.) then
           ice2=0.                  
           uu2=w2
           ice1=0.
           u1=w0
        endif                   
        PR=PRECIP 
        rspr=rs
        rslespr=rsles    
        TC=TCCC-TABS0   
        IF (TC.GT.0.) then
           TCC=0.         
        else
           TCC=TC 
        endif      
       rsles=rslespr*(1.+0.1/10./nstep*(snles+wles)*EXP(.08*TCc-
     !        21.*RSlespr)) 
        ESN=0.11/RSles-0.11   
        ESNW=ESN/(1-ESN)   
       emk=emkh*lsai 
        t0day=t0daypr-tabs0   
        
        TML=TMLL-TABS0
        TBYDC=TBYD-TABS0   
c PSN - dry layer depth at the beginning of the time step (i.e. from previous time step),cm        
C constants and parameters
C K0 - coef. of filtration in mm/time_step;       
        K0=KFILTR*24./nstep 
        k02=k0*(POR2/P)**(2*bpar+2.)                                  
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
        M=1./(1.+APS*(PSN+HPOD))              
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
               if(ml.gt.snles) ml=snles       ! snowmelt
             else
               if(abs(mL).gt.wles) ml=-wles   ! freezing of water in snow
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
c
         QBZIM=qbzim+(BZIM-VD) 
         sumqbz=sumqbz+qbzim
         bzim=vd
c now, BZIM is REAL energy consumed to snowmelt (mm/time_step)
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
         
         IF (KSI.LE.0.) KSI=0.0001 
         KSI0=KSI
C calculation of snow density
C RSPR -snow density at previous time-step
        RS=RSpr*(1.+0.1/10./nstep*(H+wlsn)*EXP(.08*TCc-21.*RSpr)) 
C calculation of snow heat conductivity  
c formula by ABELS         
c         L1=0.0068*RS**2  
c formula by Bracht (Kuchment, 1983, p. 135)          
          L1=0.0049*RS**2
C         L1=0.00122*RS
C formula by Yanson (p.60 - Gusev)
c         L1=(0.05+1.9*RS+6.*RS**4)/1000.
C calculation of snow height (cm) using snow pack (mm)
         H2=(H+wlsnpr)/10./RS 
         if (h.eq.0.) h2=0.
C calculation of 'privedennoi' snow height
         HH=H2*L2/L1+hpod*L2/LPOD 
C calculation of change in volume of liquid water in snow cover
       DV1=0.  
       ESN=0.11/RS-0.11  
       WLSN0=RS*H2*ESN*10./RB    
       IF ((TC.LT.0.).AND.(WLSNpr.GT.0.)) THEN
c freezing of water in 'wet' snow cover
c here, snowfreezing is positive
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
          WLSN=WLSNpr+VD
          IF (WLSN.LE.WLSN0) THEN
             VD=0.       
          ELSE
             VD=WLSN-WLSN0
             WLSN=WLSN0
          ENDIF             
          gwlsn=0.
       ENDIF 
C calculation of soil surface temperature    
         TP1=KSI0*t0day/(KSI0+HH)
         IF (TP1.GT.0.) TP1=0.
         IF (TP1.LE.-0.1) THEN 
            SUMTP=0.
         ELSE
            SUMTP=SUMTP+DELTAT
         END IF                 
C heat flux from unfrozen zone to the freezing depth
         istep=(jdy-1)*nstep+i
         IF (jdy.GT.182) THEN
           if (kodzero.eq.0) then
              kodzero=1
              zerot1=jdy-10.
           endif  
           QT=9./4.*KSI0**2+12.*A3*(DELTAT*(iSTEP-ZEROT1*nstep)) 
      if (qt.lt.0.) write(*,*) '111',jdy,i,ksi0,jdy,istep,zerot1*nstep,
     !        nstep,qt
         ELSE                                           
           QT=9./4.*KSI0**2+12.*A3*(DELTAT*(iSTEP+IM-ZEROT1*nstep))   
           kodzero=0
           IF (QT.LT.0.) WRITE(*,*) jdy,I,QT,ISTEP,ZEROT1*nstep,IM,
     !      (iSTEP+IM-ZEROT1*nstep),DELTAT,A3,KSI0 
         END IF
         QT=2.*L3*TZV/(SQRT(QT)-3./2.*KSI0)   
       SUMTC=SUMTC+TC  
       DQBZIM=0.
C calculation of soil freezing depth KSI
       LL=L*RB*(W1m-WH)
       IF (LL.EQ.0.) LL=1E-11 
      
       IF ((W1m-WH).le.0.01) THEN 
          KSI=0.
          ksiot=0.
       ELSE
          LZV=L*RB*(W1m-WH)+C2*ABS(TP1)/2.
          IF (LZV.EQ.0.) LZV=1E-11
          S=QT*DELTAT/LZV
          S1=(KSI0+HH)**2-(2.*L2*TCC*DELTAT/LZV)*sncover+S**2 
          IF (S1.GT.0.) GOTO 60
            KSI=1E-11
            GOTO 70
   60     KSI=-HH-S+SQRT(S1)
          IF (KSI.LE.0.) KSI=1E-11
   70     CONTINUE 
          KSI=KSI-DV1/(W1m-WH)
          if (h0.le.2000.) then
             if (ksi.gt.200.) ksi=200.     
          else
             if (ksi.gt.(h0/10.)) ksi=h0/10.
          endif
            
C calculation of soil thawing depth  KSIOT
          IF (H2.GT.0.001) THEN    
            IF (TC.GT.0.) SUMTC=SUMTC-TC    
             if (sumtc.GT.0.) then
                KSIOT=-L3/LPOD*HPOD+SQRT((L3/LPOD*HPOD)**2+
     !                 2.*L3*SUMTC*DELTAT/LL)  
             else
                ksiot=0.
             endif      
          ELSE
             IF (SUMTC.LE.0.) SUMTC=0.
          KSIOT=-L3/LPOD*HPOD+SQRT((L3/LPOD*HPOD)**2+
     !              2.*L3*SUMTC*DELTAT/LL)
          END IF                                        
       ENDIF                                          
C calculation of infiltration QVPIT and runoff STOK 
        dqbzim=0.
        IF ((KSI.LT.0.000001).OR.(KSI.LE.KSIOT)) THEN
           QVPIT=VD
           STOK=0.
           KSI=0.
           KSIOT=0.
           unn=w0
           ice=0.
        ELSE            
           QT1=8.*L2*0.1*SUMTP/(KSI-KSIOT)**2  
C calculation of unfrozen water UNN and ICE in frozen zone
           TSR=TP1
           IF (TSR.LT.-0.1) THEN 
              UNN=WH
           ELSE
              IF (QT1.LT.LL) THEN
                 UNN=WH+(W0-WH)*QT1/LL
              ELSE
                 UNN=W0
              ENDIF
           ENDIF
           ICE=(W0-UNN)*RB/RL 
           dqbzim=0.
           IF (ICE.LE.0.) then
             dqBzim= -ICE*RL*L*(KSI-KSIOT)
             ICE=0.     
           ENDIF 
C the calculation of liquid water into wet zone of frozen soil(by Ye.M.Gusev)
           IF(VD.LT.0.) WRITE(*,*) '***',II,I,ICE,VD,DELSN,DELWL,
     !           DELWL2,DELWL3
           XXUU=VD*(1.+8.*ICE)**2/K0  
           XXUU=SQRT(XXUU)
           XXUU=SQRT(XXUU)
           U=US+(P-US)*XXUU
           IF (U.LT.(P-ICE)) THEN
              QVPIT=VD
              STOK=0.
           ELSE
              QVPIT=K0*((P-ICE-US)/(P-US))**4/(1.+8.*ICE)**2
              STOK=VD-QVPIT
              U=P-ICE
           ENDIF
        ENDIF                   
C the calculation of Q2UP 
            ww=(u1*hroot+uu2*(h2soil-hroot))/h2soil
            iice=(ice1*hroot+ice2*(h2soil-hroot))/h2soil
       if (ww.le.0.) write(*,*) i,ii,ww,u1,uu2,w2,ice2,h2soil,HROOT,
     !  h2psoil
        kw2=k0*(ww/p)**(2.*bpar+3.)/(1.+8.*iice)**2
        dif=kw2*fi00/p**(-bpar)*(-bpar)*ww**(-bpar-1)/(1.+
     !      8.*iice)
        q2up=dif*(uu2-u1)/((h2soil)/2.)           
        if (q2up.lt.0.) q2up=0.
C the calculation of  EP (mm/time_step) 
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
c calculation of liquid water U1 in root zone
        ksimm=ksi*10.
        ksiotmm=ksiot*10.
        if ((ksimm.eq.0.).or.(ksiotmm.gt.hroot)) then
           ice1=0
        else    
           if (ksiotmm.eq.0.) then
             if (ksimm.gt.hroot) then 
               ice1=ice
             else
               ice1=ice*ksimm/hroot
             endif
           else            
             if (ksimm.gt.hroot) then
               ice1=ice*(hroot-ksiotmm)/hroot
             else
               ice1=ice*(ksimm-ksiotmm)/hroot
             endif
           endif        
        endif    
        
        if ((ice1*RL/RB).GT.(WK-UNN)) ICE1=(WK-UNN)*RB/RL
        if (ice1.lt.0.) ice1=0.
	  u1=wK-ice1*RL/RB  
        
C calculation of  SBR 
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
                 SBR=Wk-W3
                 Wk=W3
           ELSE
                 SBR=0.
           ENDIF         
        END IF              
        u1=wk-ice1*RL/RB
C the calculation of Q3UP, Q2DOWN, Q2SUM 
       if (kodhyd2.eq.1) then
         q3up=0.
       else 
         if (kodhyd.eq.1) then
           if (q3up.gt.0.) then
             q3up=0.
           else
             h2soil=h2psoil 
             kodhyd2=1
           endif
         else 
           if (ksimm.gt.h2soil) then
              q3up=0.
           else         
             ww=(uu2*(h2soil-hroot)+por2*(h0-h2soil))/(h0-hroot)
             iice=ice2*(h2soil-hroot)/(h0-hroot)
       IF (WW.LT.0.) WRITE(*,*) II,I,WK,U1,WW,UU2,H2SOIL,HROOT,H0,ksi
             kw3=k02*(ww/por2)**(2.*bpar+3.)/(1.+8.*iice)**2
             dif3=kw3*fi00/por2**(-bpar)*(-bpar)*ww**(-bpar-1)/
     !           (1.+8.*iice)
             fik=fi00*(ww/por2)**(-bpar)*(1.+8.*iice)
             q3up=dif3*(por2-uu2)/((h2soil-hroot)/2.)
             if (q3up.lt.0.) q3up=0.
           endif  
         endif  
       endif
        q2down=sbr*hroot
        q2sum=q2down-q2up
        DRAIN1=Q2SUM
c the calculation of UU2 and ICE2 in the 2nd layer 
        q3hgr=(por2-w2)*(h2psoil-h2soil) 
         w2=w2pr+(q2sum+q3up)/(h2soil-hroot) 
         if (w2.lt.0.) then
            q2sum=0.
            w2=w2pr+(q2sum+q3up)/(h2soil-hroot)   
         endif   
       if ((ksimm.lt.hroot).or.((ksimm.gt.h2soil).and.(ksiotmm.gt.
     !    h2soil))) then
          ice2=0.
       else   
          if (ksimm.gt.h2soil) then
             if (ksiotmm.lt.hroot) then
                ice2=ice
             else                 
                ice2=ice*(h2soil-ksiotmm)/(h2soil-hroot)
             endif
          else
             if (ksiotmm.lt.hroot) then    
                ice2=ice*(ksimm-hroot)/(h2soil-hroot)
             else   
                ice2=ice*(ksimm-ksiotmm)/(h2soil-hroot)
             end if
          endif   
       endif
       
       if ((ice2*RL/RB).GT.(W2-WH)) ICE2=(W2-WH)*RB/RL
       if(ice2.lt.0.) ice2=0.
	 uu2=w2-ice2*RL/RB
       
       sbr=0.
C calculation of SBR from the 2nd zone  
      DEL1=0.
      DEL2=0.
      if (ksimm.gt.h2soil) then
        sbr=0.
        IF ((w2-ICE2*RL/RB).gt.(por2-ice2)) then
           DEL2=(w2-(por2-ice2+ice2*RL/RB))*(H2SOIL-HROOT)
           wk=wk+DEL2/HROOT
           DRAIN1=drain1-del2
           w2=por2-ice2+ICE2*RL/RB 
           uu2=w2-ice2*RL/RB   
           IF ((wk-ICE1*RL/RB).GT.(Por2-ICE1)) THEN
              DEL1=(WK-(Por2-ICE1+ICE1*RL/RB))*hroot
              WK=Por2-ICE1+ICE1*RL/RB     
           ENDIF 
           U1=WK-ICE1*RL/RB   
        ENDIF
      else   
        IF ((Por2-ICE2).GT.NB) THEN
           if ((w2-ICE2*RL/RB).GT.NB) then
              if ((w2-ICE2*RL/RB).GT.(Por2-ICE2)) then 
                 SBR1=w2-ICE2*RL/RB-Por2+ICE2
                 W3=Por2-ICE2+ICE2*RL/RB
                 w2=ICE2*RL/RB+NB+(Por2-ICE2-NB)*EXP(-4.*K02*((nb-wzav)/
     !              (por2-wzav))**3/hroot/(por2-wzav)/(1.+8.*ICE2)**2)
                 SBR2=W3-w2
                 SBR=SBR1+SBR2       
              ELSE                  
                 W3=w2
                 w2=ICE2*RL/RB+NB+(W3-ICE2*RL/RB-NB)*EXP(-4.*K02*((nb-
     !          wzav)/(por2-wzav))**3/hroot/(por2-wzav)/(1.+8.*ICE2)**2)
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
        uu2=w2-ice2*RL/RB 
      endif  
c calculation of drainage out of the 2nd zone (with taking into account
c change in the water table        
       q3down=sbr*(h2soil-hroot)
       drain=q3down-q3up
       podzstok=drain-q3hgr  
       
      if (del1.gt.0.) then
       IF (TC.LT.0.) THEN
          SUM=SUM+DEL1
       ELSE
          IF (SUM.EQ.0.) THEN
             STOK=STOK+DEL1
          ELSE
             WLSN=WLSN+DEL1 
             IF (WLSN.GT.WLSN0) THEN
                STOK=STOK+WLSN-WLSN0
                WLSN=WLSN0
             ENDIF   
          ENDIF      
       ENDIF              
      endif 
       DELSUM=SUM-HSUM0
       DELWLSN=WLSN-WLSNPR
       H=SUM                                 
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
      subroutine hydrogr(stok,DRAIN,shydr,kodhyd,hhyd,hh1hyd)  
      
      IMPLICIT NONE
      
      REAL STOK,DRAIN,SHYDR,HHYD,HH1HYD,H,HH1,XI,lat,mn
      REAL YI,E,DY,DX,SHERN,DT,GRAD,AX,AY,H2,H3,R,H1,HH2,HH,BAL
      INTEGER KODHYD,ii,i,istep,nyear,NSTEP,jdy
      
      COMMON /TIME/ I,II,ISTEP,NYEAR,dt,NSTEP,jdy  
      common /hydro/ xi,yi,lat
c	common /manning/ SherN
      common /manning/ mn
c      
c      write(*,*) mn
c	stop
	shydr=0.
      h=hhyd
      hh1=hh1hyd                  
      e=.000001   
c      shern=0.05 
       shern=mn  
c	OPEN (44, file='d:\ISI-MIP2\res\shern.dat')   
	
c	    write(44,*) nyear,ii,jdy,i,shern 

c		 if(ii.eq.365) stop
c calculation of dx and dy for 1x1 grad grid cell   
      dy=111000./2.
      dx=111000.*cos(3.14157*lat/180.)/2. 
C hill slopes by X- and Y-axes (xi,yi) and spatial steps (dx,dy in m)
c input of the Manning coefficient (sherN)   
c calculation                                   
          grad=sqrt(xi**2.+yi**2.)
          Ax=xi/sherN/sqrt(grad)                             
          Ay=yi/sherN/sqrt(grad)
          h2=0.
c calculation of the surface runoff
      h3=h    
C  translation of mm/3hour into m/s   
        r=(stok+DRAIN)/dt/1000.
        h1=.1                                                      
        kodhyd=0
      do 11 while (abs(h1-h2).ge.e)
         hh2=hh1 
         hh1=hh2+(h-Ax*dt/dx*hh2**(5./3.)+.5*R*dt-hh2)
     *       /(1.+5./3.*Ax*dt/dx*hh2**(2./3.))
         if (hh1.lt.0.) then
              kodhyd=1
              hh1=hh2
              goto 13
         endif    
         hh=hh1
         h2=h1  
         h1=h2+(hh-Ay*dt/dy*h2**(5./3.)+.5*R*dt-h2)
     *      /(1.+5./3.*Ay*dt/dy*h2**(2./3.)) 
c         if(h1.lt.0.and.h1.gt.-e-17) h1=0.
	   if (h1.lt.0.) then
              kodhyd=1
              goto 13
         endif    
  11  continue    
       h=h1     
c balance of surface water 
       bal=r*dx*dy*dt+(h3-h)*dx*dy-Ax*hh**(5./3.)*dy*dt-
     *      Ay*h**(5./3.)*dx*dt                       
c balance in mm/3hour      
       BAL=BAL/DX/DY*1000.
c hydrograf of surface runoff (mm/3-hour)                   
       shydr=(Ax*hh**(5./3.)*dy*dt+Ay*h**(5./3.)*dx*dt)/dx/dy*1000.  
   13  hhyd=h
       hh1hyd=hh1
      return
      end
c////////////////////////////////////////////////////////////////      