c made from gswpmop.for
C              
C                   OPTIMIZATION OF PARAMETERS
C
C *******************************************************
      module dflib
      end module dflib
C 

C *************************************************************************************	

      program OPTIMIZATION
	 
       USE DFLIB
       implicit none
C
C ****************************************************************************
C      WARNING !!!
C
CCC	CHANGE in NEXT LINES  (***) amount of measured data y_i, output data a,b,c ...
CCC   and optmized parameters par_i
CCC   in accordance with experiment of optimization parameters
CCC   CHANGE NAME of main MODEL in accordance with problem too
C
C ****************************************************************************
C
	character*1 ichar1
      character*2 ichar2 
      character*5 river
	character*3 riv
	REAL, ALLOCATABLE ::   ab(:,:),S(:,:)  !   ***
	REAL, ALLOCATABLE :: RR(:),pmin(:),pmax(:),lon(:),lat(:)     ! according to number of parameters
	REAL ran,yy,yyy
	real nb,wzav,por,bpar,kf,fi0,amt,sryt,albzim,prles,sai, SherN
      real leaf,emkh,xi,yi,hroot,h0,hveg,elevat,ksiot0,a(20000),par(100)
	real long0,lati0,u0,albsnow
      integer minlat,maxlat,minlon,maxlon,i0m,i0_obs,i0m_obs,ny,nm,nd
      real run_mmday(35,12,31),run_mmon(35,12) 
      real run_obsd(35,12,31),run_obsm(35,12),y1,mn(366),n0,a_mn,b_mn,t  
	real bias,effect,sig_obs,effect_max,eff(5),bia(5),x20,bias_max 
	real sr_obs

c	integer,ALLOCATABLE ::   ncl(:)
	integer i,ii,j,i0,m, k,k0,IERR, m0,land,ico,jco,npar,ncell,iii
	integer kod_data,ik,iyr_fir,kod_exp, kod_prec,kod_cal,kod_sens
	integer iyrfir,num_year,firday,kod_calibr,kmax,j00,ncl,isen,ip

      common /opt_par/ par,kod_exp, kod_prec,albsnow,iyrfir,num_year,
     !	firday,kod_calibr,kod_cal
	common /riv/ minlat,maxlat,minlon,maxlon,long0,lati0,u0
  	common /opt/ run_mmday,run_mmon    !mm/month
c	common /manning/ SherN
      common /manning/ mn              
     
C        PARAMETER RANGE		
c  Parameters for a river:
      npar=15
c	npar=33       ! the number of optimizing parameters              (Mezen)
cccccccccccccccc      ncell=10        ! the number of grid cells       (Mezen)
c       ncell=57     ! the number of grid cells                         (Pechora)
c       ncell=3
cccccccccccccccccc      long0=46.5      ! longitude of output mouth   (Mezen)
cccccccccccccccccc	lati0=64.5      ! latitude of output mouth      (Mezen)
c          long0=52.5             ! longitude of output mouth   (Pechora)
c          lati0=67.5      ! latitude of output mouth      (Pechora)
      u0=0.36         ! effective velocity of water flow in cannals of cells, m/s
      n0=0.07         ! Manning coefficient
	albsnow=0.75 
c      i0=4932         ! number of daily time steps in 14 years (July 1982- December 1995)   (Mezen)
c      i0=4567         ! number of daily time steps in 13 years (July 1982- December 1994)
ccc      i0_obs=4748     ! number of days with observed runoff 14 years   (Mezen)
c      i0_obs=8254    ! number of days with observed runoff  (1980-2003)   (Pechora)
c      i0_obs=4348     ! number of days with observed runoff 13 years
c	i0m=162         ! number of months in 14 years (July 1982- December 1995)
	i0m_obs=156     ! number of months with observed runoff   ***********  for Mezen
c      iyr_fir=1982     ! the first year of model simulation

	m0=5                             ! number of points in a smoothing group
	k0=m0*200                       ! number of criterium calculation
c       m0=1
c       k0=m0*10
	
	write(*,*) 'Input the name of the river'
	write(*,*) '               mezen'
      write(*,*) '               pechr'
      write(*,*) '               sevdv'
      write(*,*) '               olenk'
      write(*,*) '               indig'
	write(*,*) '               kolym'
	write(*,*) '               onega'
      write(*,*) '               yanaa'
      write(*,*) '               ponoi'
      write(*,*) '               tulom'

      read (*,*) river 
	 riv=river
	 write(*,*) riv,' ',river
	
	open (21, file='c:\d\rivers\data\'//river//'\info.'//riv) 
      read(21,*) ncell,long0,lati0,i0,i0_obs,iyr_fir,iyrfir,num_year,
     !	firday,sig_obs,sr_obs
      write(*,*) ncell,long0,lati0,i0,i0_obs,iyr_fir,iyrfir,num_year,
     !	firday,sig_obs,sr_obs
      

	ALLOCATE (ab(k0,npar+3), S(k0-m0,npar+2),      !  ***
     !   RR(npar+1),pmin(npar),pmax(npar),lon(ncell),lat(ncell),
c     !   ncl(ncell),stat=ierr)                                                
     !   stat=ierr)  

c	OPEN (20, file='c:\d\rivers\res\mn.mez') 


      write(*,*) 'Which precipitation will be used: 1 - GSWP'
	write(*,*) '                                  2 - Measured'
	read(*,*) kod_prec
      
	write(*,*) kod_prec

      write(*,*) 'Which data will be used to run the model: 1 - B0-exper
     !iment'
	write(*,*) '                                          2 - P3-exper
     !iment'
      read(*,*) kod_exp
      
	kod_sens=2
      kod_cal=0 
      write(*,*) 'Do you want to calibrate the model: yes - (1)'
	write(*,*) '                               or  not  - (2) ?'
	read(*,*) kod_calibr
		
	if (kod_calibr.ne.1) then
	write(*,*) 'Do you want to calculate with a priori parameters (1)'
	write(*,*) '                 or  with calibrated parameters (2) ?'
	read (*,*) kod_cal 
      if (kod_cal.eq.2) then
	    write(*,*) 'Do you want to make a sensitivity test: yes - (1)'
      	write(*,*) '                               or  not  - (2) ?'
	    read(*,*) kod_sens
       endif
      endif
c	 stop 

      OPEN (7,file='c:\d\rivers\fixed_fields\fixed_param_'//riv//'.txt') ! for reading lat and lon
	
      do i=1,ncell
	 read (7,*) ncl,land,ico,jco,lon(i),lat(i),nb,wzav,por,bpar,kf,
     !	 fi0,amt,sryt,albzim,prles,sai,leaf,emkh,xi,yi,hroot,h0,hveg,
     !     elevat,ksiot0
      enddo
      close(7)

      minlat=lat(1)
c	minlat=int(lat(1))
	maxlat=lat(1)
	minlon=lon(1)
	maxlon=lon(1)
      do i=1,ncell
          if (lat(i).lt.minlat) minlat=lat(i)
          if (lat(i).gt.maxlat) maxlat=lat(i)
	    if (lon(i).lt.minlon) minlon=lon(i)
          if (lon(i).gt.maxlon) maxlon=lon(i)
	enddo
c	minlon=minlon-1


	if(kod_calibr.eq.2) then

	  if(kod_cal.eq.2) then
	     OPEN (9,file='c:\d\rivers\res\opt_par.csv')   ! reading of optimal parameters
! calibrated parameters
           do i=1,npar
       	     read(9,*) par(i)
           enddo
	     close(9)
	u0=par(13)   !	   RIVER VELOCITY
          
	  else
! list of a priori parameters
c       par(1)=0.       ! k_k0
c       par(2)=1.       ! k_hroot     
c       par(3)=1.       ! k_nb
c       par(4)=1.       ! k_alb_zim  
       par(5)=albsnow
c       par(6)=1.       ! k_alb_leto 
c       par(7)=1.       ! k_sw
c       par(8)=1.       ! k_lw
c       par(9)=1.       ! h0/hroot   
       par(10)=n0     ! Manning
       par(11)=0.     ! a_mn
       par(12)=0.     ! b_mn	
c       par(13)=u0     !   u0, m/s
c       par(14)=1.     ! coefficient for rainfall
c       par(15)=1.     ! coefficient for snowfall
	  endif   
	
	else
       write(*,*) 'Which data will be used for calibration: 1 - daily'
	 write(*,*) '                                         2 - monthly'
	 read(*,*) kod_data
       if (kod_data.eq.1) then
      	OPEN (11, file='c:\d\rivers\data\'//river//'\obs_run_day.'//
     ! 		riv)      
      	run_obsd=-99.
		do i=1,i0_obs
            read (11,*) ny,nm,nd,y1                           ! measured data   ***
	      if ((ny.gt.iyr_fir).or.(ny.eq.iyr_fir.and.nm.ge.7)) then
               
              ny=ny-iyr_fir+1
      	    run_obsd(ny,nm,nd)=y1
	        
            endif
		enddo
      	close (11)
       else
	  	OPEN (11, file='c:\d\rivers\data\'//river//'\obs_run_mon.'//
     !		riv)     
      	run_obsm=-99.
		do i=1,i0m_obs
            read (11,*) ny,nm,y1 
		  if (ny.ge.iyr_fir) then                          ! measured data   ***
               ny=ny-iyr_fir+1
		     run_obsm(ny,nm)=y1
            endif
	 	enddo
      	close (11) 
       endif
         pmin(1)=-8.            !k_k0
	  pmax(1)=-7.

       pmin(2)=1.     ! k_hroot
	  pmax(2)=1.3

	pmin(3)=1.                !  k_nb
	  pmax(3)=1.

	pmin(4)=1.                  ! k_alb_zim
	  pmax(4)=1.

	pmin(5)=0.8	              ! albsnow
	  pmax(5)=0.85

	pmin(6)=0.99       ! k_alb_leto    
	  pmax(6)=1.15

	pmin(7)=1.02        ! k_sw
	  pmax(7)=1.03

	pmin(8)=1.01         ! k_lw
	  pmax(8)=1.03

	pmin(9)=1.9            ! h0/hroot
	  pmax(9)=2.5

	pmin(10)=0.1        ! Manning
	  pmax(10)=0.25

	pmin(11)=0.0                  ! a_mn
	  pmax(11)=0.0

	pmin(12)=0.0                  ! b_mn
	  pmax(12)=0.0

	pmin(13)=1.2            ! u0
	  pmax(13)=2.

	pmin(14)=1.1  ! coefficient for rainfall
	  pmax(14)=1.3

	pmin(15)=0.75          ! coefficient for snowfall
	  pmax(15)=0.8


	open (20, file='c:\d\rivers\res\ab')
	endif
	
	if (kod_sens.eq.1) goto 719
      if (kod_calibr.eq.2) goto 717

C          MONTE CARLO technique to set PARAMETERS VALUES      
	CALL RANDOM_seed 
	
	Do k=1,k0
        write(*,*) 'k',k
          
                          CALL RANDOM_NUMBER (ran)   ! for the kf
               par(1)=(pmin(1)+(pmax(1)-pmin(1))*ran)
c             	write(*,*) par(i)
c	        pause
         	  
		do i=2,npar    ! for the rest parameters
             CALL RANDOM_NUMBER (ran) 	
		   par(i)=pmin(i)+(pmax(i)-pmin(i))*ran
c			write(*,*) par(i)
c	        pause		  
	    enddo	
	

  	u0=par(13)   !	   RIVER VELOCITY
  717 continue   
	  

c calculation of Manning coef. with accounting for its annual course
	n0=par(10)
	a_mn=par(11)
	b_mn= par(12)
	
        do i=61,365
	     t=-15.37*sin(2*3.1415926*(i+77.64)/365.)-0.7
            if (i.gt.60.and.i.le.180) then 
		    if(t.gt.b_mn) then 
	          mn(i)=n0
	        else
	          mn(i)=n0*exp(-a_mn*(t-b_mn))
              endif
            else
	        if(t.gt.(-b_mn)) then
	          mn(i)=n0
	        else
	          mn(i)=n0*exp(-a_mn*(t+b_mn))
              endif
            endif
         enddo

	   do i=1,60
	      mn(i)=mn(365)*(60-i)/59.+mn(61)*(i-1)/59.
         enddo
         mn(366)=mn(365)
      
                         	
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      goto 10
  
 	CALL model_plus (npar,ncell,river)                                      !      ***

	CALL model_minus (npar,ncell,river)                                         !      ***
  203     	call transform(ncell,river,i0)
c	call test(i0,iyr_fir)
c   10      write (*,*) i0
      call routing(i0,iyr_fir,river)
c        call routing(12784,iyr_fir,river)                               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         if (kod_calibr.eq.2) goto 818
       open (17, file='c:\d\rivers\res\rrr') 
	    ik=0
	   yy=0.
	   bias=0.
         do i=1,num_year
	     do j=1,12
	        if (kod_data.eq.1) then
	            do ii=1,31
	  
	              if(run_obsd(i,j,ii).gt.-99.) then
				    yy=yy+(run_mmday(i,j,ii)-run_obsd(i,j,ii))**2
	                bias=bias+run_mmday(i,j,ii)-run_obsd(i,j,ii)  
c	                write(17,400) i,j,ii,run_obsd(i,j,ii),
c     !					run_mmday(i,j,ii)
	                ik=ik+1
                    endif
	            enddo
	         else  
			   if (run_obsm(i,j).gt.-99.) then
			      yy=yy+(run_obsm(i,j)-run_mmon(i,j))**2
	              bias=bias+run_mmon(i,j)-run_obsm(i,j)
                    ik=ik+1
c				  write(17,400) i,j,ik,run_obsm(i,j),run_mmon(i,j)
                 endif
	        endif
           enddo
         enddo
c        write (*,*) 'ik ',ik,yy
         
c      stop			    
  400  format(3i5,2f10.3)	          

C             Calculation of criterium of optimization
	
C         OUTPUT of PARAMETERS and CRITERIUM	
       bias=bias/ik/sr_obs*100.
       ab(k,npar+2)=sqrt(yy/ik)
	 ab(k,npar+1)=bias
	 effect=1-(yy/ik)/sig_obs**2
       ab(k,npar+3)=effect

c		Selection of minimal effectivity
		if (k.eq.1) then
			effect_max=effect
	        bias_max=bias
			kmax=1

	    else
	        if (effect.gt.effect_max) then  
			  	effect_max=effect
				bias_max=bias
                  kmax=k
	        endif
	    endif
		Write (*,*) 'effect_max=', effect_max, 'kmax=', kmax
	Write (*,*) 'bias =', bias_max

      do i=1,npar
	   ab(k,i)=par(i)

                 
	enddo
c                              stop                         

	j00=npar+3
	write (20,'(41E12.4)') (ab(k,j), j=1,j00)      

	enddo      ! k

****************************************************
c	open (20, file='c:\d\rivers\res\ab')
	
c	do i=1,k0
c	      j00=npar+2
		  
c	write (20,'(41E12.4)') (ab(i,j), j=1,j00)
	    
c	enddo
*******************************************************

CC	Close(2)
                 WRITE  (*,*)  '     END of MODEL RUN'
	           Write  (*,*)  '       SORT ...'
C *******************************************************************
C
C             Smoothing
C
C	             SORT
 	
		do i=1,k0
		do k=1,k0-1
	   
	if (ab(k+1,npar+2).lt.ab(k,npar+2)) then
			yyy=ab(k,npar+2)
	    ab(k,npar+2)=ab(k+1,npar+2)
	      ab(k+1,npar+2)=yyy
      	do iii=1,npar+1
	       yyy=ab(k,iii)
	        ab(k,iii)=ab(k+1,iii)
	         ab(k+1,iii)=yyy
          enddo
    	endif

	enddo
	enddo	
CC	close (2)
                  
				WRITE (*,*) '     END of SORT'
	            WRITE (*,*)  '     SMOOTHING  ...'
				  	
C                 SMOOTHING
C
        do k=1,k0-m0
		 do iii=1,npar+2	
			S(k,iii)=0.
		 enddo
                         
	     do i=1,m0
	       do iii=1,npar+2
			S(k,iii)=S(k,iii)+ab(k+i-1,iii)/m0
		   enddo
	     enddo
        enddo 
                   
				   WRITE (*,*), '    END of SMOOTHING'
	               WRITE (*,*)   '  '
	               WRITE  (*,*)  '    RESULTS:'
C                  WRITING of OPTIMAL PARAMETERS VALUES
	write (*,*) (S(1,iii),iii=1,npar)

	                                   ! ****
      Write (*,*) 'min yy=',S(1,npar+2), 'bias=',S(1,npar+1) 

C                   Calculation of parameters errors
	
C	RR(1)=0.
CC	RR(2)=0.
C	RR(3)=0.
C	do k=1,m0
C	  RR(1)=RR(1)+(ab(k,1)-S(1,1))**2./(m0-1)/m0
C	  RR(2)=RR(2)+(ab(k,2)-S(1,2))**2./(m0-1)/m0
C	  RR(3)=RR(1)+(ab(k,3)-S(1,3))**2./(m0-1)/m0
C	enddo
C	rr(1)=sqrt(rr(1))
C	rr(2)=sqrt(rr(2))
C	rr(3)=sqrt(rr(3))

C      write (*,*) 'sig1=',rr(1),'  ','sig2=',rr(2),'  ','sig3=',rr(3)

c	open (3, file='c:\d\rivers\res\optima_sort')	
c	  write(3,*) 'k0 (m/s)    hroot(m)'
c	  do iii=1,npar+2
c           write (3,*) S(1,iii),S(1,iii+ncell) 
c	  enddo  
	 
c	 write(3,*) 'k_rainfall  k_snowfall'                                ! ****
c       write(3,*) S(1,npar-1),S(1,npar)

c	 write(3,*) 'b_mn     u0_m/s '                                 ! ****
c       write(3,*) S(1,npar-3),S(1,npar-2)
c       write(3,*) '  a_mn   n0  h0/hroot'  
c	 write(3,*) S(1,npar-4),S(1,npar-5),S(1,npar-6)
c	 write(3,*) ' k_lw  k_sw'
c	 write(3,*) S(1,npar-7),S(1,npar-8)
c	 write(3,*) 'kalblet      albsnow       albzim'                  ! ****
c      write(3,*) S(1,npar-9),S(1,npar-10),S(1,npar-11)
c      write(3,*) 'koeff_nb' 
c       write(3,*)  S(1, npar-12)
	 	  
c	Write (3,*) 'min yy=',S(1,npar+2)
c	Write (3,*)  'bias=', S(1,npar+1)
c	do i=1,npar
c	  write(3,*) S(1,i)
c      enddo
       
	open (4, file='c:\d\rivers\res\all')
	
	do i=1,k0
	      j00=npar+2
		  
	write (4,'(17E12.4)') (ab(i,j), j=1,j00)
	    
	enddo

c    	close (3)
      close (4)	
	close  (20)
	close (17)
      goto 818

c sensitivity tests
  719 continue
c      write(*,*) 'rrrr sensit'
	                                           
      OPEN (35, file='c:\d\rivers\res\sensit')  
	open (17, file='c:\d\rivers\res\rrr') 
      OPEN (11, file='c:\d\rivers\data\'//river//'\obs_run_day.'//
     ! 		riv) 
c        
              
          
      	run_obsd=-99.
c	write(*,*) i0_obs
	
		do i=1,i0_obs
            read (11,*) ny,nm,nd,y1                           ! measured data   ***
	      if ((ny.gt.iyr_fir).or.(ny.eq.iyr_fir.and.nm.ge.7)) then
               
              ny=ny-iyr_fir+1
      	    run_obsd(ny,nm,nd)=y1
	        
            endif
		enddo
      	close (11)
      
	do ip=1,npar
	if (ip.eq.11.or.ip.eq.12) goto 900
      
	write(*,*) '*********** parameter********', ip
	
	   x20=0.8
        if (ip.eq.1) x20=-0.2

	   eff=0.
         bia=0.
	   do isen=1,5
	write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&', ip,isen   
	      	  
		  if (ip.eq.1) then
               par(ip)=par(ip)+x20*LOG(10.)
	      else
	         par(ip)=par(ip)*x20
            endif
	Write (*,*) par
      u0=par(13)
c calculation of Manning coef. with accounting for its annual course
	n0=par(10)
	a_mn=par(11)
	b_mn= par(12)
	
        do i=61,365
	     t=-15.37*sin(2*3.1415926*(i+77.64)/365.)-0.7
            if (i.gt.60.and.i.le.180) then 
		    if(t.gt.b_mn) then 
	          mn(i)=n0
	        else
	          mn(i)=n0*exp(-a_mn*(t-b_mn))
              endif
            else
	        if(t.gt.(-b_mn)) then
	          mn(i)=n0
	        else
	          mn(i)=n0*exp(-a_mn*(t+b_mn))
              endif
            endif
         enddo

	   do i=1,60
	      mn(i)=mn(365)*(60-i)/59.+mn(61)*(i-1)/59.
         enddo
         mn(366)=mn(365)
c        write(*,*) par                 	
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c      goto 10

c       write (*,*) par



   	CALL model_plus (npar,ncell,river)                                      !      ***

	CALL model_minus (npar,ncell,river)                                         !      ***
     	call transform(ncell,river,i0)
c	call test(i0,iyr_fir)
c   10      write (*,*) i0
      call routing(i0,iyr_fir,river)
c        call routing(12784,iyr_fir,river)                               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ik=0
	 yy=0.
	 bias=0.  
         do i=1,num_year
	     do j=1,12
c	        if (kod_data.eq.1) then
	            do ii=1,31
	  
	              if(run_obsd(i,j,ii).gt.-99.) then
				    yy=yy+(run_mmday(i,j,ii)-run_obsd(i,j,ii))**2
	                bias=bias+run_mmday(i,j,ii)-run_obsd(i,j,ii)  
c	                write(17,400) i,j,ii,run_obsd(i,j,ii),
c     !					run_mmday(i,j,ii)
	                ik=ik+1
                    endif
	            enddo
c	         else  
c			   if (run_obsm(i,j).gt.-99.) then
c			      yy=yy+(run_obsm(i,j)-run_mmon(i,j))**2
c	              bias=bias+run_mmon(i,j)-run_obsm(i,j)
c                    ik=ik+1
c				  write(17,400) i,j,ik,run_obsm(i,j),run_mmon(i,j)
c                 endif
c	        endif
           enddo
         enddo

	eff(isen)=1-yy/ik/sig_obs**2
	bia(isen)=bias/ik/sr_obs*100.

	      if (ip.eq.1) then
               par(ip)=par(ip)-x20*LOG(10.)
	      else
	         par(ip)=par(ip)/x20
            endif
            x20=x20+0.1

c	write(*,*) ip,isen, eff(1),bia(1)
c        stop
       enddo  !is
c      write (*,*) ip,isen,(eff(i)*100,i=1,5),(bia(i),i=1,5)
	write (35,100) (eff(i)*100,i=1,5),(bia(i),i=1,5)
  900 continue
	enddo   !ip
  100 format(5F10.3,5f12.3)
      close(35)
	close (17)
  818	continue
  819	WRITE (*,*) '      E  N  D'
      stop
	end