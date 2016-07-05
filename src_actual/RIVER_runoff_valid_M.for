	subroutine routing(ncell,m0,iyr_fir,river,u0)
c m0 - number of daily time steps

!
!             PROGRAM OF STREAMFLOW TRANSFORMATION
!      USING RUNOFF DATE OF CELLS AND BASIN ROUTING SCHEME
!                     for grid (1o x 1o)
!   ******************************************************

		IMPLICIT NONE 


      character*5 river
	REAL, ALLOCATABLE :: Y(:)
      REAL, ALLOCATABLE :: Q(:),S(:),Runoff(:)

     	Integer  i,j,i0,i1,j0,j1,step,m,m0,m1,NN,im,iy,id
	Integer ind1,k,k0,mm,ii,ierr,ierr1,ierr2,ierr4,iyr_fir
cccc	Integer, ALLOCATABLE ::  IND2(:,:)
	Real, ALLOCATABLE ::  Din(:), Dout(:)
cccc	Real, ALLOCATABLE ::  area(:,:)

c      Real, ALLOCATABLE ::  run_day(:,:),ec_day(:,:),pr_day(:,:)
	real dt, c, Ct,  lon, lat, ara, x, x1,x2, x3, x4, x5, xx
	real long, lati, long0, lati0, areal, u0
      real run_mmday(35,12,31),run_mmon(35,12) 
	real run_day(1700,12000),ec_day(1700,12000)
	real pr_day(1700,12000)
	integer n1(12000),n2(12000),n3(12000),ilength
	integer lat_cell(1700),lon_cell(1700),icell,ncell,num_cell 
	integer IND2(720,360),ind(1700),ind_lon(1700),ind_lat(1700)
	real area(720,360),areaw(1700),sum_ec,sum_r,sum_p
	integer jlon,ilat

	Character ltt*2, lng*3, s$*12,ltt1*3, ltt2*1
      Character lng4*4, lng2*2, lng1*1, s1$*16
	
c	common /opt/ run_mmday   ! mm/day
	common /opt/ run_mmday,run_mmon    !mm/month
c	common /riv/ i0,i1,j0,j1,long0,lati0,u0
 

c	common /routing/ run_day,n1,n2,n3,lat_cell,lon_cell 
c       common /rout/ lat_cell,lon_cell,n1,n2,n3,run_day 
       common /rout/ ind,areaw,n1,n2,n3,run_day,ind_lon,ind_lat,ec_day,
     !	 pr_day  
c input of the months' length      
      
! ___________________________________________________________________________
!
!	                           WARNING !!!
!                          
!                            G U I D A N C E   !!!
!   BELLOW are the MAIN characteristics of RIVER BASIN and TIME STEPS
!  THESE  CHARACTERISTICS ARE NEEDED TO SET HERE AS INPUT DATA  BEFORE RUNING  !!!
!
!	IN ADDITION the INDEX1.TXT-file and the INPUT RUNOFF-files are needed too
!
!    1.  EXAMPLES of INPUT RUNOFF-file NAME:      _-3___46.txt
!       (they use cell latitude and longitude     -87_-123.txt
!        values and occupy 12 positions)          _34___-2.txt 
!                                                 _45__167.txt
!
!
!    2.  EXAMPLE of structure of INDEX1.TXT-file
!       (lon) (lat) (ind)(area)
!        92.5	 66.5	0	0
!        92.5	 67.5	0	0
!        92.5	 68.5	0	0
!        92.5	 69.5	3	0.436
!        93.5	 66.5	0	0
!        93.5	 67.5	0	0
!        93.5	 68.5	0	0
!        93.5	 69.5	3	0.436
!        94.5	 66.5	0	0
!        94.5	 67.5	0	0
!        94.5	 68.5	0	0
!        94.5	 69.5	2	0.437

!  
!        RESULTS (OUTPUT BASIN STREAMFLOW, mm/time_step) will be in Q-file
                      
cccccccccccccccccc	i0=62         !  
cccccccccccccccccc	i1=65         !  BOUNDARIES (altitude (-60:84), longitude (-180,180))
cccccccccccccccccc	j0=44	
c	j1=135		  !  of the land surface BOX containing the basin
cccccccccccccccccc	j1=50 
c   	j1=137       !
ccccccccccccccccc	long0=45.5   ! longitude of output mouth
ccccccccccccccccc	lati0=65.5    ! latitude of output mouth
c	long0=135.5   ! longitude of output mouth
c	lati0=48.5    ! latitude of output mouth	              !
c	m0=4932       ! total number of time_steps      
	              !
ccccccccccccccccc	u0=0.36       ! effective velocity of water flow in cannals of cells, m/s
	              ! taking into account the MEANDERING RATIO = 1.4
c      u0=.5
c       u0=0.2
     
c      write(*,*) 'routing',i0,i1,j0,j1,u0,long0,lati0
c	stop

      dt=1.         ! time_step (in days)_
!  ___________________________________________________________________
!

cccc	k0=(i1-i0+1)*(j1-j0+1)  !  suprenum of rout length (in number of cells)
c	write(*,*) 'k0',k0	

!          DETERMINATIOM of ARRAYS DIMENSIONS
                  m1=m0+1
      ALLOCATE (y(m0),Q(m0),S(m1),Runoff(m0),Din(m0),Dout(m0),stat=ierr)
c      ALLOCATE (n1(m0),n2(m0), n3(m0), stat=ierr4)
ccccccccccccc	ALLOCATE  (IND2(j0:j1,i0:i1), stat=ierr1)
ccccccccccccc	ALLOCATE  (area(j0:j1,i0:i1), stat=ierr2)
c       ALLOCATE  (run_day(ncell,m0),ec_day(ncell,m0),pr_day(ncell,m0),
c     !	  stat=ierr2)
     
      Open (5,file='d:\ISI-MIP2\res\'//'Q_day')     ! output file for streamflow from the basin
! ************************************************************************
!         preparing indicaters and areas file for current work
	IND2=0.
	area=0.   !OOOOOOOOOOOOOOOOO
      areal=0.
c       write(*,*) i0,i1,j0,j1,long0,lati0,u0

!		Read of indicater file and preparing INDICATER array		
c      Open (4, file='d:\ISI-MIP2\data\'//river//'\index11.txt')    ! input data with longitude,latitude, 
                                  	! index (or mask) and area of cells
c      read(4,*) ilength ! the number of lines in the file
c	do j=j0, j1
c		do i=i0, i1
      i0=360   !minlat
	i1=1     !maxlat
	j0=720   !minlon
	j1=1     !maxlon
c        do i=1,ilength
         do i=1,ncell
c	    	read (4,*) jlon,ilat,long, lati, ind1, ara   ! ara (cell area) in 10^-4 km^2

c  	  		     IND2(long-.5,lati-.5)=ind1
c	             area(long-.5,lati-.5)=ara
                   jlon=ind_lon(i)
				 ilat=ind_lat(i)   
                   IND2(jlon,ilat)=ind(i)
	             area(jlon,ilat)=areaw(i)
	             if (ind(i).eq.9) then
	                 lati0=ilat
                       long0=jlon
                   endif
      
          if (ilat.lt.i0) i0=ilat
          if (ilat.gt.i1) i1=ilat
          if (jlon.lt.j0) j0=jlon
	    if (jlon.gt.j1) j1=jlon
         enddo

c	close (4)
 	k0=(i1-i0+1)*(j1-j0+1)  !  suprenum of rout length (in number of cells)
c	write(*,*) 'k0',k0,i0,i1,j0,j1
c	pause
! ______________________________________________________________
!
!        START of BASIN STREAMFLOW CALCULATION (Runoff-array)
! 
c     	Open (4, file='d:\ISI-MIP2\data\'//river//'\index11.txt')    !  input INDEX1-file with cell data    
c	read(4,*) ilength ! the number of lines in the file
		
		Runoff=0.
				 
c	do j=j0, j1
c		do i=i0, i1
	
c			read (4,*) lon,lat,ind,ara
c              read (4,*) jlon,ilat,long, lati, ind1, ara 
c	                   long=lon
c					   lati=lat
c					   ind1=ind
c					                               ara=1.    ! PROBA
c	if (ind1.gt.0.) then 
       
	 do icell=1,ncell
c	  if(lat_cell(icell).eq.i.and.lon_cell(icell).eq.j) num_cell=icell
c       enddo
         long=ind_lon(icell)
	   lati=ind_lat(icell) 

	   ind1=ind(icell)
	   num_cell=icell
				step=1
	             S=0.
				
!                   Calculation of basin area				
!				areal=areal+area(long-.5,lati-.5)
	            areal=areal+areaw(icell)
				
c	  write(*,*) s$, areal*10000., 'sq.km'

!          Calculation of routing trajectory for selected cell
	       do k=1, k0
c       write(*,*) 'k0',k0,k
!              CALCULATION of CELL CHARACTERISTICS        
c				c=u0*86400./sqrt(area(long-.5,lati-.5))/100000.
              	c=u0*86400./sqrt(area(long,lati))/100000.
  		        Ct=exp(-c*dt)
		
!              Calculation of time set of RUNOFF for selected CELL				  
                 do m=1,m0
				  if(k.eq.1) then

cccc            read (3,*) n1(m), n2(m), n3(m),x1 

              Y(m)=run_day(num_cell,m)    ! Y - cell runoff in mm/time_step        

c		  Din(m)=Y(m)/1000.*area(long-.5,lati-.5)*10000.*1000000.
            Din(m)=Y(m)/1000.*area(long,lati)*10000.*1000000.
				  else 
				    Din(m)=Dout(m)   	
	              endif

 				        S(m+1)=Ct*S(m)+Din(m)/c*(1-Ct)
				        Dout(m)=Din(m)-(S(m+1)-S(m))/dt
		       enddo !m  end of calculation for one cell


!              EXIT from basin
			    	if (lati.eq.lati0.and.long.eq.long0) then 
                     		Q=Dout
						exit
					end if

!         Selection of direction in routing scheme
! схема работает для случая, когда индексация широт возрастает к югу (lati=1 на севере)
c	   write(*,*) icell,step,ind1,long,lati
c	pause
	  select case (ind1)
		case (1)
		  long=long
c	      lati=lati+1
            lati=lati-1
			case (2)		
		  long=long+1
c		  lati=lati+1   
            lati=lati-1        
		case (3)        !ok
		  long=long+1
		  lati=lati  	    
		case (4)
		  long=long+1
c		  lati=lati-1  
	      lati=lati+1  
		case (5)
		  long=long
c		  lati=lati-1  
		  lati=lati+1  
		case (6)		
 		  long=long-1
c		  lati=lati-1  
            lati=lati+1         
		case (7)         !ok
		  long=long-1
		  lati=lati  	    
		case (8)
		  long=long-1
c		  lati=lati+1  
            lati=lati-1 
          end select
     

c         ind1=IND2(long-.5,lati-.5)
          ind1=IND2(long,lati)
      
!      SUMMARIZING OF NUMBER OF STEPS in ROUTING TRAJECTORY 
!      from SELECTED CELL to Basin EXIT cell
	   step=step+1

  		end do ! k  - end of output Streamflow calculation 
                     ! for one way  

!      SUMMARIZING of  OUTPUT streamflow for EXIT CELL
		do mm=step,m0
		Runoff(mm)=Q(mm-step+1)+Runoff(mm)     !m^3/dt
		enddo !mm			

c	end if    !  end of output Streamflow calculation 
                !  for all previous CELLS  
c		close (3)

c		end do !   i - end of calculation on latitude
c	end do     !   j - end of calculation on longitude	
      enddo  !  end of calculation for all cells
	
			
!	                                TOTAL BASIN STREAMLOW mm/dt
		Do mm=1, m0
       iy=n1(mm)
	 im=n2(mm)
	 id=n3(mm)
       iy=iy-iyr_fir+1
       run_mmday(iy,im,id)=Runoff(mm)/areal/10000./1000000.*1000.       !  mm/dt
	Write (5,16)  n1(mm),im,id,run_mmday(iy,im,id)   !  mm/dt
    		enddo
	write(*,*) areal		
	close (4)
	close(5)
   16           format (3i6,e12.3)   
   
   
c             CALCULATION of MONTHLY STREAMFLOW, mm/month

	   Open (5,file='d:\ISI-MIP2\res\'//'Q_mon') 
	   	
		NN=0
		xx=0.
	
	do mm=1, m0
         iy=n1(mm)
         iy=iy-iyr_fir+1
	   im=n2(mm)
         
		if (NN.eq.im.or.mm.eq.1) then
	         xx=xx+Runoff(mm)/areal/10000./1000000.*1000. !  mm/month
	    else
			write (5, *) n1(mm-1), n2(mm-1), xx
                run_mmon((n1(mm-1)-iyr_fir+1),n2(mm-1))=xx
			xx=Runoff(mm)/areal/10000./1000000.*1000. !  mm/month
		end if
            NN=im
	end do
	    
		 run_mmon(iy,im)=xx
			write (5, *) n1(m0), n2(m0), xx
	    close (5)


c             CALCULATION of ANNUAL STREAMFLOW, mm/year

	open (5, file='d:\ISI-MIP2\res\year')
	    NN=0
		xx=0.
      do mm=1, m0
		if (NN.eq.n1(mm).or.mm.eq.1) then
	         xx=xx+Runoff(mm)/areal/10000./1000000.*1000. !  mm/y
	    else
			write (5, *) n1(mm-1),  xx
			xx=Runoff(mm)/areal/10000./1000000.*1000. !  mm/y
		end if
            NN=n1(mm)
	end do
			write (5, *) n1(m0), xx
	    close (5)
c BASIN BALANCE
      open (5, file='d:\ISI-MIP2\res\balance')
        
	  do m=1,m0
         sum_ec=0.
	   sum_r=0.
	   sum_p=0.
	     do i=1,ncell
	        sum_ec=sum_ec+ec_day(i,m)*areaw(i)/areal 
              sum_r=sum_r+run_day(i,m)*areaw(i)/areal
              sum_p=sum_p+pr_day(i,m)*areaw(i)/areal
           enddo
           write(5,111) sum_ec,sum_r,sum_p
        enddo
	close(5)
  111 format(3e14.4)   
    
	     end 