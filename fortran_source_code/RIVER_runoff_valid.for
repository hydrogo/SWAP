!
!             PROGRAM OF RIVER RUNOFF WITH TRANSFORMATION
!       USING RUNOFF DATA OF CELLS AND BASIN ROUTING SCHEME
!                     for grid (1o x 1o)
!      ******************************************************
      
	
	
      subroutine transform(ncell_tot,river,i0)
c the program for rewriting files point(land).txt into _lat_lon.txt 
c    (in so doing the day of year are re-calculating  into number of month and day of month)
		implicit none

	real dir, s, lon, lat
	character  Ncell*5,river*5 
	real long, lati, long0, lati0, areal, u0
	Character ltt*2, lng*3, s$*12,ltt1*3, ltt2*1
      Character lng4*4, lng2*2, lng1*1, s1$*16
	character f1*15, f2*12, A*35
	integer i,j,n,m, ii,nmon,im,ncell_tot,i0

c     number of output characteristics within point-files
	real A1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,MON(12),a33,a14                     
	           !!!
	real a15,a16,a17
      DATA MON/31,28,31,30,31,30,31,31,30,31,30,31/

      OPEN (1, file='c:\d\rivers\data\'//river//'\index2.txt')
	write(*,*) 'ncell_tot',ncell_tot
c	stop
c      number of cells	 
	do ii=1,ncell_tot                                                 !!!
	read (1,*) lon, lat, dir, s, Ncell

c      lat=75.5
c	lon=66.5

	i=lat-.5
	j=lon-0.5

c        TANSFORMATION of LATITUDE	
c	read(lat, '(I3)') m


		if (i.lt.-9) then 
			  write (ltt1, '(i3)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$=ltt1//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$=ltt1//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$=ltt1//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then 
			         write (lng1, '(i1)') j
                       s$=ltt1//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then
			         write (lng2, '(i2)') j
	                 s$=ltt1//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$=ltt1//'__'//lng//'.txt'
	                 goto 17
          endif

    	    if (i.ge.-9.and.i.lt.0) then
			  write (ltt, '(i2)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='_'//ltt//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$='_'//ltt//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='_'//ltt//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then
			         write (lng2, '(i2)') j
	                 s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$='_'//ltt//'__'//lng//'.txt'
	                 goto 17
          endif
						
          if (i.ge.0.and.i.lt.10) then 
		     write (ltt2, '(i1)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='__'//ltt2//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$='__'//ltt2//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$='__'//ltt2//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='__'//ltt2//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then 
			         write (lng2, '(i2)') j
	                 s$='__'//ltt2//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$='__'//ltt2//'__'//lng//'.txt'
	                 goto 17
          endif
               

	          write (ltt, '(i2)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='_'//ltt//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then 
			          write (lng, '(i3)') j
                        s$='_'//ltt//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then 
			         write (lng2, '(i2)') j
                       s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='_'//ltt//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then 
			         write (lng2, '(i2)') j
	                 s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$='_'//ltt//'__'//lng//'.txt'
	                 goto 17

  17	f1='Point'//' '//Ncell//'.txt'
      f2=s$

cccccccccccccccccccccccc       WRITE(*,*)F1
cccccccccccccccccccccccc	 WRITE(*,*)F2
c     stop

	open (2,file='c:\d\rivers\res\'//f2)
	open (3, file='c:\d\rivers\res\'//f1)

      read (3,*) A
cccccccccccccccccccccccccccccccc	write(*,*) A


c	number of time steps
	 do j=1,i0                                               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1645 
      read (3,*) A1,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17

	  	mon(2)=28
	if (abs(a1/4.-int(a1/4.)).eq.0.) MON(2)=29
      nmon=1   
	do im=1,12
          a3=a3-mon(im)
		if (a3.gt.0) then
	       nmon=nmon+1
          else
	       a2=nmon
	       a33=a3+mon(im) 
		   goto 999
	    endif
      enddo
  999 continue
	                    !!! 
	Write(2,'(3f6.0,14E14.4)') A1,a2,a33,a4,a5,a6,a7,a8,a9,a10,a11,a12
     !	,a13,a14,a15,a16,a17             !!!
	
	 enddo !i0
c	write (*,*) ii	    
	close (2)
	close (3)  
	end do   ! ii
c	write (*,*) 'end streamflow'
	end 

	subroutine routing(m0,iyr_fir,river)
c m0 - number of daily time steps

!
!             PROGRAM OF STREAMFLOW TRANSFORMATION
!      USING RUNOFF DATE OF CELLS AND BASIN ROUTING SCHEME
!                     for grid (1o x 1o)
!   ******************************************************

		IMPLICIT NONE 


      character*5 river
	REAL, ALLOCATABLE :: Y(:)
      REAL, ALLOCATABLE :: Q(:),S(:),Runoff(:),aRunoff1(:),
     !	Runoff1(:),Evap(:),Prec(:),sw(:),rdown(:),t2(:),u2(:),q2(:), 
     ! swe(:),ksi(:),ksiot(:),hsnow(:) 
	Integer  i,j,i0,i1,j0,j1,ind,step,m,m0,m1,NN,im,iy,id,iyy,imm,idd
	Integer ind1,k,k0,mm,ii,ierr,ierr1,ierr2, ierr3, ierr4,iyr_fir
	Integer, ALLOCATABLE ::  IND2(:,:)
	Real, ALLOCATABLE ::  Din(:), Dout(:)
	Real, ALLOCATABLE ::  area(:,:)
	real dt, c, Ct,  lon, lat,ara,x,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,xx
	real x11,x12,x13,x14
	integer, ALLOCATABLE :: n1(:), n2(:), n3(:)
	real long, lati, long0, lati0, areal, u0
      real run_mmday(35,12,31),run_mmon(35,12)  

	Character ltt*2, lng*3, s$*12,ltt1*3, ltt2*1
      Character lng4*4, lng2*2, lng1*1, s1$*16
	
c	common /opt/ run_mmday   ! mm/day
	common /opt/ run_mmday,run_mmon    !mm/month
	common /riv/ i0,i1,j0,j1,long0,lati0,u0
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

	k0=(i1-i0+1)*(j1-j0+1)  !  suprenum of rout length (in number of cells)
c	write(*,*) 'k0',k0	

!          DETERMINATIOM of ARRAYS DIMENSIONS
                  m1=m0+1
      ALLOCATE (y(m0),Q(m0),S(m1),Runoff(m0),Din(m0),Dout(m0),stat=ierr)
      ALLOCATE (n1(m0),n2(m0), n3(m0), stat=ierr4)
	ALLOCATE  (IND2(j0:j1,i0:i1), stat=ierr1)
	ALLOCATE  (area(j0:j1,i0:i1), stat=ierr2)
      ALLOCATE  (Runoff1(m0),aRunoff1(m0),Evap(m0),sw(m0),rdown(m0),
     !	t2(m0),u2(m0),q2(m0),swe(m0),ksi(m0),ksiot(m0),hsnow(m0), 
     !Prec(m0), stat=ierr3)   
! ************************************************************************

	    
		Open (5,file='c:\d\rivers\res\'//'Q_day')     ! output file for streamflow from the basin
c      write(*,*) 'u0',u0
c	stop 

!         preparing indicaters and areas file for current work
	IND2=0.
	area=0.   !OOOOOOOOOOOOOOOOO
      areal=0.
       write(*,*) i0,i1,j0,j1,long0,lati0,u0

!		Read of indicater file and preparing INDICATER array		
       Open (4, file='c:\d\rivers\data\'//river//'\index1.txt')    ! input data with longitude,latitude, 
                                  	! index (or mask) and area of cells
	do j=j0, j1
c	write(*,*) j
		do i=i0, i1	
c		write(*,*) i	
			read (4,*) long, lati, ind1, ara   ! ara (cell area) in 10^4 km^2
c	write(*,*) i,j,long, lati, ind1, ara
	              
c					                               ara=1.    ! PROBA				  
				  IND2(long-.5,lati-.5)=ind1
	              area(long-.5,lati-.5)=ara

		enddo
	enddo
	        close (4)

      
!     
! ______________________________________________________________
!
!        START of BASIN STREAMFLOW CALCULATION (Runoff-array)
! 
     	Open (4, file='c:\d\rivers\data\'//river//'\index1.txt')    !  input INDEX1-file with cell data    

		Runoff=0.
	
c	WRITE (*,*) 'j0, j1,i0, i1' 
						 
	do j=j0, j1
		do i=i0, i1
	
			read (4,*) lon,lat,ind,ara
	                   long=lon
					   lati=lat
					   ind1=ind
c          write(*,*) j,i,lon,lat,ind,ara
c					                               ara=1.    ! PROBA
	if (ind.gt.0.) then 


!              PREPARARION OF FILE-NAME for selected CELL
       	   
		if (i.lt.-9) then 
			  write (ltt1, '(i3)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$=ltt1//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$=ltt1//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$=ltt1//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then 
			         write (lng1, '(i1)') j
                       s$=ltt1//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then
			         write (lng2, '(i2)') j
	                 s$=ltt1//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$=ltt1//'__'//lng//'.txt'
	                 goto 17
          endif

    	    if (i.ge.-9.and.i.lt.0) then
			  write (ltt, '(i2)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='_'//ltt//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$='_'//ltt//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='_'//ltt//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then
			         write (lng2, '(i2)') j
	                 s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$='_'//ltt//'__'//lng//'.txt'
	                 goto 17
          endif
						
          if (i.ge.0.and.i.lt.10) then 
		     write (ltt2, '(i1)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='__'//ltt2//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$='__'//ltt2//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$='__'//ltt2//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='__'//ltt2//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then 
			         write (lng2, '(i2)') j
	                 s$='__'//ltt2//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$='__'//ltt2//'__'//lng//'.txt'
	                 goto 17
          endif
               

	          write (ltt, '(i2)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='_'//ltt//'_'//lng4//'.txt'
	                 goto 17
               endif 
               if(j.ge.-99.and.j.lt.-9) then 
			          write (lng, '(i3)') j
                        s$='_'//ltt//'__'//lng//'.txt'
	                  goto 17
	         endif
               if(j.ge.-9.and.j.lt.0) then 
			         write (lng2, '(i2)') j
                       s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='_'//ltt//'____'//lng1//'.txt'
	                 goto 17
	         endif
               if(j.ge.10.and.j.lt.100) then 
			         write (lng2, '(i2)') j
	                 s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 17
               endif
		             write (lng, '(i3)') j
                       s$='_'//ltt//'__'//lng//'.txt'
	                 goto 17
	                         
!		             write (ltt, '(i2)') i
!                      write (lng, '(i3)') j
!					 s$=ltt//'_'//lng//'.txt'
!                  	 write (*,*) s$
c		write (*,*) s$
c		stop
      
  17	        open (3, file='c:\d\rivers\res\'//s$)   ! OPEN of selected file with 
                                  ! input cell runoff in mm/time_step
			
c			s1$=s$//'.out'
c			open (7, file=s1$)  !OPEN of selected file with 
                                  ! otput cell streamflow in mm/time_st					        

				step=1
	             S=0.
				
!                   Calculation of basin area				
				areal=areal+area(long-.5,lati-.5)
				
c	  write(*,*) s$, areal*10000., 'sq.km'

!          Calculation of routing trajectory for selected cell
	       do k=1, k0
c       write(*,*) 'k0',k0,k
!              CALCULATION of CELL CHARACTERISTICS        
				c=u0*86400./sqrt(area(long-.5,lati-.5))/100000.
  		        Ct=exp(-c*dt)
			
			
!              Calculation of time set of RUNOFF for selected CELL				  
                 do m=1,m0
				  if(k.eq.1) then
c			        read (3,*) Y(m)    ! Y - cell runoff in mm/time_step        

c		read (3,*) x,x,x,x,x
c		read (3,*) x
      	 
c		 stop
		   
		read (3,*) n1(m), n2(m), n3(m),x,x,x,x1,x2,x,x,x,x,x,x,x,x,x 
c	read (3,*) n1(m), n3(m),x,x,x,x1,x2 
              Y(m)=(x1+x2)*3600.*24.    ! Y - cell runoff in mm/time_step        


		  Din(m)=Y(m)/1000.*area(long-.5,lati-.5)*10000.*1000000.
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
	  select case (ind1)
		case (1)
		  long=long
	      lati=lati+1
			case (2)		
		  long=long+1
		  lati=lati+1         
		case (3)
		  long=long+1
		  lati=lati  	    
		case (4)
		  long=long+1
		  lati=lati-1  		
		case (5)
		  long=long
		  lati=lati-1  		
		case (6)		
 		  long=long-1
		  lati=lati-1           
		case (7)
		  long=long-1
		  lati=lati  	    
		case (8)
		  long=long-1
		  lati=lati+1  
          end select
         ind1=IND2(long-.5,lati-.5)

!      SUMMARIZING OF NUMBER OF STEPS in ROUTING TRAJECTORY 
!      from SELECTED CELL to Basin EXIT cell
	   step=step+1

  		    end do ! k  - end of output Streamflow calculation 
                     ! for one way  
					
c	        Write (5,15) lon,lat, step, Q    !  m^3/dt
c               Write (7,15) lon,lat, step, Q    !  m^3/dt

 15           format (2f8.1,i4,4932e12.3)


!      SUMMARIZING of  OUTPUT streamflow for EXIT CELL
		do mm=step,m0
		Runoff(mm)=Q(mm-step+1)+Runoff(mm)     !m^3/dt
		enddo !mm			


	end if    !  end of output Streamflow calculation 
                !  for all previous CELLS  
		close (3)

c		                        close(7)

		end do !   i - end of calculation on latitude
	end do     !   j - end of calculation on longitude	
		
!	                                TOTAL BASIN STREAMLOW mm/dt
c			write (5,*) long0,lati0, 1
		Do mm=1, m0
       iy=n1(mm)
	 im=n2(mm)
	 id=n3(mm)
       iy=iy-iyr_fir+1
       run_mmday(iy,im,id)=Runoff(mm)/areal/10000./1000000.*1000.       !  mm/dt
    
	 Write (5,16)  n1(mm),im,id,run_mmday(iy,im,id)   !  mm/dt
c	 Write (5,16)  n1(mm),n2(mm),n3(mm),Runoff(mm)/areal/10000.
c     !/1000000.*1000. !  mm/dt
        
 16           format (3i6,e12.3)                       !4932
		enddo
				
	close (4)
	close(5)
cccccccccccccccccccccccccccccccccc	               write(*,*)  '      areal=', areal*10000.,'sq km'



c       START of BASIN WATER BALANCE
! 
     	Open (4, file='c:\d\rivers\data\'//river//'\index1.txt')    !  input INDEX1-file with cell data    

		Runoff1=0.
	    Evap=0.
	sw=0.
	rdown=0.
	t2=0.
	u2=0.
	q2=0.
	    Prec=0.
		k=0
						 
	do j=j0, j1
		do i=i0, i1
			read (4,*) lon,lat,ind,ara
	                   long=lon
					   lati=lat
					   ind1=ind

	if (ind.gt.0.) then 


!              PREPARARION OF FILE-NAME for selected CELL
       	   
		if (i.lt.-9) then 
			  write (ltt1, '(i3)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$=ltt1//'_'//lng4//'.txt'
	                 goto 19
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$=ltt1//'__'//lng//'.txt'
	                  goto 19
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$=ltt1//'___'//lng2//'.txt'
	                 goto 19 
	         endif
               if(j.ge.0.and.j.lt.10) then 
			         write (lng1, '(i1)') j
                       s$=ltt1//'____'//lng1//'.txt'
	                 goto 19
	         endif
               if(j.ge.10.and.j.lt.100) then
			         write (lng2, '(i2)') j
	                 s$=ltt1//'___'//lng2//'.txt'
	                 goto 19
               endif
		             write (lng, '(i3)') j
                       s$=ltt1//'__'//lng//'.txt'
	                 goto 19
          endif

    	    if (i.ge.-9.and.i.lt.0) then
			  write (ltt, '(i2)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='_'//ltt//'_'//lng4//'.txt'
	                 goto 19
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$='_'//ltt//'__'//lng//'.txt'
	                  goto 19
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 19 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='_'//ltt//'____'//lng1//'.txt'
	                 goto 19
	         endif
               if(j.ge.10.and.j.lt.100) then
			         write (lng2, '(i2)') j
	                 s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 19
               endif
		             write (lng, '(i3)') j
                       s$='_'//ltt//'__'//lng//'.txt'
	                 goto 19
          endif
						
          if (i.ge.0.and.i.lt.10) then 
		     write (ltt2, '(i1)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='__'//ltt2//'_'//lng4//'.txt'
	                 goto 19
               endif 
               if(j.ge.-99.and.j.lt.-9) then
			          write (lng, '(i3)') j
                        s$='__'//ltt2//'__'//lng//'.txt'
	                  goto 19
	         endif
               if(j.ge.-9.and.j.lt.0) then
			         write (lng2, '(i2)') j
                       s$='__'//ltt2//'___'//lng2//'.txt'
	                 goto 19 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='__'//ltt2//'____'//lng1//'.txt'
	                 goto 19
	         endif
               if(j.ge.10.and.j.lt.100) then 
			         write (lng2, '(i2)') j
	                 s$='__'//ltt2//'___'//lng2//'.txt'
	                 goto 19
               endif
		             write (lng, '(i3)') j
                       s$='__'//ltt2//'__'//lng//'.txt'
	                 goto 19
          endif
               

	          write (ltt, '(i2)') i
	          if(j.lt.-99) then 
			         write (lng4, '(i4)') j
                       s$='_'//ltt//'_'//lng4//'.txt'
	                 goto 19
               endif 
               if(j.ge.-99.and.j.lt.-9) then 
			          write (lng, '(i3)') j
                        s$='_'//ltt//'__'//lng//'.txt'
	                  goto 19
	         endif
               if(j.ge.-9.and.j.lt.0) then 
			         write (lng2, '(i2)') j
                       s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 19 
	         endif
               if(j.ge.0.and.j.lt.10) then
			         write (lng1, '(i1)') j
                       s$='_'//ltt//'____'//lng1//'.txt'
	                 goto 19
	         endif
               if(j.ge.10.and.j.lt.100) then 
			         write (lng2, '(i2)') j
	                 s$='_'//ltt//'___'//lng2//'.txt'
	                 goto 19
               endif
		             write (lng, '(i3)') j
                       s$='_'//ltt//'__'//lng//'.txt'
	                 goto 19
	                         
      
  19	        open (3, file='c:\d\rivers\res\'//s$)   ! OPEN of selected file with 
               k=k+1                   ! input cell water balance componennts in mm/time_step
			
ccccccccccccccccccccccccccccccccc			WRITE(*,*) s$
		        do m=1,m0	       
	read (3,*) x,x,x,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14 
c		read (3,*) x,x,x1,x2,x3,x4,x5 
		    prec(m)=prec(m)+(x2+x3)*area(long-.5,lati-.5)/areal
	        evap(m)=evap(m)+x1*area(long-.5,lati-.5)/areal
	sw(m)=sw(m)+x6*area(long-.5,lati-.5)/areal
	rdown(m)=rdown(m)+x7*area(long-.5,lati-.5)/areal
	t2(m)=t2(m)+x8*area(long-.5,lati-.5)/areal
	u2(m)=u2(m)+x9*area(long-.5,lati-.5)/areal
	q2(m)=q2(m)+x10*area(long-.5,lati-.5)/areal
	swe(m)=swe(m)+x11*area(long-.5,lati-.5)/areal
	ksi(m)=ksi(m)+x12*area(long-.5,lati-.5)/areal
	ksiot(m)=ksiot(m)+x13*area(long-.5,lati-.5)/areal
	hsnow(m)=hsnow(m)+x14*area(long-.5,lati-.5)/areal
	        
			Runoff1(m)=Runoff1(m)+(x4+x5)*area(long-.5,lati-.5)/areal
				end do

	end if
	close(3)
	  enddo
	enddo
		close(4)

			write(*,*) k

		open (9, file='c:\d\rivers\res\balance')
	write (9,*)' PRECIP  EVAP  RUNOFF(no_transf) SW    
     !RDOWN    T2   U2   Q2'
	Do mm=1, m0
	       arunoff1(mm)=Runoff1(mm)*86400.
	  Write (9,20) n1(mm),n2(mm),n3(mm),  
     !	            Prec(mm)*86400.,
     !                Evap(mm)*86400., 
     !	            Runoff1(mm)*86400.,
     !        sw(mm),rdown(mm),t2(mm),u2(mm),q2(mm),
     !         (Runoff(mm)/areal/10000./1000000.*1000.) 
     !,swe(mm),
     !ksi(mm),ksiot(mm),hsnow(mm)
	enddo

20          format (3i5,13e14.3)

ccccccccccccccccccccccccccccccccccccccccccc	close(9)
c               END of WATER BALANCE


c             CALCULATION of MONTHLY STREAMFLOW, mm/month

	   open (5, file='c:\d\rivers\res\Q_mon')
	
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

	open (5, file='c:\d\rivers\res\year')
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

	write  (*,*) '       E  N  D'
	
	     end 