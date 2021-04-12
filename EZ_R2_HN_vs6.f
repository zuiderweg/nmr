! program to calculate R2 relaxation of amide protons (HN or H ) from a pdb file
! allows another correlation time for HN-methyl relaxation vector motion
! Erik R.P. Zuiderweg, The Uiversity of Michigan, Ann Arbor MI USA, Radboud Univesrity, Nijmegen, The Netherlands
! March 2021
! Donot distribute without the author's consent.



	implicit none

	integer:: i,j, k,l, m, n
	integer:: ios



	! latest PDB format 
	character(len=4):: atom 
	character(len=4):: atname, resname 
	character(len=1):: chain, altloc,  icode 
	character(len=2):: ID,  charge 
	integer::  resnumber,  atnumber 	     
	real:: x,y,z, occ, Bfact
	

	character(len=40):: pdbname, file_out
	character(len=200):: line
	real:: cutoff, distance, x_ref, y_ref, z_ref, R2
	integer:: resnumber_ref
	character(len=4):: atname_ref, answer
	real:: tc, spec, wh, B0, tc_methyl


!!! start program




	
 6	format(a4,2x,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,10x,2a2) ! latest pdb format
				
 21	print 1
 1	format(' pdb file  > ', $)
	read(*,2) pdbname
 2	format(a)	

 	open(unit=10,file=pdbname,status='old',form='formatted', iostat=ios)
	if(ios.gt.0) then
	print 20
 20	format(' file not found. try again (y/n) > ', $)
 	read(*,*) answer
 	if(answer(1:1).eq.'n'.or.answer(1:1).eq.'N') stop
 	goto 21
 	endif

	open(unit=40, file='scratch', form='formatted')   ! make internal copy

 100	read(unit=10,fmt=2,err=101, end=101 )line
	write(unit=40,fmt=2 )line
	goto 100

 101	continue
 	rewind(10)
 	rewind(40)


	print 7
 7	format(' distance cutoff  (A) > ', $)
	read(*,*) cutoff

	print 11
 11	format(' enter tc in seconds ',$)
	read(*,*) tc

	print 12
 12	format(' enter methyl tc in seconds ',$)
	read(*,*) tc_methyl

	
	print 13
 13	format(' enter spectrometer 500, 600 etc ', $)
 	read(*,*) spec
	B0=spec/500.*11.76
	
	wh=2.6752d8*B0

 23	print 14
 14	format(' output file name > ', $)
	read(*,2) file_out
	open(unit=11,file=file_out,status='new',form='formatted', iostat=ios)

	if(ios.gt.0) then
	print 22
 22	format(' file exists. overwrite?  (y/n) > ', $)
 	read(*,*) answer
 	if(answer(1:1).eq.'n'.or.answer(1:1).eq.'N') 	goto 23
 	open(unit=11,file=file_out,status='unknown',form='formatted')
 	endif	






	
	write(11,*) ' resnumber atname Lw(Hz)'


!find all the protons of interest one by one and create a sphere around them 
	rewind(40)  ! not really necessary
! MASTERLOOP over all protons in the pdb file

 1000	read(unit=40,fmt=2,err=1002, end=1002 )line
	if(line(1:4).ne.'ATOM') goto 1000
	backspace(40)
	read(unit=40,fmt=6)
 	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x,y,z,
	1	occ, Bfact, ID, charge
 
 	if(atname(2:3).eq.'H '.or.atname(2:3).eq.'HN') then
	write(*,*)' working on ', resnumber, atname(1:4)

  	resnumber_ref=resnumber
	atname_ref=atname(1:4)
	x_ref=x
	y_ref=y
	z_ref=z





 	goto 1001
	endif
	goto 1000

 1001	continue		



	R2=0.
	rewind(10)   ! necessary
 2000	read(unit=10,fmt=2,err=2001, end=2001 )line
	if(line(1:4).ne.'ATOM') goto 2000
		backspace(10)
 		read(unit=10,fmt=6)
	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x,y,z,
	1	occ, Bfact, ID, charge
 
	if(atname(1:1).eq.'H'.or.atname(2:2).eq.'H' ) then
!	
	distance=sqrt((x_ref-x)**2+(y_ref-y)**2+(z_ref-z)**2)
	if(distance.gt.cutoff) goto 2000
	if(distance.lt.1.e-3) goto 2000   ! skip the reference proton 

	if(resname.eq.'ALA') then
	if(atname(2:3).eq.'HB') then
	R2=R2+5.60e11/20.*(5*tc_methyl + 9*tc_methyl/(1.+(tc_methyl*wh)**2) + 6*tc_methyl/(1.+(tc_methyl*wh*2.)**2))/(distance**6)  ! unlike , fast motion
	goto 2000
	endif
	endif

	if(resname.eq.'VAL') then
	if(atname(1:2).eq.'HG') then
	R2=R2+5.60e11/20.*(5*tc_methyl + 9*tc_methyl/(1.+(tc_methyl*wh)**2) + 6*tc_methyl/(1.+(tc_methyl*wh*2.)**2))/(distance**6)  ! unlike , fast motion
	goto 2000
	endif
	endif

	if(resname.eq.'LEU') then
	if(atname(1:2).eq.'HD') then
	R2=R2+5.60e11/20.*(5*tc_methyl + 9*tc_methyl/(1.+(tc_methyl*wh)**2) + 6*tc_methyl/(1.+(tc_methyl*wh*2.)**2))/(distance**6)  ! unlike , fast motion
	goto 2000
	endif
	endif


	if(resname.eq.'MET') then
	if(atname(2:3).eq.'HE') then
	R2=R2+5.60e11/20.*(5*tc_methyl + 9*tc_methyl/(1.+(tc_methyl*wh)**2) + 6*tc_methyl/(1.+(tc_methyl*wh*2.)**2))/(distance**6)  ! unlike , fast motion
	goto 2000
	endif
	endif

	if(resname.eq.'ILE') then
	if(atname(1:2).eq.'HD') then
	R2=R2+5.60e11/20.*(5*tc_methyl + 9*tc_methyl/(1.+(tc_methyl*wh)**2) + 6*tc_methyl/(1.+(tc_methyl*wh*2.)**2))/(distance**6)  ! unlike , fast motion
	goto 2000
	endif
	if(atname(1:3).eq.'HG2') then
	R2=R2+5.60e11/20.*(5*tc_methyl + 9*tc_methyl/(1.+(tc_methyl*wh)**2) + 6*tc_methyl/(1.+(tc_methyl*wh*2.)**2))/(distance**6)  ! unlike , fast motion
	goto 2000
	endif
	endif

	R2=R2+5.60e11/20.*(5*tc + 9*tc/(1.+(tc*wh)**2) + 6*tc/(1.+(tc*wh*2.)**2))/(distance**6)  ! unlike everything else
	goto 2000
	endif

	goto 2000
 2001	continue

 !	hbar=1.054d-34
!	mu=1.d-7
	
!	yh=2.6752d8
!	yn=-2.712d7
!	yc= 6.728d7
	
!	wh=yh*B0
!	wn=yn*B0

! parsing out the constants we get
! mu^2 yh^4 hbar^2 =5.60e-49
! angstrom correction 1e60
! R2 = 5.60e11 / 20 * ( 9tc + 15 tc/(1+wh2tc^2) + 6 tc/(1+4w^2tc^2))   ! LIKE
! R2 = 5.60e11 / 20 * ( 5tc + 9 tc/(1+wh2tc^2) + 6 tc/(1+4w^2tc^2))   ! UNLIKE


 	write(11,*)  resnumber_ref, atname_ref , R2/3.1416


 	goto 1000

 1002 stop

 	end