! program to calculate the sum dipolar field of surrounding hydrogens on a center amide hydrogen in a pdb file
! taking dipolar R2 cross-correlated relaxation (interference) into account
! by permuting the diplar sign distribution of the 10 closest protons
! the program creates a powder distribution of the total dipolar field by rotating the magnetic field in the PDB frame according to
! a spherical distribution (http://corysimon.github.io/articles/uniformdistn-on-sphere/).
! The solution R2 is calculated as the second moment of the powder pattern.
! Written by Erik R.P. Zuiderweg, The University of Michigan, Ann Arbor, MI, USA
! and Radboud University, Nijmegen, The Netherlands
! December 2020 (zuiderwe@umich.edu or erikzuiderweg@gmail.com)
! Donot distribute without the author's consent.


	implicit none

	integer:: i,j, k,l, m, n, mm

	integer::k_count, l_count, ios
	character(len=40):: pdbname, pdbname1, output_name, fullfilename, moment_name, moment_name_full, concise
	character(len=4):: flag
	character(len=200):: line
	character(len=50)::glued
	character(len=3):: answer


	! latest PDB format 
	character(len=4):: atom 
	character(len=4):: atname, resname 
	character(len=1):: chain, altloc,  icode 
	character(len=2):: ID,  charge 
	integer::  resnumber,  atnumber 	     
	real:: x,y,z, occ, Bfact

	integer::resnumber_save, swap
	character(len=4)::resname_save
	character(len=4)::atname_save
	real::x_save, y_save, z_save, distance_save, sign_save



	type:: pdb
		integer::resnumberPDB
		character(len=4)::resnamePDB
		character(len=4)::atnamePDB
		real::xPDB
		real::yPDB
		real::zPDB
		real::distPDB
		real::signPDB
	end type pdb
	type(pdb), dimension(400)::shell


	type::dist
		real::xDIST
		real::yDIST
		real::zDIST
	end type dist
	type(dist), dimension(5000)::distribution

	real:: cutoff, distance, x_ref, y_ref, z_ref, x_shell, y_shell, z_shell, x_new, y_new, z_new, fraction
	real:: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dis,dihed
	integer:: resnumber_ref
	character(len=4):: atname_ref, resname_ref
	integer:: nuc_count, self_id
	real::  dip_total(5000), dip, distance_to_origin, dip_min, dip_max, sign, pi, theta, phi, cosdef
	integer:: icount

	integer::k1,k2,k3,k4, k5, k6, k7, k8, k9, k10


	real::super_sum, second, average, second_J, second_noJ
	real:: J_HNHA(1000)
	real::second_stat(1024), second_stat_average_no_J, second_stat_average_with_J, inverse_no_J, inverse_with_J
	real::deuter, tc, B0, spec, wh, ffactor



!!! start program




	
 6	format(a4,2x,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,10x,2a2) ! latest pdb format
				
 21	print 1
 1	format(' pdb file > ', $)
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


 23	print 5
 5	format('  output name of complete R2 documentation > ', $)
 	read(*,2) moment_name
	open(unit=90,file=moment_name,status='new',form='formatted', iostat=ios)
	if(ios.gt.0) then
	print 22
 22	format(' file exists. overwrite?  (y/n) > ', $)
 	read(*,*) answer
 	if(answer(1:1).eq.'y'.or.answer(1:1).eq.'Y') 	then
	open(unit=90,file=moment_name,status='unknown',form='formatted')
 	goto 24
 	endif	
 	goto 23
 	endif

 24	print 7
 7	format('  output name of concise R2 documentation  > ', $)
 	read(*,2) concise
 	open(unit=80,file=concise,status='new',form='formatted', iostat=ios)
	if(ios.gt.0) then
	print 25
 25	format(' file exists. overwrite?  (y/n) > ', $)
 	read(*,*) answer
 	if(answer(1:1).eq.'y'.or.answer(1:1).eq.'Y') 	then
	open(unit=80,file=concise,status='unknown',form='formatted')
 	goto 26
 	endif	
 	goto 24
 	endif


 26	write(80,*) ' resnumber_ref  resname_ref atname_ref JHNHA',  
	1	' dip_min dip_max linewdth_no_J linewidth_with_J, peakheight_no_J, peakheight_with_J '
 	print 8
 8	format(' distance cutoff  (A) > ', $)
	read(*,*) cutoff



	print 11
 11	format(' enter tc in seconds ',$)
	read(*,*) tc
	
	print 12
 12	format(' enter spectrometer 500, 600 etc ', $)
 	read(*,*) spec
	B0=spec/500.*11.76
	
	wh=2.6752e8*B0
	
	
	write(*,*) 'initializing '
	

	open(unit=30, file='SELECTION', form='formatted', status='unknown')
	
	pi=3.1415926536
	
	open(unit=11,status='scratch',form='formatted')
	open(unit=12,status='scratch',form='formatted')
	open(unit=13,status='scratch',form='formatted')
	open(unit=40,status='scratch',form='formatted')

 50	read(unit=10,fmt=2, err=51, end=51) line  ! make internal copies of the PDB file
 	write(11,2)line
 	write(12,2) line
 	write(13,2) line
 	write(40,2) line
 	go to 50
 51	continue

 	rewind(10)
	rewind(11)
	rewind(12)
	rewind(13)
	rewind(40)



 	! calculate the 3JHNHa couplings


 100	read(unit=10,fmt=2, err=1111, end=1111) line
	if(line(1:4).eq.'ATOM') then
		backspace(unit=10)
		goto 110
	endif
	goto 100

 110	read(unit=10,fmt=6)
	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x1,y1,z1,
	1	occ, Bfact, ID, charge

	
	if(atname(2:3).ne.'C ') go to 100       ! 


	rewind(11)
  	rewind(12)
  	rewind(13)


 200	read(unit=11,fmt=2, err=100,end=100) line
	if(line(1:4).eq.'ATOM') then
		backspace(unit=11)
		goto 210
	endif
	goto 200

 210	read(unit=11,fmt=6, err=100,end=100)
	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x2,y2,z2,
	1	occ, Bfact, ID, charge

	
	if(atname(2:3).ne.'N ') go to 200         ! find N

	dis=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
	if(dis.lt.0.1.or.dis.gt.2.0)   goto 200   !  find close N

	
 300	read(unit=12,fmt=2, err=100,end=100) line
	if(line(1:4).eq.'ATOM') then
		backspace(unit=12)
		goto 310
	endif
	goto 300

 310	read(unit=12,fmt=6,err=100,end=100)
	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x3,y3,z3,
	1	occ, Bfact, ID, charge

	if(atname(2:3).ne.'CA') go to 300         ! eg find CA
	dis=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
	if(dis.lt.0.1.or.dis.gt.2.0)   goto 300   !  find close CA

	

 400	read(unit=13,fmt=2, err=100,end=100) line
	if(line(1:4).eq.'ATOM') then
		backspace(unit=13)
		goto 410
	endif
	goto 400


 410	read(unit=13,fmt=6,err=100,end=100)
	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x4,y4,z4,
	1	occ, Bfact, ID, charge
	if(atname(2:3).ne.'C ') go to 400           ! eg find C
	dis=sqrt((x3-x4)**2+(y3-y4)**2+(z3-z4)**2)
	if(dis.lt.0.1.or.dis.gt.2.0) go to 400   ! find close C


	call  dihedral(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dihed)

	dihed=dihed/pi*180.
	J_HNHA(resnumber)=7.97*(cos((dihed-60)*pi/180.))**2-1.26*(cos((dihed-60)*pi/180.))+ 0.63  ! Hz

	write(77,*) resnumber, 'dihed', dihed, J_HNHA(resnumber)


	rewind(12)
	rewind(11)
	rewind(13)
 	goto 100	

 1111	rewind(10)




!find all the protons of interest one by one and create a sphere around them 
	rewind(40)  ! not really necessary
! MASTERLOOP over all protons in the pdb file

 1000	read(unit=40,fmt=2,err=1003, end=1002 )line
	if(line(1:4).ne.'ATOM') goto 1000
	backspace(40)
	read(unit=40,fmt=6)
 	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x,y,z,
	1	occ, Bfact, ID, charge
 
 	if(atname(2:3).eq.'H '.or.atname(2:3).eq.'HN') then
	write(*,*)' working on ', resnumber, atname(1:4)
	write(30,*)' working on ', resnumber, atname(1:4)

  	resnumber_ref=resnumber
	atname_ref=atname(1:4)
	resname_ref=resname
	x_ref=x
	y_ref=y
	z_ref=z





 	goto 1001
	endif
	goto 1000

 1001	continue		





! now create new pdb file sphere and store in the SHELL structure
	write(*,*) ' reading the relevant sphere of proton coordinates into memory '

	nuc_count=0
	rewind(10)   ! necessary
 2000	read(unit=10,fmt=2,err=2001, end=2001 )line
	if(line(1:4).ne.'ATOM') goto 2000
		backspace(10)
 		read(unit=10,fmt=6)
	1	atom,atnumber,atname,altloc,resname,chain,resnumber,icode,x,y,z,
	1	occ, Bfact, ID, charge
 
	if(atname(1:1).eq.'H'.or.atname(2:2).eq.'H') then




	distance=sqrt((x_ref-x)**2+(y_ref-y)**2+(z_ref-z)**2)

	if(distance.lt.1.e-3) goto 2000   ! skip the reference proton ! self_id=nuc_count
	
	if(distance.le.cutoff) then  ! collect the protons within cutoff only and reset the origin at the proton of interest


	write(*,*) resname,chain,resnumber, atname, distance
	write(30,*) resname,chain,resnumber, atname, distance

	nuc_count=nuc_count+1
	if(nuc_count.gt.399) goto 1500

	shell(nuc_count)%resnumberPDB=resnumber
	shell(nuc_count)%resnamePDB=resname
	shell(nuc_count)%atnamepdb=atname
	shell(nuc_count)%xPDB=x-x_ref
	shell(nuc_count)%yPDB=y-y_ref	
	shell(nuc_count)%zPDB=z-z_ref
	shell(nuc_count)%distPDB=distance
	shell(nuc_count)%signPDB=1.

	endif
	endif
	goto 2000
 2001	continue

 	write(*,*) " The shell contains "  , nuc_count , "protons around the center proton "
	write(30,*) " The shell contains "  , nuc_count , "protons around the center proton "

 	if(nuc_count.eq.0) goto 1000

! sort this SHELL STRUCTURE small to large
	swap=1

 2500 if(swap.eq.0) goto 2501
	do  i=1,nuc_count
	do  j=0,nuc_count-i
	swap=0
	if(shell(i)%distPDB.gt.shell(i+j)%distPDB) then

	swap=1
	resnumber_save=shell(i)%resnumberPDB
	resname_save=shell(i)%resnamePDB
	atname_save=shell(i)%atnamepdb
	x_save=shell(i)%xPDB
	y_save=shell(i)%yPDB	
	z_save=shell(i)%zPDB
	distance_save=shell(i)%distPDB
	sign_save=shell(i)%signPDB

	shell(i)%resnumberPDB=shell(i+j)%resnumberPDB
	shell(i)%resnamePDB=shell(i+j)%resnamePDB
	shell(i)%atnamePDB=shell(i+j)%atnamePDB
	shell(i)%xPDB=shell(i+j)%xPDB
	shell(i)%yPDB=shell(i+j)%yPDB
	shell(i)%zPDB=shell(i+j)%zPDB
  	shell(i)%distPDB=shell(i+j)%distPDB
	shell(i)%signPDB=shell(i+j)%signPDB

	shell(i+j)%resnumberPDB=resnumber_save
	shell(i+j)%resnamePDB=resname_save
	shell(i+j)%atnamePDB=atname_save
	shell(i+j)%xPDB=x_save
	shell(i+j)%yPDB=y_save
	shell(i+j)%zPDB=z_save
  	shell(i+j)%distPDB=distance_save
	shell(i+j)%signPDB=sign_save

	endif
	end do
	end do


	go to 2500

 2501	continue


	write(*,*) (shell(k)%distPDB, k=1, nuc_count)

 	open(unit=60,status='scratch',form='formatted')

 		do 710 k10=1,1024

! start by attributing random signs to all dipoles in the shell to simulate cross correlations
	
	do j=1,400
	shell(j)%signPDB=100.
	end do

	do i=1, nuc_count
	call random_number(fraction)
	if(fraction.ge.0.5) then 
	shell(i)%signPDB=1.
	else
	shell(i)%signPDB=-1.
	endif
	end do


	

! finished with permuting signs of PDB records

! now generate a random but isotropic spherical distribution needed below

	do i=1,5000
			call random_number(fraction)
			theta=2.*pi*fraction
			call random_number(fraction)
			phi=acos(1.-2.*fraction)
			distribution(i)%xDIST=sin(phi)*cos(theta)
			distribution(i)%yDIST=sin(phi)*sin(theta)
			distribution(i)%zDIST=cos(phi)
	end do

! finally start computing dipolar interaction
! start by taking the PDB frame as the magnetic frame 
! calculate and add the dipolar field of the surrounding protons on location 0,0,0, z component needed only. 
! then use the ditribution above for new z-orientations as the molecule lies in different directions
	
	dip_min=+1.d12
	dip_max=-1.d12

	do 600, icount=1,5000		! number of orientations in powder pattern
	dip_total(icount)=0.0d0

	do 700 n=1,nuc_count          ! number of protons in shell


	x_shell=shell(n)%xPDB
	y_shell=shell(n)%yPDB
	z_shell=shell(n)%zPDB
	sign=shell(n)%signPDB
	distance_to_origin=shell(n)%distPDB


	cosdef=(shell(n)%xPDB*distribution(icount)%xDIST+shell(n)%yPDB*distribution(icount)%yDIST+shell(n)%zPDB*distribution(icount)%zDIST)/shell(n)%distPDB/1.
	dip=2.81e5*sign*(3.*(cosdef)**2-1.)/(shell(n)%distPDB**3)      ! solid state with flip-flop term the dipole of a certain proton pair with certain permutation in a given orientation
	dip_total(icount)=dip_total(icount)+dip   ! the sum dipole of a all proton pairs with certain permutation in a given orientation

  700	continue   ! nulcei in shell

 	if(dip_total(icount).ge.dip_max) dip_max=dip_total(icount)
	if(dip_total(icount).le.dip_min) dip_min=dip_total(icount)



 !	write(99,*) resnumber_ref,  resname_ref, atname_ref,  shell(1)%signPDB, shell(2)%signPDB, shell(3)%signPDB, shell(4)%signPDB, shell(5)%signPDB, shell(6)%signPDB,shell(7)%signPDB, shell(8)%signPDB,
 !	1	distribution(icount)%xDIST, distribution(icount)%yDIST, distribution(icount)%zDIST, dip_total(icount), dip_min, dip_max  !(shell(l)%signPDB, l=5,100)
  600 continue         ! orientations

  	super_sum=0.      ! now calculate proton linewidth using second moment of powder pattern from dip_total
  	second=0.
  	do icount =1,5000
  	super_sum=super_sum+dip_total(icount)
  	end do
  	average=super_sum/5000.
  	do icount =1,5000
  	second= second+(dip_total(icount)-average)**2
  	end do
  	second=second/5000. ! The "D2 term'  in radians
 
	write(60,*) resnumber_ref,  resname_ref, atname_ref,second !! we are writing this away for 1024 different permutations
	write(90,*) resnumber_ref,  resname_ref, atname_ref,  shell(1)%signPDB, shell(2)%signPDB, shell(3)%signPDB, shell(4)%signPDB, shell(5)%signPDB, shell(6)%signPDB, shell(7)%signPDB, shell(8)%signPDB,shell(9)%signPDB,shell(10)%signPDB, second, second_J

 710		continue
	
	write(90,*)	'  '
	close(99)

	dip_min=+1.d12
	dip_max=-1.d12
	second_stat_average_no_J=0.
	second_stat_average_with_J=0.
	inverse_no_J=0.0
	inverse_with_J=0.0
	rewind(60)


	
! here we collect the results for a single HN ...!!!
	do l=1,1024
	read(60,*) resnumber_ref,  resname_ref, atname_ref, second_stat(l)    ! second_stat(l) is the "D2" term for each permutation
	if(second_stat(l).le.dip_min) dip_min=second_stat(l)
	if(second_stat(l).ge.dip_max) dip_max=second_stat(l)
	second_stat_average_no_J=second_stat_average_no_J+second_stat(l)/pi*4.*tc  
	second_stat_average_with_J=second_stat_average_with_J+second_stat(l)/pi*4.*tc  + J_HNHA(resnumber_ref)
	inverse_no_J=inverse_no_J+1.0/ (second_stat(l)/pi*4.*tc)
	inverse_with_J=inverse_with_J+1.0/(second_stat(l)/pi*4.*tc+J_HNHA(resnumber_ref)) 

	end do

	close(60)

	write(80,*) resnumber_ref,  resname_ref, atname_ref,  J_HNHA(resnumber_ref),
	1	dip_min, dip_max, second_stat_average_no_J/1024., second_stat_average_with_J/1024., inverse_no_J, inverse_with_J

 	goto 1000   !! next HN....

 1002 close(90)
 		stop

 1003 write(*,*) ' PDB read error'
  	stop

 1500	write(*,*) "too many shell protons around residue ", resnumber_ref, "please reduce cutoff "
 	stop

	end



	



	subroutine dihedral(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dihed)
	implicit none
	integer:: i,j, k,l, m, n
	real::x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dihed
	real::ax,ay,az,bx,by,bz,cx,cy,cz
	real::qx,qy,qz,rx,ry,rz
	real::rq,rr,qq, result,sx,sy,sz, test

	ax=x1-x2	! vect1
	ay=y1-y2
	az=z1-z2
	bx=x2-x3	! vect2
 	by=y2-y3
	bz=z2-z3
	cx=x3-x4	! vect3
	cy=y3-y4
	cz=z3-z4
	qx=ay*bz-az*by	! vext1*vect2 normal to plane of 1 and 2
	qy=az*bx-ax*bz
	qz=ax*by-ay*bx

	rx=by*cz-bz*cy	! vext2*vect3 normal to plane of 2 and 3
	ry=bz*cx-bx*cz
	rz=bx*cy-by*cx

	rq=rx*qx+ry*qy+rz*qz  ! dot product between the two normals
	rr=rx*rx+ry*ry+rz*rz	! modulus
	qq=qx*qx+qy*qy+qz*qz

	result=rq/(sqrt(rr*qq))
	if(result.ge.1.0)       result=0.99999
	if(result.le.-1.0)       result=-0.99999
	dihed=acos(result)  ! radians ! angle in radians between the normals

	sx=qy*rz-qz*ry	!! construct form normals  vector parallel to vector 2
	sy=qz*rx-qx*rz
	sz=qx*ry-qy*rx

	test=bx*sx+by*sy+bz*sz !! paralel or antiparalel?
	if(test.gt.0.) dihed=-dihed
	
	return
	end

