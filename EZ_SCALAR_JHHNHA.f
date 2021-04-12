! program to calculate 3JHNHA couplings from a pdb file, no protons needed
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
	real:: pi
	real:: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dihed, J_HNHA, dis
	integer:: resnumber1, resnumber2, resnumber3, resnumber4
	character(len=4) atname1, atname2, atname3, atname4, answer


!!! start program







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

 23	print 14
 14	format(' output file name > ', $)
	read(*,2) file_out
	open(unit=77,file=file_out,status='new',form='formatted', iostat=ios)

	if(ios.gt.0) then
	print 22
 22	format(' file exists. overwrite?  (y/n) > ', $)
 	read(*,*) answer
 	if(answer(1:1).eq.'n'.or.answer(1:1).eq.'N') 	goto 23
 	open(unit=77,file=file_out,status='unknown',form='formatted')
 	endif	



	
	
	write(*,*) 'initializing '
	
	write(77,*) ' resnumber atname  resnumber  atname resnumber atname resnumber atname  dihed J_HNHA '
	
	pi=3.1415926536
	
	open(unit=11,status='scratch',form='formatted')
	open(unit=12,status='scratch',form='formatted')
	open(unit=13,status='scratch',form='formatted')
	

 50	read(unit=10,fmt=2, err=51, end=51) line  ! make internal copies of the PDB file
 	write(11,2)line
 	write(12,2) line
 	write(13,2) line
 	go to 50
 51	continue

 	rewind(10)
	rewind(11)
	rewind(12)
	rewind(13)




 	! calculate the 3JHNHa couplings


 100	read(unit=10,fmt=2, err=1111, end=1111) line
	if(line(1:4).eq.'ATOM') then
		backspace(unit=10)
		goto 110
	endif
	goto 100

 110	read(unit=10,fmt=6)
	1	atom,atnumber,atname1,altloc,resname,chain,resnumber1,icode,x1,y1,z1,
	1	occ, Bfact, ID, charge

	
	if(atname1(2:3).ne.'C ') go to 100       ! 


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
	1	atom,atnumber,atname2,altloc,resname,chain,resnumber2,icode,x2,y2,z2,
	1	occ, Bfact, ID, charge

	
	if(atname2(2:3).ne.'N ') go to 200         ! find N

	dis=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
	if(dis.lt.0.1.or.dis.gt.2.0)   goto 200   !  find close N

	
 300	read(unit=12,fmt=2, err=100,end=100) line
	if(line(1:4).eq.'ATOM') then
		backspace(unit=12)
		goto 310
	endif
	goto 300

 310	read(unit=12,fmt=6,err=100,end=100)
	1	atom,atnumber,atname3,altloc,resname,chain,resnumber3,icode,x3,y3,z3,
	1	occ, Bfact, ID, charge

	if(atname3(2:3).ne.'CA') go to 300         ! eg find CA
	dis=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
	if(dis.lt.0.1.or.dis.gt.2.0)   goto 300   !  find close CA

	

 400	read(unit=13,fmt=2, err=100,end=100) line
	if(line(1:4).eq.'ATOM') then
		backspace(unit=13)
		goto 410
	endif
	goto 400


 410	read(unit=13,fmt=6,err=100,end=100)
	1	atom,atnumber,atname4,altloc,resname,chain,resnumber4,icode,x4,y4,z4,
	1	occ, Bfact, ID, charge
	if(atname4(2:3).ne.'C ') go to 400           ! eg find C
	dis=sqrt((x3-x4)**2+(y3-y4)**2+(z3-z4)**2)
	if(dis.lt.0.1.or.dis.gt.2.0) go to 400   ! find close C


	call  dihedral(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dihed)

	dihed=dihed/pi*180.
	J_HNHA=7.97*(cos((dihed-60)*pi/180.))**2-1.26*(cos((dihed-60)*pi/180.))+ 0.63  ! Hz

	write(77,*) resnumber1, atname1, resnumber2, atname2, resnumber3, atname3, resnumber4, atname4,  dihed, J_HNHA


	rewind(12)
	rewind(11)
	rewind(13)
 	goto 100	

 1111	rewind(10)



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

