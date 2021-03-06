



	MODULE GammaMatrices
	!
	! This module contains
	!	- Definitions of the gamma matrices
	!	- Definitions of different objects constructed as combinations
	! 	  of gamma matrices.
	!	- Functions for multiplying a matrix with a particular combination
	!	  of gamma matrices.
	!
	! Note: There might be some reundancy in the definitions of some objets.
	!       Some objets might dupicated. For example, all the gamma matrices
	!	are defined individually as g0, g1, g2 and g3, but then they are
	! 	collected again within a vector of gamma matrices called gVec.
	!	This will be left as it is, since some functions later on deppend
	!	on the individual gamma matrices while others deppend on the whole vector.
	!
	USE Parameters
	
	implicit none

	!Identity, gamma matrices and projectors
	integer,parameter :: nCor=4, spin=4	
	complex(kind=8),dimension(4,4) :: II, Pplus, Pminus	
	complex(kind=8),dimension(4,4) :: Qplus, Qminus		!XSF projectors
	complex(kind=8),dimension(4,4,4) :: Q5,Q5p		!XSF boundary gamma structure	
	complex(kind=8),dimension(4,4,3,4) :: Qk, Qkp			
	complex(kind=8),dimension(4,4,3) :: Pkplus, Pkminus			
	complex(kind=8),dimension(4,4) :: g0, g1, g2, g3, g5
	complex(kind=8),dimension(4,4) :: gA, gP, gV, gS
	complex(kind=8),dimension(4,4) :: g5_Pplus, g5_Pminus
	complex(kind=8),dimension(4,4,4) :: gX		!Vector with gamma matrices in fermion bilinears X
	complex(kind=8),dimension(4,4,16) :: gBil	!Vector with gamma matrices for all bilinears
	complex(kind=8),dimension(4,4,42) :: g4f	!Vector with gamma matrices the 4-fermion operators
	complex(kind=8),dimension(4,4,4,3) :: gY	!Vector with gamma matrices in fermion bilinears Yk
	complex(kind=8),dimension(4,4,0:3) ::  gVec	!Vector with gamma matrices g0, g1, g2 and g3
	complex(kind=8),dimension(4,4,0:3,0:3) :: sigma		!sigma tensor: 4-spin 4-spin 4-vector 4-vector

	CONTAINS

	SUBROUTINE GammaStructureConstruction
	!
	! Subroutine to initialise all the gamma matrices an everything 
	! constructed as a combination of gamma matrices.
	!
	implicit none
	integer :: i,j,k

	!We set everything to 0. 
	g0(:,:) = dcmplx(0.0d0,0.0d0)		
	g1(:,:) = dcmplx(0.0d0,0.0d0)
	g2(:,:) = dcmplx(0.0d0,0.0d0)
	g3(:,:) = dcmplx(0.0d0,0.0d0)
	g5(:,:) = dcmplx(0.0d0,0.0d0)
	gA(:,:) = dcmplx(0.0d0,0.0d0)
	gP(:,:) = dcmplx(0.0d0,0.0d0)
	gV(:,:) = dcmplx(0.0d0,0.0d0)
	gS(:,:) = dcmplx(0.0d0,0.0d0)
	Pplus(:,:) = dcmplx(0.0d0,0.0d0)	
	Pminus(:,:) = dcmplx(0.0d0,0.0d0)
	Qplus(:,:) = dcmplx(0.0d0,0.0d0)	
	Qminus(:,:) = dcmplx(0.0d0,0.0d0)
	II(:,:) = dcmplx(0.0d0,0.0d0)
	Q5(:,:,:) = dcmplx(0.0d0,0.0d0)
	Q5p(:,:,:) = dcmplx(0.0d0,0.0d0)
	Qk(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	Qkp(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	Pkplus(:,:,:) = dcmplx(0.0d0,0.0d0)
	Pkminus(:,:,:) = dcmplx(0.0d0,0.0d0)
	gX(:,:,:) = dcmplx(0.0d0,0.0d0)
	gBil(:,:,:) = dcmplx(0.0d0,0.0d0)
	g4f(:,:,:) = dcmplx(0.0d0,0.0d0)
	gY(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	g5_Pplus(:,:) = dcmplx(0.0d0,0.0d0)	
	g5_Pminus(:,:) = dcmplx(0.0d0,0.0d0)

	!We define explicitelly the non-zero matrix elements
	!Identity
	II(1,1) = dcmplx(1.0d0,0.0d0)
	II(2,2) = dcmplx(1.0d0,0.0d0)
	II(3,3) = dcmplx(1.0d0,0.0d0)
	II(4,4) = dcmplx(1.0d0,0.0d0)

	
	!Gamma0
	g0(1,1) = dcmplx(1.0d0,0.0d0)
	g0(2,2) = dcmplx(1.0d0,0.0d0)
	g0(3,3) = dcmplx(-1.0d0,0.0d0)
	g0(4,4) = dcmplx(-1.0d0,0.0d0)

	!Gamma1
	g1(4,1) = dcmplx(0.0d0,1.0d0)
	g1(3,2) = dcmplx(0.0d0,1.0d0)
	g1(2,3) = dcmplx(0.0d0,-1.0d0)
	g1(1,4) = dcmplx(0.0d0,-1.0d0)
	
	!Gamma2
	g2(4,1) = dcmplx(-1.0d0,0.0d0)
	g2(3,2) = dcmplx(1.0d0,0.0d0)
	g2(2,3) = dcmplx(1.0d0,0.0d0)
	g2(1,4) = dcmplx(-1.0d0,0.0d0)

	!Gamma3
	g3(3,1) = dcmplx(0.0d0,1.0d0)
	g3(4,2) = dcmplx(0.0d0,-1.0d0)
	g3(1,3) = dcmplx(0.0d0,-1.0d0)
	g3(2,4) = dcmplx(0.0d0,1.0d0)

	!Gamma5
!	g5(3,1) = dcmplx(-1.0d0,0.0d0)
!	g5(4,2) = dcmplx(-1.0d0,0.0d0)
!	g5(1,3) = dcmplx(-1.0d0,0.0d0)
!	g5(2,4) = dcmplx(-1.0d0,0.0d0)
	g5 = matmul(matmul(g0,g1),matmul(g2,g3))


	!Bilinears
	gA = matmul(g0,g5)
	gP = g5
	gV = g0
	gS = II

	!P+ projector	
!	Pplus(1,1) = dcmplx(1.0d0,0.0d0)
!	Pplus(2,2) = dcmplx(1.0d0,0.0d0)
	Pplus = 0.5d0*(II + g0)	

	!P- projector
!	Pminus(3,3) = dcmplx(1.0d0,0.0d0)
!	Pminus(4,4) = dcmplx(1.0d0,0.0d0)
	Pminus = 0.5d0*(II - g0)


	!Qplus
	Qplus = 0.5d0*(II + imag*matmul(g0,g5))

	!Qminus
	Qminus = 0.5d0*(II - imag*matmul(g0,g5))

	! Boundary at x0=0: 	
	Q5(:,:,1) = matmul(g0,matmul(g5,Qminus))	!uu
	Q5(:,:,2) = matmul(g5,Qplus)			!ud
	Q5(:,:,3) = matmul(g5,Qminus)			!du
	Q5(:,:,4) = matmul(g0,matmul(g5,Qplus))	!dd

	! Boundary at x0=T: 	
	Q5p(:,:,1) = -matmul(g0,matmul(g5,Qplus))	!uu
	Q5p(:,:,2) = matmul(g5,Qminus)			!ud
	Q5p(:,:,3) = matmul(g5,Qplus)			!du
	Q5p(:,:,4) = -matmul(g0,matmul(g5,Qminus))	!dd

	

	!A vector with gamma matrices: 
	gVec(:,:,0) = g0(:,:)
	gVec(:,:,1) = g1(:,:)
	gVec(:,:,2) = g2(:,:)
	gVec(:,:,3) = g3(:,:)


	do k=1,3

	! Boundary at x0=0: (correct projectors)		
	Qk(:,:,k,1) = matmul(gVec(:,:,k),Qminus)		!uu
	Qk(:,:,k,2) = matmul(g0,matmul(gVec(:,:,k),Qplus))	!ud
	Qk(:,:,k,3) = matmul(g0,matmul(gVec(:,:,k),Qminus))	!du
	Qk(:,:,k,4) = matmul(gVec(:,:,k),Qplus)		!dd

	! Boundary at x0=T: (correct projectors)		
	Qkp(:,:,k,1) = matmul(gVec(:,:,k),Qplus)		!uu
	Qkp(:,:,k,2) = -matmul(g0,matmul(gVec(:,:,k),Qminus))	!ud
	Qkp(:,:,k,3) = -matmul(g0,matmul(gVec(:,:,k),Qplus))	!du
	Qkp(:,:,k,4) = matmul(gVec(:,:,k),Qminus)		!dd

	Pkplus(:,:,k) = (II(:,:)+gVec(:,:,k))/2.0d0
	Pkminus(:,:,k) = (II(:,:)-gVec(:,:,k))/2.0d0

	enddo


	!the sigma tensor
	do i=0,3
	do j=0,3
	  sigma(:,:,i,j) = imag*(matmul(gVec(:,:,i),gVec(:,:,j)) - matmul(gVec(:,:,j),gVec(:,:,i)))/2.0d0
	enddo
	enddo


	!The vector with the gamma structure for the X bilinear fermion operators
	gX(:,:,1) = gA
	gX(:,:,2) = gP
	gX(:,:,3) = gV
	gX(:,:,4) = gS

	

	do k=1,3
	!The vector with the gamma structure for the X bilinear fermion operators
	gY(:,:,1,k) = matmul(gVec(:,:,k),g5)	!A_k
	gY(:,:,2,k) = gVec(:,:,k)			!V_K
	gY(:,:,3,k) = imag*sigma(:,:,k,0)			!T_k0
	gY(:,:,4,k) = imag*matmul(g5,sigma(:,:,k,0))	!Tt_k0
	
	enddo

	gBil(:,:,1) = gA
	gBil(:,:,2) = gP
	gBil(:,:,3) = gV
	gBil(:,:,4) = gS
	gBil(:,:,5) = gY(:,:,1,1)   !A_1
	gBil(:,:,6) = gY(:,:,1,2)   !A_2
	gBil(:,:,7) = gY(:,:,1,3)   !A_3
	gBil(:,:,8) = gY(:,:,2,1)   !V_1
	gBil(:,:,9) = gY(:,:,2,2)   !V_2
	gBil(:,:,10) = gY(:,:,2,3)   !V_3
	gBil(:,:,11) = gY(:,:,3,1)   !T_10
	gBil(:,:,12) = gY(:,:,3,2)   !T_20
	gBil(:,:,13) = gY(:,:,3,3)   !T_30
	gBil(:,:,14) = gY(:,:,4,1)   !Tt_10
	gBil(:,:,15) = gY(:,:,4,2)   !Tt_20
	gBil(:,:,16) = gY(:,:,4,3)   !Tt_30


	g4f(:,:,1) = II(:,:)
	g4f(:,:,2) = g5(:,:)
	g4f(:,:,3) = g0(:,:)
	g4f(:,:,4) = g1(:,:)
	g4f(:,:,5) = g2(:,:)
	g4f(:,:,6) = g3(:,:)
	g4f(:,:,7) = matmul(g0(:,:),g5(:,:))
	g4f(:,:,8) = matmul(g1(:,:),g5(:,:))
	g4f(:,:,9) = matmul(g2(:,:),g5(:,:))
	g4f(:,:,10) = matmul(g3(:,:),g5(:,:))
	g4f(:,:,11) = imag*sigma(:,:,0,0)
	g4f(:,:,12) = imag*sigma(:,:,0,1)
	g4f(:,:,13) = imag*sigma(:,:,0,2)
	g4f(:,:,14) = imag*sigma(:,:,0,3)
	g4f(:,:,15) = imag*sigma(:,:,1,0)
	g4f(:,:,16) = imag*sigma(:,:,1,1)
	g4f(:,:,17) = imag*sigma(:,:,1,2)
	g4f(:,:,18) = imag*sigma(:,:,1,3)
	g4f(:,:,19) = imag*sigma(:,:,2,0)
	g4f(:,:,20) = imag*sigma(:,:,2,1)
	g4f(:,:,21) = imag*sigma(:,:,2,2)
	g4f(:,:,22) = imag*sigma(:,:,2,3)
	g4f(:,:,23) = imag*sigma(:,:,3,0)
	g4f(:,:,24) = imag*sigma(:,:,3,1)
	g4f(:,:,25) = imag*sigma(:,:,3,2)
	g4f(:,:,26) = imag*sigma(:,:,3,3)
	g4f(:,:,27) = imag*matmul(g5(:,:),sigma(:,:,0,0))
	g4f(:,:,28) = imag*matmul(g5(:,:),sigma(:,:,0,1))
	g4f(:,:,29) = imag*matmul(g5(:,:),sigma(:,:,0,2))
	g4f(:,:,30) = imag*matmul(g5(:,:),sigma(:,:,0,3))
	g4f(:,:,31) = imag*matmul(g5(:,:),sigma(:,:,1,0))
	g4f(:,:,32) = imag*matmul(g5(:,:),sigma(:,:,1,1))
	g4f(:,:,33) = imag*matmul(g5(:,:),sigma(:,:,1,2))
	g4f(:,:,34) = imag*matmul(g5(:,:),sigma(:,:,1,3))
	g4f(:,:,35) = imag*matmul(g5(:,:),sigma(:,:,2,0))
	g4f(:,:,36) = imag*matmul(g5(:,:),sigma(:,:,2,1))
	g4f(:,:,37) = imag*matmul(g5(:,:),sigma(:,:,2,2))
	g4f(:,:,38) = imag*matmul(g5(:,:),sigma(:,:,2,3))
	g4f(:,:,39) = imag*matmul(g5(:,:),sigma(:,:,3,0))
	g4f(:,:,40) = imag*matmul(g5(:,:),sigma(:,:,3,1))
	g4f(:,:,41) = imag*matmul(g5(:,:),sigma(:,:,3,2))
	g4f(:,:,42) = imag*matmul(g5(:,:),sigma(:,:,3,3))



	!Useful for the fermion propagator
	g5_Pplus = matmul(g5,Pplus)
	g5_Pminus = matmul(g5,Pminus)

	END SUBROUTINE GammaStructureConstruction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	G-shuffling functions
!

	!Functions to compute products with gamma matrices

	FUNCTION gX_shuf(i,A)
	!
	! Multiply gX A
	!
	integer :: i
	complex(kind=8),dimension(4,4) :: A
	complex(kind=8),dimension(4,4) :: gX_shuf

	if(i.eq.1)then	!if gA
	  gX_shuf = gA_shuf(A)
	else if(i.eq.2)then !if gP
	  gX_shuf = gP_shuf(A)
	else if(i.eq.3)then !if gV
	  gX_shuf = gV_shuf(A)
	else if(i.eq.4)then !if gS
	  gX_shuf = A
	endif

	END FUNCTION gX_shuf



	FUNCTION gA_shuf(A)
	
	complex(kind=8),dimension(4,4) :: A
	complex(kind=8),dimension(4,4) :: gA_shuf
	!
	! Multiply gA A
	!
	gA_shuf(1,:) = - A(3,:)
	gA_shuf(2,:) = - A(4,:)
	gA_shuf(3,:) =  A(1,:)
	gA_shuf(4,:) =  A(2,:)

	END FUNCTION gA_shuf

	FUNCTION gP_shuf(A)
	!
	! Multiply gP A
	!	
	complex(kind=8),dimension(4,4) :: A
	complex(kind=8),dimension(4,4) :: gP_shuf

	gP_shuf(1,:) = - A(3,:)
	gP_shuf(2,:) = - A(4,:)
	gP_shuf(3,:) = - A(1,:)
	gP_shuf(4,:) = - A(2,:)

	END FUNCTION gP_shuf

	FUNCTION gV_shuf(A)
	!
	! Multiply gV A
	!	
	complex(kind=8),dimension(4,4) :: A
	complex(kind=8),dimension(4,4) :: gV_shuf

	gV_shuf(1,:) =  A(1,:)
	gV_shuf(2,:) =  A(2,:)
	gV_shuf(3,:) = - A(3,:)
	gV_shuf(4,:) = - A(4,:)

	END FUNCTION gV_shuf


	FUNCTION g5Ag5_shuf(A)
	!
	! Multiply g5 A g5
	!
	complex(kind=8),dimension(4,4) :: A
	complex(kind=8),dimension(4,4) :: g5Ag5_shuf

	g5Ag5_shuf(1,1) =  A(3,3)
	g5Ag5_shuf(1,2) =  A(3,4)
	g5Ag5_shuf(1,3) =  A(3,1)
	g5Ag5_shuf(1,4) =  A(3,2)

	g5Ag5_shuf(2,1) =  A(4,3)
	g5Ag5_shuf(2,2) =  A(4,4)
	g5Ag5_shuf(2,3) =  A(4,1)
	g5Ag5_shuf(2,4) =  A(4,2)

	g5Ag5_shuf(3,1) =  A(1,3)
	g5Ag5_shuf(3,2) =  A(1,4)
	g5Ag5_shuf(3,3) =  A(1,1)
	g5Ag5_shuf(3,4) =  A(1,2)

	g5Ag5_shuf(4,1) =  A(2,3)
	g5Ag5_shuf(4,2) =  A(2,4)
	g5Ag5_shuf(4,3) =  A(2,1)
	g5Ag5_shuf(4,4) =  A(2,2)

	END FUNCTION g5Ag5_shuf

		





	END MODULE GammaMatrices
