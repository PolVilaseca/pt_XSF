
	MODULE OneLoop_4f_Gk2
	!
	! Contains the set of routines for the calculaton of the 1-loop
	! diagrams contributing to the Gk2_4f correlation functions (including
	! the contribution of counterterms.
	!
	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators

	implicit none

	CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Diagrams for G2:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_1a(l,g_i1,g_i2,flav,Gk24f_1a)
	!
	! Subroutine to compute diagram Gk2_4f_1a 
	!
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_1a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: Gk24f_0
	real(kind=8) :: multi	


	!initialization
	Gk24f_1a(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	Gk24f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call Gk2_4f_0(l,g_i1,g_i2,flav,Gk24f_0)


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   loop = loop + D_prop(0,0,0,0)*multi

	enddo
	enddo
	enddo

	Gk24f_1a(:) = -0.5d0*C2_R*loop*Gk24f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE Gk2_4f_1c(l,g_i1,g_i2,flav,Gk24f_1c)
	!
	! Subroutine to compute diagram Gk2_4f_1c
	!
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_1c
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: Gk24f_0
	real(kind=8) :: multi	


	!initialization
	Gk24f_1c(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	Gk24f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call Gk2_4f_0(l,g_i1,g_i2,flav,Gk24f_0)


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   loop = loop + D_prop(0,0,l-1,l-1)*multi

	enddo
	enddo
	enddo

	Gk24f_1c(:) = -0.5d0*C2_R*loop*Gk24f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_1c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_2a(l,g_i1,g_i2,flav,Gk24f_2a)
	!
	! Subroutine to compute diagram Gk2_4f_2a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_2a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0, S_prop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,x0,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	Gk24f_2a(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

	f = f1f2(flav(4),flav(3))
	fp = f1f2(flav(2),flav(1))
        call Fermion_Propagator(l,nvec_2,S_Prop_0)
	
!	do k=1,3
	loop(:) = dcmplx(0.0d0,0.0d0)



	!Construct the trace in the loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	
	   call Fermion_Propagator(l,nvec_1+qloop,S_Prop)
  	   call Gluon_propagator(l,qloop,D_prop)	


	   do x0=c0,c1

	     mat1 = matmul(Qk(:,:,kdir,f),S_prop(:,:,1,x0,flav(3)))
	     mat2 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))
	     mat4 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,0,0)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	Gk24f_2a(:) = Gk24f_2a(:) + loop(:)

!	enddo !k

	Gk24f_2a(:) = -C2_R*Gk24f_2a(:)/(3.0d0*(l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_2a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_2b(l,g_i1,g_i2,flav,Gk24f_2b)
	!
	! Subroutine to compute diagram Gk2_4f_2b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_2b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp, x0,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	Gk24f_2b(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)


	f = f1f2(flav(4),flav(3))
	fp = f1f2(flav(2),flav(1))
	call Fermion_Propagator(l,nvec_1,S_Prop_0)

!	do k=1,3
	loop(:) = dcmplx(0.0d0,0.0d0)


	!Construct the trace in the loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	
	   
	   call Fermion_Propagator(l,nvec_2+qloop,S_Prop)
  	   call Gluon_propagator(l,qloop,D_prop)	


	   do x0=c0,c1

	     mat1 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))
	     mat2 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Qkp(:,:,kdir,fp),S_prop(:,:,l-1,x0,flav(1)))
	     mat4 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,l-1,l-1)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	Gk24f_2b(:) = Gk24f_2b(:) + loop(:)

!	enddo !k
	Gk24f_2b(:) = -C2_R*Gk24f_2b(:)/(3.0d0*(l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_2b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_3a(l,g_i1,g_i2,flav,Gk24f_3a)
	!
	! Subroutine to compute diagram Gk2_4f_3a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_3a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_3a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)


	do x0 = c0,c1
	do s0=0,l
		
	    mat1 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))
	    mat2 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat3 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))
	    mat4 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(4)))

	    S_G(:,:,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	Gk24f_3a(:) = Gk24f_3a(:) + g_loop(:)

!	enddo !k


	Gk24f_3a(:) = C2_R*dimrep*Gk24f_3a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_3a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_3b(l,g_i1,g_i2,flav,Gk24f_3b)
	!
	! Subroutine to compute diagram Gk2_4f_3b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_3b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_3b(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)


	do t0=0,l
	do x0 = c0,c1
		
	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(3)),g4f(:,:,g_i2))
	    mat2 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat3 = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	    mat4 = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))

	    S_G(:,:,t0,x0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_prop(:,:,1,s0,flav(3)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	Gk24f_3b(:) = Gk24f_3b(:) +g_loop(:)

!	enddo !k

	Gk24f_3b(:) = -C2_R*dimrep*Gk24f_3b(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_3b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_3c(l,g_i1,g_i2,flav,Gk24f_3c)
	!
	! Subroutine to compute diagram Gk2_4f_3c 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_3c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_3c(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(1)),g4f(:,:,g_i1))
	    mat2 = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))		
	    mat3 = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	    mat4 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))

	    S_G(:,:,t0,x0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_prop(:,:,l-1,s0,flav(1)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	Gk24f_3c(:) = Gk24f_3c(:) + g_loop(:)

!	enddo	!k
	Gk24f_3c(:) = C2_R*dimrep*Gk24f_3c(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_3c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_3d(l,g_i1,g_i2,flav,Gk24f_3d)
	!
	! Subroutine to compute diagram Gk2_4f_3d
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_3d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_3d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))
	    mat2 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))		
	    mat3 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))
	    mat4 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(2)))

	    S_G(:,:,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,l-1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	Gk24f_3d(:) = Gk24f_3d(:) + g_loop(:)

!	enddo!k

	Gk24f_3d(:) = -C2_R*dimrep*Gk24f_3d(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_3d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_4a(l,g_i1,g_i2,flav,Gk24f_4a)
	!
	! Subroutine to compute diagram Gk2_4f_4a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_4a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_4a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1
		
	    mat1 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat2 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),g4f(:,:,g_i1))

	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,1,flav(4)),Qk(:,:,kdir,f),S_prop_0(:,:,1,s0,flav(3)),V1mu(:,:),S_prop(:,:,t0,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_4a(:) = Gk24f_4a(:) + g_loop(:)

!	enddo !k

	Gk24f_4a(:) = C2_R*dimrep*Gk24f_4a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_4a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_4b(l,g_i1,g_i2,flav,Gk24f_4b)
	!
	! Subroutine to compute diagram Gk2_4f_4b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_4b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_4b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1

	    mat1 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat2 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))		
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),g4f(:,:,g_i1))

	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,s0,flav(4)),V1mu(:,:),S_prop_0(:,:,t0,1,flav(4)),Qk(:,:,kdir,f),S_prop(:,:,1,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	Gk24f_4b(:) = Gk24f_4b(:) +g_loop(:)

!	enddo !k

	Gk24f_4b(:) = -C2_R*dimrep*Gk24f_4b/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_4b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_4c(l,g_i1,g_i2,flav,Gk24f_4c)
	!
	! Subroutine to compute diagram Gk2_4f_4c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_4c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_4c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1
		
	    mat1 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))
	    mat2 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),g4f(:,:,g_i2))

	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,s0,flav(2)),V1mu(:,:),S_prop_0(:,:,t0,l-1,flav(2)),Qkp(:,:,kdir,fp),S_prop(:,:,l-1,x0,flav(1)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	Gk24f_4c(:) = Gk24f_4c(:) +g_loop(:)	

!	enddo!

	Gk24f_4c(:) = C2_R*dimrep*Gk24f_4c(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_4c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_4d(l,g_i1,g_i2,flav,Gk24f_4d)
	!
	! Subroutine to compute diagram Gk2_4f_4d
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_4d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_4d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1

	    mat1 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))
	    mat2 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))		
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),g4f(:,:,g_i2))

	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,s0,flav(1)),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	Gk24f_4d(:) = Gk24f_4d(:) +g_loop(:)
!	enddo !k

	Gk24f_4d(:) = -C2_R*dimrep*Gk24f_4d(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_4d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_5a(l,g_i1,g_i2,flav,Gk24f_5a)
	!
	! Subroutine to compute diagram Gk2_4f_5a
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_5a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_5a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Qk(:,:,kdir,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))		
	    mat3 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat4 = matmul(matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1)),S_prop_0(:,:,x0,s0,flav(4)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1
	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0p,x0,s0),V1mu(:,:),S_prop(:,:,t0,s0p,flav(4)),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)



	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO
	g_loop(x0) = f_temp

	enddo	!x0

	Gk24f_5a(:) = Gk24f_5a(:) + g_loop(:)

!	enddo	!k	

	Gk24f_5a(:) = C2_R*dimrep*Gk24f_5a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_5a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_5b(l,g_i1,g_i2,flav,Gk24f_5b)
	!
	! Subroutine to compute diagram Gk2_4f_5b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_5b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_5b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(3)),g4f(:,:,g_i2))
	    mat2 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))		
	    mat3 = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f)),S_prop_0(:,:,1,s0,flav(3)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0p,x0,s0),V1mu(:,:),S_prop(:,:,t0,s0p,flav(3)),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)



	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO
	g_loop(x0) = f_temp

	enddo	!x0

	Gk24f_5b(:) = Gk24f_5b(:) + g_loop(:)

!	enddo !k

	Gk24f_5b(:) = C2_R*dimrep*Gk24f_5b(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_5b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_5c(l,g_i1,g_i2,flav,Gk24f_5c)
	!
	! Subroutine to compute diagram Gk2_4f_5c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_5c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_5c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)


!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(1)),g4f(:,:,g_i1))
	    mat2 = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))		
	    mat3 = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp)),S_prop_0(:,:,l-1,s0,flav(1)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0p,x0,s0),V1mu(:,:),S_prop(:,:,t0,s0p,flav(1)),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)



	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	!$OMP END PARALLEL DO

	g_loop(x0) = f_temp

	enddo	!x0

	Gk24f_5c(:) = Gk24f_5c(:) + g_loop(:)

!	enddo !k

	Gk24f_5c(:) = C2_R*dimrep*Gk24f_5c(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_5c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_5d(l,g_i1,g_i2,flav,Gk24f_5d)
	!
	! Subroutine to compute diagram Gk2_4f_5d
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_5d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_5d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))		
	    mat3 = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))
	    mat4 = matmul(matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2)),S_prop_0(:,:,x0,s0,flav(2)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0p,x0,s0),V1mu(:,:),S_prop(:,:,t0,s0p,flav(2)),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)



	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO
	g_loop(x0) = f_temp

	enddo	!x0


	Gk24f_5d(:) = Gk24f_5d(:) + g_loop(:)

!	enddo!k


	Gk24f_5d(:) = C2_R*dimrep*Gk24f_5d(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_5d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_6a(l,g_i1,g_i2,flav,Gk24f_6a)
	!
	! Subroutine to compute diagram Gk2_4f_6a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_6a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_6a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)


!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)	

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Qk(:,:,kdir,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))		
	    mat3 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat4 = matmul(matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1)),S_prop_0(:,:,x0,s0,flav(4)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)
  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b	

	  do mu=0,3	

	  call V2(l,nvec_1,-nvec_1,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo
	!$OMP END PARALLEL DO
	g_loop(x0) = f_temp

	enddo	!x0

	Gk24f_6a(:) = Gk24f_6a(:) +g_loop(:)


!	enddo !k


	  Gk24f_6a(:) = -C2_R*dimrep*Gk24f_6a(:)/(3.0d0*(l**3))

	END SUBROUTINE Gk2_4f_6a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_6b(l,g_i1,g_i2,flav,Gk24f_6b)
	!
	! Subroutine to compute diagram Gk2_4f_6b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_6b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_6b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(3)),g4f(:,:,g_i2))
	    mat2 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))		
	    mat3 = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f)),S_prop_0(:,:,1,s0,flav(3)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	

	   trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b	

	  do mu=0,3	

	  call V2(l,nvec_1,-nvec_1,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	!$OMP END PARALLEL DO
	g_loop(x0) = f_temp

	enddo	!x0

	Gk24f_6b(:) = Gk24f_6b(:) + g_loop(:)

!	enddo !k

	  Gk24f_6b(:) = -C2_R*dimrep*Gk24f_6b(:)/(3.0d0*(l**3))

	END SUBROUTINE Gk2_4f_6b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_6c(l,g_i1,g_i2,flav,Gk24f_6c)
	!
	! Subroutine to compute diagram Gk2_4f_6c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_6c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_6c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)	

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(1)),g4f(:,:,g_i1))
	    mat2 = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))		
	    mat3 = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp)),S_prop_0(:,:,l-1,s0,flav(1)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_2,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1		
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	

	   trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b	

	  do mu=0,3	

	  call V2(l,nvec_2,-nvec_2,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	!$OMP END PARALLEL DO
	g_loop(x0) = f_temp

	enddo	!x0

	Gk24f_6c(:) = Gk24f_6c(:) +g_loop(:)

!	enddo !k

	  Gk24f_6c(:) = -C2_R*dimrep*Gk24f_6c(:)/(3.0d0*(l**3))

	END SUBROUTINE Gk2_4f_6c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_6d(l,g_i1,g_i2,flav,Gk24f_6d)
	!
	! Subroutine to compute diagram Gk2_4f_6d
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_6d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	Gk24f_6d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

!	do k=1,3
	g_loop = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))		
	    mat3 = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))
	    mat4 = matmul(matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2)),S_prop_0(:,:,x0,s0,flav(2)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_2,chunk,f,fp,flav,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1		
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b	

	  do mu=0,3	

	  call V2(l,nvec_2,-nvec_2,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	!$OMP END PARALLEL DO
	g_loop(x0) = f_temp


	enddo	!x0

	Gk24f_6d(:) = Gk24f_6d(:) + g_loop(:)

!	enddo !k


	  Gk24f_6d(:) = -C2_R*dimrep*Gk24f_6d(:)/(3.0d0*(l**3))

	END SUBROUTINE Gk2_4f_6d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_7a(l,g_i1,g_i2,flav,Gk24f_7a)
	!
	! Subroutine to compute diagram Gk2_4f_7a 
	!

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_7a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4,c0:c1) :: G_S_O_S_G
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k


	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_7a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G_S_O_S_G(:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)
	do x0 = c0,c1
	   mat1(:,:) = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	   mat2(:,:) = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	   G_S_O_S_G(:,:,x0) = matmul(mat1(:,:),matmul(Qkp(:,:,kdir,fp),mat2(:,:)))
	enddo

	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,1,flav(4)),matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,s0,flav(3))))
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1
	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_S,G_S_O_S_G,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(V1mu(:,:),S_O5_S(:,:,t0,s0p),V1nu(:,:),S_prop(:,:,t0p,x0,flav(3)),G_S_O_S_G(:,:,x0),S_prop(:,:,x0,s0,flav(4)))
	  f_temp = f_temp + multi*trace*D_prop(mu,nu,u0,u0p)

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	g_loop(x0) = f_temp


	enddo	!x0

	Gk24f_7a(:) = Gk24f_7a(:) + g_loop(:)

!	enddo 	!k


	  Gk24f_7a(:) = C2_R*dimrep*Gk24f_7a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_7a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_7b(l,g_i1,g_i2,flav,Gk24f_7b)
	!
	! Subroutine to compute diagram Gk2_4f_7b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_7b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4,c0:c1) :: G_S_O_S_G
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_7b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G_S_O_S_G(:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)
	do x0 = c0,c1
	   mat1(:,:) = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))
	   mat2(:,:) = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	   G_S_O_S_G(:,:,x0) = matmul(mat1(:,:),matmul(Qk(:,:,kdir,f),mat2(:,:)))
	enddo

	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,s0,flav(1))))
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_S,G_S_O_S_G,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(G_S_O_S_G(:,:,x0),S_prop(:,:,x0,s0p,flav(2)),V1nu(:,:),S_O5_S(:,:,t0p,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo
!$OMP END PARALLEL DO

	g_loop(x0) = f_temp
	enddo	!x0


	Gk24f_7b(:) = Gk24f_7b(:) + g_loop(:)

!	enddo !k

	  Gk24f_7b(:) = C2_R*dimrep*Gk24f_7b(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_7b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_8a(l,g_i1,g_i2,flav,Gk24f_8a)
	!
	! Subroutine to compute diagram Gk2_4f_8a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_8a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0, S_prop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,x0,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	Gk24f_8a(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f = f1f2(flav(4),flav(3))
	fp = f1f2(flav(2),flav(1))
    call Fermion_Propagator(l,nvec_2,S_Prop_0)

!	do k=1,3
	loop(:) = dcmplx(0.0d0,0.0d0)	

	!Construct the trace in the loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	
	   call Fermion_Propagator(l,nvec_1+qloop,S_Prop)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))
	     mat2 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Qkp(:,:,kdir,fp),S_prop(:,:,l-1,x0,flav(1)))
	     mat4 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + D_prop(0,0,0,l-1)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	Gk24f_8a(:) = Gk24f_8a(:) + loop(:)

!	enddo !k

	Gk24f_8a(:) = C2_R*dimrep*Gk24f_8a(:)/(3.0d0*(l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_8a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_8b(l,g_i1,g_i2,flav,Gk24f_8b)
	!
	! Subroutine to compute diagram Gk2_4f_8b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_8b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp, x0,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	Gk24f_8b(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f = f1f2(flav(4),flav(3))
	fp = f1f2(flav(2),flav(1))
	call Fermion_Propagator(l,nvec_1,S_Prop_0)

!	do k=1,3
	loop(:) =dcmplx(0.0d0,0.0d0)

	!Construct the trace in the loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	   call Fermion_Propagator(l,nvec_2+qloop,S_Prop)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Qk(:,:,kdir,f),S_prop(:,:,1,x0,flav(3)))
	     mat2 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))
	     mat4 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + D_prop(0,0,l-1,0)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo


	Gk24f_8b(:) = Gk24f_8b(:) +loop(:)

!	enddo !k

	Gk24f_8b(:) = C2_R*dimrep*Gk24f_8b(:)/(3.0d0*(l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_8b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_9a(l,g_i1,g_i2,flav,Gk24f_9a)
	!
	! Subroutine to compute diagram Gk2_4f_9a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_9a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0, S_prop_q, S_prop_mq
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,x0,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	Gk24f_9a(:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f = f1f2(flav(4),flav(3))
	fp = f1f2(flav(2),flav(1))
    call Fermion_Propagator(l,nvec_2,S_Prop_0)
	
!	do k=1,3
	loop(:) = dcmplx(0.0d0,0.0d0)

	!Construct the trace in the loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))
	     mat2 = matmul(g4f(:,:,g_i2),S_prop_mq(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))
	     mat4 = matmul(g4f(:,:,g_i1),S_prop_q(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + D_prop(0,0,0,l-1)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	Gk24f_9a(:) = Gk24f_9a(:) + loop(:)

!	enddo !k

	Gk24f_9a(:) = -C2_R*dimrep*Gk24f_9a(:)/(3.0d0*(l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_9a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_9b(l,g_i1,g_i2,flav,Gk24f_9b)
	!
	! Subroutine to compute diagram Gk2_4f_9b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_9b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop,mqloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0, S_prop_q, S_prop_mq
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp, x0,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	Gk24f_9b(:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f = f1f2(flav(4),flav(3))
	fp = f1f2(flav(2),flav(1))
	call Fermion_Propagator(l,nvec_1,S_Prop_0)

!	do k=1,3
	loop(:) = dcmplx(0.0d0,0.0d0)

	!Construct the trace in the loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	   call Fermion_Propagator(l,nvec_1+qloop,S_Prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_Prop_mq)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Qk(:,:,kdir,f),S_prop_mq(:,:,1,x0,flav(3)))
	     mat2 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Qkp(:,:,kdir,fp),S_prop_q(:,:,l-1,x0,flav(1)))
	     mat4 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + D_prop(0,0,0,l-1)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo


	Gk24f_9b(:) = Gk24f_9b(:) + loop(:)

!	enddo !k

	Gk24f_9b(:) = -C2_R*dimrep*Gk24f_9b(:)/(3.0d0*(l*1.0d0)**3)

	END SUBROUTINE Gk2_4f_9b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_10a(l,g_i1,g_i2,flav,Gk24f_10a)
	!
	! Subroutine to compute diagram Gk2_4f_10a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_10a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_10a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,x0,flav(3)))
	    mat2 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat3 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,s0,flav(1)))

	    OSGSOS(:,:,x0,s0) = matmul(matmul(mat1,mat2),mat3)

	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_10a(:) = Gk24f_10a(:) + g_loop(:)
!	enddo !k

	Gk24f_10a(:) = C2_R*dimrep*Gk24f_10a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_10a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_10b(l,g_i1,g_i2,flav,Gk24f_10b)
	!
	! Subroutine to compute diagram Gk2_4f_10b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_10b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_10b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	    mat3 = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))

	    OSGSOS(:,:,t0,x0) = matmul(matmul(mat1,mat2),mat3)

	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,t0,x0),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,s0,flav(2)),V1mu(:,:))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_10b(:) = Gk24f_10b(:) + g_loop(:)

!	enddo !k

	Gk24f_10b(:) = -C2_R*dimrep*Gk24f_10b(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_10b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_10c(l,g_i1,g_i2,flav,Gk24f_10c)
	!
	! Subroutine to compute diagram Gk2_4f_10c 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_10c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_10c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)


	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Qk(:,:,kdir,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	    mat3 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))

	    OSGSOS(:,:,t0,x0) = matmul(matmul(mat1,mat2),mat3)

	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,t0,x0),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,s0,flav(4)),V1mu(:,:))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_10c(:) = Gk24f_10c(:) +g_loop(:)

!	enddo !k

	Gk24f_10c(:) = C2_R*dimrep*Gk24f_10c(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_10c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_10d(l,g_i1,g_i2,flav,Gk24f_10d)
	!
	! Subroutine to compute diagram Gk2_4f_10d
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_10d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_10d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,x0,flav(1)))
	    mat2 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))
	    mat3 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,s0,flav(3)))

	    OSGSOS(:,:,x0,s0) = matmul(matmul(mat1,mat2),mat3)

	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_10d(:) = Gk24f_10d(:) +g_loop(:)

!	enddo !k

	Gk24f_10d(:) = -C2_R*dimrep*Gk24f_10d(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_10d




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_11a(l,g_i1,g_i2,flav,Gk24f_11a)
	!
	! Subroutine to compute diagram Gk2_4f_11a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_11a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_11a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	    SOSG(:,:,t0,x0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(Qk(:,:,kdir,f),matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1-qloop,-nvec_1,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_mq(:,:,x0,s0,flav(2)),V1mu(:,:),SOSG(:,:,t0,x0),S_prop_q(:,:,x0,1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo



	Gk24f_11a(:) = Gk24f_11a(:) + g_loop(:)

!	enddo !k


	Gk24f_11a(:) = C2_R*dimrep*Gk24f_11a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_11a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_11b(l,g_i1,g_i2,flav,Gk24f_11b)
	!
	! Subroutine to compute diagram Gk2_4f_11b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_11b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_11b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


!	do k=1,3 
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat2 = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,s0,flav(1)))
	    SOSG(:,:,x0,s0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(g4f(:,:,g_i1),matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_mq(:,:,1,x0,flav(3)),SOSG(:,:,x0,s0),V1mu(:,:),S_prop_q(:,:,t0,x0,flav(1)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_11b(:) = Gk24f_11b(:) + g_loop(:)

!	enddo !k
	Gk24f_11b(:) = -C2_R*dimrep*Gk24f_11b(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_11b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_11c(l,g_i1,g_i2,flav,Gk24f_11c)
	!
	! Subroutine to compute diagram Gk2_4f_11c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_11c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_11c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)


	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))
	    mat2 = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,s0,flav(3)))
	    SOSG(:,:,x0,s0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(g4f(:,:,g_i2),matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1+qloop,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_q(:,:,l-1,x0,flav(1)),SOSG(:,:,x0,s0),V1mu(:,:),S_prop_mq(:,:,t0,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_11c(:) = Gk24f_11c(:) + g_loop(:)

!	enddo !k


	Gk24f_11c(:) = C2_R*dimrep*Gk24f_11c(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_11c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_11d(l,g_i1,g_i2,flav,Gk24f_11d)
	!
	! Subroutine to compute diagram Gk2_4f_11d 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_11d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_11d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Qk(:,:,kdir,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	    SOSG(:,:,t0,x0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(Qkp(:,:,kdir,fp),matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_q(:,:,x0,s0,flav(4)),V1mu(:,:),SOSG(:,:,t0,x0),S_prop_mq(:,:,x0,l-1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo


	Gk24f_11d(:) = Gk24f_11d(:) + g_loop(:)

!	enddo !k

	Gk24f_11d(:) = -C2_R*dimrep*Gk24f_11d(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_11d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_12a(l,g_i1,g_i2,flav,Gk24f_12a)
	!
	! Subroutine to compute diagram Gk2_4f_12a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_12a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_12a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3

	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	do t0=0,l
	do x0=c0,c1

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,1,flav(4)),Qk(:,:,kdir,f))
	    mat2(:,:) = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	    mat3(:,:) = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Qkp(:,:,kdir,fp))

	    S_O5_S(:,:,t0,x0,s0) = matmul(mat1(:,:),matmul(mat2(:,:),matmul(mat3(:,:),S_prop_0(:,:,l-1,s0,flav(1)))))

	enddo
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_S,g_i1,g4f,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0p,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,s0p,flav(4)),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo
!$OMP END PARALLEL DO

	g_loop(x0) = f_temp

	enddo	!x0


	Gk24f_12a(:) = Gk24f_12a(:) + g_loop(:)

!	enddo !k

	  Gk24f_12a(:) = C2_R*dimrep*Gk24f_12a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_12a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Gk2_4f_12b(l,g_i1,g_i2,flav,Gk24f_12b)
	!
	! Subroutine to compute diagram Gk2_4f_12b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_12b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_12b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)


	do s0=0,l
	do t0=0,l
	do x0=c0,c1

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat2(:,:) = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	    mat3(:,:) = matmul(S_prop_0(:,:,x0,1,flav(4)),Qk(:,:,kdir,f))

	    S_O5_S(:,:,t0,x0,s0) = matmul(mat1(:,:),matmul(mat2(:,:),matmul(mat3(:,:),S_prop_0(:,:,1,s0,flav(3)))))

	enddo
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,S_O5_S,g_i2,g4f,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)



	do n1=0,l-1		
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0p,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,s0p,flav(2)),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	g_loop(x0) = f_temp


	enddo	!x0


	Gk24f_12b(:) = Gk24f_12b(:) + g_loop(:)

!	enddo !k


	  Gk24f_12b(:) = C2_R*dimrep*Gk24f_12b(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_12b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_13a(l,g_i1,g_i2,flav,Gk24f_13a)
	!
	! Subroutine to compute diagram Gk2_4f_13a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_13a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q, S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: SOSG1, SOSG2
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop,mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_13a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG1(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SOSG2(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)

	do t0=0,l
	do x0=c0,c1

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,1,flav(4)),Qk(:,:,kdir,f))
	    mat2(:,:) = matmul(S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2))
	    SOSG1(:,:,t0,x0) = matmul(mat1(:,:),mat2(:,:))

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Qkp(:,:,kdir,fp))
	    mat2(:,:) = matmul(S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1))
	    SOSG2(:,:,t0,x0) = matmul(mat1(:,:),mat2(:,:))

	enddo
	enddo	

	!Construct the gluon loop
	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,SOSG1,SOSG2,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1-qloop,-nvec_1,qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(SOSG1(:,:,t0,x0),S_prop_mq(:,:,x0,s0p,flav(2)),V1nu(:,:),SOSG2(:,:,t0p,x0),S_prop_q(:,:,x0,s0,flav(4)),V1mu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(mu,nu,u0,u0p)

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	g_loop(x0) = f_temp


	enddo	!x0


	Gk24f_13a(:) = Gk24f_13a(:) + g_loop(:)

!	enddo !k

	  Gk24f_13a(:) = C2_R*dimrep*Gk24f_13a(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_13a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_13b(l,g_i1,g_i2,flav,Gk24f_13b)
	!
	! Subroutine to compute diagram Gk2_4f_13b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f_13b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q, S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: SOSG1, SOSG2
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop,mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f_13b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG1(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SOSG2(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(3))
	fp=f1f2(flav(2),flav(1))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)

!	do k=1,3
	g_loop(:) = dcmplx(0.0d0,0.0d0)


	do x0=c0,c1
	do s0=0,l

	    mat1(:,:) = matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,1,flav(4)))
	    mat2(:,:) = matmul(Qk(:,:,kdir,f),S_prop_0(:,:,1,s0,flav(3)))
	    SOSG1(:,:,x0,s0) = matmul(mat1(:,:),mat2(:,:))

	    mat1(:,:) = matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat2(:,:) = matmul(Qkp(:,:,kdir,fp),S_prop_0(:,:,l-1,s0,flav(1)))
	    SOSG2(:,:,x0,s0) = matmul(mat1(:,:),mat2(:,:))
	enddo
	enddo	

	!Construct the gluon loop
	do x0 = c0,c1


	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,nvec_1,chunk,f,fp,flav,x0,SOSG1,SOSG2,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)


	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1,-nvec_1+qloop,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(SOSG1(:,:,x0,s0p),V1nu(:,:),S_prop_mq(:,:,t0p,x0,flav(3)),SOSG2(:,:,x0,s0),V1mu(:,:),S_prop_q(:,:,t0,x0,flav(1)))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	g_loop(x0) = f_temp

	enddo	!x0


	Gk24f_13b(:) = Gk24f_13b(:) +g_loop(:)

!	enddo !k

	  Gk24f_13b(:) = C2_R*dimrep*Gk24f_13b(:)/(3.0d0*l**3)

	END SUBROUTINE Gk2_4f_13b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_dm(l,g_i1,g_i2,flav,Gk24fdm)
	!
	! Subroutine to compute diagram Gk2_4f_dm
	!
	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24fdm
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: Gk24fdm_1, Gk24fdm_2, Gk24fdm_3, Gk24fdm_4
	complex(kind=8) :: part1, part2
	complex(kind=8),dimension(4,4) :: mat1,mat2,mat3,mat4
	integer, dimension(3) :: nvec
	integer :: f1,f2,x0,s0,k

	nvec(:) = 0

	!initialization
	Gk24fdm(:) = dcmplx(0.0d0,0.0d0)
	Gk24fdm_1(:) = dcmplx(0.0d0,0.0d0)
	Gk24fdm_2(:) = dcmplx(0.0d0,0.0d0)
	Gk24fdm_3(:) = dcmplx(0.0d0,0.0d0)
	Gk24fdm_4(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(2),flav(1))
	f2 = f1f2(flav(4),flav(3))

	call Fermion_Propagator(l,nvec,S_prop)

!	do k=1,3

	Gk24fdm_1(:) = dcmplx(0.0d0,0.0d0)
	Gk24fdm_2(:) = dcmplx(0.0d0,0.0d0)
	Gk24fdm_3(:) = dcmplx(0.0d0,0.0d0)
	Gk24fdm_4(:) = dcmplx(0.0d0,0.0d0)


	do x0 = c0,c1

	   do s0=0,l	
 
	   mat1 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat2 = matmul(g4f(:,:,g_i1),matmul(S_prop(:,:,x0,s0,flav(4)),S_prop(:,:,s0,1,flav(4))))
	   mat3 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat4 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))

	   Gk24fdm_1(x0) = Gk24fdm_1(x0) + Trace4(mat1,mat2,mat3,mat4)

	   enddo

	   do s0=0,l	
 
	   mat1 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat2 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))
	   mat3 = matmul(Qk(:,:,kdir,f2),matmul(S_prop(:,:,1,s0,flav(3)),S_prop(:,:,s0,x0,flav(3))))
	   mat4 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))

	   Gk24fdm_2(x0) = Gk24fdm_2(x0) + Trace4(mat1,mat2,mat3,mat4)

	   enddo

	   do s0=0,l	
 
	   mat1 = matmul(Qkp(:,:,kdir,f1),matmul(S_prop(:,:,l-1,s0,flav(1)),S_prop(:,:,s0,x0,flav(1))))
	   mat2 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))
	   mat3 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat4 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))

	   Gk24fdm_3(x0) = Gk24fdm_3(x0) + Trace4(mat1,mat2,mat3,mat4)

	   enddo

	   do s0=0,l	
 
	   mat1 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat2 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))
	   mat3 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat4 = matmul(g4f(:,:,g_i2),matmul(S_prop(:,:,x0,s0,flav(2)),S_prop(:,:,s0,l-1,flav(2))))

	   Gk24fdm_4(x0) = Gk24fdm_4(x0) + Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo


	Gk24fdm(:) = Gk24fdm(:) + dimrep*(Gk24fdm_1(:)+Gk24fdm_2(:)+Gk24fdm_3(:)+Gk24fdm_4(:))

!	enddo !k

	Gk24fdm(:) = Gk24fdm(:)/3.0d0


	END SUBROUTINE Gk2_4f_dm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_zf(l,g_i1,g_i2,flav,Gk24fzf)
	!
	! Subroutine to compute diagram Gk2_4f_zf 
	!
	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24fzf
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: Gk24fzf_1, Gk24fzf_2, Gk24fzf_3, Gk24fzf_4
	complex(kind=8) :: part1, part2
	complex(kind=8),dimension(4,4) :: mat1,mat2,mat3,mat4,mat5
	integer, dimension(3) :: nvec
	integer :: f1,f2,x0,k

	nvec(:) = 0

	!initialization
	Gk24fzf(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf_1(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf_2(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf_3(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf_4(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(2),flav(1))
	f2 = f1f2(flav(4),flav(3))

	call Fermion_Propagator(l,nvec,S_prop)

!	do k=1,3
	Gk24fzf_1(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf_2(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf_3(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf_4(:) = dcmplx(0.0d0,0.0d0)


	do x0 = c0,c1
 
	   mat1 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat2 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	   mat3 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat4 = matmul(g4f(:,:,g_i1),matmul(S_prop(:,:,x0,0,flav(4)),S_prop(:,:,0,1,flav(4))))
	   mat5 = matmul(g4f(:,:,g_i1),matmul(S_prop(:,:,x0,l,flav(4)),S_prop(:,:,l,1,flav(4))))

	   Gk24fzf_1(x0) = Gk24fzf_1(x0) + Trace4(mat1,mat2,mat3,mat4+mat5)


	   mat1 = matmul(Qk(:,:,kdir,f2),matmul(S_prop(:,:,1,0,flav(3)),S_prop(:,:,0,x0,flav(3))))
	   mat2 = matmul(Qk(:,:,kdir,f2),matmul(S_prop(:,:,1,l,flav(3)),S_prop(:,:,l,x0,flav(3))))
	   mat3 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	   mat4 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat5 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	   Gk24fzf_2(x0) = Gk24fzf_2(x0) + Trace4(mat1+mat2,mat3,mat4,mat5)


	   mat1 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat2 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	   mat3 = matmul(Qkp(:,:,kdir,f1),matmul(S_prop(:,:,l-1,0,flav(1)),S_prop(:,:,0,x0,flav(1))))
	   mat4 = matmul(Qkp(:,:,kdir,f1),matmul(S_prop(:,:,l-1,l,flav(1)),S_prop(:,:,l,x0,flav(1))))
	   mat5 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	   Gk24fzf_3(x0) = Gk24fzf_3(x0) + Trace4(mat1,mat2,mat3+mat4,mat5)


	   mat1 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat2 = matmul(g4f(:,:,g_i2),matmul(S_prop(:,:,x0,l,flav(2)),S_prop(:,:,l,l-1,flav(2))))
	   mat3 = matmul(g4f(:,:,g_i2),matmul(S_prop(:,:,x0,0,flav(2)),S_prop(:,:,0,l-1,flav(2))))
	   mat4 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat5 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	   Gk24fzf_4(x0) = Gk24fzf_4(x0) + Trace4(mat1,mat2+mat3,mat4,mat5)
	


	enddo


	Gk24fzf(:) = Gk24fzf(:) + dimrep*(Gk24fzf_1(:)+Gk24fzf_2(:)+Gk24fzf_3(:)+Gk24fzf_4(:))

!	enddo !k

	Gk24fzf(:) = Gk24fzf(:)/3.0d0

	END SUBROUTINE Gk2_4f_zf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_ds(l,g_i1,g_i2,flav,Gk24fds)
	!
	! Subroutine to compute diagram Gk2_4f_ds
	!
	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24fds
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: Gk24fds_1, Gk24fds_2, Gk24fds_3, Gk24fds_4
	complex(kind=8) :: part1, part2
	complex(kind=8),dimension(4,4) :: mat1,mat2,mat3,mat4,mat5
	integer, dimension(3) :: nvec
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: Amat
	integer :: f1,f2,x0,k

	nvec(:) = 0

	!initialization
	Gk24fds(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds_1(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds_2(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds_3(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds_4(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(2),flav(1))
	f2 = f1f2(flav(4),flav(3))

	pe_plus(:) = nvec(:)*(2.0d0*pi)/(l*1.0d0) + theta/l
	ptwit(:) = sin(pe_plus(:))
	phat(:) = 2.0d0*sin(pe_plus(:)/2.0d0)

	Amat = imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)
	Amat = Amat + 0.5d0*II*( phat(1)**2 + phat(2)**2 + phat(3)**2 )

	call Fermion_Propagator(l,nvec,S_prop)

!	do k=1,3

	Gk24fds_1(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds_2(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds_3(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds_4(:) = dcmplx(0.0d0,0.0d0)

	do x0 = c0,c1
 
	   mat1 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat2 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	   mat3 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat4 = matmul(g4f(:,:,g_i1),matmul(S_prop(:,:,x0,0,flav(4)),matmul(Amat(:,:),S_prop(:,:,0,1,flav(4)))))
	   mat5 = matmul(g4f(:,:,g_i1),matmul(S_prop(:,:,x0,l,flav(4)),matmul(Amat(:,:),S_prop(:,:,l,1,flav(4)))))

	   Gk24fds_1(x0) = Gk24fds_1(x0) + Trace4(mat1,mat2,mat3,mat4+mat5)


	   mat1 = matmul(Qk(:,:,kdir,f2),matmul(S_prop(:,:,1,0,flav(3)),matmul(Amat(:,:),S_prop(:,:,0,x0,flav(3)))))
	   mat2 = matmul(Qk(:,:,kdir,f2),matmul(S_prop(:,:,1,l,flav(3)),matmul(Amat(:,:),S_prop(:,:,l,x0,flav(3)))))
	   mat3 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	   mat4 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat5 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	   Gk24fds_2(x0) = Gk24fds_2(x0) + Trace4(mat1+mat2,mat3,mat4,mat5)


	   mat1 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat2 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))
	   mat3 = matmul(Qkp(:,:,kdir,f1),matmul(S_prop(:,:,l-1,0,flav(1)),matmul(Amat(:,:),S_prop(:,:,0,x0,flav(1)))))
	   mat4 = matmul(Qkp(:,:,kdir,f1),matmul(S_prop(:,:,l-1,l,flav(1)),matmul(Amat(:,:),S_prop(:,:,l,x0,flav(1)))))
	   mat5 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	   Gk24fds_3(x0) = Gk24fds_3(x0) + Trace4(mat1,mat2,mat3+mat4,mat5)


	   mat1 = matmul(Qk(:,:,kdir,f2),S_prop(:,:,1,x0,flav(3)))
	   mat2 = matmul(g4f(:,:,g_i2),matmul(S_prop(:,:,x0,l,flav(2)),matmul(Amat(:,:),S_prop(:,:,l,l-1,flav(2)))))
	   mat3 = matmul(g4f(:,:,g_i2),matmul(S_prop(:,:,x0,0,flav(2)),matmul(Amat(:,:),S_prop(:,:,0,l-1,flav(2)))))
	   mat4 = matmul(Qkp(:,:,kdir,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat5 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))

	   Gk24fds_4(x0) = Gk24fds_4(x0) + Trace4(mat1,mat2+mat3,mat4,mat5)
	


	enddo

	Gk24fds(:) = Gk24fds(:) + dimrep*(Gk24fds_1(:)+Gk24fds_2(:)+Gk24fds_3(:)+Gk24fds_4(:))

!	enddo !k

	Gk24fds(:) = Gk24fds(:)/3.0d0


	END SUBROUTINE Gk2_4f_ds



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE OneLoop_4f_Gk2
