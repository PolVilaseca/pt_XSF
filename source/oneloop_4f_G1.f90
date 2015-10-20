




	MODULE OneLoop_4f_G1
	!
	! Contains the set of routines for the calculaton of the 1-loop
	! diagrams contributing to the G1_4f correlation functions (including
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
!!!
!!!	Subroutines for the G1_4f 1loop diagrams
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE G1_4f_1a(l,g_i1,g_i2,flav,G14f_1a)
	!
	! Subroutine to compute diagram G1_4f_1a 
	!
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_1a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(c0:c1) :: G14f_0
	real(kind=8) :: multi	


	!initialization
	G14f_1a(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G14f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call G1_4f_0(l,g_i1,g_i2,flav,G14f_0)

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

	G14f_1a(:) = -0.5d0*C2_R*loop*G14f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE G1_4f_1c(l,g_i1,g_i2,flav,G14f_1c)
	!
	! Subroutine to compute diagram G1_4f_1c 
	!
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_1c
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(c0:c1) :: G14f_0
	real(kind=8) :: multi	


	!initialization
	G14f_1c(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G14f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call G1_4f_0(l,g_i1,g_i2,flav,G14f_0)


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

	G14f_1c(:) = -0.5d0*C2_R*loop*G14f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_1c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G1_4f_2a(l,g_i1,g_i2,flav,G14f_2a)
	!
	! Subroutine to compute diagram G1_4f_2a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_2a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8) :: part1, part2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G14f_2a(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))	

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
  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,qloop,S_prop)

	   do x0=c0,c1

	     part1 = Trace4(Q5(:,:,fp),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))
	     loop(x0) = loop(x0) + D_prop(0,0,0,0)*multi*part1

	  enddo !x0

	enddo
	enddo
	enddo


	!construct the tree-level part outside the loop:
	
	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   part2 = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	   G14f_2a(x0) = loop(x0)*part2

	enddo

	G14f_2a(:) = C2_R*dimrep*dimrep*G14f_2a(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_2a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_2b(l,g_i1,g_i2,flav,G14f_2b)
	!
	! Subroutine to compute diagram G1_4f_2b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_2b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8) :: part1, part2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G14f_2b(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))	

	!Construct the trace in the loop
	do n1=0,l-1			
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,qloop,S_prop)

	   do x0=c0,c1

	     part1 = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))
	     loop(x0) = loop(x0) + D_prop(0,0,l-1,l-1)*multi*part1

	  enddo !x0

	enddo
	enddo
	enddo


	!construct the tree-level part outside the loop:
	
	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   part2 = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	   G14f_2b(x0) = part2*loop(x0)

	enddo

	G14f_2b(:) = C2_R*dimrep*dimrep*G14f_2b(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_2b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3a(l,g_i1,g_i2,flav,G14f_3a)
	!
	! Subroutine to compute diagram G1_4f_3a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
		
	    S_G(:,:,x0,s0) = matmul(S_prop_0(:,:,1,x0,flav(3)),matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(4))))

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
	  trace = Trace4(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,1,flav(4)),Q5(:,:,fp))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo


	  G14f_3a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3b(l,g_i1,g_i2,flav,G14f_3b)
	!
	! Subroutine to compute diagram G1_4f_3b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0 = 0,l
	do x0 = c0,c1		
	    S_G(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0,flav(3)),matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4))))
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

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0=c0,c1
	  trace = Trace4(Q5(:,:,fp),S_prop(:,:,1,s0,flav(3)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo


	  G14f_3b(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3c(l,g_i1,g_i2,flav,G14f_3c)
	!
	! Subroutine to compute diagram G1_4f_3c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0 = 0,l
	do x0 = c0,c1		
	    S_G(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0,flav(1)),matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2))))
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
	call Fermion_Propagator(l,nvec_2+qloop,S_prop)

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_2+qloop,-nvec_2,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0=c0,c1
	  trace = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,s0,flav(1)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_3c(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3d(l,g_i1,g_i2,flav,G14f_3d)
	!
	! Subroutine to compute diagram G1_4f_3d 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
		
	    S_G(:,:,x0,s0) = matmul(S_prop_0(:,:,l-1,x0,flav(1)),matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(2))))

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
	call Fermion_Propagator(l,nvec_2+qloop,S_prop)


	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_2,-nvec_2-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace4(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,l-1,flav(2)),Q5p(:,:,f))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_3d(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4a(l,g_i1,g_i2,flav,G14f_4a)
	!
	! Subroutine to compute diagram G1_4f_4a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

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

	   !Build S_G
	   do t0=0,l
	   do x0=c0,c1
	      S_G(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0,flav(3)),matmul(g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4))))
	   enddo
	   enddo

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)

	do x0=c0,c1
	  trace = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,s0,flav(3)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo

	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo


	  G14f_4a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4b(l,g_i1,g_i2,flav,G14f_4b)
	!
	! Subroutine to compute diagram G1_4f_4b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

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

	   !Build S_G
	   
	   do x0 = c0,c1
	   do s0=0,l
	      S_G(:,:,x0,s0) = matmul(S_prop(:,:,1,x0,flav(3)),matmul(g4f(:,:,g_i2),S_prop(:,:,x0,s0,flav(4))))
	   enddo
	   enddo


	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)

	do x0 = c0,c1
	  trace = Trace4(Q5(:,:,fp),S_G(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo





	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo


	  G14f_4b(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4c(l,g_i1,g_i2,flav,G14f_4c)
	!
	! Subroutine to compute diagram G1_4f_4c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_2,S_prop_0)

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
	   call Fermion_Propagator(l,nvec_2+qloop,S_prop)

	   !Build S_G
	   
	   do x0 = c0,c1
	   do s0=0,l
	      S_G(:,:,x0,s0) = matmul(S_prop(:,:,l-1,x0,flav(1)),matmul(g4f(:,:,g_i1),S_prop(:,:,x0,s0,flav(2))))
	   enddo
	   enddo


	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_2+qloop,-nvec_2,-qloop,s0,t0,u0,mu,V1mu)

	do x0 = c0,c1
	  trace = Trace4(Q5p(:,:,f),S_G(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,l-1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo





	!Prepare the right trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_4c(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4d(l,g_i1,g_i2,flav,G14f_4d)
	!
	! Subroutine to compute diagram G1_4f_4d 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_2,S_prop_0)

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
	   call Fermion_Propagator(l,nvec_2+qloop,S_prop)

	   !Build S_G
	   do t0=0,l
	   do x0=c0,c1
	      S_G(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0,flav(1)),matmul(g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2))))
	   enddo
	   enddo

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_2,-nvec_2-qloop,qloop,s0,t0,u0,mu,V1mu)

	do x0=c0,c1
	  trace = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,flav(1)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo





	!Prepare the right trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_4d(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5a(l,g_i1,g_i2,flav,G14f_5a)
	!
	! Subroutine to compute diagram G1_4f_5a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
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
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,flav(4)),matmul(matmul(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3))),matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(4)))))
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






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo


	  G14f_5a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5b(l,g_i1,g_i2,flav,G14f_5b)
	!
	! Subroutine to compute diagram G1_4f_5b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
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
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(3)),matmul(matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4))),matmul(Q5(:,:,fp),S_prop_0(:,:,1,s0,flav(3)))))
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

	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo



	  G14f_5b(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5c(l,g_i1,g_i2,flav,G14f_5c)
	!
	! Subroutine to compute diagram G1_4f_5c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
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
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace


	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(1)),matmul(matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2))),matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,flav(1)))))
	enddo	
	enddo
	enddo


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
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_2+qloop,S_prop)

	trace = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,nvec_2,-nvec_2-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_2+qloop,-nvec_2,-qloop,s0p,t0p,u0p,nu,V1nu)

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




	!Prepare the left trace
	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_5c(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5d(l,g_i1,g_i2,flav,G14f_5d)
	!
	! Subroutine to compute diagram G1_4f_5d 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
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
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),matmul(matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1))),matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(2)))))
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






	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_5d(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6a(l,g_i1,g_i2,flav,G14f_6a)
	!
	! Subroutine to compute diagram G1_4f_6a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6a

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
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,flav(4)),matmul(matmul(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3))),matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(4)))))
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






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo


	  G14f_6a(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6b(l,g_i1,g_i2,flav,G14f_6b)
	!
	! Subroutine to compute diagram G1_4f_6b
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6b

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
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(3)),matmul(matmul(g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4))),matmul(Q5(:,:,fp),S_prop_0(:,:,1,s0,flav(3)))))
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






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo


	  G14f_6b(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6c(l,g_i1,g_i2,flav,G14f_6c)
	!
	! Subroutine to compute diagram G1_4f_6c
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6c

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
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)



	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(1)),matmul(matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2))),matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,flav(1)))))
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






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_6c(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6d(l,g_i1,g_i2,flav,G14f_6d)
	!
	! Subroutine to compute diagram G1_4f_6d
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6d

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
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)



	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),matmul(matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1))),matmul(g4f(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(2)))))
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






	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_6d(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_7a(l,g_i1,g_i2,flav,G14f_7a)
	!
	! Subroutine to compute diagram G1_4f_7a 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_7a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
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
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_7a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,1,flav(4)),matmul(Q5(:,:,fp),S_prop_0(:,:,1,s0,flav(3))))
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

	  trace = Trace6(S_O5_S(:,:,t0p,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,s0p,flav(4)),V1nu(:,:))
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






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))

	enddo

	  G14f_7a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_7a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_7b(l,g_i1,g_i2,flav,G14f_7b)
	!
	! Subroutine to compute diagram G1_4f_7b 
	!
	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_7b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
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
	integer :: f,fp

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_7b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_S
	call Fermion_Propagator(l,nvec_2,S_prop_0)


	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,flav(1))))
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

	  trace = Trace6(S_O5_S(:,:,t0p,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,s0p,flav(2)),V1nu(:,:))
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






	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,fp),S_prop_0(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))

	enddo


	  G14f_7b(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_7b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_dm(l,g_i1,g_i2,flav,G14fdm)
	!
	! Subroutine to compute diagram G1_4f_dm
	!
	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14fdm
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: G14fdm_1, G14fdm_2, G14fdm_3, G14fdm_4
	complex(kind=8) :: part1, part2
	complex(kind=8),dimension(4,4) :: mat1
	integer, dimension(3) :: nvec
	integer :: f,x0,s0

	nvec(:) = 0

	!initialization
	G14fdm(:) = dcmplx(0.0d0,0.0d0)
	G14fdm_1(:) = dcmplx(0.0d0,0.0d0)
	G14fdm_2(:) = dcmplx(0.0d0,0.0d0)
	G14fdm_3(:) = dcmplx(0.0d0,0.0d0)
	G14fdm_4(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,nvec,S_prop)

	do x0 = c0,c1

	  f = f1f2(flav(2),flav(1))
	  part1 = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))

	  f = f1f2(flav(4),flav(3))
	  part2 = dcmplx(0.0d0,0.0d0)
	  do s0=0,l
	  part2 = part2 + Trace5(Q5(:,:,f),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,s0,flav(4)),S_prop(:,:,s0,1,flav(4)))
	  enddo
	  G14fdm_1(x0) = part1*part2

	  part2 = dcmplx(0.0d0,0.0d0)
	  do s0=0,l
	  part2 = part2 + Trace5(Q5(:,:,f),S_prop(:,:,1,s0,flav(3)),S_prop(:,:,s0,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))
	  enddo
	  G14fdm_2(x0) = part1*part2



	  f = f1f2(flav(4),flav(3))
	  part1 = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))

	  f = f1f2(flav(2),flav(1))
	  part2 = dcmplx(0.0d0,0.0d0)
	  do s0=0,l
	  part2=part2+Trace5(Q5p(:,:,f),S_prop(:,:,l-1,s0,flav(1)),S_prop(:,:,s0,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))
	  enddo
	  G14fdm_3(x0) = part1*part2


	  part2 = dcmplx(0.0d0,0.0d0)
	  do s0=0,l
	  part2=part2+Trace5(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,s0,flav(2)),S_prop(:,:,s0,l-1,flav(2)))
	  enddo
	  G14fdm_4(x0) = part1*part2

	enddo


	G14fdm(:) = -(dimrep**2)*(G14fdm_1(:)+G14fdm_2(:)+G14fdm_3(:)+G14fdm_4(:))

	END SUBROUTINE G1_4f_dm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G1_4f_zf(l,g_i1,g_i2,flav,G14fzf)
	!
	! Subroutine to compute diagram G1_4f_zf
	!
	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14fzf
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: G14fzf_1, G14fzf_2, G14fzf_3, G14fzf_4
	complex(kind=8) :: part1, part2
	complex(kind=8),dimension(4,4) :: mat1
	integer, dimension(3) :: nvec
	integer :: f,x0

	nvec(:) = 0

	!initialization
	G14fzf(:) = dcmplx(0.0d0,0.0d0)
	G14fzf_1(:) = dcmplx(0.0d0,0.0d0)
	G14fzf_2(:) = dcmplx(0.0d0,0.0d0)
	G14fzf_3(:) = dcmplx(0.0d0,0.0d0)
	G14fzf_4(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,nvec,S_prop)

	do x0 = c0,c1

	  f = f1f2(flav(2),flav(1))
	  part1 = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))

	  f = f1f2(flav(4),flav(3))
	  mat1(:,:) = matmul(S_prop(:,:,x0,0,flav(4)),S_prop(:,:,0,1,flav(4)))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,x0,l,flav(4)),S_prop(:,:,l,1,flav(4)))
	  part2 = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),mat1(:,:))

	  G14fzf_1(x0) = part1*part2


	  mat1(:,:) = matmul(S_prop(:,:,1,0,flav(3)),S_prop(:,:,0,x0,flav(3)))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,1,l,flav(3)),S_prop(:,:,l,x0,flav(3)))
	  part2 = Trace4(Q5(:,:,f),mat1(:,:),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))

	  G14fzf_2(x0) = part1*part2



	  f = f1f2(flav(4),flav(3))
	  part1 = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))

	  f = f1f2(flav(2),flav(1))
	  mat1(:,:) = matmul(S_prop(:,:,l-1,0,flav(1)),S_prop(:,:,0,x0,flav(1)))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,l-1,l,flav(1)),S_prop(:,:,l,x0,flav(1)))
	  part2 = Trace4(Q5p(:,:,f),mat1(:,:),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))

	  G14fzf_3(x0) = part1*part2


	  mat1(:,:) = matmul(S_prop(:,:,x0,0,flav(2)),S_prop(:,:,0,l-1,flav(2)))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,x0,l,flav(2)),S_prop(:,:,l,l-1,flav(2)))
	  part2 = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),mat1(:,:))

	  G14fzf_4(x0) = part1*part2

	enddo


	G14fzf(:) = -(dimrep**2)*(G14fzf_1(:)+G14fzf_2(:)+G14fzf_3(:)+G14fzf_4(:))

	END SUBROUTINE G1_4f_zf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G1_4f_ds(l,g_i1,g_i2,flav,G14fds)
	!
	! Subroutine to compute diagram G1_4f_ds
	!
	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14fds
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: G14fds_1, G14fds_2, G14fds_3, G14fds_4
	complex(kind=8) :: part1, part2
	integer, dimension(3) :: nvec
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: mat1
	complex(kind=8),dimension(4,4) :: Amat
	integer :: f,x0

	nvec(:) = 0

	!initialization
	G14fds(:) = dcmplx(0.0d0,0.0d0)
	G14fds_1(:) = dcmplx(0.0d0,0.0d0)
	G14fds_2(:) = dcmplx(0.0d0,0.0d0)
	G14fds_3(:) = dcmplx(0.0d0,0.0d0)
	G14fds_4(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	pe_plus(:) = nvec(:)*(2.0d0*pi)/(l*1.0d0) + theta/l
	ptwit(:) = sin(pe_plus(:))
	phat(:) = 2.0d0*sin(pe_plus(:)/2.0d0)

	Amat = imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)
	Amat = Amat + 0.5d0*II*( phat(1)**2 + phat(2)**2 + phat(3)**2 )


	call Fermion_Propagator(l,nvec,S_prop)

	do x0 = c0,c1

	  f = f1f2(flav(2),flav(1))
	  part1 = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))

	  f = f1f2(flav(4),flav(3))
	  mat1(:,:) = matmul(S_prop(:,:,x0,0,flav(4)),matmul(Amat(:,:),S_prop(:,:,0,1,flav(4))))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,x0,l,flav(4)),matmul(Amat(:,:),S_prop(:,:,l,1,flav(4))))
	  part2 = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),mat1(:,:))

	  G14fds_1(x0) = part1*part2


	  mat1(:,:) = matmul(S_prop(:,:,1,0,flav(3)),matmul(Amat(:,:),S_prop(:,:,0,x0,flav(3))))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,1,l,flav(3)),matmul(Amat(:,:),S_prop(:,:,l,x0,flav(3))))
	  part2 = Trace4(Q5(:,:,f),mat1(:,:),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))

	  G14fds_2(x0) = part1*part2



	  f = f1f2(flav(4),flav(3))
	  part1 = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))

	  f = f1f2(flav(2),flav(1))
	  mat1(:,:) = matmul(S_prop(:,:,l-1,0,flav(1)),matmul(Amat(:,:),S_prop(:,:,0,x0,flav(1))))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,l-1,l,flav(1)),matmul(Amat(:,:),S_prop(:,:,l,x0,flav(1))))
	  part2 = Trace4(Q5p(:,:,f),mat1(:,:),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))

	  G14fds_3(x0) = part1*part2


	  mat1(:,:) = matmul(S_prop(:,:,x0,0,flav(2)),matmul(Amat(:,:),S_prop(:,:,0,l-1,flav(2))))
	  mat1(:,:) = mat1(:,:) + matmul(S_prop(:,:,x0,l,flav(2)),matmul(Amat(:,:),S_prop(:,:,l,l-1,flav(2))))
	  part2 = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),mat1(:,:))

	  G14fds_4(x0) = part1*part2

	enddo


	G14fds(:) = -(dimrep**2)*(G14fds_1(:)+G14fds_2(:)+G14fds_3(:)+G14fds_4(:))

	END SUBROUTINE G1_4f_ds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





	END MODULE OneLoop_4f_G1