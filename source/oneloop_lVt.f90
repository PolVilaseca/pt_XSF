




	MODULE OneLoop_lVt
	!
	! Contains the set of routines for the calculaton of the 1-loop
	! diagrams contributing to the lVt correlation functions (including
	! the contribution of counterterms.
	!
	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators
	USE TreeLevel_lVt
	USE omp_lib

	implicit none

	CONTAINS


	SUBROUTINE lVt1a(l,f1,f2,lVt_1a)
	!
	! Subroutine to compute diagram lVt_1a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_1a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: lVt_0
	real(kind=8) :: multi	


	!initialization
	lVt_1a(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	lVt_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call lVt0(l,f1,f2,lVt_0)

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


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

	lVt_1a(:) = -0.5d0*C2_R*loop*lVt_0(:)/((l*1.0d0)**3)

	END SUBROUTINE lVt1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE lVt2(l,f1,f2,lVt_2)
	!
	! Subroutine to compute diagram lVt_2
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_2
	
	integer :: x0,n1,n2,n3,i
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1,3) :: loop_1, loop_2
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop	
	complex(kind=8) :: trace_1, trace_2
	complex(kind=8),dimension(3) :: q

	real(kind=8) :: multi	
	integer :: f,k

	!initialization
	lVt_2(:) = dcmplx(0.0d0,0.0d0)
	loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	f=f1f2(f2,f1)
	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,qloop,S_prop)

	   !Construction of the tree level f with momentum q: 

	     do k=1,3
	     do x0=c0,c1

	   trace_1 = 0.5d0*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,1,f2))*exp(-imag*q(k))
	   trace_2 = 0.5d0*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,1,f2))*exp(imag*q(k))

	       loop_1(x0,k) = loop_1(x0,k) + D_prop(0,0,0,0)*trace_1*multi
	       loop_2(x0,k) = loop_2(x0,k) + D_prop(0,0,0,0)*trace_2*multi

	     enddo
	     enddo

	enddo
	enddo
	enddo

	lVt_2(:) = C2_R*(loop_1(:,1)+loop_1(:,2)+loop_1(:,3) - loop_2(:,1)-loop_2(:,2)-loop_2(:,3))/(3.0d0*(l)**3)


	END SUBROUTINE lVt2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE lVt3a(l,f1,f2,lVt_3a)
	!
	! Subroutine to compute diagram lVt_3a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_3a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l,3) :: S_Gk_1, S_Gk_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f,k

	lVt_3a(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_1(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_2(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)


	do k=1,3
	do x0 = c0,c1
	do s0=0,l
		
	    S_Gk_1(:,:,x0,s0,k) = matmul(S_prop_0(:,:,1,x0,f1),matmul(Pkplus(:,:,k),S_prop_0(:,:,x0,s0,f2)))
	    S_Gk_2(:,:,x0,s0,k) = matmul(S_prop_0(:,:,1,x0,f1),matmul(Pkminus(:,:,k),S_prop_0(:,:,x0,s0,f2)))

	enddo
	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)


	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do k=1,3
	do x0 = c0,c1
	  trace_1 = Trace4(S_Gk_1(:,:,x0,s0,k),V1mu(:,:),S_prop(:,:,t0,1,f2),Qk(:,:,k,f))
	  trace_2 = Trace4(S_Gk_2(:,:,x0,s0,k),V1mu(:,:),S_prop(:,:,t0,1,f2),Qk(:,:,k,f))
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(0,mu,0,u0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(0,mu,0,u0)*trace_2*multi
	enddo
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	lVt_3a(:) = -C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3) - l_loop_2(:,1)-l_loop_2(:,2)-l_loop_2(:,3))/(6.0d0*l**3)


	END SUBROUTINE lVt3a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVt3b(l,f1,f2,lVt_3b)
	!
	! Subroutine to compute diagram lVt_3b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_3b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,3) :: S_Gk_1, S_Gk_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2 
	integer :: f,k

	lVt_3b(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_1(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_2(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_2 = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do k=1,3
	do t0 = 0,l
	do x0 = c0,c1		
	    S_Gk_1(:,:,t0,x0,k) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(Pkplus(:,:,k),S_prop_0(:,:,x0,1,f2)))
	    S_Gk_2(:,:,t0,x0,k) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(Pkminus(:,:,k),S_prop_0(:,:,x0,1,f2)))
	enddo
	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)
	
	do k=1,3
	do x0=c0,c1
	  trace_1 = Trace4(Qk(:,:,k,f),S_prop(:,:,1,s0,f1),V1mu(:,:),S_Gk_1(:,:,t0,x0,k))
	  trace_2 = Trace4(Qk(:,:,k,f),S_prop(:,:,1,s0,f1),V1mu(:,:),S_Gk_2(:,:,t0,x0,k))
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(mu,0,u0,0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(mu,0,u0,0)*trace_2*multi
	enddo
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  lVt_3b(:) = C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3) - l_loop_2(:,1)-l_loop_2(:,2)-l_loop_2(:,3))/(6.0d0*l**3)

	END SUBROUTINE lVt3b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE lVt4a(l,f1,f2,lVt_4a)
	!
	! Subroutine to compute diagram lVt_4a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_4a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop,S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,3) :: S_Gk_1, S_Gk_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu, i
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi 
	complex(kind=8) :: trace_1, trace_2
	integer :: f,k
	complex(kind=8),dimension(3) :: q

	lVt_4a(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_1(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_2(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!fermion propagator 0 momenta
	call Fermion_Propagator(l,p_ext,S_prop_0)

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,p_ext+qloop,S_prop)

	   !Build S_G
	   do k=1,3
	   do t0=0,l
	   do x0=c0,c1
	      S_Gk_1(:,:,t0,x0,k) = matmul(S_prop(:,:,t0,x0,f1),matmul(Pkplus(:,:,k),S_prop(:,:,x0,1,f2)))*exp(-imag*q(k))
	      S_Gk_2(:,:,t0,x0,k) = matmul(S_prop(:,:,t0,x0,f1),matmul(Pkminus(:,:,k),S_prop(:,:,x0,1,f2)))*exp(imag*q(k))
	   enddo
	   enddo
	   enddo

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	do k=1,3
	do x0=c0,c1
	  trace_1 = Trace4(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1),V1mu(:,:),S_Gk_1(:,:,t0,x0,k))
	  trace_2 = Trace4(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1),V1mu(:,:),S_Gk_2(:,:,t0,x0,k))
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(0,mu,0,u0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(0,mu,0,u0)*trace_2*multi
	enddo	!x0
	enddo


	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo

	lVt_4a(:) = -C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3) - l_loop_2(:,1)-l_loop_2(:,2)-l_loop_2(:,3))/(6.0d0*l**3)

	END SUBROUTINE lVt4a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVt4b(l,f1,f2,lVt_4b)
	!
	! Subroutine to compute diagram lVt_4b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_4b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l,3) :: S_Gk_1, S_Gk_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu, i
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f,k
	complex(kind=8),dimension(3) :: q

	lVt_4b(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_1(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_2(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!fermion propagator 0 momenta
	call Fermion_Propagator(l,p_ext,S_prop_0)

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,p_ext+qloop,S_prop)

	   !Build S_G
	   
	   do k=1,3
	   do x0 = c0,c1
	   do s0=0,l
	      S_Gk_1(:,:,x0,s0,k) = matmul(S_prop(:,:,1,x0,f1),matmul(Pkplus(:,:,k),S_prop(:,:,x0,s0,f2)))*exp(-imag*q(k))
	      S_Gk_2(:,:,x0,s0,k) = matmul(S_prop(:,:,1,x0,f1),matmul(Pkminus(:,:,k),S_prop(:,:,x0,s0,f2)))*exp(imag*q(k))
	   enddo
	   enddo
	   enddo


	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)

	do k=1,3
	do x0 = c0,c1
	  trace_1 = Trace4(Qk(:,:,k,f),S_Gk_1(:,:,x0,s0,k),V1mu(:,:),S_prop_0(:,:,t0,1,f2))
	  trace_2 = Trace4(Qk(:,:,k,f),S_Gk_2(:,:,x0,s0,k),V1mu(:,:),S_prop_0(:,:,t0,1,f2))
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(mu,0,u0,0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(mu,0,u0,0)*trace_2*multi
	enddo	!x0
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo

	lVt_4b(:) = C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3) - l_loop_2(:,1)-l_loop_2(:,2)-l_loop_2(:,3))/(6.0d0*l**3)

	END SUBROUTINE lVt4b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVt5a(l,f1,f2,lVt_5a)
	!
	! Subroutine to compute diagram lVt_5a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_5a
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l,3) :: S_Ok_SGS_1, S_Ok_SGS_2	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1, trace_2
	complex(kind=8),dimension(c0:c1,3) :: l_loop
	real(kind=8) :: multi
	integer :: f,k
	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Ok_SGS_1(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Ok_SGS_2(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	lVt_5a(:) = dcmplx(0.0d0,0.0d0)
	l_loop(:,:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do k=1,3
	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_Ok_SGS_1(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Qk(:,:,k,f),S_prop_0(:,:,1,x0,f1)),matmul(Pkplus(:,:,k),S_prop_0(:,:,x0,s0,f2))))
	    S_Ok_SGS_2(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Qk(:,:,k,f),S_prop_0(:,:,1,x0,f1)),matmul(Pkminus(:,:,k),S_prop_0(:,:,x0,s0,f2))))
	enddo	
	enddo
	enddo
	enddo


	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_Ok_SGS_1,S_Ok_SGS_2,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)



	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	
	call Fermion_Propagator(l,p+q,S_prop)

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  do k=1,3
	  trace_1 = Trace4(S_Ok_SGS_1(:,:,t0,x0,s0,k),V1mu(:,:),S_prop(:,:,t0p,s0p,f2),V1nu(:,:))
	  trace_2 = Trace4(S_Ok_SGS_2(:,:,t0,x0,s0,k),V1mu(:,:),S_prop(:,:,t0p,s0p,f2),V1nu(:,:))

	  f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(nu,mu,u0p,u0)

	  enddo


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

	lVt_5a(x0) = f_temp
	
	enddo	!x0

 	lVt_5a(:) = -dimrep*C2_R*lVt_5a(:)/((l**3)*6.0d0)

	
	END SUBROUTINE lVt5a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE lVt5b(l,f1,f2,lVt_5b)
	!
	! Subroutine to compute diagram lVt_5b 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_5b
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l,3) :: SGS_Ok_S_1, SGS_Ok_S_2	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1, trace_2
	complex(kind=8),dimension(c0:c1,3) :: l_loop
	real(kind=8) :: multi
	integer :: f,k
	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)



	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_Ok_S_1(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_Ok_S_2(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	lVt_5b(:) = dcmplx(0.0d0,0.0d0)
	l_loop(:,:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do k=1,3
	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    SGS_Ok_S_1(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(Pkplus(:,:,k),S_prop_0(:,:,x0,1,f2)),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1))))
	    SGS_Ok_S_2(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(Pkminus(:,:,k),S_prop_0(:,:,x0,1,f2)),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1))))
	enddo	
	enddo
	enddo
	enddo



	do x0 = c0,c1
	f_temp = dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,SGS_Ok_S_1,SGS_Ok_S_2,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	
	call Fermion_Propagator(l,p+q,S_prop)

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  do k=1,3
	  trace_1 = Trace4(SGS_Ok_S_1(:,:,t0,x0,s0,k),V1mu(:,:),S_prop(:,:,t0p,s0p,f1),V1nu(:,:))
	  trace_2 = Trace4(SGS_Ok_S_2(:,:,t0,x0,s0,k),V1mu(:,:),S_prop(:,:,t0p,s0p,f1),V1nu(:,:))
	  f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(nu,mu,u0p,u0)
	  enddo

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

	lVt_5b(x0) = f_temp
	
	enddo	!x0

 	lVt_5b(:) = -dimrep*C2_R*(lVt_5b(:))/((l**3)*6.0d0)

	
	END SUBROUTINE lVt5b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE lVt6a(l,f1,f2,lVt_6a)
	!
	! Subroutine to compute diagram lVt_6a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_6a
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l,3) :: S_Ok_SGkS_1, S_Ok_SGkS_2 	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1, trace_2
	complex(kind=8),dimension(c0:c1,3) :: l_loop
	real(kind=8) :: multi
	integer :: f,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Ok_SGkS_1(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Ok_SGkS_2(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	lVt_6a(:) = dcmplx(0.0d0,0.0d0)
	l_loop(:,:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do k=1,3
	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_Ok_SGkS_1(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Qk(:,:,k,f),S_prop_0(:,:,1,x0,f1)),matmul(Pkplus(:,:,k),S_prop_0(:,:,x0,s0,f2))))
	    S_Ok_SGkS_2(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Qk(:,:,k,f),S_prop_0(:,:,1,x0,f1)),matmul(Pkminus(:,:,k),S_prop_0(:,:,x0,s0,f2))))

	enddo	
	enddo
	enddo
	enddo


	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_Ok_SGkS_1,S_Ok_SGkS_2,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 

	  do mu=0,3	
!	  do nu=0,3

	  call V2(l,p,-p,q,-q,s0,t0,u0,u0,mu,V2mu)	!Remember that this guy has a delta(mu,nu) inside

	  do k=1,3
	  trace_1 = Trace2(S_Ok_SGkS_1(:,:,t0,x0,s0,k),V2mu(:,:))
	  trace_2 = Trace2(S_Ok_SGkS_2(:,:,t0,x0,s0,k),V2mu(:,:))

	  f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(mu,mu,u0,u0)	!There is a delta(mu,nu) in the vertex
	  enddo

!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	lVt_6a(x0) = f_temp
	enddo	!x0

 	lVt_6a(:) = dimrep*C2_R*lVt_6a(:)/((l**3)*6.0d0)

	
	END SUBROUTINE lVt6a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE lVt6b(l,f1,f2,lVt_6b)
	!
	! Subroutine to compute diagram lVt_6b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_6b
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l,3) :: SGkS_Ok_S_1, SGkS_Ok_S_2 	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1 ,trace_2
	real(kind=8) :: multi
	complex(kind=8),dimension(c0:c1,3) :: l_loop
	integer :: f,k
	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGkS_Ok_S_1(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGkS_Ok_S_2(:,:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	lVt_6b(:) = dcmplx(0.0d0,0.0d0)
	l_loop(:,:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do k=1,3
	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	      SGkS_Ok_S_1(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(Pkplus(:,:,k),S_prop_0(:,:,x0,1,f2)),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1))))
	      SGkS_Ok_S_2(:,:,t0,x0,s0,k) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(Pkminus(:,:,k),S_prop_0(:,:,x0,1,f2)),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1))))
	enddo	
	enddo
	enddo
	enddo


	do x0 = c0,c1
	f_temp=dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,SGkS_Ok_S_1,SGkS_Ok_S_2,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 

	  do mu=0,3	
!	  do nu=0,3

	  call V2(l,p,-p,q,-q,s0,t0,u0,u0,mu,V2mu)	!Remember that this guy has a delta(mu,nu) inside

	  do k=1,3
	  trace_1 = Trace2(SGkS_Ok_S_1(:,:,t0,x0,s0,k),V2mu(:,:))
	  trace_2 = Trace2(SGkS_Ok_S_2(:,:,t0,x0,s0,k),V2mu(:,:))
	  f_temp = f_temp + multi*(trace_1 - trace_2)*D_prop(mu,mu,u0,u0)	
	  enddo

!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	lVt_6b(x0) = f_temp
	enddo	!x0

 	lVt_6b(:) = dimrep*C2_R*(lVt_6b(:))/((l**3)*6.0d0)

	
	END SUBROUTINE lVt6b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE lVt7(l,f1,f2,lVt_7)
	!
	! Subroutine to compute diagram lVt_7 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_7
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,3) :: S_Ok_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	complex(kind=8),dimension(c0:c1,3) :: l_loop
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,qloop
	real(kind=8), dimension(3) :: q
	complex(kind=8) :: trace_1,trace_2
	real(kind=8) :: multi
	integer :: f,k

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Ok_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	l_loop(:,:) = dcmplx(0.0d0,0.0d0)
	lVt_7(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_S

	call Fermion_Propagator(l,p,S_prop_0)

	do k=1,3
	do s0=0,l
	do t0=0,l
	    S_Ok_S(:,:,t0,s0,k) = matmul(S_prop_0(:,:,t0,1,f2),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1)))
	enddo	
	enddo
	enddo



	do x0 = c0,c1
	f_temp = dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_Ok_S,Pkplus,Pkminus,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p+qloop,S_prop)


	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,p,-p-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+qloop,-p,-qloop,s0p,t0,u0p,nu,V1nu)

	  do k=1,3

	  trace_1 = Trace6(S_Ok_S(:,:,t0,s0,k),V1mu(:,:),S_prop(:,:,t0p,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,s0p,f2),V1nu(:,:))*exp(-imag*q(k))
	  trace_2 = Trace6(S_Ok_S(:,:,t0,s0,k),V1mu(:,:),S_prop(:,:,t0p,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,s0p,f2),V1nu(:,:))*exp(imag*q(k))
	  f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(nu,mu,u0p,u0)

	  enddo

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

	lVt_7(x0) = f_temp
	enddo	!x0


 	lVt_7(:) = -(1.0d0/6.0d0)*dimrep*C2_R*lVt_7(:)/((l**3)*1.0d0)

	
	END SUBROUTINE lVt7


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVt8(l,f1,f2,lVt_8)
	!
	! Subroutine to compute diagram lVt_8 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_8
	
	integer :: n1,n2,n3,x0,k,f
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1,3) :: loop_1,loop_2
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop	
	complex(kind=8) :: trace_1, trace_2
	integer,dimension(3) :: nvec
	real(kind=8) :: multi	

	f=f1f2(f2,f1)
	nvec(:) = 0

	!initialization
	lVt_8(:) = dcmplx(0.0d0,0.0d0)
	loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,nvec,S_Prop)
	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)
  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)

	   do x0 = c0,c1
	   do k=1,3

	   trace_1 = 0.5*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,1,f2))
	   trace_2 = 0.5*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,1,f2))
	   loop_1(x0,k) = loop_1(x0,k) + D_prop(k,k,x0,x0)*trace_1*multi
	   loop_2(x0,k) = loop_2(x0,k) + D_prop(k,k,x0,x0)*trace_2*multi

	   enddo
	   enddo
	
	enddo
	enddo
	enddo

	lVt_8(:) = -C2_R*(loop_1(:,1)+loop_1(:,2)+loop_1(:,3)-loop_2(:,1)-loop_2(:,2)-loop_2(:,3))/(6.0d0*(l)**3)


	END SUBROUTINE lVt8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE lVt9a(l,f1,f2,lVt_9a)
	!
	! Subroutine to compute diagram lVt_9a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_9a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,3) :: S_Gk_1, S_Gk_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f,k
	complex(kind=8),dimension(3) :: q

	lVt_9a(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do x0 = c0,c1
	do k=1,3			
 
	    S_Gk_1(:,:,x0,k) = matmul(Qk(:,:,k,f),matmul(S_prop_0(:,:,1,x0,f1),Pkplus(:,:,k)))
	    S_Gk_2(:,:,x0,k) = matmul(Qk(:,:,k,f),matmul(S_prop_0(:,:,1,x0,f1),Pkminus(:,:,k)))

	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)
	
	do k=1,3
	do x0 = c0,c1
	  trace_1 = Trace2(S_Gk_1(:,:,x0,k),S_prop(:,:,x0,1,f2))*exp(-imag*q(k)/2.0d0)
	  trace_2 = Trace2(S_Gk_2(:,:,x0,k),S_prop(:,:,x0,1,f2))*exp(imag*q(k)/2.0d0)
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(0,k,0,x0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(0,k,0,x0)*trace_2*multi
	enddo
	enddo

	enddo	!momenta
	enddo
	enddo

	lVt_9a(:) = -C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3)+l_loop_2(:,1)+l_loop_2(:,2)+l_loop_2(:,3))/(6.0d0*l**3)


	END SUBROUTINE lVt9a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVt9b(l,f1,f2,lVt_9b)
	!
	! Subroutine to compute diagram lVt_9b 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_9b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,3) :: S_Gk_1, S_Gk_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2 
	integer :: f,k
	complex(kind=8),dimension(3) :: q

	lVt_9b(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_Gk_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_2 = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do k=1,3
	do x0 = c0,c1		
	    S_Gk_1(:,:,x0,k) = matmul(matmul(Pkplus(:,:,k),S_prop_0(:,:,x0,1,f2)),Qk(:,:,k,f))
	    S_Gk_2(:,:,x0,k) = matmul(matmul(Pkminus(:,:,k),S_prop_0(:,:,x0,1,f2)),Qk(:,:,k,f))
	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)

	do k=1,3
	do x0=c0,c1
	  trace_1 = Trace2(S_prop(:,:,1,x0,f1),S_Gk_1(:,:,x0,k))*exp(-imag*q(k)/2.0d0)
	  trace_2 = Trace2(S_prop(:,:,1,x0,f1),S_Gk_2(:,:,x0,k))*exp(imag*q(k)/2.0d0)
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(k,0,x0,0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(k,0,x0,0)*trace_2*multi
	enddo
	enddo

	enddo	!momenta
	enddo
	enddo

	lVt_9b(:) = C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3)+l_loop_2(:,1)+l_loop_2(:,2)+l_loop_2(:,3))/(6.0d0*l**3)

	END SUBROUTINE lVt9b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVt10a(l,f1,f2,lVt_10a)
	!
	! Subroutine to compute diagram lVt_10a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_10a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,3) :: SQSG_1, SQSG_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f,k
	complex(kind=8),dimension(c0:c1) :: l_1, l_2
	real(kind=8),dimension(3) :: q

	lVt_10a(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	l_1(:) = dcmplx(0.0d0,0.0d0)
	l_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SQSG_1(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SQSG_2(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do k=1,3
	
	    SQSG_1(:,:,t0,x0,k) = matmul(S_prop_0(:,:,t0,1,f2),matmul(Qk(:,:,k,f),matmul(S_prop_0(:,:,1,x0,f1),Pkplus(:,:,k))))
	    SQSG_2(:,:,t0,x0,k) = matmul(S_prop_0(:,:,t0,1,f2),matmul(Qk(:,:,k,f),matmul(S_prop_0(:,:,1,x0,f1),Pkminus(:,:,k))))

	enddo
	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)


	
        do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	do k=1,3

	  trace_1 = Trace3(SQSG_1(:,:,t0,x0,k),S_prop(:,:,x0,s0,f2),V1mu(:,:))*exp(-imag*q(k)/2.0d0)
	  trace_2 = Trace3(SQSG_2(:,:,t0,x0,k),S_prop(:,:,x0,s0,f2),V1mu(:,:))*exp(imag*q(k)/2.0d0)
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(mu,k,u0,x0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(mu,k,u0,x0)*trace_2*multi

	enddo
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo
	
	l_1(:) = -C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3))/(6.0d0*l**3)
	l_2(:) = C2_R*dimrep*(l_loop_2(:,1)+l_loop_2(:,2)+l_loop_2(:,3))/(6.0d0*l**3)


	lVt_10a(:) = l_1(:) - l_2(:)


	END SUBROUTINE lVt10a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVt10b(l,f1,f2,lVt_10b)
	!
	! Subroutine to compute diagram lVt_10b 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_10b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l,3) :: GSQS_1, GSQS_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1,3) :: l_loop_1, l_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2 
	integer :: f,k
	complex(kind=8),dimension(3) :: q

	lVt_10b(:) = dcmplx(0.0d0,0.0d0)
	l_loop_1(:,:) = dcmplx(0.0d0,0.0d0)
	l_loop_2(:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	GSQS_1(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	GSQS_2(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_2 = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)


	do k=1,3
	do x0 = c0,c1		
	do s0 = 0,l
	    GSQS_1(:,:,x0,s0,k) = matmul(Pkplus(:,:,k),matmul(S_prop_0(:,:,x0,1,f2),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1))))
	    GSQS_2(:,:,x0,s0,k) = matmul(Pkminus(:,:,k),matmul(S_prop_0(:,:,x0,1,f2),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f1))))
	enddo
	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	 q(1)=(2.0d0*pi*n1)/(l*1.0d0)				
	 q(2)=(2.0d0*pi*n2)/(l*1.0d0)
	 q(3)=(2.0d0*pi*n3)/(l*1.0d0)

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)

	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do k=1,3
	do x0=c0,c1
	  trace_1 = Trace3(GSQS_1(:,:,x0,s0,k),V1mu(:,:),S_prop(:,:,t0,x0,f1))*exp(-imag*q(k)/2.0d0)
	  trace_2 = Trace3(GSQS_2(:,:,x0,s0,k),V1mu(:,:),S_prop(:,:,t0,x0,f1))*exp(imag*q(k)/2.0d0)
	  l_loop_1(x0,k) = l_loop_1(x0,k) + D_prop(k,mu,x0,u0)*trace_1*multi
	  l_loop_2(x0,k) = l_loop_2(x0,k) + D_prop(k,mu,x0,u0)*trace_2*multi
	enddo
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  lVt_10b(:) = -C2_R*dimrep*(l_loop_1(:,1)+l_loop_1(:,2)+l_loop_1(:,3)+l_loop_2(:,1)+l_loop_2(:,2)+l_loop_2(:,3))/(6.0d0*l**3)

	END SUBROUTINE lVt10b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVtdm(l,f1,f2,lVt_dm)
	!
	! Subroutine to compute diagram lVt_dm
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_dm
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1,3) :: lVt_dm_1, lVt_dm_2
	integer, dimension(3) :: p_ext
	integer :: f,s0,x0,k

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	lVt_dm(:) = dcmplx(0.0d0,0.0d0)
	lVt_dm_1(:,:) = dcmplx(0.0d0,0.0d0)
	lVt_dm_2(:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do k=1,3
	do x0 = c0,c1
	do s0 = 0,l

	lVt_dm_1(x0,k) = lVt_dm_1(x0,k)-dimrep*Trace5(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,s0,f2),S_prop(:,:,s0,1,f2))

	lVt_dm_2(x0,k) = lVt_dm_2(x0,k)-dimrep*Trace5(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,s0,f2),S_prop(:,:,s0,1,f2))

	enddo
	enddo
	enddo

	lVt_dm(:) = (lVt_dm_1(:,1) + lVt_dm_1(:,2) + lVt_dm_1(:,3) - lVt_dm_2(:,1)- lVt_dm_2(:,2)- lVt_dm_2(:,3))/3.0d0


	END SUBROUTINE lVtdm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVtzf(l,f1,f2,lVt_zf)
	!
	! Subroutine to compute diagram lVt_zf
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_zf
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: lVt1_zf_1, lVt1_zf_2
	complex(kind=8),dimension(c0:c1) :: lVt2_zf_1, lVt2_zf_2
	integer, dimension(3) :: p_ext
	integer :: f,x0,k

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	lVt_zf(:) = dcmplx(0.0d0,0.0d0)
	lVt1_zf_1(:) = dcmplx(0.0d0,0.0d0)
	lVt1_zf_2(:) = dcmplx(0.0d0,0.0d0)
	lVt2_zf_1(:) = dcmplx(0.0d0,0.0d0)
	lVt2_zf_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do k=1,3
	do x0 = c0,c1

	lVt1_zf_1(x0) = lVt1_zf_1(x0)+Trace5(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,0,f2),S_prop(:,:,0,1,f2))
	lVt1_zf_2(x0) = lVt1_zf_2(x0)+Trace5(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,l,f2),S_prop(:,:,l,1,f2))

	lVt2_zf_1(x0) = lVt2_zf_1(x0)+Trace5(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,0,f2),S_prop(:,:,0,1,f2))
	lVt2_zf_2(x0) = lVt2_zf_2(x0)+Trace5(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,l,f2),S_prop(:,:,l,1,f2))

	enddo
	enddo


	lVt_zf(:) = -dimrep*(lVt1_zf_1(:)+lVt1_zf_2(:)-lVt2_zf_1(:)-lVt2_zf_2(:))/3.0d0

	END SUBROUTINE lVtzf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lVtds(l,f1,f2,lVt_ds)
	!
	! Subroutine to compute diagram lVt_ds 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_ds
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: lVt1_ds_1, lVt1_ds_2
	complex(kind=8),dimension(c0:c1) :: lVt2_ds_1, lVt2_ds_2
	integer, dimension(3) :: p_ext
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: Amat
	integer :: f,x0,k

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	lVt_ds(:) = dcmplx(0.0d0,0.0d0)
	lVt1_ds_1(:) = dcmplx(0.0d0,0.0d0)
	lVt1_ds_2(:) = dcmplx(0.0d0,0.0d0)
	lVt2_ds_1(:) = dcmplx(0.0d0,0.0d0)
	lVt2_ds_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	pe_plus(:) = p_ext(:)*(2.0d0*pi)/(l*1.0d0) + theta/l
	ptwit(:) = sin(pe_plus(:))
	phat(:) = 2.0d0*sin(pe_plus(:)/2.0d0)

	Amat = imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)
	Amat = Amat + 0.5d0*II*( phat(1)**2 + phat(2)**2 + phat(3)**2 )



	call Fermion_Propagator(l,p_ext,S_prop)

	do k=1,3
	do x0 = c0,c1

	lVt1_ds_1(x0) = lVt1_ds_1(x0)+Trace6(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,0,f2),Amat,S_prop(:,:,0,1,f2))
	lVt1_ds_2(x0) = lVt1_ds_2(x0)+Trace6(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,l,f2),Amat,S_prop(:,:,l,1,f2))

	lVt2_ds_1(x0) = lVt2_ds_1(x0)+Trace6(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,0,f2),Amat,S_prop(:,:,0,1,f2))
	lVt2_ds_2(x0) = lVt2_ds_2(x0)+Trace6(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,l,f2),Amat,S_prop(:,:,l,1,f2))

	enddo
	enddo

	lVt_ds(:) = -dimrep*(lVt1_ds_1(:)+lVt1_ds_2(:)-lVt2_ds_1(:)-lVt2_ds_2(:))/3.0d0

	END SUBROUTINE lVtds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	END MODULE OneLoop_lVt


