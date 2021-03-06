



	MODULE Vertices
	!
	! Module containing subroutines for the calculation of the different vertices
	! appearing at 1-loop order.
	!
	USE GammaMatrices
	USE Input
	Use Parameters

	implicit none
	CONTAINS





      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE V_1_mu_scalar(s,x0,y0,z0,mu,V1mu)			!4-spinor, 4-spinor
!	-For the XSF.
!	-Computes, independently from the color structure, the fermion-fermion-
!	 gluon vertex for the Wilson fermion action, at a given time indeces,
!	 for a given momenta and for a given vector component.
!
!	-Input: -s: momenta
!		-x0,y0,z0: temporal indeces
!		-mu: space-time vectorial component
!		
!	-Output: -V1mu: the vertex, 4x4 spinor matrix
	implicit none

	integer,intent(in) :: x0,y0,z0,mu
	real(kind=8),dimension(3),intent(in) :: s
	complex(kind=8),dimension(4,4),intent(out) :: V1mu	!4-spinor, 4-spinor
	real(kind=8),dimension(3) :: s_plus
	
	!initialize:
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	if(mu.eq.0)then

		V1mu(:,:) = -Pminus*delta(x0+1,y0)*delta(x0,z0)*(1-delta(x0,lat_size)) + Pplus*delta(x0-1,y0)*delta(y0,z0)*(1-delta(x0,0))

	else	
	!k components

		s_plus(mu) = s(mu) + theta/lat_size
		V1mu(:,:) = delta(x0,y0)*delta(x0,z0)*(gVec(:,:,mu)*cos(s_plus(mu))-imag*II*sin(s_plus(mu)))
		V1mu(:,:) = V1mu(:,:)*(1-delta(x0,0)-delta(x0,lat_size))

	endif


	END SUBROUTINE V_1_mu_scalar







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE S_1_mu_scalar(q,x0,z0,mu,S1mu)
!	-For the XSF.
!	-Computes, independently from the color structure, the fermion-fermion-
!	 gluon vertex for the Clover part of the fermion action, 
!	 at a given time indeces,
!	 for a given momenta and for a given vector component.
!
!	-Input: -s: momenta
!		-x0,y0,z0: temporal indeces
!		-mu: space-time vectorial component
!		
!	-Output: -S1mu: the vertex, 4x4 spinor matrix

	implicit none

	integer,intent(in) :: x0,z0,mu
	real(kind=8),dimension(3),intent(in) :: q
	complex(kind=8),dimension(4,4),intent(out) :: S1mu	!4-spinor, 4-spinor
	real(kind=8),dimension(3) :: qdot

	integer :: k
	complex(kind=8) :: piece1
	complex(kind=8),dimension(4,4) :: piece2, piece3

	!initialize
	S1mu(:,:) = dcmplx(0.0d0,0.0d0)
	piece1 = dcmplx(0.0d0,0.0d0)
	piece2(:,:) = dcmplx(0.0d0,0.0d0)
	piece3(:,:) = dcmplx(0.0d0,0.0d0)

	!build the qdot
	do k=1,3
	  qdot(k)=sin(q(k))
	enddo

	!0 component
	if(mu.eq.0)then

		S1mu(:,:) = (delta(x0,z0)+delta(x0-1,z0))*(sigma(:,:,1,0)*qdot(1)+sigma(:,:,2,0)*qdot(2)+sigma(:,:,3,0)*qdot(3))/4.0d0
	!	S1mu(:,:) = S1mu(:,:)*(1 - delta(x0,T) - delta(x0,0))	!Removing from the boundaries

	else
	!k components
	  	piece1 = 0.5d0*cos(q(mu)/2.0d0)
	  	piece2(:,:) = imag*(delta(x0+1,z0)-delta(x0-1,z0))*sigma(:,:,0,mu)/2.0d0
	  	piece3(:,:) = delta(x0,z0)*(sigma(:,:,1,mu)*qdot(1)+sigma(:,:,2,mu)*qdot(2)+sigma(:,:,3,mu)*qdot(3))
	  	S1mu(:,:) = piece1*(piece2(:,:) + piece3(:,:))!*(1 - delta(x0,T) - delta(x0,0))	!Removing from the boundaries

	endif

	END SUBROUTINE S_1_mu_scalar 


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE V1(l,nvec_pp,nvec_p,nvec_q,s0,t0,u0,mu,VS1mu)	!
!	-For the XSF.
!	-Computes, independently from the color structure, the fermion-fermion-
!	 gluon vertex for the full fermion action, at a given time indeces,
!	 for given momenta and for a given vector component.
!	-For the computation it calls V_1_mu_scalar and S_1_mu_scalar and
!	 recombines its results
!
!	-Input: -l: lattice size
!		-nvec_pp,nvec_p,nvec_q: momenta ineces
!		-s0,t0,u0: temporal indeces
!		-mu: space-time vectorial component
!		
!	-Output: -Vs1mu: the vertex, 4x4 spinor matrix

	implicit none

	integer,intent(in) :: l
	integer,intent(in) :: s0,t0,u0,mu
	integer,dimension(3),intent(in) :: nvec_p,nvec_pp,nvec_q
	complex(kind=8),dimension(4,4),intent(out) :: VS1mu 	!4-spinor, 4-spinor
	complex(kind=8),dimension(4,4) :: V1mu, S1mu		!4-spinor, 4-spinor
	
	integer :: i
	real(kind=8),dimension(3) :: pp, p, q
	!initialize
	VS1mu(:,:) = dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	S1mu(:,:) = dcmplx(0.0d0,0.0d0)
	

	

	 pp(1)=(2.0d0*pi*nvec_pp(1))/(l*1.0d0)				
	 pp(2)=(2.0d0*pi*nvec_pp(2))/(l*1.0d0)
	 pp(3)=(2.0d0*pi*nvec_pp(3))/(l*1.0d0)

	 p(1)=(2.0d0*pi*nvec_p(1))/(l*1.0d0)				
	 p(2)=(2.0d0*pi*nvec_p(2))/(l*1.0d0)
	 p(3)=(2.0d0*pi*nvec_p(3))/(l*1.0d0)

	 q(1)=(2.0d0*pi*nvec_q(1))/(l*1.0d0)				
	 q(2)=(2.0d0*pi*nvec_q(2))/(l*1.0d0)
	 q(3)=(2.0d0*pi*nvec_q(3))/(l*1.0d0)


	call V_1_mu_scalar((pp-p)/2.0d0,s0,t0,u0,mu,V1mu)
	call S_1_mu_scalar(q,s0,u0,mu,S1mu)

	VS1mu = V1mu + csw*S1mu*delta(s0,t0)*(1 - delta(s0,l) - delta(s0,0))	!Removing from the boundaries


	END SUBROUTINE V1



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE V_2_mu_scalar(s,x0,y0,z0,mu,V2mu)			!4-spinor, 4-spinor
!	-For the XSF.
!	-Computes, independently from the color structure, the fermion-fermion-
!	 gluon vertex for the Wilson fermion action, at a given time indeces,
!	 for a given momenta and for a given vector component.
!
!
!	-Input: -s: momenta
!		-x0,y0,z0: temporal indeces
!		-mu: space-time vectorial component (we only feed in a Loentz
!		     component due to the delta_munu in the definition of the 
!		     vertex
!		
!	-Output: -V2mu: the vertex, 4x4 spinor matrix
	implicit none

	integer,intent(in) :: x0,y0,z0,mu
	real(kind=8),dimension(3),intent(in) :: s
	complex(kind=8),dimension(4,4),intent(out) :: V2mu	!4 spinor 1 vector
	real(kind=8),dimension(3) :: s_plus
		
	!initialize:
	V2mu(:,:) = dcmplx(0.0d0,0.0d0)

	
	!0 component
	if(mu.eq.0)then
	
		V2mu(:,:) = -Pminus*delta(x0+1,y0)*delta(x0,z0)*(1-delta(x0,lat_size)) - Pplus*delta(x0-1,y0)*delta(y0,z0)*(1-delta(x0,0))

	else
		  s_plus(mu) = s(mu) + theta/lat_size
		  V2mu(:,:) = delta(x0,y0)*delta(x0,z0)*(imag*gVec(:,:,mu)*sin(s_plus(mu)) - II*cos(s_plus(mu)))
		  V2mu(:,:) = V2mu(:,:)*(1-delta(x0,0)-delta(x0,lat_size))
	endif

	END SUBROUTINE V_2_mu_scalar

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	SUBROUTINE V2(l,nvec_pp,nvec_p,nvec_q,nvec_qp,s0,t0,u0,u0p,mu,VS2mn)
!	-For the XSF.
!	-Computes, independently from the color structure, the fermion-fermion-
!	 gluon vertex for the full fermion action, at a given time indeces,
!	 for given momenta and for a given vector component.
!	-For the computation it calls V_2_mu_scalar.
!	-Note that for correlation functions of fermion bilinears at 1-loop, 
!	 this vertex only appears in tadpole diagrams. 
!	 Due to the symmetries under gluon exchange, 
! 	 the term with the conmutator vanishes
!	 and there is no contribution from the Clover Term.  
!	 Other diagrams would indeed receive
!	 contribution from the clover term.
!
!	-Input: -l: lattice size
!		-nvec_pp,nvec_p,nvec_q,nvec_qp: momenta ineces
!		-s0,t0,u0: temporal indeces
!		-mu: space-time vectorial component
!		
!	-Output: -VS2mu: the vertex, 4x4 spinor matrix


	implicit none

	integer,intent(in) :: l
	integer,intent(in) :: s0,t0,u0,u0p,mu
	integer,dimension(3),intent(in) :: nvec_p,nvec_pp,nvec_q,nvec_qp
	complex(kind=8),dimension(4,4),intent(out) :: VS2mn	!4-spinor, 4-spinor
	complex(kind=8),dimension(4,4) :: V2mu 			!4-spinor, 4-spinor
	real(kind=8),dimension(3) :: p,pp,q,qp
	integer :: i
	!initialize
	VS2mn(:,:) = dcmplx(0.0d0,0.0d0)
	V2mu(:,:) = dcmplx(0.0d0,0.0d0)
	
	!Construct the momenta
	 pp(1)=(2.0d0*pi*nvec_pp(1))/(l*1.0d0)				
	 pp(2)=(2.0d0*pi*nvec_pp(2))/(l*1.0d0)
	 pp(3)=(2.0d0*pi*nvec_pp(3))/(l*1.0d0)

	 p(1)=(2.0d0*pi*nvec_p(1))/(l*1.0d0)				
	 p(2)=(2.0d0*pi*nvec_p(2))/(l*1.0d0)
	 p(3)=(2.0d0*pi*nvec_p(3))/(l*1.0d0)

	 q(1)=(2.0d0*pi*nvec_q(1))/(l*1.0d0)				
	 q(2)=(2.0d0*pi*nvec_q(2))/(l*1.0d0)
	 q(3)=(2.0d0*pi*nvec_q(3))/(l*1.0d0)

	 qp(1)=(2.0d0*pi*nvec_qp(1))/(l*1.0d0)				
	 qp(2)=(2.0d0*pi*nvec_qp(2))/(l*1.0d0)
	 qp(3)=(2.0d0*pi*nvec_qp(3))/(l*1.0d0)



	call V_2_mu_scalar((pp-p)/2.0d0,s0,t0,u0,mu,V2mu)
	
				!Note that this already contains the delta_mu_nu
	  VS2mn(:,:) = 0.5d0*V2mu(:,:)*delta(u0,u0p)
	
	
	END SUBROUTINE V2	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE Vertices
