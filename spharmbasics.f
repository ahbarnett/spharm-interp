C     Fortran library for sph harm projection and eval on grids & circles
C
C     Stuff from Leslie/Zydrunas FMMLIB3D, and localexp3d.
C     pulled by Barnett 7/31/15
C     modified by Greengard 8/1/15
C     Edited down, Barnett 8/25/22. Depends on yrecursion.f
c
c     Note there is also a "p" version of this lib with packed cnm storage.
C
C***********************************************************************
C     calling sequences changed from Barnett's spharmrouts.f 
C***********************************************************************
c
C***********************************************************************
      subroutine projloc3d(nterms,ldl,nquad,nquadm,xnodes,wts,
     1           phival,local)
C***********************************************************************
C     Usage:
C
C           compute spherical harmonic expansion on unit sphere
C           of function tabulated at nquad*nquad grid points.
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl    = dimension parameter for local expansion
C           nquad  = number of quadrature nodes in theta direction.
C           nquadm = number of quadrature nodes in phi direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           NOTE:    We assume phi_i = (i-1)*2*pi/nquadm, as do the 
C                    routines in projection.f. However, we permit
C                    different numbers of nodes in theta and phi.
C***********************************************************************
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C
C     NOTE:
C
C     yrecursion.f produces Ynm with a nonstandard scaling:
C     (without the 1/sqrt(4*pi)). Thus the orthogonality relation
C     is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C     In the first loop below, you see
C
Cccc        marray(jj,m) = sum*2*pi/nquadm
C           marray(jj,m) = sum/(2*nquadm)
C
C     The latter has incorporated the 1/(4*pi) normalization factor
C     into the azimuthal quadrature weight (2*pi/nquadm).
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm,ldl
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 phival(nquadm,nquad)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      complex *16, allocatable :: marray(:,:)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(marray(nquad,-nterms:nterms))
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
            local(l,m) = 0.0d0
         enddo
      enddo

c
c     create marray (intermediate array)
c
      do m=-nterms,nterms
         emul = cdexp(imag*m*2*pi/nquadm)
         do jj=1,nquad
            sum = 0
            ephi = 1.0d0
            do kk = 1,nquadm
ccc               ephi = cdexp(imag*m*(kk-1)*2*pi/nquadm)
               sum = sum + phival(kk,jj)*dconjg(ephi)
               ephi = ephi*emul
            enddo
ccc         marray(jj,m) = sum*2*pi/nquadm
            marray(jj,m) = sum/(2*nquadm)
         enddo
      enddo
c
c     get local exp
c
      do jj=1,nquad
         cthetaj = xnodes(jj)
         call ylgndr(nterms,cthetaj,ynm)
         do m=-nterms,nterms
            zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local(l,m) = local(l,m) + zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo

      return
      end



C     from localexp3d/h3dwrappers.f =====================================

c**********************************************************************
      subroutine flattenlocexpz(nterms,in,out,stride)
c     unpacks a complex *16 2d local expansion array (in) into a 1d list (out).
c     Output array is written to using given stride.
c     Alex Barnett 3/28/12
      implicit none
      integer i,n,m,nterms,stride
      complex *16 in(0:nterms,-nterms:nterms), out(*)
      
      i = 1
      do n=0,nterms
         do m=-n,n
            out(i) = in(n,m)
            i = i+stride
         enddo
      enddo
      end

c**********************************************************************
      subroutine stacklocexpz(nterms,in,out,stride)
c     packs a complex *16 1d local expansion list (in) into a 2d array (out).
c     Input array is read using given stride.
c     Alex Barnett 2/4/13
      implicit none
      integer i,n,m,nterms,stride
      complex *16 out(0:nterms,-nterms:nterms), in((nterms+1)**2)
      
      i = 1
      do n=0,nterms
         do m=-n,n
            out(n,m) = in(i)
            i = i+stride
         enddo
      enddo
      end
C
C
C
C***********************************************************************
      subroutine shevalsphere(local,phival,
     1           nterms,lmp,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion on an 
C     nquad x nquadm grid on the unit sphere. 
C
C---------------------------------------------------------------------
C     INPUT:
C
C     local    : coefficients of spherical harmonic exp.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension parameter for local
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C     xnodes   : Legendre nodes in theta (x_j = cos theta_j).
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : function value on tensor product
C                mesh on target sphere. phi is the fast variable, theta slow
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquadm,nquad)
      complex *16, allocatable :: phitemp(:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(phitemp(-nterms:nterms,nquad))
c
c$OMP PARALLEL DO
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(m,jj) = 0.0d0
      enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL PRIVATE(ynm)
c     note each thread gets own allocation of ynm output from ylgndr
      allocate(ynm(0:nterms,0:nterms))
C$OMP DO PRIVATE(ctheta,stheta,m,mabs,ix,n)
      do jj=1,nquad
         ctheta = xnodes(jj)
         stheta = dsqrt(1.0d0 - ctheta**2)
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               phitemp(m,jj) = phitemp(m,jj) +
     1                local(n,m)*ynm(n,mabs)
            enddo
         enddo
      enddo
C$OMP END DO
c$OMP END PARALLEL
c
c$OMP PARALLEL DO PRIVATE(kk,ephik,ephi,m)
      do jj = 1,nquad
      do kk = 1,nquadm
         phival(kk,jj) = 0.0d0
         ephik = cdexp(2*pi*(kk-1)*imag/nquadm)
         ephi = ephik**(-nterms)
         do m = -nterms,nterms
            phival(kk,jj) = phival(kk,jj) + phitemp(m,jj)*ephi
            ephi = ephi*ephik
         enddo
      enddo
      enddo
c$OMP END PARALLEL DO
      return
      end
C
C
C
C
c-----------------------------------------------------------------------
c
c      sheval3d: computes potential due to SH expansion 
c                f_j =  A * {M_nm}  f_j is value at point j given
c                       by pts (cos theta_j,phi_j)
c
c      shproj3d: adjoint of sheval3d wth diagonal prefactor
c                 a collection of charges  
c                 N_nm = A^H D f_j
c
c**********************************************************************
      subroutine sheval3d(mpole,nterms,ctheta,phi,pot)
c**********************************************************************
c
c     This subroutine evaluates the spherical harmonic expansion
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi) 
c             n   m
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ctheta,phi  :    target location: ctheta=cos(theta)=z; phi azim.
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    value at target
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ier
      integer lpp,ipp,lephi,lused
      real *8 ctheta,phi
      real *8, allocatable :: w(:)
      complex *16 pot
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: wephi(:)
c
      ier=0
c
c     Carve up workspace for Ynm 
c
      lpp=(nterms+1)**2+5
      ipp=1
      lused = ipp+lpp
      allocate(w(lused))
ccc      iephi = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
ccc      lephi=2*(2*nterms+3)+5 
      lephi=(2*nterms+3)+5 
c
      lused=1+lephi
      allocate(wephi(lused))
c
      call sheval3d0(mpole,nterms,ctheta,phi,pot,w(ipp),wephi)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine sheval3d0(mpole,nterms,ctheta,phi,pot,ynm,ephi)
c**********************************************************************
c
c     See sheval3d for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms
      integer i,n,m
      real *8 ctheta,phi
      real *8 ynm(0:nterms,0:nterms)
      complex *16 pot,ephi1
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1)
c
      real *8 done,cphi,sphi,stheta
      complex *16 eye
      complex *16 ztmp2
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=done
      ephi(1)=ephi1
      cphi = dreal(ephi1)
      sphi = dimag(ephi1)
      ephi(-1)=dconjg(ephi1)
      do i=2,nterms+1
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
c
c     get the associated Legendre functions
c     and scale by 1/sqrt(2l+1)
c
      call ylgndr(nterms,ctheta,ynm)
c
      pot=mpole(0,0)
c
      do n=1,nterms
         pot=pot+mpole(n,0)*ynm(n,0)
         do m=1,n
            ztmp2 = mpole(n,m)*ephi(m) + mpole(n,-m)*ephi(-m)
            pot=pot+ynm(n,m)*ztmp2
         enddo
      enddo
      return
      end
c
c
c
c
