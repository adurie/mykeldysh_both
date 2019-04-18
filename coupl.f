C     Modified so only left-hand lead potentials change


c     *****************************************************************
c     PROGRAM jmbcocu.f
c
c     This program calculates the coupling for an arbitrary lattice
c     and growth direction, using the 'torque formula' at finite
c     temperature.
c     The program is written for nspin tight-binding bands, with
c     upto numnn.th nearest neighbour interactions.
c     The routine param depends on atomic composition, nspin and numnn.
c     The routine eint depends on nspin (and numnn ???).
c     ******************************************************************
c     
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)

      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
      complex*16 ec
      common /fermi/ ef
      common/lbl/ pi, sq2, sq3
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/geom/a1(3),a2(3),a3(3),xn(3),numnn 
      dimension b1(3),b2(3)
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax

      dimension vcuu(0:nlayx),vcdd(0:nlayx),vcud(0:nlayx),vcdu(0:nlayx)
      common/parts1/vruu(0:nlayx),vrdd(0:nlayx),vrud(0:nlayx),
     $  vrdu(0:nlayx)


c     ------------------------------------------------------------
c
c     Note that nmax must be greater than 2 and nmin greater or equal
c     to 2. nmin and nmax are the minimum and maximum number of Cu planes
c     in the spacer.
c     temp  = kbT (effective temperature in Ryd)
c     iwmax = number of Matsubara frequencies used
c     nq    = set up partition in the two dimensional IBZ
c     ef    = Fermi energy (chemical potential)
c

      pi    = 4.0d0*datan(1.0d0)
      sq2   = dsqrt(2.0d0)
      sq3   = dsqrt(3.0d0)
c     ------------------------------------------------------------
c     DATA

      nspin=1
      nsl=1
c     order of highest nearest neighbour interaction
      numnn=1
 
c     fraction of an atomic layer to be adlayered
      frat=1.d0/1.d0
 
      ifrac=int(1.d-10+(1.d0/frat))
      nmin=0*nsl
      nmax=20
      ndiff=(nmax-nmin)*ifrac
      nmat=nsl*nspin

      tfac     = 8.617d-05/13.6058d0
      temp  = 300.d0*tfac
      ef    = 0.d0
      iwmax = 15
      nq    = 1


c     3 primitive lattice vectors (units of alat)
c     FCC
      a1(1)=1.0
      a1(2)=0.0
      a1(3)=0.0
      a2(1)=0.0
      a2(2)=1.0
      a2(3)=0.0
      a3(1)=0.0
      a3(2)=0.0
      a3(3)=1.0
 
c     lattice constant
      alat=1.0d0             !! redundant
      do k=1,3
        a1(k)=a1(k)*alat
        a2(k)=a2(k)*alat
        a3(k)=a3(k)*alat
      enddo

c     perpendicular direction n
      xn(1)=0
      xn(2)=0
      xn(3)=1
 

      if(ndiff.gt.nlayx)then
        write(*,*)' ERROR : ndiff > nlayx ',ndiff,nlayx
        stop
      endif
c
c     -----------------------------------------------------------------
c
      write(*,*)' Co-Cu-Co'
      write(*,*)
      write(*,*)' FCC latice'
      write(*,*)
      write(*,*)' lattice vectors'
      write(*,*)(a1(i),i=1,3)
      write(*,*)(a2(i),i=1,3)
      write(*,*)(a3(i),i=1,3)
      write(*,*)
      write(*,*)' growth direction'
      write(*,*)(xn(i),i=1,3)
      write(*,*)
      write(*,*)' temp =',temp
      write(*,*)' iwmax =',iwmax
      write(*,*)' nq =',nq
      write(*,*)' ef =',ef
c
c     -----------------------------------------------------------------
c
c     determine electronic and geometric structure
      call struct(a1,a2,a3,xn,numnn)
      call param

c     determine the reciprocal lattice structure perp'r to xn
      call recip(xn,b1,b2,irecip)

c     the next line checks that the hamiltonian can be written
c     in semi-block U-T diagonal form with U and T of dimension
c     (nspin*nsl)X(nspin*nsl)
      if(ipmax.gt.nsl)then
        write(*,*)' ERROR - HAMIL : ipmax > nsl !! '
        stop
      endif
c     -----------------------------------------------------------------
      write(*,*)
      write(*,*)' perpr reciprocal lattice vectors'
      write(*,*)(b1(i),i=1,3)
      write(*,*)(b2(i),i=1,3)
      write(*,*)
c     -----------------------------------------------------------------

      do in=0,ndiff
        vcuu(in) = 0.d0
        vcud(in) = 0.d0
        vcdu(in) = 0.d0
        vcdd(in) = 0.d0
      enddo
c
      fact = (pi + pi)*temp
      do iw=1,iwmax
        wm = pi*dfloat(2*iw-1)*temp
        ec = dcmplx(ef,wm)
        call sumk(irecip,ec,nq,b1,b2)
CCC     call kcon(irecip,ec,fact,b1,b2)
        do in=0,ndiff
          vcuu(in) = vcuu(in) - fact*vruu(in)
          vcud(in) = vcud(in) - fact*vrud(in)
          vcdu(in) = vcdu(in) - fact*vrdu(in)
          vcdd(in) = vcdd(in) - fact*vrdd(in)
        enddo         
        write(*,*)' iw =',iw,'   complete'
        call flush(6)
      enddo  
c
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)'          COUPLING '
      do in=0,ndiff
        write(*,*)nmin+in*frat,(vcuu(in)+vcdd(in)-vcud(in)-vcdu(in))
      enddo
      stop
      end
c
c     *****************************************************************
c
      subroutine kcon(irecip,zener,fact,b1,b2)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's

      common/lbl/ pi, sq2, sq3
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/geom/a1(3),a2(3),a3(3),xn(3),numnn
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
      dimension b1(3),b2(3),d1(3),d2(3),xk(3)
      dimension zresu(0:nlayx),zresd(0:nlayx)
      dimension zresud(0:nlayx),zresdu(0:nlayx)
      dimension zconu(0:nlayx),zcond(0:nlayx)
      dimension zconud(0:nlayx),zcondu(0:nlayx)
      dimension xc(0:nlayx),xcold(0:nlayx)
      common/parts1/vruu(0:nlayx),vrdd(0:nlayx),vrud(0:nlayx),
     $  vrdu(0:nlayx)


c     ------------------------------------------------------------
      do in=0,ndiff
        xc(in) = 0.d0
        xcold(in) = 0.d0
      enddo
c     ------------------------------------------------------------
      xcon=5.d-6

      do kk=1,10
        nk=1+(2**(kk+2))

c       calculate the function
        call sumk(irecip,zener,nk,b1,b2)
        do in=0,ndiff
          xc(in)=fact*
     $      (vruu(in)+vrdd(in)-vrud(in)-vrdu(in))
        enddo


c       check for convergence
        diff=0.d0
        iamcon=1
        if(kk.eq.1)iamcon=0
        do in=0,ndiff
          diffc=abs(xc(in)-xcold(in))
          diff=max(diff,diffc)
          if((diff.gt.xcon).or.(kk.le.1))iamcon=0
        enddo
        if(iamcon.eq.0)then
          write(*,'(" NOT converged for nk =",i5," ;  xcon =",f11.8,
     $" ;  error =",f11.8)')nk,xcon,diff
          if(nk.gt.150)then
            write(*,*)' FORGETTING THIS ENERGY POINT !!!!!!'
            write(*,*)' K-POINT CONVERGENCE FAILED !!!!!!',diff
            goto 123
          endif

          do in=0,ndiff
            xcold(in)=xc(in)
          enddo
        endif
        if(iamcon.eq.1)then
          write(*,'(" converged for nk =",i5," ;  xcon =",f11.8,
     $" ;  error =",f11.8)')nk,xcon,diff
          goto 123
        endif
      enddo
123   continue
      return
      end
c
c     *****************************************************************
c
      subroutine sumk(irecip,ec,nk,b1,b2)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's
      common/lbl/ pi, sq2, sq3
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/geom/a1(3),a2(3),a3(3),xn(3),numnn 
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
      dimension b1(3),b2(3),d1(3),d2(3),xk(3)
      common/parts/vquu(0:nlayx),vqdd(0:nlayx),vqud(0:nlayx),
     $  vqdu(0:nlayx)
      common/parts1/vruu(0:nlayx),vrdd(0:nlayx),vrud(0:nlayx),
     $  vrdu(0:nlayx)
      complex*16 ec

c     -----------------------------------------------------------------
      if(nk.ne.1)dq=1.d0/dfloat(nk-1)
      if(nk.eq.1)dq=1.d0
      max=nk-1
      nktot=0
      do in=0,ndiff
        vruu(in) = 0.d0
        vrud(in) = 0.d0
        vrdu(in) = 0.d0
        vrdd(in) = 0.d0
      enddo
c     -----------------------------------------------------------------
      if(irecip.eq.0)then
        write(*,*)' SUMK : irecip = 0 ---- NOT CODED '
        stop
      endif
      if(irecip.eq.1)goto 10
      if(irecip.eq.2)goto 20
      if(irecip.eq.3)goto 30
      if(irecip.eq.4)goto 40
c     -----------------------------------------------------------------
c     CUBIC
10    continue
      do k=1,3
        d1(k)=b1(k)/2.d0
        d2(k)=b2(k)/2.d0
      enddo
      do i=0,max
        do j=0,i
          y=dq*j
C     ?????????????????????????????????????????????????????????????????
C     ?????????????????????????????????????????????????????????????????
        x=2.d0+dq*i
C     ?????????????????????????????????????????????????????????????????
C     ?????????????????????????????????????????????????????????????????
          do k=1,3
            xk(k)=x*d1(k)+y*d2(k)
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
          iwght=8
          if(i.eq.j)iwght=4
          if(j.eq.0)iwght=4
          if(i.eq.max)iwght=4
          if(i.eq.0)iwght=1
          if(j.eq.max)iwght=1
          if((i.eq.max).and.(j.eq.0))iwght=2
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          nktot=nktot+iwght
          call acopl(ec,xk)
c
          do in=0,ndiff
            vruu(in) = vruu(in) + vquu(in)*iwght
            vrud(in) = vrud(in) + vqud(in)*iwght
            vrdu(in) = vrdu(in) + vqdu(in)*iwght
            vrdd(in) = vrdd(in) + vqdd(in)*iwght
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        enddo
      enddo
      goto 5000
c     -----------------------------------------------------------------
c     PRIMITIVE-RECTANGULAR
20    continue
      do k=1,3
        d1(k)=b1(k)/2.d0
        d2(k)=b2(k)/2.d0
      enddo
      do i=0,max
        x=dq*i
        do j=0,max
          y=dq*j
          do k=1,3
            xk(k)=x*d1(k)+y*d2(k)
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
          iwght=4
          if(i.eq.0)iwght=2
          if(j.eq.0)iwght=2
          if(i.eq.max)iwght=2
          if(j.eq.max)iwght=2
          if((i.eq.0).and.(j.eq.0))iwght=1
          if((i.eq.0).and.(j.eq.max))iwght=1
          if((i.eq.max).and.(j.eq.0))iwght=1
          if((i.eq.max).and.(j.eq.max))iwght=1
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          nktot=nktot+iwght
          call acopl(ec,xk)
c
          do in=0,ndiff
            vruu(in) = vruu(in) + vquu(in)*iwght
            vrud(in) = vrud(in) + vqud(in)*iwght
            vrdu(in) = vrdu(in) + vqdu(in)*iwght
            vrdd(in) = vrdd(in) + vqdd(in)*iwght
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        enddo
      enddo
      goto 5000
c     -----------------------------------------------------------------
c     HEXAGONAL
30    continue
      do k=1,3
        d1(k)=b1(k)/2.d0
        d2(k)=(b1(k)+b2(k))/3.d0 - b1(k)/2.d0
      enddo
      do i=0,max
        x=dq*i
        do j=0,i
          y=dq*j
          do k=1,3
            xk(k)=x*d1(k)+y*d2(k)
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
          iwght=12
          if(i.eq.j)iwght=6
          if(j.eq.0)iwght=6
          if(i.eq.max)iwght=6
          if(i.eq.0)iwght=1
          if(j.eq.max)iwght=2
          if((i.eq.max).and.(j.eq.0))iwght=3
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          nktot=nktot+iwght
          call acopl(ec,xk)
c
          do in=0,ndiff
            vruu(in) = vruu(in) + vquu(in)*iwght
            vrud(in) = vrud(in) + vqud(in)*iwght
            vrdu(in) = vrdu(in) + vqdu(in)*iwght
            vrdd(in) = vrdd(in) + vqdd(in)*iwght
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        enddo
      enddo
      goto 5000
c     -----------------------------------------------------------------
c     CENTRED-RECTANGULAR
40    continue
c     ?????????????????????????????????????????????????????????????????
      write(*,*)'CENTRED-RECTANGULAR NOT CODED YET '
      stop
c     ?????????????????????????????????????????????????????????????????
C     b1b1=dot(b1,b1,3)
C     b1b2=dot(b1,b2,3)
C     alpha=b1b1/(2.d0*(b1b1+b1b2))
C     do k=1,3
C       d1(k)=(b1(k)-b2(k))/2.d0
C       s(k)=(0.5d0-alpha)*(b1(k)+b2(k))
C       t(k)=(2.d0*alpha-0.5d0)*(b1(k)+b2(k))
C     enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sum over d1 x s  square
C     do i=0,max
C       x=dq*i
C       do j=0,max
C         y=dq*j
C         do k=1,3
C           xk(k)=x*d1(k)+y*s(k)
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
C         iwght=4
C         if(i.eq.0)iwght=2
C         if(j.eq.0)iwght=2
C         if(i.eq.max)iwght=2
C         if(j.eq.max)iwght=2
C         if((i.eq.0).and.(j.eq.0))iwght=1
C         if((i.eq.0).and.(j.eq.max))iwght=1
C         if((i.eq.max).and.(j.eq.0))iwght=1
C         if((i.eq.max).and.(j.eq.max))iwght=1
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C         nktot=nktot+iwght
C         call acopl(ec,xk)
c
C         do in=0,ndiff
C           vruu(in) = vruu(in) + vquu(in)*iwght
C           vrud(in) = vrud(in) + vqud(in)*iwght
C           vrdu(in) = vrdu(in) + vqdu(in)*iwght
C           vrdd(in) = vrdd(in) + vqdd(in)*iwght
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C       enddo
C     enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sum over d1 x t  triangle
C     do i=0,max
C       x=dq*i
C       do j=0,max
C         y=dq*j
C         do k=1,3
C           xk(k)=x*d1(k)+y*s(k)
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
C         iwght=4
C         if(i.eq.0)iwght=2
C         if(j.eq.0)iwght=2
C         if(i.eq.max)iwght=2
C         if(j.eq.max)iwght=2
C         if((i.eq.0).and.(j.eq.0))iwght=1
C         if((i.eq.0).and.(j.eq.max))iwght=1
C         if((i.eq.max).and.(j.eq.0))iwght=1
C         if((i.eq.max).and.(j.eq.max))iwght=1
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C         nktot=nktot+iwght
C         call acopl(ec,xk)
C
C         do in=0,ndiff
C           vruu(in) = vruu(in) + vquu(in)*iwght
C           vrud(in) = vrud(in) + vqud(in)*iwght
C           vrdu(in) = vrdu(in) + vqdu(in)*iwght
C           vrdd(in) = vrdd(in) + vqdd(in)*iwght
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     goto 5000
c     -----------------------------------------------------------------
c     =================================================================
5000  continue
      do in=0,ndiff
        vruu(in) = vruu(in)/dfloat(nktot)
        vrud(in) = vrud(in)/dfloat(nktot)
        vrdu(in) = vrdu(in)/dfloat(nktot)
        vrdd(in) = vrdd(in)/dfloat(nktot)
      enddo
      return
      end
c
c     *****************************************************************
c
      subroutine acopl(ec,xk)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
      complex*16 ec
      dimension zglu(nmatx,nmatx),zgld(nmatx,nmatx)
      dimension zgru(nmatx,nmatx),zgrd(nmatx,nmatx)
      dimension zflu(nmatx,nmatx),zfld(nmatx,nmatx)
      dimension zt(nmatx,nmatx),ztdag(nmatx,nmatx)
      dimension zruu(nmatx,nmatx),zrdd(nmatx,nmatx)
      dimension zrud(nmatx,nmatx),zrdu(nmatx,nmatx)

      complex*16 garu(nmatx,nmatx,0:(nlayx+1))
      complex*16 gard(nmatx,nmatx,0:(nlayx+1))
      dimension wkspc(nmatx,nmatx),wkspc2(nmatx)
      dimension xk(3),indx(nmatx)
      complex*16 wkspc
      common/lbl/ pi, sq2, sq3
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/geom/a1(3),a2(3),a3(3),xn(3),numnn 
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
      common/parts/vquu(0:nlayx),vqdd(0:nlayx),vqud(0:nlayx),
     $  vqdu(0:nlayx)
c
c
c     This subroutine calculates the integrand
c
c     -----------------------------------------------------------------
c     first calculate Glu(0-->nlay), Gld(0-->nlay), Gru(nsl), Grd(nsl),
c     and the nsl-->(nlay-nsl) hopping zt, ztdag
      call dos(ec,xk,zglu,zgld,garu,gard,zt)
c     -----------------------------------------------------------------
      do ir=1,nmat
        do is=1,nmat
          ztdag(ir,is)=dconjg(zt(is,ir))
        enddo
      enddo
      call multiply(zglu,zt,zflu,nmat,nmatx)
      call multiply(zgld,zt,zfld,nmat,nmatx)
      call multiply(ztdag,zflu,zglu,nmat,nmatx)
      call multiply(ztdag,zfld,zgld,nmat,nmatx)
c     -----------------------------------------------------------------
c
c
      do 10 in=0,ndiff     !!! we want DOS from nmin --> nmax 
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        do ir=1,nmat
          do is=1,nmat
            zgru(ir,is)=garu(ir,is,in)
            zgrd(ir,is)=gard(ir,is,in)
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
C       EVALUATE MY FORMULA
        call multiply(zgru,zglu,zruu,nmat,nmatx)
        call multiply(zgru,zgld,zrud,nmat,nmatx)
        call multiply(zgrd,zglu,zrdu,nmat,nmatx)
        call multiply(zgrd,zgld,zrdd,nmat,nmatx)
        do ir=1,nmat
          zruu(ir,ir)=zruu(ir,ir)-1.d0
          zrud(ir,ir)=zrud(ir,ir)-1.d0
          zrdu(ir,ir)=zrdu(ir,ir)-1.d0
          zrdd(ir,ir)=zrdd(ir,ir)-1.d0
        enddo


c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ifail = 0
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       ?????????????????????????????????????????????????????????????????
c       ?????????????????????????????????????????????????????????????????
        call f03ahf(nmat,zruu,nmatx,detr,deti,id,wkspc,indx)
        xcon=2.d0**id
        detr = detr*xcon
        deti = deti*xcon
        vquu(in) = dreal(cdlog(dcmplx(detr,deti)))/pi

        call f03ahf(nmat,zrud,nmatx,detr,deti,id,wkspc,indx)
        xcon=2.d0**id
        detr = detr*xcon
        deti = deti*xcon
        vqud(in) = dreal(cdlog(dcmplx(detr,deti)))/pi

        call f03ahf(nmat,zrdu,nmatx,detr,deti,id,wkspc,indx)
        xcon=2.d0**id
        detr = detr*xcon
        deti = deti*xcon
        vqdu(in) = dreal(cdlog(dcmplx(detr,deti)))/pi

        call f03ahf(nmat,zrdd,nmatx,detr,deti,id,wkspc,indx)
        xcon=2.d0**id
        detr = detr*xcon
        deti = deti*xcon
        vqdd(in) = dreal(cdlog(dcmplx(detr,deti)))/pi
c       ?????????????????????????????????????????????????????????????????
c       ?????????????????????????????????????????????????????????????????
c       
10    continue     
      return
      end
c
c     *****************************************************************
c
      subroutine dos(zener,xk,zglu2,zgld2,zgru,zgrd,zt1)
c     this program produces the DOS for Co-(Cu-up) trilayer system
c     evaluate zgl(1-->npl*nsl)
c     Co-Cu(down) may be produced by setting ind=2
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c     max No of planes, principal layers & spin components
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
c     zgl and zgr are the left and right GF diagonal elements
      dimension zgl(nmatx,nmatx,0:(nlayx+1))
      dimension zgr(nmatx,nmatx,0:(nlayx+1))
      dimension zu0(nmatx,nmatx),zt0(nmatx,nmatx)
      dimension zu0a(nmatx,nmatx),zt0a(nmatx,nmatx)
      dimension zu1(nmatx,nmatx),zt1(nmatx,nmatx)
      dimension xk(3)
      dimension zgru(nmatx,nmatx,0:(nlayx+1))
      dimension zgrd(nmatx,nmatx,0:(nlayx+1))
      dimension zglu2(nmatx,nmatx)
      dimension zgld2(nmatx,nmatx)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/geom/a1(3),a2(3),a3(3),xn(3),numnn 
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
      common/lbl/pi,sq2,sq3
c     *****************************************************************
c     -----------------------------------------------------------------
c     define spacer hamiltonian elements  -  u1, t1
      ind=3
      call hamil(zu1,zt1,xk,ind)
c     -----------------------------------------------------------------
c     FIND GF for UP BANDS
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     define bulk hamiltonian elements  - u0, t0 for LH lead
      ind=1
      call hamil(zu0,zt0,xk,ind)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     evaluate zgr(0-->npl*nsl) and zgl(nsl*i) for i=0,1,2,...
c     corresponding to the LH surface GF for all atomic layers
c     and the RH surface GF for all principal layers
      call green(zu0,zt0,zu1,zt1,zener,xk,ind,zglu2,zgru)
c
c     -----------------------------------------------------------------
c     FIND GF for DOWN BANDS
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     define bulk hamiltonian elements  - u0, t0 for LH lead
      ind=2
      call hamil(zu0,zt0,xk,ind)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     evaluate zgl(0-->npl*nsl) and zgr(nsl*i) for i=0,1,2,...
c     corresponding to the LH surface GF for all atomic layers
c     and the RH surface GF for all principal layers
      call green(zu0,zt0,zu1,zt1,zener,xk,ind,zgld2,zgrd)
c     -----------------------------------------------------------------
      return
      end
c
c     *****************************************************************
c
      subroutine green(zu0,zt0,zu1,zt1,zener,xk,ind,zsurfl,
     $zgr)
c     this routine evaluates zgl(1) and
c     zgr(nsl*i) for i=0,1,2,...
c     corresponding to the LH surface GF for all atomic layers
c     and the RH surface GF for all principal layers
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (n2x=2*nmatx)
      dimension zgr(nmatx,nmatx,0:(nlayx+1))
      dimension zu0(nmatx,nmatx),zt0(nmatx,nmatx)
      dimension zu0a(nmatx,nmatx),zt0a(nmatx,nmatx)
      dimension zu1(nmatx,nmatx),zt1(nmatx,nmatx)
      dimension zu01(nmatx,nmatx),zt01(nmatx,nmatx),zt01dag(nmatx,nmatx)
      dimension zsurfl(nmatx,nmatx),zsurfr(nmatx,nmatx)
      dimension zfoo(nmatx,nmatx)
      dimension xk(3)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/lbl/pi,sq2,sq3
      common/saveus/zevc(n2x,n2x),zevcinv(n2x,n2x),zevl(n2x),iamdiag
c     *****************************************************************
c     -----------------------------------------------------------------
c     evaluate left and right surface GF's
      call surfacenew(zu0,zt0,zener,zsurfl,zsurfr)
c     -----------------------------------------------------------------

      do ir=1,nmat
        do is=1,nmat
          zgr(ir,is,0)=zsurfr(ir,is)
        enddo
      enddo

c     =================================================================
c     load RH GF for 1 atomic sub-adlayer
c
C  NB : zgr has different dimensions to zgpl
      iamdiag=0
      xmin=nmin-2*0
      call adlaysave(zgr,zu1,zt1,zener,1,1,xmin)   !! -- gr(nmin)
      do ir=1,nmat
        do is=1,nmat
          zgr(ir,is,0)=zgr(ir,is,1)
        enddo
      enddo
      iamdiag=1
      call adlaysave(zgr,zu1,zt1,zener,1,ndiff,frat)   !! -- gr(nmin-->nlay)
c     -----------------------------------------------------------------

      return
      end
c
c     *****************************************************************
c
      subroutine surface(zu0,zt0,zgl,zgr,zener,err)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c     max No of planes, principal layers & spin components
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
c     zgl,zgr and zgn are the left, right and total GF diagonal elements
      dimension zgl(nmatx,nmatx,0:(nplx+1)),zgr(nmatx,nmatx,0:(nplx+1))
      dimension zu0(nmatx,nmatx),zt0(nmatx,nmatx)
      dimension zgl1(nmatx,nmatx),zgr1(nmatx,nmatx)
      dimension zgl2(nmatx,nmatx),zgr2(nmatx,nmatx)
      dimension zaib(nmatx,nmatx),zbia(nmatx,nmatx)
      dimension zinv(nmatx,nmatx),zinva(nmatx,nmatx),zinvb(nmatx,nmatx)
      dimension zunit(nmatx,nmatx)
      dimension zalpha1(nmatx,nmatx),zalpha2(nmatx,nmatx)
      dimension zbeta1(nmatx,nmatx),zbeta2(nmatx,nmatx)
      dimension zeps1(nmatx,nmatx),zeps2(nmatx,nmatx)
      dimension zepsl1(nmatx,nmatx),zepsl2(nmatx,nmatx)
      dimension zepsr1(nmatx,nmatx),zepsr2(nmatx,nmatx)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff

      do ir=1,nmat
        do is=1,nmat
          zunit(ir,is)=0.d0
        enddo
        zunit(ir,ir)=1.d0
      enddo

c     -----------------------------------------------------------------
      maxit=100
c     -----------------------------------------------------------------
c     zeroth iteration
      do ir=1,nmat
        do is=1,nmat
          zeps1(ir,is)=zu0(ir,is)
          zepsl1(ir,is)=zu0(ir,is)
          zepsr1(ir,is)=zu0(ir,is)
          zalpha1(ir,is)=zt0(ir,is)
          zbeta1(is,ir)=dconjg(zt0(ir,is))
        enddo
      enddo
c     determine zgl1 := (E - H00)**-1
      do ir=1,nmat
        do is=1,nmat
          zgl1(ir,is)=zener*zunit(ir,is)-zu0(ir,is)
        enddo
      enddo
c     find inverse to zgl1
      call invers(zgl1,nmat,nmatx)
      do ir=1,nmat
        do is=1,nmat
          zgr1(ir,is)=zgl1(ir,is)
        enddo
      enddo
c     -----------------------------------------------------------------
c     calculate gl(0) and gr(0) iteratively
      do its=1,maxit
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       first calculate zinva = [(E - zeps1)**-1]*zalpha1
c                and    zinvb = [(E - zeps1)**-1]*zbeta1
        do ir=1,nmat
          do is=1,nmat
            zinv(ir,is)=zener*zunit(ir,is)-zeps1(ir,is)
          enddo
        enddo
        call invers(zinv,nmat,nmatx)
        call multiply(zinv,zalpha1,zinva,nmat,nmatx)
        call multiply(zinv,zbeta1,zinvb,nmat,nmatx)
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       now determine zalpha2, zbeta2, zeps2, zepsl2
        call multiply(zalpha1,zinva,zalpha2,nmat,nmatx)
        call multiply(zbeta1,zinvb,zbeta2,nmat,nmatx)
        call multiply(zalpha1,zinvb,zaib,nmat,nmatx)
        call multiply(zbeta1,zinva,zbia,nmat,nmatx)
        do ir=1,nmat
          do is=1,nmat
            zeps2(ir,is)=zeps1(ir,is)+zaib(ir,is)+zbia(ir,is)
            zepsl2(ir,is)=zepsl1(ir,is)+zaib(ir,is)
            zepsr2(ir,is)=zepsr1(ir,is)+zbia(ir,is)
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       now determine G(0,0):= zgl2 and G(n,n) := zgr2
        do ir=1,nmat
          do is=1,nmat
            zgl2(ir,is)=zener*zunit(ir,is)-zepsl2(ir,is)
            zgr2(ir,is)=zener*zunit(ir,is)-zepsr2(ir,is)
          enddo
        enddo
        call invers(zgl2,nmat,nmatx)
        call invers(zgr2,nmat,nmatx)
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       now test for convergence
        do ir=1,nmat
          do is=1,nmat
            difl=cdabs(zgl1(ir,is)-zgl2(ir,is))
            difr=cdabs(zgr1(ir,is)-zgr2(ir,is))
            if(difl.gt.err)goto 50
            if(difr.gt.err)goto 50
          enddo
        enddo
        goto 100   !! convergence achieved
50      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       now load up alpha1, beta1, eps1, epsl1, zgl1
        do ir=1,nmat
CCC       write(*,*)(zgl1(ir,is),is=1,nmat)
          do is=1,nmat
            zeps1(ir,is)=zeps2(ir,is)
            zepsl1(ir,is)=zepsl2(ir,is)
            zepsr1(ir,is)=zepsr2(ir,is)
            zalpha1(ir,is)=zalpha2(ir,is)
            zbeta1(ir,is)=zbeta2(ir,is)
            zgl1(is,ir)=zgl2(is,ir)
            zgr1(is,ir)=zgr2(is,ir)
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      enddo
c     -----------------------------------------------------------------
      write(*,*)' SURFACE : failed to converge for maxit=',maxit
      stop

100   continue
c     -----------------------------------------------------------------
c     now load zgl and zgr (!! NB: zgr1-->zgl and zgl1-->zgr)
      do ir=1,nmat
        do is=1,nmat
          zgr(ir,is,0)=zgl1(ir,is)
          zgl(ir,is,0)=zgr1(ir,is)
        enddo
      enddo
c     -----------------------------------------------------------------
      return
      end
c
c     *****************************************************************
c
      subroutine adlayer1(zgl,zu1,zt1,zener)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c     max No of planes, principal layers & spin components
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
c     zgl,zgr and zgn are the left, right and total GF diagonal elements
      dimension zgl(nmatx,nmatx)
      dimension zu1(nmatx,nmatx),zt1(nmatx,nmatx)
      dimension zgt(nmatx,nmatx),ztgt(nmatx,nmatx)
      dimension zt1dag(nmatx,nmatx)
      dimension zunit(nmatx,nmatx)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff

      do ir=1,nmat
        do is=1,nmat
          zunit(ir,is)=0.d0
          zt1dag(is,ir)=dconjg(zt1(ir,is))
        enddo
        zunit(ir,ir)=1.d0
      enddo
c     -----------------------------------------------------------------
c     adlayer ontop of gl
      call multiply(zgl,zt1,zgt,nmat,nmatx)
      call multiply(zt1dag,zgt,ztgt,nmat,nmatx)
      do ir=1,nmat
        do is=1,nmat
          zgl(ir,is)=zener*zunit(ir,is)-zu1(ir,is)-ztgt(ir,is)
        enddo
      enddo
c     find inverse to zginv
      call invers(zgl,nmat,nmatx)
c     -----------------------------------------------------------------

      return
      end
c
c     *****************************************************************
c     *****************************************************************
C     HAMILTONIAN ROUTINES -- GENERAL STRUCTURE
c     *****************************************************************
c     *****************************************************************
      subroutine hamil(zu,zt,xk,ind)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)

      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
      common/lbl/ pi, sq2, sq3
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/geom/a1(3),a2(3),a3(3),xn(3),numnn
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
      dimension xk(3)
      dimension zu(nmatx,nmatx),zt(nmatx,nmatx)
      dimension zh(nspinx,nspinx)
      common/par1/ s0(5),p0(5),d0e(5),d0t(5)
      common/par2/ sss(5,2),sps(5,2),pps(5,2),ppp(5,2),sds(5,2)
      common/par3/ pds(5,2),pdp(5,2),dds(5,2),ddp(5,2),ddd(5,2)

c     -----------------------------------------------------------------
c     compute the hamiltonian zu and zt
      do i=1,nsl
        do j=1,nsl
          iu=j-i
          call helement(iu,xk,zh,ind)
          do isp=1,nspin
            do jsp=1,nspin
              iii=(i-1)*nspin+isp
              jjj=(j-1)*nspin+jsp
              zu(iii,jjj)=zh(isp,jsp)
            enddo
          enddo
        enddo
      enddo
      do i=1,nsl
        do j=1,i
          it=nsl+j-i
          call helement(it,xk,zh,ind)
          do isp=1,nspin
            do jsp=1,nspin
              iii=(i-1)*nspin+isp
              jjj=(j-1)*nspin+jsp
              zt(iii,jjj)=zh(isp,jsp)
            enddo
          enddo
        enddo
      enddo
c     -----------------------------------------------------------------

      return
      end
c     *****************************************************************
      subroutine helement(ipl,xk,zh,ind)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
      common/lbl/ pi, sq2, sq3
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
      dimension d(3),c(3),xk(3)
      dimension zh(nspinx,nspinx),zt(nspinx,nspinx),rt(nspinx,nspinx)
      common/par1/ s0(5),p0(5),d0e(5),d0t(5)
      common/par2/ sss(5,2),sps(5,2),pps(5,2),ppp(5,2),sds(5,2)
      common/par3/ pds(5,2),pdp(5,2),dds(5,2),ddp(5,2),ddd(5,2)

      do ir=1,nspin
        do is=1,nspin
          zh(ir,is)=0.d0
          zt(ir,is)=0.d0
          rt(ir,is)=0.d0
        enddo
      enddo
c     for each point in the ipl.th plane
      do i=1,np(ipl)
        nn=nnp(ipl,i)   !!! this point is an nn.th N.N to 0
        do k=1,3
          d(k)=rp(ipl,i,k)   !!! it has position vector d
          c(k)=cosp(ipl,i,k) !!! it has direction cosines
        enddo
c     ?????????????????????????????????????????????????????????????????
c     NB : units of pi NOT 2pi -- since lattice const = 2
CCCCC   dk=pi*dot(xk,d,3)
        dk=dot(xk,d,3)
c     ?????????????????????????????????????????????????????????????????
CCC     zdk=dcmplx(0.d0,dk)
CCC     zexdk=cdexp(zdk)
        zexdk=dcmplx(cos(dk),sin(dk))
c       find the hopping for this nn.th N.N in direction d
c     ?????????????????????????????????????????????????????????????????
c     ?????????????????????????????????????????????????????????????????
c       ensure that this routine gives U for d = 0
        call eint11(nn,ind,zt)
c     ?????????????????????????????????????????????????????????????????
c     ?????????????????????????????????????????????????????????????????
C        dd=dot(d,d,3)
C        if(dd.eq.0)then
C          rt(1,1)=s0(ind)
C          rt(2,2)=p0(ind)
C          rt(3,3)=p0(ind)
C          rt(4,4)=p0(ind)
C          do ir=5,7
C            rt(ir,ir)=d0t(ind)
C          enddo
C          do ir=8,9
C            rt(ir,ir)=d0e(ind)
C          enddo
C        else
C          g1=sss(ind,nn)
C          g2=sps(ind,nn)
C          g3=pps(ind,nn)
C          g4=ppp(ind,nn)
C          g5=sds(ind,nn)
C          g6=pds(ind,nn)
C          g7=pdp(ind,nn)
C          g8=dds(ind,nn)
C          g9=ddp(ind,nn)
C          g10=ddd(ind,nn)
C          call eint1(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,c,rt)
C        endif
C        do ir=1,nspin
C          do is=1,nspin
C            zt(ir,is)=rt(ir,is)
C          enddo
C        enddo
c     ?????????????????????????????????????????????????????????????????
        do ir=1,nspin
          do is=1,nspin
            zh(ir,is)=zh(ir,is)+zt(ir,is)*zexdk
          enddo
        enddo
      enddo
      return
      end
c     *****************************************************************
      subroutine struct(a1,a2,a3,xn,numnn)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nnx=2,npmax=100)   !!! upto nnx th N.N
c                                    npmax - max No of points per plane
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (lattx=(2*nnx+1)*(2*nnx+1)*(2*nnx+1))
c     construct a lattice from the primitive lattice vectors {a1,a2,a3}
c     and order according to distance from origin
c     then sort into those perp'r and not-perp'r to a direction xn
      dimension r(lattx,3),a1(3),a2(3),a3(3),v(3),d(lattx),iorder(lattx)
      dimension xn(3),rtmp(lattx,3),dtmp(lattx)
      dimension rxn(lattx),rxntmp(lattx),nnn(lattx),nnntmp(lattx)
      dimension ixn(lattx),x1(3),x2(3),x3(3)
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
c     -----------------------------------------------------------------
c     INPUT : {a1,a2,a3} - primitive lattice vectors
c             xn         - growth direction
c             numnn      - upto order numnn.th N.N interactions
c
c     OUTPUT : ipmax     - the fathest plane (i.e. there are numnn.th
c                          N.N's on the ipmax.th and -ipmax.th plane
c                          in the +xn direction.
c              np        - np(i) ; i=-ipmax,...,ipmax is the number of
c                          points on the ith plane which are either 
c                          1st, 2nd, ... ,numnn.th N.N's.
c              rp        - rp(i,j,k) ; 
c                          i=-ipmax,...,ipmax ; j=1,...,np(i) ; k=1,2,3
c                          the location of the jth point on the ith plane
c              cosp      - cosp(i,j,k) ;
c                          i=-ipmax,...,ipmax ; j=1,...,np(i) ; k=1,2,3
c                          the (x,y,z) direction cosines of the jth point
c                          on the ith plane
c              nnp       - nnp(i,j) ; i=-ipmax,...,ipmax ; j=1,...,np(i)
c                          tells you if the jth point on the ith plane
c                          is a 1st, 2nd, ... N.N.
c     -----------------------------------------------------------------
c     normalise xn
      xnorm=dot(xn,xn,3)
      do i=1,3
        xn(i)=xn(i)/sqrt(xnorm)
      enddo
c     -----------------------------------------------------------------
c     now construct lattice
      icnt=0
      do i=-numnn,numnn,1
        do j=-numnn,numnn,1
          do k=-numnn,numnn,1
            icnt=icnt+1
            iorder(icnt)=icnt
            do ll=1,3
              v(ll)=i*a1(ll)+j*a2(ll)+k*a3(ll)
              r(icnt,ll)=v(ll)
            enddo
              d(icnt)=sqrt(v(1)**2 + v(2)**2 + v(3)**2)
          enddo
        enddo
      enddo
      ntot=icnt
      if(ntot.gt.lattx)then
        write(*,*)' ERROR - STRUCT : ntot > lattx !! '
        stop
      endif
c     -----------------------------------------------------------------
c     sort into distance from origin
      call sort(d,iorder,ntot)
      do i=1,ntot
        do j=1,3
          rtmp(i,j)=r(iorder(i),j)
        enddo
        dtmp(i)=d(iorder(i))
      enddo
c     -----------------------------------------------------------------
c     count up number of points with same distance from origin
      d0=dtmp(1)
      id=1       !!! id th most distant set from origin
      nnntmp(1)=id-1
      icnt=1
      do i=2,ntot
        if(abs(d0-dtmp(i)).gt.1.d-8)then
          d0=dtmp(i)
          id=id+1
        endif
        nnntmp(i)=id-1
        if(nnntmp(i).gt.numnn)goto 100
        icnt=i
      enddo
100   continue
      ntot=icnt   !!! total number of points upto numnn.th N.N
      if(ntot.gt.((2*nslx+1)*npmax))then
        write(*,*)' ERROR - STRUCT : ntot > (2*nslx+1)*npmax !! '
        stop
      endif
c     -----------------------------------------------------------------
c     now order into distance along the xn direction
      do i=1,ntot
        iorder(i)=i
        do k=1,3
          v(k)=rtmp(i,k)
        enddo
        rxntmp(i)=dot(v,xn,3)
      enddo
      call sort(rxntmp,iorder,ntot)
      do i=1,ntot
        rxn(i)=rxntmp(iorder(i))
        nnn(i)=nnntmp(iorder(i))
        d(i)=dtmp(iorder(i))
        do k=1,3
          r(i,k)=rtmp(iorder(i),k)
        enddo
      enddo
c     -----------------------------------------------------------------
c     normalise rxn
      rmin=abs(rxn(1))
      do i=2,ntot
        ri=rxn(i)
        if((abs(ri).lt.rmin).and.(abs(ri).gt.1.d-8))rmin=abs(ri)
      enddo
      do i=1,ntot
        rxn(i)=rxn(i)/rmin
        if(rxn(i).gt.0)ixn(i)=int(rxn(i)+1.d-8)
        if(rxn(i).lt.0)ixn(i)=int(rxn(i)-1.d-8)
      enddo
c     -----------------------------------------------------------------
c     find out how many points on the ith plane -- np(i)
      ipmax=ixn(ntot)   !!! is the fathest plane
      if(ipmax.gt.nslx)then
        write(*,*)' ERROR - STRUCT : ipmax > nslx !! '
        stop
      endif
      ip0=ixn(1)
      np(ip0)=1
      do i=2,ntot
        if(abs(ip0-ixn(i)).lt.1.d-8)then  ! r(i) on same plane as r(i-1)
          np(ip0)=np(ip0)+1
          if(np(ip0).gt.npmax)then
            write(*,*)' ERROR - STRUCT : np(ip0) > npmax !! '
            stop
          endif
        else                              ! r(i) on next plane as r(i-1)
          ip0=ixn(i)
          np(ip0)=1
        endif
      enddo
c     load points in ith plane into rp(i)
      icnt=0
      do i=-ipmax,ipmax,1
        do j=1,np(i)
          icnt=icnt+1
          nnp(i,j)=nnn(icnt)
          do k=1,3
            rp(i,j,k)=r(icnt,k)
          enddo
        enddo
      enddo
c     find the direction cosines cosp(i)
      x1(1)=1.d0
      x1(2)=0.d0
      x1(3)=0.d0
      x2(1)=0.d0
      x2(2)=1.d0
      x2(3)=0.d0
      x3(1)=0.d0
      x3(2)=0.d0
      x3(3)=1.d0
      do i=-ipmax,ipmax,1
        do j=1,np(i)
          do k=1,3
            v(k)=rp(i,j,k)
          enddo
          rpnorm=sqrt(v(1)**2 + v(2)**2 + v(3)**2)
          if(rpnorm.eq.0)then
            cosp(i,j,1)=0.d0
            cosp(i,j,2)=0.d0
            cosp(i,j,3)=0.d0
          else
            cosp(i,j,1)=dot(x1,v,3)/rpnorm
            cosp(i,j,2)=dot(x2,v,3)/rpnorm
            cosp(i,j,3)=dot(x3,v,3)/rpnorm
          endif
        enddo
      enddo
c     -----------------------------------------------------------------
C     icnt=0
C     do i=-ipmax,ipmax,1
C       do j=1,np(i)
C         icnt=icnt+1
C         write(*,'(i5,3x,3f10.6,5x,f12.8,5x,i3,5x,i3)')icnt,
C    $      (rp(i,j,k),k=1,3),d(icnt),nnn(icnt),ixn(icnt)
C       enddo
C     enddo
c     -----------------------------------------------------------------
      return
      end
c
c     *****************************************************************
c
      subroutine recip(xn,b1,b2,irecip)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=2,npmax=100) !! upto nnx.th N.N's
      common/lbl/ pi, sq2, sq3
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      common/structure/rp(-nslx:nslx,npmax,3),cosp(-nslx:nslx,npmax,3),
     $  np(-nslx:nslx),nnp(-nslx:nslx,npmax),ipmax
      dimension an1(3),an2(3),d(3),v(3),b1(3),b2(3),xn(3)

c     now construct reciprocal lattice perp'r to xn
c     check nnp(0) ordered  (it should be !!)
      do i=2,np(0)
        if(nnp(0,i-1).gt.nnp(0,i))then
          write(*,*)' ERROR : NPP not ordered !!!!'
          stop
        endif
      enddo
c     now determine the 2 (lin'y independent) nearest points in this
c     plane (with smallest angle between them)
      an1(1)=rp(0,2,1)
      an1(2)=rp(0,2,2)
      an1(3)=rp(0,2,3)
      do i=3,np(0)
        do k=1,3
          d(k)=rp(0,i,k)
        enddo
        call cross(an1,d,v,3)
        vv=dot(v,v,3)
        if(vv.gt.1.d-8)then
          an2(1)=d(1)
          an2(2)=d(2)
          an2(3)=d(3)
          goto 300
        endif
      enddo
300   continue
      an1an1=dot(an1,an1,3)
      an2an2=dot(an2,an2,3)
      if((an1an1.lt.1.d-8).or.(an2an2.lt.1.d-8))then
        write(*,*)' ERROR : RECIP : an1 or an2 = 0 !!',an1an1,an2an2
        write(*,*)' Try increasing numnn1 !! '
      endif
c     an1 and an2 are the primitive lattice vectors of this plane
c     -----------------------------------------------------------------
c     construct the third (normalised) lattice vector
c     N.B by construction d.an1 = d.an2 = 0
      xnorm=dot(xn,xn,3)
      do i=1,3
        d(i)=xn(i)/sqrt(xnorm)
      enddo
c     -----------------------------------------------------------------
c     now determine the reciprocal lattice vectors
      call cross(an2,d,b1,3)
      call cross(an1,d,b2,3)
      an1b1=dot(an1,b1,3)
      an2b2=dot(an2,b2,3)
      do k=1,3
        b1(k)=2*pi*b1(k)/an1b1
        b2(k)=2*pi*b2(k)/an2b2
      enddo
c     now determine whether b2 or -b2 has the smallest angle with b1
      b1b1=dot(b1,b1,3)
      b2b2=dot(b2,b2,3)
      b1b2=dot(b1,b2,3)
      if(b1b2.lt.0)then
        do k=1,3
          b2(k)=-b2(k)
        enddo
      endif
      b1b2=-b1b2
      cos12=acos(b1b2/sqrt(b1b1*b2b2))
c     -----------------------------------------------------------------
c     now classify the reciprocal lattice
      irecip=0
      if((abs(b1b2).lt.1.d-8).and.(abs(b1b1-b2b2).lt.1.d-8))then
        write(*,*)' reciprocal lattice is cubic'
        irecip=1
        goto 400
      elseif((abs(b1b2).lt.1.d-8).and.(abs(b1b1-b2b2).ge.1.d-8))then
        write(*,*)' reciprocal lattice is primitive-rectangular'
        irecip=2
        goto 400
      elseif((abs(cos12-pi/3.d0).lt.1.d-8).and.
     $  (abs(b1b1-b2b2).lt.1.d-8))then
        write(*,*)' reciprocal lattice is hexagonal'
        irecip=3
        goto 400
      elseif((abs(b1b2).gt.1.d-8).and.(abs(b1b1-b2b2).lt.1.d-8))then
        write(*,*)' reciprocal lattice is centred-rectangular'
        irecip=4
        goto 400
      endif
400   continue
      write(*,*)
      write(*,*)
      return
      end
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
C     MATHEMATICAL ROUTINES
c     *****************************************************************
c     *****************************************************************
c
      double precision function dot(x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(n),y(n)
      xy=0.d0
      do i=1,n
        xy=xy+x(i)*y(i)
      enddo
      dot=xy
      return
      end
c
c     *****************************************************************
c
      subroutine cross(x,y,xy,n)
      implicit double precision (a-h,o-z)
      dimension x(n),y(n),xy(n)
      xy(1)=x(2)*y(3)-x(3)*y(2)
      xy(2)=x(3)*y(1)-x(1)*y(3)
      xy(3)=x(1)*y(2)-x(2)*y(1)
      return
      end
c
c     *****************************************************************
c
      subroutine multiply(zx,zy,zxy,nmat,nmatx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zx(nmatx,nmatx),zy(nmatx,nmatx),zxy(nmatx,nmatx)
      do ir=1,nmat
        do is=1,nmat
          zxy(ir,is)=0.d0
          do k=1,nmat
            zxy(ir,is)=zxy(ir,is)+zx(ir,k)*zy(k,is)
          enddo
        enddo
      enddo
      return
      end
c
c     *****************************************************************
c
      subroutine invers(array,n,nx) 
      complex*16 array(nx,nx)
      complex*16 amax, save
      parameter (iworkx=100)
      integer ik(iworkx),jk(iworkx) 
      if(iworkx.lt.n)then
        write(*,*)' ERROR INVERS : iworkx < n'
        stop
      endif
11    do 100 k=1,n 
        amax=0.d0 
21      do 30 i=k,n 
          do 30 j=k,n 
23          if(cdabs(amax)-cdabs(array(i,j)))24,24,30 
24          amax=array(i,j) 
            ik(k)=i 
            jk(k)=j 
30      continue 
41      i=ik(k) 
        if(i-k)21,51,43 
43      do 50 j=1,n 
          save=array(k,j) 
          array(k,j)=array(i,j) 
50        array(i,j)=-save 
51      j=jk(k) 
        if(j-k)21,61,53 
53      do 60 i=1,n 
          save=array(i,k) 
          array(i,k)=array(i,j) 
60        array(i,j)=-save 
61      do 70 i=1,n 
          if(i-k)63,70,63 
63        array(i,k)=-array(i,k)/amax 
70      continue 
71      do 80 i=1,n 
          do 80 j=1,n 
            if(i-k)74,80,74 
74          if(j-k)75,80,75 
75          array(i,j)=array(i,j)+array(i,k)*array(k,j) 
80      continue 
81      do 90 j=1,n 
          if(j-k)83,90,83 
83        array(k,j)=array(k,j)/amax 
90      continue 
        array(k,k)=1.d0/amax 
100   continue 
101   do 130 l=1,n 
        k=n-l+1 
        j=ik(k) 
        if(j-k)111,111,105 
105     do 110 i=1,n 
          save=array(i,k) 
          array(i,k)=-array(i,j) 
110       array(i,j)=save 
111     i=jk(k) 
        if(i-k)130,130,1013
1013    do 120 j=1,n 
          save=array(k,j) 
          array(k,j)=-array(i,j) 
120       array(i,j)=save 
130   continue 
140   return 
      end
c
c     *****************************************************************
c
      subroutine cdet(nn,a,ns,detr,deti,au,indx)
      implicit double precision (a-h,o-z)
      integer nn,ns,i,j
      integer indx(ns)
      double precision d,detr,deti
      complex*16 det
      complex*16 a(ns,ns),au(ns,ns)
c
c     calcula o determinante de uma matriz complexa
c
c     input:  a    - matriz de dimensao (ns,ns)
c             ns   - dimensao da matriz
c     output: indx - array de trabalho de dimensao nsize
c             detr - parte real do determinante
c             deti - parte imaginaria do determinante
c
      do 10 i=1,ns
      do 10 j=1,ns
         au(i,j) = a(i,j)
 10   continue
      call ludcmp(au,nn,ns,indx,d)
      det = dcmplx(1.d0,0.d0)
      do 11 j=1,nn
         det = det*dcmplx(d,0.d0)*au(j,j)
 11   continue
      detr = dreal(det)
      deti = dimag(det)
      return
      end
c
c     *****************************************************************
c
      subroutine ludcmp(a,n,np,indx,d)
      implicit double precision (a-h,o-z)
      integer  nmatx,n,np,i,j,k,imax
      double precision   tiny,d,aamax,dum
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      integer    indx(np)
      double precision     vv(nmatx)
      complex*16 a(np,np)
      complex*16 czero,sum,zdum
c
      tiny=1.0d-20
      d = 1.d0
      czero = (0.d0,0.d0)
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (cdabs(a(i,j)).gt.aamax) aamax=cdabs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'Singular matrix.'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*cdabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            zdum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=zdum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.czero)a(j,j)=dcmplx(tiny,tiny)
          zdum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*zdum
18        continue
        endif
19    continue
      if(a(n,n).eq.czero)a(n,n)=dcmplx(tiny,tiny)
      return
      end   
c
c     *****************************************************************
c
      subroutine sort(x,ka,l)
      implicit double precision (a-h,o-z)
c
c     subroutine sorts the array x
c     the indices (ascending order) are given in array ka
c
      parameter (nvecx=10000)
      dimension x(l),ka(l),kb(nvecx)
      if(l.gt.nvecx)then
        write(*,*)' ERROR SORT.F : l > nvecx ',l,nvecx
        stop
      endif
      do 13 j=1,l
        kb(j)=j
13      ka(j)=j
      l1=1
200   j=1
20    limi1=l1*2*(j-1)+1
      if(limi1.gt.l) go to 28
      limi2=limi1+l1
      lims1=min(limi2-1,l)
      if(limi2.gt.l) go to 28
      lims2=min(limi2+l1-1,l)
      ir=limi1
      i1=ir
      i2=limi2
21    k1=ka(i1)
      k2=ka(i2)
      if(x(k1).le.x(k2)) go to 22
      kb(ir)=k2
      i2=i2+1
      ir=ir+1
      if(i2.le.lims2) go to 21
      go to 23
22    kb(ir)=k1
      i1 = i1+1
      ir=ir+1
      if(i1.le.lims1) go to 21
24    k2=ka(i2)
      kb(ir) = k2
      i2 = i2+1
      ir = ir+1
      if(i2.le.lims2) go to 24
      go to 25
23    k1=ka(i1)
      kb(ir) = k1
      i1=i1+1
      ir=ir+1
      if(i1.le.lims1) go to 23
25    j=j+1
      go to 20
28    do 280 k=1,l
280     ka(k)=kb(k)
      l1 = l1*2
      if(l1.lt.2*l) go to 200
      return
      end
c
c     *****************************************************************
c     *****************************************************************
c     ATOM DEPENDENT ROUTINES
c     *****************************************************************
c     *****************************************************************
c
      subroutine param
c     THIS ROUTINE IS ATOM DEPENDENT :-
c     it is written for Cu and Co-up and Co-down, with s,p,d bands.
      implicit double precision (a-h,o-z)
      common/par1/ s0(5),p0(5),d0e(5),d0t(5)
      common/par2/ sss(5,2),sps(5,2),pps(5,2),ppp(5,2),sds(5,2)
      common/par3/ pds(5,2),pdp(5,2),dds(5,2),ddp(5,2),ddd(5,2)
c
c     two centre integrals; slater-koster parameters for spin polarized
c     bcc Iron was obtained starting from the paramagnetic parametrization
c     of R H Victora "Magnetic and Electronic Properties of Transition
c     Metals and Overlayers". Self consistency is performed shifting
c     the centre of the d bands only (Udd=1.eV Usp=0)
c     bcc Cr bulk spin bands:
c     Also from R H Victora's paper (Udd=Usp=0 eV)
c
c
c     -----------------------------------------------------------------
c     The first index in the tight binding parameter arrays refers to
c     1: Co up
c     2: Co down
c     3: Cu
c
c     -----------------------------------------------------------------
c     Co up:
c
      cshift = .575530d0 - .715751d0
      delta = 1.113608939931278d-001
      vex = delta*0.5d0
      s0(1) =  1.12946d0 + cshift
      p0(1) =  1.75262d0 + cshift
      d0t(1) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift - vex
      d0e(1) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift - vex
c
c     first n.n.
c
      sss(1,1) = -0.09043d0
      sps(1,1) =  0.13649d0
      pps(1,1) =  0.23748d0
      ppp(1,1) = -0.00142d0
      sds(1,1) = -0.03806d0
      pds(1,1) = -0.04069d0
      pdp(1,1) =  0.02797d0
      dds(1,1) = -0.04213d0
      ddp(1,1) =  0.02976d0
      ddd(1,1) = -0.00684d0
c
c     second n.n.
c
      sss(1,2) = -0.00337d0
      sps(1,2) =  0.00135d0
      pps(1,2) =  0.02849d0
      ppp(1,2) =  0.01099d0
      sds(1,2) = -0.01119d0
      pds(1,2) = -0.01061d0
      pdp(1,2) =  0.01134d0
      dds(1,2) = -0.00759d0
      ddp(1,2) =  0.00495d0
      ddd(1,2) = -0.00016d0
c
c     -----------------------------------------------------------------
c     Co down:
c
      cshift = .575530d0 - .715751d0
      delta = 1.113608939931278d-001
      vex = delta*0.5d0
      s0(2) =  1.12946d0 + cshift
      p0(2) =  1.75262d0 + cshift
      d0t(2) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift + vex
      d0e(2) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift + vex
c
c     first n.n.
c
      sss(2,1) = -0.09043d0
      sps(2,1) =  0.13649d0
      pps(2,1) =  0.23748d0
      ppp(2,1) = -0.00142d0
      sds(2,1) = -0.03806d0
      pds(2,1) = -0.04069d0
      pdp(2,1) =  0.02797d0
      dds(2,1) = -0.04213d0
      ddp(2,1) =  0.02976d0
      ddd(2,1) = -0.00684d0
c
c     second n.n.
c
      sss(2,2) = -0.00337d0
      sps(2,2) =  0.00135d0
      pps(2,2) =  0.02849d0
      ppp(2,2) =  0.01099d0
      sds(2,2) = -0.01119d0
      pds(2,2) = -0.01061d0
      pdp(2,2) =  0.01134d0
      dds(2,2) = -0.00759d0
      ddp(2,2) =  0.00495d0
      ddd(2,2) = -0.00016d0
c
c     -----------------------------------------------------------------
c     Cu:
c
      s0(3) =  0.79466d0
      p0(3) =  1.35351d0
      d0t(3) =  0.5d0*(0.37307d0 + 0.37180d0)
      d0e(3) =  0.5d0*(0.37307d0 + 0.37180d0)
c
c     first n.n.
c
      sss(3,1) = -0.07518d0
      sps(3,1) =  0.11571d0
      pps(3,1) =  0.19669d0
      ppp(3,1) =  0.01940d0
      sds(3,1) = -0.03107d0
      pds(3,1) = -0.03289d0
      pdp(3,1) =  0.01753d0
      dds(3,1) = -0.02566d0
      ddp(3,1) =  0.01800d0
      ddd(3,1) = -0.00408d0
c
c     second n.n.
c
      sss(3,2) = -0.00092d0
      sps(3,2) =  0.01221d0
      pps(3,2) =  0.05389d0
      ppp(3,2) =  0.00846d0
      sds(3,2) = -0.00852d0
      pds(3,2) = -0.00536d0
      pdp(3,2) =  0.00321d0
      dds(3,2) = -0.00451d0
      ddp(3,2) =  0.00241d0
      ddd(3,2) = -0.00029d0
c
c     -----------------------------------------------------------------
c     Co(UP)/Cu Interface
c
      s0(4) =  0.79466d0
      p0(4) =  1.35351d0
      d0t(4) =  0.5d0*(0.37307d0 + 0.37180d0)
      d0e(4) =  0.5d0*(0.37307d0 + 0.37180d0)
c
c     first n.n.
c
      sss(4,1) = -0.082453
      sps(4,1) =  0.125671
      pps(4,1) =  0.216125
      ppp(4,1) =  0.005249
      sds(4,1) = -0.034388
      pds(4,1) = -0.036583
      pdp(4,1) =  0.022143
      dds(4,1) = -0.032879
      ddp(4,1) =  0.023145
      ddd(4,1) = -0.005283
c
c     second n.n.
c
      sss(4,2) = -0.001761
      sps(4,2) =  0.004060
      pps(4,2) =  0.039183
      ppp(4,2) =  0.009642
      sds(4,2) = -0.009764
      pds(4,2) = -0.007541
      pdp(4,2) =  0.006033
      dds(4,2) = -0.005851
      ddp(4,2) =  0.003454
      ddd(4,2) = -0.000215
c
c     -----------------------------------------------------------------
c     Co(UP)/Cu Interface
c
      s0(5) =  0.79466d0
      p0(5) =  1.35351d0
      d0t(5) =  0.5d0*(0.37307d0 + 0.37180d0)
      d0e(5) =  0.5d0*(0.37307d0 + 0.37180d0)
c
c     first n.n.
c
      sss(5,1) = -0.082453
      sps(5,1) =  0.125671
      pps(5,1) =  0.216125
      ppp(5,1) =  0.005249
      sds(5,1) = -0.034388
      pds(5,1) = -0.036583
      pdp(5,1) =  0.022143
      dds(5,1) = -0.032879
      ddp(5,1) =  0.023145
      ddd(5,1) = -0.005283
c
c     second n.n.
c
      sss(5,2) = -0.001761
      sps(5,2) =  0.004060
      pps(5,2) =  0.039183
      ppp(5,2) =  0.009642
      sds(5,2) = -0.009764
      pds(5,2) = -0.007541
      pdp(5,2) =  0.006033
      dds(5,2) = -0.005851
      ddp(5,2) =  0.003454
      ddd(5,2) = -0.000215
c
c     
      return
      end
c
c     *****************************************************************
c     *****************************************************************
c
      subroutine eint1(sss,sps,pps,ppp,ss,ps,pp,ds,dp,dd,w,b)
C     THIS ROUTINE IS SPIN DEPENDENT :-
C     IT IS WRITTEN FOR s,p and d BANDS ONLY
      implicit double precision (a-h,o-z)      
      dimension b(9,9)
      dimension w(3)
      x=w(1)
      y=w(2)
      z=w(3)
      xx=x*x
      xy=x*y
      yy=y*y
      yz=y*z
      zz=z*z
      zx=z*x
      xxyy=xx*yy
      yyzz=yy*zz
      zzxx=zz*xx
      aux=pps-ppp
      r3=dsqrt(3.d0)
      aux1=r3*ss
      f8=3.d0*zz-1.d0
      f1=xx+yy
      f2=xx-yy
      f3=zz-.5d0*f1
      g1=1.5d0*f2*ds
      g2=r3*f3*ds
      b(1,1)=sss
      b(1,2)=x*sps
      b(1,3)=y*sps
      b(1,4)=z*sps
      b(2,1)=-b(1,2)
      b(2,2)=xx*pps+(1.d0-xx)*ppp
      b(2,3)=xy*aux
      b(2,4)=zx*aux
      b(3,1)=-b(1,3)
      b(3,2)=b(2,3)
      b(3,3)=yy*pps+(1.d0-yy)*ppp
      b(3,4)=yz*aux
      b(4,1)=-b(1,4)
      b(4,2)=b(2,4)
      b(4,3)=b(3,4)
      b(4,4)=zz*pps+(1.d0-zz)*ppp
      b(1,5)=xy*aux1
      b(1,6)=yz*aux1
      b(1,7)=zx*aux1
      b(1,8)=.5d0*f2*aux1
      b(1,9)=.5d0*f8*ss
      b(5,1)=b(1,5)
      b(6,1)=b(1,6)
      b(7,1)=b(1,7)
      b(8,1)=b(1,8)
      b(9,1)=b(1,9)
      f4=.5d0*r3*f2*ps
      f5=.5d0*f8*ps
      aux2=r3*xx*ps+(1.d0-2.d0*xx)*pp
      b(2,5)=aux2*y
      b(2,6)=(r3*ps-2.d0*pp)*xy*z
      b(2,7)=aux2*z
      b(2,8)=(f4+(1.d0-f2)*pp)*x
      b(2,9)=(f5-r3*zz*pp)*x
      aux3=(r3*yy*ps+(1.d0-2.d0*yy)*pp)
      b(3,5)=aux3*x
      b(3,6)=aux3*z
      b(3,7)=b(2,6)
      b(3,8)=(f4-(1.d0+f2)*pp)*y
      b(3,9)=(f5-r3*zz*pp)*y
      aux4=r3*zz*ps+(1.d0-2.d0*zz)*pp
      b(4,5)=b(2,6)
      b(4,6)=aux4*y
      b(4,7)=aux4*x
      b(4,8)=(f4-f2*pp)*z
      b(4,9)=(f5+r3*f1*pp)*z
      b(5,2)=-b(2,5)
      b(6,2)=-b(2,6)
      b(7,2)=-b(2,7)
      b(8,2)=-b(2,8)
      b(9,2)=-b(2,9)
      b(5,3)=-b(3,5)
      b(6,3)=-b(3,6)
      b(7,3)=-b(3,7)
      b(8,3)=-b(3,8)
      b(9,3)=-b(3,9)
      b(5,4)=-b(4,5)
      b(6,4)=-b(4,6)
      b(7,4)=-b(4,7)
      b(8,4)=-b(4,8)
      b(9,4)=-b(4,9)
      b(5,5)=3.d0*xxyy*ds+(f1-4.d0*xxyy)*dp+(zz+xxyy)*dd
      b(5,6)=(3.d0*yy*ds+(1.d0-4.d0*yy)*dp+(yy-1.d0)*dd)*zx
      b(5,7)=(3.d0*xx*ds+(1.d0-4.d0*xx)*dp+(xx-1.d0)*dd)*yz
      b(5,8)=(g1-2.d0*f2*dp+.5d0*f2*dd)*xy
      b(5,9)=(g2-2.d0*r3*zz*dp+.5d0*r3*(1.d0+zz)*dd)*xy
      b(6,5)=b(5,6)
      b(6,6)=3.d0*yyzz*ds+(yy+zz-4.d0*yyzz)*dp+(xx+yyzz)*dd
      b(6,7)=(3.d0*zz*ds+(1.d0-4.d0*zz)*dp+(zz-1.d0)*dd)*xy
      b(6,8)=(g1-(1.d0+2.d0*f2)*dp+(1.d0+.5d0*f2)*dd)*yz
      b(6,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*yz
      b(7,5)=b(5,7)
      b(7,6)=b(6,7)
      b(7,7)=3.d0*zzxx*ds+(zz+xx-4.d0*zzxx)*dp+(yy+zzxx)*dd
      b(7,8)=(g1+(1.d0-2.d0*f2)*dp-(1.d0-.5d0*f2)*dd)*zx
      b(7,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*zx
      b(8,5)=b(5,8)
      b(8,6)=b(6,8)
      b(8,7)=b(7,8)
      b(8,8)=.75d0*f2*f2*ds+(f1-f2*f2)*dp+(zz+.25d0*f2*f2)*dd
      b(8,9)=.5d0*f2*g2-r3*zz*f2*dp+.25d0*r3*(1.d0+zz)*f2*dd
      b(9,5)=b(5,9)
      b(9,6)=b(6,9)
      b(9,7)=b(7,9)
      b(9,8)=b(8,9)
      b(9,9)=f3*f3*ds+3.d0*zz*f1*dp+.75d0*f1*f1*dd
      return
      end
c
c     *****************************************************************
      subroutine eint11(nn,ind,zt)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=250,nslx=1,nspinx=1)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
      dimension zt(nspinx,nspinx)
      dimension d(3),c(3)
c     ?????????????????????????????????????????????????????????????????
c     ensure that this routine gives U for d = 0
c     ?????????????????????????????????????????????????????????????????
c     ?????????????????????????????????????????????????????????????????
      do ir=1,nspin
        do is=1,nspin
          zt(ir,is)=0.0d0
        enddo
      enddo
      if(nn.eq.0)then
        if(ind.eq.1)zui=-2.9d0
        if(ind.eq.2)zui=-2.55d0
        if(ind.eq.3)zui=-2.80d0
        if(ind.eq.4)zui=7.00d0
      elseif(nn.eq.1)then
        zui=0.5d0
      elseif(nn.eq.2)then
        zui=0.25d0
      endif
      do ir=1,nspin
        zt(ir,ir)=zui
      enddo
      return
      end
c
c     *****************************************************************
c
      subroutine surfacenew(zu,zt,zener,zsurfl,zsurfr)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (n2x=2*nmatx)
      dimension zgamma(nmatx,nmatx)
      dimension zunit(nmatx,nmatx)
      dimension ztmp1(nmatx,nmatx),ztmp2(nmatx,nmatx),ztmp3(nmatx,nmatx)
      dimension zu(nmatx,nmatx),zt(nmatx,nmatx)
      dimension ztinv(nmatx,nmatx),zsinv(nmatx,nmatx)
      dimension zs(nmatx,nmatx)
      dimension zp(n2x,n2x)
      dimension zsurfl(nmatx,nmatx),zsurfr(nmatx,nmatx)
      dimension isort(n2x)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      dimension wk1(n2x),wk2(n2x),wk3(n2x)
      dimension rr(n2x),ri(n2x),evlab(n2x),zevl(n2x)
      dimension ar(n2x,n2x),ai(n2x,n2x)
      dimension vr(n2x,n2x),vi(n2x,n2x)
c     -----------------------------------------------------------------

      nmat2=2*nmat
      do i=1,nmat
        do j=1,nmat
          zunit(i,j)=0.d0
        enddo
        zunit(i,i)=1.d0
      enddo

c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     now calculate GF's from closed form
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     define zp
c     redefine gamma and delta
      do i=1,nmat
        do j=1,nmat
          ztinv(i,j)=zt(i,j)
          zs(i,j)=dconjg(zt(j,i))
          zsinv(i,j)=zs(i,j)
          ztmp1(i,j)=zener*zunit(i,j)-zu(i,j)
        enddo
      enddo
      call invers(ztinv,nmat,nmatx)
      call invers(zsinv,nmat,nmatx)
      call multiply(ztmp1,ztinv,zgamma,nmat,nmatx)
      do ir=1,nmat
        do is=1,nmat
          zp(ir,is)=0.d0
          zp(ir,is+nmat)=ztinv(ir,is)
          zp(ir+nmat,is)=-zs(ir,is)
          zp(ir+nmat,is+nmat)=zgamma(ir,is)
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     diagonalise zp
      do ir=1,nmat2
        do is=1,nmat2
          ar(ir,is)=dreal(zp(ir,is))
          ai(ir,is)=dimag(zp(ir,is))
        enddo
      enddo 
      idoevec=1
      ifail=0
      call cg(n2x,nmat2,ar,ai,rr,ri,idoevec,vr,vi,wk1,wk2,wk3,ifail)
      if(ifail.ne.0)then
        write(*,*)'ADLAYER : ifail =',ifail
        stop 
      endif
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sort the |evals| into increasing order
c     and check eval(nmat+1)>eval(nmat)
c
      do ir=1,nmat2
        evlab(ir)=sqrt((rr(ir)**2)+(ri(ir)**2))
        isort(ir)=ir
      enddo
      call sort(evlab,isort,nmat2)
      evln1=evlab(isort(nmat+1))
      evln=evlab(isort(nmat))
      if(abs(evln1-evln).lt.1.d-10)then
        write(*,*)' ERROR SURFACE : degenerate eigenstates '
        stop
      endif
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     load the highest nmat eigenvectors
      do ir=1,nmat
        do is=1,nmat
          ik=isort(is+nmat)
          ztmp1(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
          ztmp2(ir,is)=dcmplx(vr(ir+nmat,ik),vi(ir+nmat,ik))
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate L.H. zsurf
      call invers(ztmp2,nmat,nmatx)
      call multiply(ztmp1,ztmp2,zsurfl,nmat,nmatx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate R.H. zsurf
      do ir=1,nmat2
        ik=isort(ir)
        zevl(ir)=dcmplx(rr(ik),ri(ik))
      enddo
      do ir=1,nmat
        do is=1,nmat
          ik=isort(is)
          ztmp1(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
          ztmp2(ir,is)=ztmp1(ir,is)
        enddo
      enddo
      call invers(ztmp2,nmat,nmatx)
      do ir=1,nmat
        do is=1,nmat
          ztmp1(ir,is)=ztmp1(ir,is)*zevl(is)
        enddo
      enddo
      call multiply(ztmp1,ztmp2,ztmp3,nmat,nmatx)
      call multiply(ztmp3,zsinv,zsurfr,nmat,nmatx)
c
c     -----------------------------------------------------------------
      return
      end
c
c     *****************************************************************
c
      subroutine adlaynew(zfn,zubig,ztbig,zener,n1,n2,frac)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (n2x=2*nmatx)
      dimension ztmp1(nmatx,nmatx)
      dimension zubig(nmatx,nmatx),ztbig(nmatx,nmatx)
      dimension zf0(nmatx,nmatx)
      dimension zfn(nmatx,nmatx,0:(nlayx+1))
      dimension zg0(nmatx,nmatx),zgn(nmatx,nmatx)
      dimension zp(n2x,n2x),zevc(n2x,n2x),zevcinv(n2x,n2x)
      dimension isort(n2x)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      dimension wk1(n2x),wk2(n2x),wk3(n2x)
      dimension rr(n2x),ri(n2x),evlab(n2x),evlphi(n2x)
      dimension zroot(n2x)
      common/evals/zevl(n2x)
      dimension ar(n2x,n2x),ai(n2x,n2x)
      dimension vr(n2x,n2x),vi(n2x,n2x)
c     -----------------------------------------------------------------
      dimension zt1(nspinx,nspinx,nslx),zu(nspinx,nspinx)
      dimension zt1dag(nspinx,nspinx,nslx)
      dimension ztndin(nspinx,nspinx)
      dimension zunit(nspinx,nspinx)
      dimension zdum1(nspinx,nspinx),zdum2(nspinx,nspinx)
      dimension ztmpa(nmatx,nmatx),ztmpb(nmatx,nmatx)
      dimension ztmpc(nmatx,nmatx),ztmpd(nmatx,nmatx)
c     -----------------------------------------------------------------
      pi=acos(-1.d0)

      nmat2=2*nmat
      do ir=1,nmat
        do is=1,nmat
          zf0(ir,is)=zfn(ir,is,0)
        enddo
      enddo
      do i=1,nspin
        do j=1,nspin
          zunit(i,j)=0.d0
          zu(i,j)=zubig(i,j)
        enddo
        zunit(i,i)=1.d0
      enddo
      do il=0,nsl-1,1
        do i=1,nspin
          do j=1,nspin
            zt1(i,j,nsl-il)=ztbig(i+il*nspin,j)
            zt1dag(j,i,nsl-il)=dconjg(ztbig(i+il*nspin,j))
          enddo
        enddo
      enddo
      do i=1,nspin
        do j=1,nspin
          ztndin(i,j)=zt1dag(i,j,nsl)
        enddo
      enddo
      call invers(ztndin,nspin,nspinx)

c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     DEFINE ZP
c     -----------------------------------------------------------------
c     FIRST DEFINE SUB-ELEMENTS A,B,C,D

c     top RH element
      do i=1,nspin
        do j=1,nspin
          zdum1(i,j)=zener*zunit(i,j)-zu(i,j)
        enddo
      enddo
      call multiply(zdum1,ztndin,zdum2,nspin,nspinx)
      do i=1,nspin
        do j=1,nspin
          ztmpa(i,nmat-nspin+j)=0.d0
          ztmpb(i,nmat-nspin+j)=ztndin(i,j)
          ztmpc(i,nmat-nspin+j)=-zt1(i,j,nsl)
          ztmpd(i,nmat-nspin+j)=zdum2(i,j)
        enddo
      enddo

      if(nsl.ge.2)then
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       top row
        do jl=0,nsl-2
          do j=1,nspin
            jj=jl*nspin+j
            do i=1,nspin
              ztmpa(i,jj)=0.d0
              ztmpb(i,jj)=0.d0
              ztmpc(i,jj)=-zt1(i,j,jl+1)
              ztmpd(i,jj)=0.d0
            enddo
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       RH column
        do il=1,nsl-1
          do i=1,nspin
            do j=1,nspin
              zdum2(i,j)=zt1dag(i,j,il)
            enddo
          enddo
          call multiply(zdum2,ztndin,zdum1,nspin,nspinx)

          do i=1,nspin
            ii=il*nspin+i
            do j=1,nspin
              jj=nmat-nspin+j
              ztmpa(ii,jj)=0.d0
              ztmpb(ii,jj)=0.d0
              ztmpc(ii,jj)=0.d0
              ztmpd(ii,jj)=-zdum1(i,j)
            enddo
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       bottom LH block
        do i=1,nmat-nspin
          ii=nspin+i
          do j=1,nmat-nspin
            ztmpa(ii,j)=0.d0
            ztmpb(ii,j)=0.d0
            ztmpc(ii,j)=0.d0
            ztmpd(ii,j)=0.d0
          enddo
          ztmpa(ii,i)=1.d0
          ztmpd(ii,i)=1.d0
        enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      endif

c     NOW LOAD ZP
      do ir=1,nmat
        do is=1,nmat
          zp(ir,is)=ztmpa(ir,is)
          zp(ir,is+nmat)=ztmpb(ir,is)
          zp(ir+nmat,is)=ztmpc(ir,is)
          zp(ir+nmat,is+nmat)=ztmpd(ir,is)
        enddo
      enddo
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------

c     diagonalise zp
      do ir=1,nmat2
        do is=1,nmat2
          ar(ir,is)=dreal(zp(ir,is))
          ai(ir,is)=dimag(zp(ir,is))
        enddo
      enddo 
      idoevec=1
      ifail=0
      call cg(n2x,nmat2,ar,ai,rr,ri,idoevec,vr,vi,wk1,wk2,wk3,ifail)
      if(ifail.ne.0)then
        write(*,*)'ADLAYER : ifail =',ifail
        stop 
      endif
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sort the |evals| into increasing order
c
      do ir=1,nmat2
        evlab(ir)=sqrt((rr(ir)**2)+(ri(ir)**2))
        isort(ir)=ir
      enddo
      call sort(evlab,isort,nmat2)
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     load the eigenvectors and calculate zg0
      do ir=1,nmat2
        do is=1,nmat2
          ik=isort(is)
          zevc(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
          zevcinv(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
        enddo
        zevl(ir)=dcmplx(rr(ir),ri(ir))
      enddo
      do ir=1,nmat2
        ik=isort(ir)
        zevl(ir)=dcmplx(rr(ik),ri(ik))
        evlab(ir)=cdabs(zevl(ir))
      enddo
      call invers(zevcinv,nmat2,n2x)

      call mobius(zevcinv,zf0,zg0,nmat)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     determine a  root of zevc for a fractional atomic sublayer
CCC   write(*,*)' EIGENVALUES '
c     define -pi < phi < pi
      do ir=1,nmat2
        xr=dreal(zevl(ir))
        xi=dimag(zevl(ir))
        if((xr.ge.0).and.(xi.ge.0))evlphi(ir)=atan(xi/xr)
        if((xr.ge.0).and.(xi.le.0))evlphi(ir)=atan(xi/xr)
        if((xr.le.0).and.(xi.ge.0))evlphi(ir)=atan(xi/xr)+pi
        if((xr.le.0).and.(xi.le.0))evlphi(ir)=atan(xi/xr)-pi
CCC     write(*,*)ir,evlab(ir),evlphi(ir)
        arg=evlphi(ir)
        ztemp=dcmplx(cos(arg),sin(arg))*evlab(ir)
        if(cdabs(ztemp-zevl(ir)).gt.1.d-10)then
          write(*,*)'ERROR in finding phase',ztemp,zevl(ir)
          stop
        endif
        xrt=abs(evlab(ir)**frac)
        zroot(ir)=dcmplx(cos(arg*frac),sin(arg*frac))*xrt
      enddo
CCC   write(*,*)
CCC   write(*,*)

CCC   do ir=1,nmat2
CCC     zroot(ir)=zevl(ir)**frac
CCC   enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate zfn
      do its=n1,n2,1

c       calculate zgn
        xpow=dfloat(its)
        do ir=1,nmat
          do is=1,nmat
            zgn(ir,is)=((zroot(ir)/zroot(is+nmat))**xpow)*zg0(ir,is)
          enddo
        enddo
        call mobius(zevc,zgn,ztmp1,nmat)
        do ir=1,nmat
          do is=1,nmat
            zfn(ir,is,its)=ztmp1(ir,is)
          enddo
        enddo

      enddo
c
c     -----------------------------------------------------------------
      return
      end
c
c     *****************************************************************
c
      subroutine adlaysave(zfn,zubig,ztbig,zener,n1,n2,frac)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (nlayx=nplx*nslx)
      parameter (n2x=2*nmatx)
      dimension ztmp1(nmatx,nmatx)
      dimension zubig(nmatx,nmatx),ztbig(nmatx,nmatx)
      dimension zf0(nmatx,nmatx)
      dimension zfn(nmatx,nmatx,0:(nlayx+1))
      dimension zg0(nmatx,nmatx),zgn(nmatx,nmatx)
      dimension zp(n2x,n2x)
      dimension isort(n2x)
      common/data/frat,ifrac,nmat,nspin,nsl
      common/layer/nmin,ndiff
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      dimension wk1(n2x),wk2(n2x),wk3(n2x)
      dimension rr(n2x),ri(n2x),evlab(n2x),evlphi(n2x)
      dimension zroot(n2x)
      dimension ar(n2x,n2x),ai(n2x,n2x)
      dimension vr(n2x,n2x),vi(n2x,n2x)
      common/saveus/zevc(n2x,n2x),zevcinv(n2x,n2x),zevl(n2x),iamdiag
c     -----------------------------------------------------------------
      dimension zt1(nspinx,nspinx,nslx),zu(nspinx,nspinx)
      dimension zt1dag(nspinx,nspinx,nslx)
      dimension ztndin(nspinx,nspinx)
      dimension zunit(nspinx,nspinx)
      dimension zdum1(nspinx,nspinx),zdum2(nspinx,nspinx)
      dimension ztmpa(nmatx,nmatx),ztmpb(nmatx,nmatx)
      dimension ztmpc(nmatx,nmatx),ztmpd(nmatx,nmatx)
c     -----------------------------------------------------------------
      pi=acos(-1.d0)

      nmat2=2*nmat
      do ir=1,nmat
        do is=1,nmat
          zf0(ir,is)=zfn(ir,is,0)
        enddo
      enddo
      do i=1,nspin
        do j=1,nspin
          zunit(i,j)=0.d0
          zu(i,j)=zubig(i,j)
        enddo
        zunit(i,i)=1.d0
      enddo
      do il=0,nsl-1,1
        do i=1,nspin
          do j=1,nspin
            zt1(i,j,nsl-il)=ztbig(i+il*nspin,j)
            zt1dag(j,i,nsl-il)=dconjg(ztbig(i+il*nspin,j))
          enddo
        enddo
      enddo
      do i=1,nspin
        do j=1,nspin
          ztndin(i,j)=zt1dag(i,j,nsl)
        enddo
      enddo
      call invers(ztndin,nspin,nspinx)

c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     DEFINE ZP --- IF IAMDIAG = 0 ONLY

      if(iamdiag.eq.1)goto 1234

c     -----------------------------------------------------------------
c     FIRST DEFINE SUB-ELEMENTS A,B,C,D

c     top RH element
      do i=1,nspin
        do j=1,nspin
          zdum1(i,j)=zener*zunit(i,j)-zu(i,j)
        enddo
      enddo
      call multiply(zdum1,ztndin,zdum2,nspin,nspinx)
      do i=1,nspin
        do j=1,nspin
          ztmpa(i,nmat-nspin+j)=0.d0
          ztmpb(i,nmat-nspin+j)=ztndin(i,j)
          ztmpc(i,nmat-nspin+j)=-zt1(i,j,nsl)
          ztmpd(i,nmat-nspin+j)=zdum2(i,j)
        enddo
      enddo

      if(nsl.ge.2)then
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       top row
        do jl=0,nsl-2
          do j=1,nspin
            jj=jl*nspin+j
            do i=1,nspin
              ztmpa(i,jj)=0.d0
              ztmpb(i,jj)=0.d0
              ztmpc(i,jj)=-zt1(i,j,jl+1)
              ztmpd(i,jj)=0.d0
            enddo
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       RH column
        do il=1,nsl-1
          do i=1,nspin
            do j=1,nspin
              zdum2(i,j)=zt1dag(i,j,il)
            enddo
          enddo
          call multiply(zdum2,ztndin,zdum1,nspin,nspinx)

          do i=1,nspin
            ii=il*nspin+i
            do j=1,nspin
              jj=nmat-nspin+j
              ztmpa(ii,jj)=0.d0
              ztmpb(ii,jj)=0.d0
              ztmpc(ii,jj)=0.d0
              ztmpd(ii,jj)=-zdum1(i,j)
            enddo
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       bottom LH block
        do i=1,nmat-nspin
          ii=nspin+i
          do j=1,nmat-nspin
            ztmpa(ii,j)=0.d0
            ztmpb(ii,j)=0.d0
            ztmpc(ii,j)=0.d0
            ztmpd(ii,j)=0.d0
          enddo
          ztmpa(ii,i)=1.d0
          ztmpd(ii,i)=1.d0
        enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      endif

c     NOW LOAD ZP
      do ir=1,nmat
        do is=1,nmat
          zp(ir,is)=ztmpa(ir,is)
          zp(ir,is+nmat)=ztmpb(ir,is)
          zp(ir+nmat,is)=ztmpc(ir,is)
          zp(ir+nmat,is+nmat)=ztmpd(ir,is)
        enddo
      enddo
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------

c     diagonalise zp
      do ir=1,nmat2
        do is=1,nmat2
          ar(ir,is)=dreal(zp(ir,is))
          ai(ir,is)=dimag(zp(ir,is))
        enddo
      enddo 
      idoevec=1
      ifail=0
      call cg(n2x,nmat2,ar,ai,rr,ri,idoevec,vr,vi,wk1,wk2,wk3,ifail)
      if(ifail.ne.0)then
        write(*,*)'ADLAYER : ifail =',ifail
        stop 
      endif
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sort the |evals| into increasing order
c
      do ir=1,nmat2
        evlab(ir)=sqrt((rr(ir)**2)+(ri(ir)**2))
        isort(ir)=ir
      enddo
      call sort(evlab,isort,nmat2)
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     load the eigenvectors and calculate zg0
      do ir=1,nmat2
        do is=1,nmat2
          ik=isort(is)
          zevc(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
          zevcinv(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
        enddo
        zevl(ir)=dcmplx(rr(ir),ri(ir))
      enddo
      do ir=1,nmat2
        ik=isort(ir)
        zevl(ir)=dcmplx(rr(ik),ri(ik))
      enddo
      call invers(zevcinv,nmat2,n2x)

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1234  continue
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mobius(zevcinv,zf0,zg0,nmat)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     determine a  root of zevc for a fractional atomic sublayer
CCC   write(*,*)' EIGENVALUES '
c     define -pi < phi < pi
      do ir=1,nmat2
        xr=dreal(zevl(ir))
        xi=dimag(zevl(ir))
        evlab(ir)=cdabs(zevl(ir))
        if((xr.ge.0).and.(xi.ge.0))evlphi(ir)=atan(xi/xr)
        if((xr.ge.0).and.(xi.le.0))evlphi(ir)=atan(xi/xr)
        if((xr.le.0).and.(xi.ge.0))evlphi(ir)=atan(xi/xr)+pi
        if((xr.le.0).and.(xi.le.0))evlphi(ir)=atan(xi/xr)-pi
CCC     write(*,*)ir,evlab(ir),evlphi(ir)
        arg=evlphi(ir)
        ztemp=dcmplx(cos(arg),sin(arg))*evlab(ir)
        if(cdabs(ztemp-zevl(ir)).gt.1.d-10)then
          write(*,*)'ERROR in finding phase',ztemp,zevl(ir)
          stop
        endif
        xrt=abs(evlab(ir)**frac)
        zroot(ir)=dcmplx(cos(arg*frac),sin(arg*frac))*xrt
      enddo
CCC   write(*,*)
CCC   write(*,*)

CCC   do ir=1,nmat2
CCC     zroot(ir)=zevl(ir)**frac
CCC   enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate zfn
      do its=n1,n2,1

c       calculate zgn
        xpow=dfloat(its)
        do ir=1,nmat
          do is=1,nmat
            zgn(ir,is)=((zroot(ir)/zroot(is+nmat))**xpow)*zg0(ir,is)
          enddo
        enddo
        call mobius(zevc,zgn,ztmp1,nmat)
        do ir=1,nmat
          do is=1,nmat
            zfn(ir,is,its)=ztmp1(ir,is)
          enddo
        enddo

      enddo
c
c     -----------------------------------------------------------------
      return
      end
c
c     *****************************************************************
c
      subroutine mobius(zmat,zf,zres,nmat)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c
c     calculates the n dimensional (right-hand) mobius
c     transformation (zres) of zmat on zf.
c
      parameter (nplx=250,nslx=1,nspinx=1)
      parameter (nmatx=nslx*nspinx)
      parameter (n2x=2*nmatx)
c
      dimension zmat(n2x,n2x),zf(nmatx,nmatx),zres(nmatx,nmatx)
      dimension zm11(nmatx,nmatx),zm12(nmatx,nmatx)
      dimension zm21(nmatx,nmatx),zm22(nmatx,nmatx)
      dimension ztmp1(nmatx,nmatx),ztmp2(nmatx,nmatx)
c
c     -----------------------------------------------------------------
c     load the sub-matrices zmij
      do ir=1,nmat
        do is=1,nmat
          zm11(ir,is)=zmat(ir,is)
          zm12(ir,is)=zmat(ir,is+nmat)
          zm21(ir,is)=zmat(ir+nmat,is)
          zm22(ir,is)=zmat(ir+nmat,is+nmat)
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate the numerator
      call multiply(zm11,zf,ztmp1,nmat,nmatx)
      do ir=1,nmat
        do is=1,nmat
          ztmp1(ir,is)=ztmp1(ir,is)+zm12(ir,is)
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate the denominator
      call multiply(zm21,zf,ztmp2,nmat,nmatx)
      do ir=1,nmat
        do is=1,nmat
          ztmp2(ir,is)=ztmp2(ir,is)+zm22(ir,is)
        enddo
      enddo
      call invers(ztmp2,nmat,nmatx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call multiply(ztmp1,ztmp2,zres,nmat,nmatx)
      return
      end
c
c     *****************************************************************
c
c
c     *****************************************************************
c
      subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
c
      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
c
c     *****************************************************************
c
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
c
c     *****************************************************************
c
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
c
c     *****************************************************************
c
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      double precision s,ars,ais,brs,bis
      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
c
c     *****************************************************************
c
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
c     *****************************************************************
c
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0d0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
c     *****************************************************************
c
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
c
c     *****************************************************************
c
      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi
c
c     (yr,yi) = complex dsqrt(xr,xi) 
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      double precision s,tr,ti,pythag
      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end
c
c     *****************************************************************
c
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
c     *****************************************************************
c
