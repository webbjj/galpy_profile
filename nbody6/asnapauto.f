ccc  This makes 'snap.dat.##' from OUT3 for every model
      program snap
      parameter (maxp=100000,maxm=50000,pi=3.141592)
      real    zmb0,rb0,tsc0,vs0,trh0
      real    rg(maxm),xg(3,maxm),vg(3,maxm),vgc(3),xp(3),vp(3),vrat(2)
      real    time,rnpairs,rbar,zmbar,rtide,tidal4,rdens(3),timetcr
      real    tscale,vstar,rc,rnc,vc,rhom,cmax,rscale,rsmin,dmin1
      real    body(maxp),x(3,maxp),v(3,maxp),dummy(20),cmrdot(3)
      real    rhos(maxp),xns(maxp),phis(maxp)
      integer ntot,model,nrun,nk,name(maxp)
      character outfile*22

c      open(1,file='init.unit')
c      read(1,*) zmb0, rb0, tsc0, vs0, trh0, beta, str
c      close(1)
      zmb0=1
      rb0=1
      tsc0=1
      vs0=1
      trh0=1
      beta=1
      str=1
c     maxm:particle
c     maxp:time
c      open(2,file='clorbit.dat')
c      do l = 1, maxm
c             read(2,*,end=101) rg(l),(xg(k,l),k=1,2),(vg(k,l),k=1,3)
c             xg(3,l) = 0.0
c      enddo
 101   modf = maxm
c      print*, "**** the final model number is",modf
c      close(2)

      radinv = 180.0/pi

      open(3,file='conf.3',form='unformatted')
      rewind(3)
      do 100 l = 1, modf
        read(3,end=100) ntot,model,nrun,nk
        print*,l,ntot,nk
        read(3) time,rnpairs,rbar,zmbar,rtide,tidal4,(rdens(k),k=1,3),
     +         timetcr,tscale,vstar,rc,rnc,vc,rhom,cmax,rscale,rsmin,
     +         dmin1,(body(j),j=1,ntot),(rhos(j),j=1,ntot),
     +         (xns(j),j=1,ntot),((x(k,j),k=1,3),j=1,ntot),
     +         ((v(k,j),k=1,3),j=1,ntot),(phis(j),j=1,ntot),
     +         (name(j),j=1,ntot)
         ntot = ntot - int(rnpairs+0.001)
         time = time!*tsc0  ! [step]
         rgc = rg(l)
         xgc = xg(1,l)              ! gc : guiding center
         ygc = xg(2,l)
         xdc = (rdens(1)*rb0+xg(1,l))   ! dc : density center
         ydc = (rdens(2)*rb0+xg(2,l))

         tdy = 0.          ! dynamical time
         veld = 0.         ! velocity dispersion
         velo = 0.     

         do 50 j =1,ntot 
            velo = velo +body(j)*(v(1,j)+v(2,j)+v(3,j))**2.
 50             continue
         veld = sqrt(velo/zmbar)
         tdy = (2.*rtide)/veld

c       transform the inertial coordinate into rotating coordinate
         theta = atan(ygc/xgc)
         cth = cos(theta)
         sth = sin(theta)
         txgc = xgc*cth+ygc*sth
         tygc = -xgc*sth+ygc*cth

         outfile(1:18)='snapshot/snap.dat.'
         write(outfile(19:22),'(i4.4)') model
         open(4,file=outfile,status='new',err=100)
c         write(4,*) '    time      rgc      xgc      ygc   theta',
c     +              '    txgc     tygc'
c         write(4,'(7f9.2)') time,rgc,xgc,ygc,theta*radinv,txgc,tygc
         do 60 j=1,ntot
            do 55 k = 1,3
               xp(k) = x(k,j) ! in galactic center coordinate [pc]
               vp(k) = v(k,j)      ! (V_tot - V_sys) [km/sec]
 55                   continue
            txp = (x(1,j)*cth+x(2,j)*sth)*rb0 + txgc
            typ = (-x(1,j)*sth+x(2,j)*cth)*rb0 + tygc
            write(4,'(7f15.8,I15)') body(j),(xp(k),k=1,3),
c           write(4,'(7f15.8,I15)') body(j),(xp(k),k=1,3),
c    +                         (vp(k),k=1,3),name(j) 
     +                         (vp(k),k=1,3),name(j)
 60             continue

         open(8,file='modin.dat',access='append')
         if(model.eq.1) write(8,*)'    #     time   rbar  tscale',
     +                  '  vstar    rc    rtid    tdy     vc    veld'
         write(8,99) model,time,rbar,tscale,vstar,rc*rb0,rtide*rb0,
     +               tdy*tsc0,vc*vs0,veld*vs0
 99          format(i6,0pf9.2,3f7.3,f9.4,2f8.2,f8.1,f8.3)
         open(11,file='rdens.dat',access='append')
         if(model.eq.1) write(11,*)'    #     time   rdens(x,y,z)'
         write(11,98) model,time,rdens(1),rdens(2),rdens(3)
 98          format(i6,0pf9.2,3f10.5)
         close(11)
         close(4)
         close(8)
      
 100      continue

 999       stop 
      end

