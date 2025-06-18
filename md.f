      PROGRAM md
       
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (nmax=648, nam=4, carb=216, nli=484)
      DIMENSION rat(3,nam,648), fat(3,nam,648)
      DIMENSION vat(3,nam,648), aat(3,nam,648)
      DIMENSION amat(nam,648), fnet(3,648)
      DIMENSION rij(3)
      DIMENSION nat(648)
      DIMENSION q(nam,648), sig(nam,648), eps(nam,648)
      CHARACTER*3 atom
      
      OPEN(20,file='energy.dat',status='unknown')
      OPEN(30,file='temperature.dat',status='unknown')
      OPEN(40,file='positions_equil1200.dat',status='unknown')
      OPEN(50,file='velocities_equil1200.dat',status='old')
      OPEN(60,file='trajectoryp.dat',status='unknown')
      OPEN(70,file='trajectoryv.dat',status='unknown')
      OPEN(100,file='mdpositions.dat',status='old')
      OPEN(300,file='potential.dat',status='unknown')
      OPEN(310,file='pressure.dat',status='unknown')
      OPEN(320,file='volume.dat',status='unknown')

      a=0.2*1.0e10
      req=1.16*1.0e-10
      dt=05.0*1.0e-16
      amli=7.0*1.66*1.0e-27
      amc=12.0*1.66*1.0e-27
      amo=16.0*1.66*1.0e-27
      amk=39.1*1.66*1.0e-27
      e=1.602*1.0e-19
      qc=1.04085*e
      qo=-0.89429*e
      qli=0.82101*e
      qk=qli
      cke=8.9875*1.0e9
      p=0.22*1.0e-10
      boxl=24.58901398*1.0e-10
      rc=boxl/2.0
      pi=3.14159
      dvx=0
      dvy=0
      dvz=0
      dxi=0
      dyi=0
      dzi=0
      tave=900.0
      boltz=1.38064852*1.0d-23
      bo=4.3369*1.0e-15
      p=0.22*1.0e-10
      forcek=1015.969777
	  
      DO j=1,carb
         nat(j)=4
      END DO
      DO j=carb+1,nmax
         nat(j)=1
      END DO
	  
      DO i=1,carb
         amat(1,i)=amc
         q(1,i)=qc
         sig(1,i)=1.70*1.0e-10
         eps(1,i)=SQRT(5.97515775*1.0e-22)
         DO iat=2,4
            amat(iat,i)=amo
            q(iat,i)=qo
            sig(iat,i)=1.48*1.0e-10
            eps(iat,i)=SQRT(1.45905015*1.0e-21)
         END DO
      END DO
      DO i=carb+1,nli
         amat(1,i)=amli
         q(1,i)=qli
         sig(1,i)=1.013*1.0e-10
         eps(1,i)=SQRT(1.271458*1.0e-22)
      END DO
      DO i=nli+1,nmax
         amat(1,i)=amk
         q(1,i)=qk
         sig(1,i)=2.368*1.0e-10
         eps(1,i)=SQRT(2.27889738*1.0e-24)
      END DO
	  
	  DO i=1,nmax
	   DO iat=1,nat(i)
	    READ(100,*) rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
	    READ(50,*) vat(1,iat,i), vat(2,iat,i), vat(3,iat,i)
	   END DO
	  END DO
	  
          DO istep=1,100000
             t=t+dt
             IF(MOD(istep,1000).EQ.0)PRINT*, istep
         
c     zero out some stuff with loops
      
      DO i=1,nmax
         DO iat=1,nam
            fat(1,iat,i)=0.0
            fat(2,iat,i)=0.0
            fat(3,iat,i)=0.0
            fnet(1,i)=0.0
            fnet(2,i)=0.0
            fnet(3,i)=0.0
            v=0
            vb=0
            vbo=0
            v2=0
            vintra=0
            vljtot=0
            vctot=0
            sumk=0
            pvir=0
            palt=0
         END DO
      END DO

c     step loop // keep as much as possible out of here
      
      DO i=1,nmax
         DO j=i+1,nmax
            DO iat=1,nat(i)
               DO jat=1,nat(j)
                  rij(1)=rat(1,jat,j)-rat(1,iat,i)
                  rij(2)=rat(2,jat,j)-rat(2,iat,i)
                  rij(3)=rat(3,jat,j)-rat(3,iat,i)
                  rij(1)=rij(1)-IDNINT(rij(1)/boxl)*boxl
                  rij(2)=rij(2)-IDNINT(rij(2)/boxl)*boxl
                  rij(3)=rij(3)-IDNINT(rij(3)/boxl)*boxl 
                  r=DSQRT((rij(1))**2+(rij(2))**2+(rij(3))**2)
                  bracket=derfc(a*r)/(r**2)+2.0*a*DEXP(-a**2*r**2)/(r*DS
     *                 QRT(pi))-(b1+2.0*a*DEXP(-a**2*rc**2)/(rc*DSQRT(pi
     *                 )))
                  IF(r.LE.rc)THEN  
                     charge=cke*q(iat,i)*q(jat,j)
                     vijc1=derfc(a*r)
                     vijc2=derfc(a*rc)
                     vijc3=2.0*a*dexp(-a**2*rc**2)
                     vijc4=2.0*a*dexp(-a**2*r**2)
                     vijc=charge*(vijc1/r-vijc2/rc+(vijc2/rc**2
     *                    +vijc3/rc*DSQRT(pi))*(r-rc))
                     ep=eps(iat,i)*eps(jat,j)
                     sigma=sig(iat,i)+sig(jat,j)
                     sr=sigma/r
                     sr6=sr**6
                     sr12=sr6*sr6
                     vlj=4.0*ep*(sr12-sr6)
                     vljtot=vljtot+vlj
                     vctot=vctot+vijc
                     v=v+vijc+vlj
                     flj=4.0*ep*(-12.*sr12/r+6.*sr6/r)
                     fc=-charge*((vijc1/r**2+vijc4/r/DSQRT(pi))-
     *                    (vijc2/rc**2+vijc3/rc/DSQRT(pi)))
                     fat(1,iat,i)=fat(1,iat,i)+(flj+fc)*rij(1)/r
                     fat(2,iat,i)=fat(2,iat,i)+(flj+fc)*rij(2)/r
                     fat(3,iat,i)=fat(3,iat,i)+(flj+fc)*rij(3)/r
                     fat(1,jat,j)=fat(1,jat,j)-(flj+fc)*rij(1)/r
                     fat(2,jat,j)=fat(2,jat,j)-(flj+fc)*rij(2)/r
                     fat(3,jat,j)=fat(3,jat,j)-(flj+fc)*rij(3)/r
                     pvir=pvir-(flj+fc)*r
                  END IF
               END DO
            END DO
         END DO
      END DO

c     intramolecular potential
      
      DO i=1,carb
         DO jat=2,4
            rij(1)=rat(1,jat,i)-rat(1,1,i)             
            rij(2)=rat(2,jat,i)-rat(2,1,i)          
            rij(3)=rat(3,jat,i)-rat(3,1,i)
            r=DSQRT((rij(1))**2+(rij(2))**2+(rij(3))**2)
            vijb=0.5*forcek*(r-req)**2
            vb=vb+vijb
            fat(1,1,i)=fat(1,1,i)+forcek*(r-req)*(1/r)*rij(1)
            fat(2,1,i)=fat(2,1,i)+forcek*(r-req)*(1/r)*rij(2)
            fat(3,1,i)=fat(3,1,i)+forcek*(r-req)*(1/r)*rij(3)
            fat(1,jat,i)=fat(1,jat,i)-forcek*(r-req)*(1/r)*rij(1)
            fat(2,jat,i)=fat(2,jat,i)-forcek*(r-req)*(1/r)*rij(2)
            fat(3,jat,i)=fat(3,jat,i)-forcek*(r-req)*(1/r)*rij(3)
         END DO
      END DO
      
      DO i=1,carb
         DO iat=2,3
            DO jat=iat+1,4
               rij(1)=rat(1,jat,i)-rat(1,iat,i)
               rij(2)=rat(2,jat,i)-rat(2,iat,i)
               rij(3)=rat(3,jat,i)-rat(3,iat,i)
               r=DSQRT((rij(1))**2+(rij(2))**2+(rij(3))**2)
               vijbo=bo*dexp(-r/p)
               vbo=vbo+vijbo
               fat(1,iat,i)=fat(1,iat,i)-((bo/p)*dexp(-r/p))*rij(1)/r
               fat(2,iat,i)=fat(2,iat,i)-((bo/p)*dexp(-r/p))*rij(2)/r
               fat(3,iat,i)=fat(3,iat,i)-((bo/p)*dexp(-r/p))*rij(3)/r
               fat(1,jat,i)=fat(1,jat,i)+((bo/p)*dexp(-r/p))*rij(1)/r
               fat(2,jat,i)=fat(2,jat,i)+((bo/p)*dexp(-r/p))*rij(2)/r
               fat(3,jat,i)=fat(3,jat,i)+((bo/p)*dexp(-r/p))*rij(3)/r
            END DO
         END DO
      END DO
	  vintra=vb+vbo
      
c     calculation of total force
       
      DO i=1,nmax
         DO iat=1,nat(i)
            fnet(1,i)=fnet(1,i)+fat(1,iat,i)
            fnet(2,i)=fnet(2,i)+fat(2,iat,i)
            fnet(3,i)=fnet(3,i)+fat(3,iat,i)
         END DO
      END DO

c     calculate pressure with alternate expression
       
      DO i=1,nmax
         DO idim=1,3
            palt=palt+fnet(idim,i)*rat(idim,1,i)
         END DO
      END DO
       
c     calculation of du/dl

       v2=0
       dl=1.0d-14
       boxl2=boxl+dl
       DO i=1,nmax-1
          DO j=i+1,nmax
             DO iat=1,nat(i)
                DO jat=1,nat(j)
                   rij(1)=rat(1,jat,j)-rat(1,iat,i)
                   rij(2)=rat(2,jat,j)-rat(2,iat,i)
                   rij(3)=rat(3,jat,j)-rat(3,iat,i)
                   rij(1)=rij(1)-IDNINT(rij(1)/boxl2)*boxl2
                   rij(2)=rij(2)-IDNINT(rij(2)/boxl2)*boxl2
                   rij(3)=rij(3)-IDNINT(rij(3)/boxl2)*boxl2
                   r=DSQRT((rij(1))**2+(rij(2))**2+(rij(3))**2)
                   IF(r.LE.rc)THEN
                      charge=cke*q(iat,i)*q(jat,j)
                      vijc1=derfc(a*r)
                      vijc2=derfc(a*rc)
                      vijc3=2.0*a*dexp(-a**2*rc**2)
                      vijc4=2.0*a*dexp(-a**2*r**2)
                      vijc=charge*(vijc1/r-vijc2/rc+(vijc2/rc**2
     *                     +vijc3/rc*DSQRT(pi))*(r-rc))
                      ep=eps(iat,i)*eps(jat,j)
                      sigma=(sig(iat,i)+sig(jat,j))
                      sr=sigma/r
                      sr6=sr**6
                      sr12=sr6*sr6
                      vlj=4.0*ep*(sr12-sr6)
                      v2=v2+vijc+vlj 
                   END IF
                END DO
             END DO
          END DO
       END DO
       deriv=(v2-v)/dl                     
       
c     step forward in time

       DO i=1,nmax
          DO iat=1,nat(i)
             aat(1,iat,i)=fat(1,iat,i)/amat(iat,i)
             aat(2,iat,i)=fat(2,iat,i)/amat(iat,i)
             aat(3,iat,i)=fat(3,iat,i)/amat(iat,i)
             dvx=aat(1,iat,i)*dt/2.0
             dvy=aat(2,iat,i)*dt/2.0
             dvz=aat(3,iat,i)*dt/2.0
             vat(1,iat,i)=vat(1,iat,i)+dvx
             vat(2,iat,i)=vat(2,iat,i)+dvy
             vat(3,iat,i)=vat(3,iat,i)+dvz
             sumk=sumk+0.5*amat(iat,i)*(vat(1,iat,i)**2+vat(2,iat,i)*
     *            *2+vat(3,iat,i)**2)            
             vat(1,iat,i)=vat(1,iat,i)+dvx
             vat(2,iat,i)=vat(2,iat,i)+dvy
             vat(3,iat,i)=vat(3,iat,i)+dvz
             dxi=vat(1,iat,i)*dt
             dyi=vat(2,iat,i)*dt
             dzi=vat(3,iat,i)*dt
             rat(1,iat,i)=rat(1,iat,i)+dxi
             rat(2,iat,i)=rat(2,iat,i)+dyi
             rat(3,iat,i)=rat(3,iat,i)+dzi              
          END DO
       END DO
       
       vol=boxl**3
       pvir=pvir/3.0/vol+nmax*boltz*tave/vol
       palt=palt/3.0/vol-deriv/3.0/boxl/boxl+nmax*boltz*tave/vol

       ipressure=1
       IF(ipressure.EQ.1)THEN
       gamma=5.0d-14
       pdes=1.0e9
       pscalevir=(1+gamma*(pvir-pdes))**(1.0/3.0)
       boxl=boxl*pscalevir
       
       DO i=1,nmax
          DO iat=1,nat(i)
             rat(1,iat,i)=rat(1,iat,i)*pscalevir
             rat(2,iat,i)=rat(2,iat,i)*pscalevir
             rat(3,iat,i)=rat(3,iat,i)*pscalevir
          END DO
       END DO
      END IF
       
       WRITE(310,*) REAL(t),REAL(pvir),REAL(palt)
	   
	   totenergy=v+vintra+sumk
       
       WRITE(20,*) REAL(t), REAL(v+vintra), REAL(sumk),REAL(totenergy)
          temp=sumk/(1.5*boltz*(6*carb-1))
          WRITE(30,*) REAL(t), temp
          WRITE(300,*) REAL(t), REAL(v+vintra)
          WRITE(320,*) REAL(t), REAL(boxl**3)
       IF(istep.GT.0)THEN
          DO i=1,nmax
             DO iat=1,nat(i)
                vat(1,iat,i)=vat(1,iat,i)*DSQRT(tave/temp)
                vat(2,iat,i)=vat(2,iat,i)*DSQRT(tave/temp)
                vat(3,iat,i)=vat(3,iat,i)*DSQRT(tave/temp)
             END DO
          END DO
       END IF
       IF(MOD(istep,5).EQ.0)THEN
      DO i=1,nmax
         DO iat=1,nat(i)
            WRITE(60,*) REAL(rat(1,iat,i)), REAL(rat(2,iat,i)), REAL(rat
     *           (3,iat,i))
            WRITE(70,*) REAL(vat(1,iat,i)), REAL(vat(2,iat,i)), REAL(vat
     *           (3,iat,i))
         END DO
      END DO
      END IF
      END DO
      OPEN(40,file='positions.dat',status='unknown')
      OPEN(50,file='velocities.dat',status='unknown')
      DO i=1,nmax
         DO iat=1,nat(i)
            WRITE(40,*) rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
            WRITE(50,*) vat(1,iat,i), vat(2,iat,i), vat(3,iat,i)
         END DO
      END DO
      STOP
      END
	
