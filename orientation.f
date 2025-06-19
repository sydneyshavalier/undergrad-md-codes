      PROGRAM orientation

      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER nmax=648, nam=4, carb=216, ibinmax=100, nstep=100000
      DIMENSION rat(3,nam,nmax), nat(nmax)
      DIMENSION histxl(ibinmax), histxm1(ibinmax), histxm2(ibinmax)
      DIMENSION histxm3(ibinmax), histxu(ibinmax)
      DIMENSION histyl(ibinmax), histym1(ibinmax), histym2(ibinmax)
      DIMENSION histym3(ibinmax), histyu(ibinmax)
      DIMENSION histzl(ibinmax), histzm1(ibinmax), histzm2(ibinmax)
      DIMENSION histzm3(ibinmax), histzu(ibinmax)

      OPEN(10,file='trajectoryp.dat',status='old')
      OPEN(20,file='orientationl.dat',status='unknown')
      OPEN(30,file='orientationm1.dat',status='unknown')
      OPEN(40,file='orientationm2.dat',status='unknown')
      OPEN(50,file='orientationm3.dat',status='unknown')
      OPEN(60,file='orientationu.dat',status='unknown')
      OPEN(70,file='totorientation.dat',status='unknown')

      d=1.0/ibinmax

      DO i=1,carb
         nat(i)=4
      END DO
      DO i=carb+1,nmax
         nat(i)=1
      END DO

      DO ibin=1,ibinmax
         histxl(ibin)=0
         histxm1(ibin)=0
         histxm2(ibin)=0
         histxm3(ibin)=0
         histxu(ibin)=0
         histyl(ibin)=0
         histym1(ibin)=0
         histym2(ibin)=0
         histym3(ibin)=0
         histyu(ibin)=0
         histzl(ibin)=0
         histzm1(ibin)=0
         histzm2(ibin)=0
         histzm3(ibin)=0
         histzu(ibin)=0
      END DO

      DO istep=1,nstep
         DO i=1,carb
            DO iat=1,nat(i)
               READ(10,*) rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
            END DO
         END DO
      
         DO i=1,carb
            co2x=rat(1,1,i)-rat(1,2,i)
            co2y=rat(2,1,i)-rat(2,2,i)
            co2z=rat(3,1,i)-rat(3,2,i)
            co2length=DSQRT((co2x**2)+(co2y**2)+(co2z**2))
            unito2x=co2x/co2length
            unito2y=co2y/co2length
            unito2z=co2z/co2length
            co3x=rat(1,1,i)-rat(1,3,i)
            co3y=rat(2,1,i)-rat(2,3,i)
            co3z=rat(3,1,i)-rat(3,3,i)
            co3length=DSQRT((co3x**2)+(co3y**2)+(co3z**2))
            unito3x=co3x/co3length
            unito3y=co3y/co3length
            unito3z=co3z/co3length
            tlength=DSQRT(((unito2y*unito3z)-(unito3y*unito2z))**2
     *           +((unito2z*unito3x)-(unito3z*unito2x))**2
     *           +((unito2x*unito3y)-(unito3x*unito2y))**2)
            crossx=DABS(((unito2y*unito3z)-(unito3y*unito2z))/tlength)
            crossy=DABS(((unito2z*unito3x)-(unito3z*unito2x))/tlength)
            crossz=DABS(((unito2x*unito3y)-(unito3x*unito2y))/tlength)
            ibinx=INT(crossx/d)+1.0
            ibiny=INT(crossy/d)+1.0
            ibinz=INT(crossz/d)+1.0
            IF(rat(3,1,i).LE.40*1.0e-10)THEN
               histxl(ibinx)=histxl(ibinx)+1.0
               histyl(ibiny)=histyl(ibiny)+1.0
               histzl(ibinz)=histzl(ibinz)+1.0
            ELSE IF(rat(3,1,i).GT.40*1.0e-10.AND.rat(3,1,i)
     *              .LE.46.7*1.0e-10)THEN
               histxm1(ibinx)=histxm1(ibinx)+1.0
               histym1(ibiny)=histym1(ibiny)+1.0
               histzm1(ibinz)=histzm1(ibinz)+1.0
            ELSE IF(rat(3,1,i).GT.46.7*1.0e-10.AND.rat(3,1,i)
     *              .LE.53.3*1.0e-10)THEN
               histxm2(ibinx)=histxm2(ibinx)+1.0
               histym2(ibiny)=histym2(ibiny)+1.0
               histzm2(ibinz)=histzm2(ibinz)+1.0
            ELSE IF(rat(3,1,i).GT.53.3*1.0e-10.AND.rat(3,1,i)
     *              .LE.60*1.0e-10)THEN
               histxm3(ibinx)=histxm3(ibinx)+1.0
               histym3(ibiny)=histym3(ibiny)+1.0
               histzm3(ibinz)=histzm3(ibinz)+1.0
            ELSE IF(rat(3,1,i).GT.60*1.0e-10)THEN
               histxu(ibinx)=histxu(ibinx)+1.0
               histyu(ibiny)=histyu(ibiny)+1.0
               histzu(ibinz)=histzu(ibinz)+1.0
            END IF
         END DO
      END DO

      DO ibin=1,ibinmax

         WRITE(20,*) REAL((ibin+0.5)*d), REAL(histxl(ibin)/carb/nstep),
     *REAL(histyl(ibin)/carb/nstep), REAL(histzl(ibin)/carb/nstep)
         WRITE(30,*) REAL((ibin+0.5)*d), REAL(histxm1(ibin)/carb/nstep),
     *REAL(histym1(ibin)/carb/nstep), REAL(histzm1(ibin)/carb/nstep)
         WRITE(40,*) REAL((ibin+0.5)*d),REAL(histxm2(ibin)/carb/nstep),
     *REAL(histym2(ibin)/carb/nstep), REAL(histzm2(ibin)/carb/nstep)
         WRITE(50,*) REAL((ibin+0.5)*d),REAL(histxm3(ibin)/carb/nstep),
     *REAL(histym3(ibin)/carb/nstep), REAL(histzm3(ibin)/carb/nstep)
         WRITE(60,*) REAL((ibin+0.5)*d),REAL(histxu(ibin)/carb/nstep),
     *REAL(histyu(ibin)/carb/nstep), REAL(histzu(ibin)/carb/nstep)

          WRITE(70,*) REAL((ibin+0.5)*d), REAL(histxl(ibin)/carb/nstep),
     *REAL(histyl(ibin)/carb/nstep), REAL(histzl(ibin)/carb/nstep)
         WRITE(70,*) REAL((ibin+0.5)*d), REAL(histxm1(ibin)/carb/nstep),
     *REAL(histym1(ibin)/carb/nstep), REAL(histzm1(ibin)/carb/nstep)
         WRITE(70,*) REAL((ibin+0.5)*d),REAL(histxm2(ibin)/carb/nstep),
     *REAL(histym2(ibin)/carb/nstep), REAL(histzm2(ibin)/carb/nstep)
         WRITE(70,*) REAL((ibin+0.5)*d),REAL(histxm3(ibin)/carb/nstep),
     *REAL(histym3(ibin)/carb/nstep), REAL(histzm3(ibin)/carb/nstep)
         WRITE(70,*) REAL((ibin+0.5)*d),REAL(histxu(ibin)/carb/nstep),
     *REAL(histyu(ibin)/carb/nstep), REAL(histzu(ibin)/carb/nstep)
      END DO
      STOP
      END
      
               

            
            
            


      
      
