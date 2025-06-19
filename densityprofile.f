      PROGRAM densityprofile

      PARAMETER nmax=648, nam=4, carb=216, ibinmax=500, nstep=100000
      DIMENSION rat(3,nam,nmax), nat(nmax)
      DIMENSION histc(ibinmax), histli(ibinmax), rij(3)
      OPEN(10,file='trajectoryp.dat',status='old')
      OPEN(20,file='densityprofile.dat',status='unknown')

      boxl=24.50329302*1.0e-10
      dz=boxl/50.0

      DO i=1,carb
         nat(i)=4
      END DO
      DO i=carb+1,nmax
         nat(i)=1
      END DO
      
      DO ibin=1,ibinmax
         histc(ibin)=0
         histli(ibin)=0
      END DO
      
         DO istep=1,nstep            
            DO i=1,nmax
               DO iat=1,nat(i)
                  READ(10,*) rat(3,iat,i)
               END DO
            END DO
            
            DO i=1,nmax
               rij(3)=rat(3,1,i)
               z=rij(3)
 5             IF(z.GT.boxl)THEN
                  z=z-boxl
                  GOTO 5
               END IF
 6             IF(z.LT.0)THEN
                  z=z+boxl
                  GOTO 6
               END IF
               ibin=INT(z/dz)+1.0
               IF(i.LE.carb)THEN
                  histc(ibin)=histc(ibin)+1.0
               ELSE
                  histli(ibin)=histli(ibin)+1.0
               END IF
            END DO
         END DO
         
         DO ibin=1,ibinmax
            vslab=dz*(boxl**2)
            WRITE(20,*) (ibin+0.5)*dz, histc(ibin)/vslab/nstep, histli(i
     *           bin)/vslab/nstep
         END DO
         STOP
         END
      
      
