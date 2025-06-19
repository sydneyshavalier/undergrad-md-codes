      PROGRAM diffusionli

      PARAMETER nmax=648, nam=4, carb=216, nstep=100000, nli=484
      PARAMETER jmax=100
      DIMENSION rat(3,nam,nmax)
      DIMENSION rij(3), rij0(3,nmax), nat(nmax)
      DIMENSION rsum(jmax)
      
      OPEN(10,file='diffusionk2.dat',status='unknown')

      dt=5.0*1.0e-16
      
      DO i=1,carb
         nat(i)=4
      END DO
      DO i=carb+1,nmax
         nat(i)=1
      END DO

      DO jstep=1,jmax
         rsum(jstep)=0
      END DO
      
      DO istep=1,nstep-jmax,1000
         OPEN(20,file='trajectorypac.dat',status='old')
         DO iskip=1,istep-1
            DO i=1,nmax
               DO iat=1,nat(i)
                  READ(20,*) rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
               END DO
            END DO
         END DO
         DO jstep=1,jmax
            DO i=1,nmax
               DO iat=1,nat(i)
                  READ(20,*) rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
               END DO
            END DO
            
            DO i=nli+1,nmax
               IF(jstep.EQ.1)THEN
                  rij0(1,i)=rat(1,1,i)
                  rij0(2,i)=rat(2,1,i)
                  rij0(3,i)=rat(3,1,i)
                  r=0
               ELSE
                  rij(1)=rat(1,1,i)-rij0(1,i)
                  rij(2)=rat(2,1,i)-rij0(2,i)
                  rij(3)=rat(3,1,i)-rij0(3,i)
                  r=rij(1)**2+rij(2)**2+rij(3)**2
               END IF
               rsum(jstep)=rsum(jstep)+r
            END DO
         END DO
         CLOSE(20)
      END DO
      DO jstep=1,jmax
         ravg=rsum(jstep)/(163.0*(nstep-jmax)/1000.0)
         t=(jstep-1)*5.0*dt
         WRITE(10,*) REAL(t), ravg
      END DO
      STOP
      END
      
