      PROGRAM vcf

      PARAMETER nmax=648, nam=4, carb=216, nstep=10000, nli=432
      PARAMETER jmax=1000
      DIMENSION vat(3,nam,nmax)
      DIMENSION vij(3), vij0(3,nmax), nat(nmax)
      DIMENSION vsum(jmax)
      
      OPEN(10,file='vcf10000.dat',status='unknown')

      dt=5.0*1.0e-16
      
      DO i=1,carb
         nat(i)=4
      END DO
      DO i=carb+1,nmax
         nat(i)=1
      END DO

      DO jstep=1,jmax
         vsum(jstep)=0
      END DO
      
      DO istep=1,nstep-jmax,100
         OPEN(20,file='trajectoryv.dat',status='old')
         DO iskip=1,istep-1
            DO i=1,nmax
               DO iat=1,nat(i)
                  READ(20,*) vat(1,iat,i), vat(2,iat,i), vat(3,iat,i)
               END DO
            END DO
         END DO
         DO jstep=1,jmax
            DO i=1,nmax
               DO iat=1,nat(i)
                  READ(20,*) vat(1,iat,i), vat(2,iat,i), vat(3,iat,i)
               END DO
            END DO
            
            DO i=carb+1,nmax
               IF(jstep.EQ.1)THEN
                  vij0(1,i)=vat(1,1,i)
                  vij0(2,i)=vat(2,1,i)
                  vij0(3,i)=vat(3,1,i)
                  v=vij0(1,i)**2+vij0(2,i)**2+vij0(3,i)**2
               ELSE
                  vij(1)=vat(1,1,i)*vij0(1,i)
                  vij(2)=vat(2,1,i)*vij0(2,i)
                  vij(3)=vat(3,1,i)*vij0(3,i)
                  v=vij(1)+vij(2)+vij(3)
               END IF
               vsum(jstep)=vsum(jstep)+v
            END DO
         END DO
         CLOSE(20)
      END DO
      DO jstep=1,jmax
         vavg=vsum(jstep)/(nli*(nstep-jmax))
         t=(jstep-1)*5.0*dt
         WRITE(10,*) REAL(t), vavg
      END DO
      STOP
      END
