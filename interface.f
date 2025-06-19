      PROGRAM interface

      PARAMETER nmax=648, nam=4, carb=216, nstep=5000
      DIMENSION rat(3,nam,nmax), vat(3,nam,nmax), nat(nmax), rij(3)
      OPEN(10,file='trajectoryp.dat',status='unknown')
      OPEN(20,file='interface.xyz',status='unknown')
      WRITE(20,*) '1296'
      WRITE(20,*) 'interface.xyz'

      DO i=1,carb
         nat(i)=4
      END DO
      DO i=carb+1,nmax
         nat(i)=1
      END DO

      DO istep=1,nstep
         DO i=1,nmax
            DO iat=1,nat(i)
               READ(10,*) rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
               vat(1,iat,i)=0
               vat(2,iat,i)=0
               vat(3,iat,i)=0
            END DO
         END DO
      END DO

      boxlxy=24.50329302*1.0e-10
      boxlz=100.0*1.0e-10
      DO i=1,nmax
         DO idim=1,3
 5          IF(rat(idim,1,i).GT.boxlxy)THEN
               DO iat=1,nat(i)
                  rat(idim,iat,i)=rat(idim,iat,i)-boxlxy
               END DO
               GOTO 5
            END IF
 6          IF(rat(idim,1,i).LT.0)THEN
               DO iat=1,nat(i)
                  rat(idim,iat,i)=rat(idim,iat,i)+boxlxy
               END DO
               GOTO 6
            END IF
         END DO
      END DO

      DO i=1,nmax
         DO iat=1,nat(i)
            rat(3,iat,i)=rat(3,iat,i)+((boxlz/2.0)-(boxlxy/2.0))
         END DO
      END DO           
  
      DO i=1,nmax
         DO iat=1,nat(i)
            IF(iat.EQ.1.AND.i.GT.carb) WRITE(20,*) 'Li',rat(1,iat,
     *           i)*1.0e10,rat(2,iat,i)*1.0e10,rat(3,iat,i)*1.0e10
            IF(iat.EQ.1.AND.i.LE.carb) WRITE(20,*) 'C',rat(1,iat,i
     *           )*1.0e10,rat(2,iat,i)*1.0e10,rat(3,iat,i)*1.0e10
            IF(iat.NE.1) WRITE(20,*) 'O',rat(1,iat,i)*1.0e10,rat(2
     *           ,iat,i)*1.0e10,rat(3,iat,i)*1.0e10
         END DO
      END DO
      
      OPEN(30,file='positions_equil.dat',status='unknown')
      OPEN(40,file='velocities_equil.dat',status='unknown')
      DO i=1,nmax
         DO iat=1,nat(i)
            WRITE(30,*)  rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
            WRITE(40,*) vat(1,iat,i), vat(2,iat,i), vat(3,iat,i)
         END DO
      END DO
      STOP
      END
      
      
      
