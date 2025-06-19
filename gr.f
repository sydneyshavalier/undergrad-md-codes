      PROGRAM gr

      PARAMETER nmax=648, nam=4, carb=216, ibinmax=2500, nstep=20000
      DIMENSION rat(3,nam,nmax), nat(nmax)
      DIMENSION hist(ibinmax), rij(3)
      OPEN(10,file='trajectoryp.dat',status='old')
      OPEN(20,file='gr.dat',status='unknown')

      boxl=24.50329302*1.0e-10
      dr=boxl/2500.0
      pi=3.14159
      p=(nmax-carb)/(boxl**3)

      DO i=1,carb
         nat(i)=4
      END DO
      DO i=carb+1,nmax
         nat(i)=1
      END DO
      
      DO ibin=1,ibinmax
         hist(ibin)=0
      END DO
      
         DO istep=1,nstep            
            DO i=1,nmax
               DO iat=1,nat(i)
                  READ(10,*) rat(1,iat,i), rat(2,iat,i), rat(3,iat,i)
               END DO
            END DO            
            DO i=carb+1,nmax-1
               DO j=i+1,nmax
                  rij(1)=rat(1,1,j)-rat(1,1,i)
                  rij(2)=rat(2,1,j)-rat(2,1,i)
                  rij(3)=rat(3,1,j)-rat(3,1,i)
                  rij(1)=rij(1)-NINT(rij(1)/boxl)*boxl
                  rij(2)=rij(2)-NINT(rij(2)/boxl)*boxl
                  rij(3)=rij(3)-NINT(rij(3)/boxl)*boxl
                  r=SQRT((rij(1))**2+(rij(2))**2+(rij(3))**2)
                  ibin=INT(r/dr)+1.0
                  hist(ibin)=hist(ibin)+1.0
               END DO
            END DO
         END DO
            
         DO ibin=1,ibinmax
            ri=(ibin-1.0)*dr
            vshell=4.0*pi*(dr*ri**2+ri*dr**2+(1/3)*dr**3)
            IF(ibin.EQ.1)vshell=1
            g=hist(ibin)/(vshell*p)
            IF(ri.LT.boxl/2) WRITE(20,*) (ibin+0.5)*dr, g/(216*nstep)
         END DO
         STOP
         END
