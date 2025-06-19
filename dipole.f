	PROGRAM dipole
	
	IMPLICIT REAL*8 (a-h,o-z)
	PARAMETER nmax=512, nstep=100000
	DIMENSION x(nmax), y(nmax), z(nmax)
	CHARACTER*2 atype(nmax)
	
	OPEN(10,file='postraj.dat',status='old')
	OPEN(20,file='dipole.dat',status='unknown')
	
	dt=5.0*1.0e-15            
	e=1.602*1.0e-19
	qcl=-1.0*e
	qna=1.0*e
	boltz=1.38065*1.0e-23
	temp=2000.0
	vol=9.0882216*1.0e-26
	
	DO istep=1,nstep
	   DO i=1,nmax
	      READ(10,*) atype(i), x(i), y(i), z(i)
	   END DO
	   
	   dipx=0.0
	   dipy=0.0
	   dipz=0.0
	   
	   DO i=1,nmax
	      IF(atype(i).EQ.'Cl')THEN
		 dipx=dipx+x(i)*qcl
		 dipy=dipy+y(i)*qcl
		 dipz=dipz+z(i)*qcl
	      ELSE IF(atype(i).EQ.'Na')THEN
		 dipx=dipx+x(i)*qna
		 dipy=dipy+y(i)*qna
		 dipz=dipz+z(i)*qna
	      END IF
	   END DO
	   WRITE (20,*) dipx,dipy,dipz
	END DO
	STOP
	END

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
			
				
				
	
