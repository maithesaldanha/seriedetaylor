      program serietaylor
c
      implicit real*8 (a-h,o-z)
      real*8 h
      external func
c
200   write(*,*) "Determine um valor para X0 inferior a 1"
      !X0 deve ser menor que um para manter uma precis∆o razo†vel
      !Contudo, essa restriá∆o pode variar
      read(*,*) x0
      if (x0.GE.1) goto 200
c
      n=100
      xmax=1
c
      deltax=(xmax-x0)/n

      do 1 i=1,n
      x=x0+deltax*(i-1)
c
      !Derivada segunda simetrizada
      dr2s=(dr1s(x0+h)-dr1s(x0-h))/(2*h)
c
      !SÇrie de Taylor
      st=func(x0)+dr1s(x0)*(x-x0)+(dr2s*(x-x0)**2)/2
c
      !SÇrie de Taylor (resultado anal°tico)
      sta=1+x+(x**2/2)
c
      open(unit=1, file='expo.dat')
      open(unit=2, file='st.dat')
      open(unit=3, file='sta.dat')
c
      !Comparaá∆o entre os resultados
      write(*,*)"x0",x0,"x=",x,"func(x)",func(x),"st=",st,"sta=",sta
      write(1,*)x,func(x)
      write(2,*)x,st
      write(3,*)x,sta
c
 1    enddo
c
      stop
      end
c
      real*8 function func(x)
      implicit real*8 (a-h,o-z)
      func=exp(x)
      return 
      end
c
      !Derivada primeira simetrizada
      real*8 function dr1s(x)
      implicit real*8 (a-h,o-z)
      h=0.001
      dr1s=(func(x+h)-func(x-h))/(2*h)
      return
      end
