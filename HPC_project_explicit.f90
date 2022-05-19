Program ExplicitProject
    implicit none
    real*8 :: rho,c,k,a,l,pi,b,dt,h,g,t
    real*8,allocatable :: u(:),uold(:),f(:),x(:),m(:),analy(:),error3(:)
    integer :: n,i,it,iter
    real*8 :: leftend, rightend,dx,error,error2

    n=100
    !!!! allocate array
    allocate( u(0:n))
    allocate( uold(0:n))
    allocate( x(0:n))
    allocate( f(0:n))
    allocate( m(0:n))
    allocate( analy(0:n))

    !!! define parameter
    pi=4.D0*DATAN(1.D0)
    leftend=0d0
    rightend=1d0
    dx=(rightend-leftend)/(n)
    dt=5.0e-5
    k=1d0
    rho=1d0
    c=1d0
    l=1d0
    g=1d0
    t=0d0
    h=0.001d0
    iter=30000
    a=k/(rho*c)
    b=a*dt/dx**2  !!!! CFL number
    write(*,*) "b=",b

    allocate( error3(1:iter))

    
    !!!!! initialize value
    do i=0,n
        x(i)=leftend+(i)*dx
        u(i)=exp(x(i))
        f(i)=sin(l*pi*x(i))
        m(i)=dt*f(i)/(rho*c)
        analy(i)=f(i)/k/pi/pi  !!!! analytic solution,exact value
    end do

    u(0)=0d0        ! left BC
    u(n)=0d0      ! right BC
        

    !!! main-loop
    do it=1,iter

        u(1)=0d0        ! left BC
        u(n+1)=0d0      !! right BC

        uold=u

        do i=1,n-1
            u(i)=(1-2*b)*uold(i)+b*u(i-1)+b*u(i+1)+m(i)   
        end do 
        
        t=dt+t

        do i=0,n
            error=error+abs(analy(i)-u(i))
        end do 

        error=error/(n+1)

        error3(it)=error

    end do

    !!! main loop end 

    write(*,*) "done"


    !!!  plot output part

    OPEN(3,FILE="result.plt", STATUS='UNKNOWN') 
    write(3,*) "x","y","m","f"
    do i=1,n+1
        write(3,*) x(i),u(i),f(i),m(i)
    end do


    !!! error 
    error=0d0;error2=0d0
    
    do i=1,n+1
        error=error+abs(analy(i)-u(i))
    end do 
    error=error/(n+1)
    do i=1,n+1
        error2=max(error2,(abs(analy(i)-u(i))))
    end do 
    write(*,*) "ave-error=",error
    write(*,*) "max-error=",error2

    !!! error-t plot

        OPEN(5,FILE="resultt.plt", STATUS='UNKNOWN') 
        write(5,*) "x","y"
        
        do i=1,iter
            write(5,*) i,error3(i)
        end do


    close(3)
    close(4)
    close(5)
    deallocate(u)
    deallocate(uold)
    deallocate(x)
    deallocate(f)
    deallocate(m)

end Program