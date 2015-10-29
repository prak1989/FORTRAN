program blasius
	integer, parameter :: N = 100
	integer :: i,j
	real :: z_older, z_old, z_new, del_x
	real :: eta = 10.0, g1_old_inf, g1_older_inf
	real, dimension(N+1) :: f, g1, g2, x
	real, dimension(N+1) :: k1_f, k2_f, k3_f, k4_f
	real, dimension(N+1) :: k1_g1, k2_g1, k3_g1, k4_g1
	real, dimension(N+1) :: k1_g2, k2_g2, k3_g2, k4_g2
	
	open(10, file = 'Blasius.dat')
	
	del_x = eta/N
	
	do i = 1,N
	  x(i) = i*del_x
	enddo

	z_older = 0.20
	call Runge_Kutta(z_older)
	g1_older_inf = g1(N);

	z_old = 0.25;
	call Runge_Kutta(z_old);
	g1_old_inf = g1(N);

	z_new = z_old - (((g1_old_inf - 1.0)*(z_old - z_older))/(g1_old_inf - g1_older_inf)); !**** Newton - Raphson method ****!

	if ((z_new - z_old)>1E-15) then
		z_older = z_old;
		z_old = z_new;		
	 
		call Runge_Kutta(z_older);
		g1_older_inf = g1(N);
	
		call Runge_Kutta(z_old);
		g1_old_inf = g1(N);

		z_new = z_old - (((g1_old_inf - 1.0)*(z_old - z_older))/(g1_old_inf - g1_older_inf)); !*** Newton-Raphson Method ***!
	end if
	print*, znew
	do i = 1,N
		write(10,*) x(i), f(i), g1(i), g2(i)
	enddo
end program blasius


subroutine Runge_Kutta(z)
	integer, parameter :: N = 10
	real :: z
	real, dimension(N+1) :: f, g1, g2
	real, dimension(N+1) :: k1_f, k2_f, k3_f, k4_f
	real, dimension(N+1) :: k1_g1, k2_g1, k3_g1, k4_g1
	real, dimension(N+1) :: k1_g2, k2_g2, k3_g2, k4_g2
	f(1) = 0.0
	g1(1) = 0.0
	g2(1) = z
	
	del_x = eta/N
	do i=1,N-1
		k1_f(i) = del_x*g1(i);
		k1_g1(i) = del_x*g2(i);
		k1_g2(i) = del_x*(-0.5*f(i)*g2(i))

		k2_f(i) = del_x*(g1(i) + (0.5*k1_g1(i)))
		k2_g1(i) = del_x*(g2(i) + (0.5*k1_g2(i)))
		k2_g2(i) = del_x*(-0.5*(f(i) + (0.5*k1_f(i)))*(g2(i) + (0.5*k1_g2(i))))

		k3_f(i) = del_x*(g1(i) + (0.5*k2_g1(i)))
		k3_g1(i) = del_x*(g2(i) + (0.5*k2_g2(i)))
		k3_g2(i) = del_x*(-0.5*(f(i) + (0.5*k2_f(i)))*(g2(i) + (0.5*k2_g2(i))))
		
		k4_f(i) = del_x*(g1(i) + k3_g1(i))
		k4_g1(i) = del_x*(g2(i) + k3_g2(i))
		k4_g2(i) = del_x*(-0.5*(f(i) + k3_f(i))*(g2(i) + k3_g2(i)))

		f(i+1) = f(i) + ((1.0/6.0)*(k1_f(i) + 2.0*k2_f(i) + 2.0*k3_f(i) + k4_f(i)))
		g1(i+1) = g1(i) + ((1.0/6.0)*(k1_g1(i) + 2.0*k2_g1(i) + 2.0*k3_g1(i) + k4_g1(i)))
		g2(i+1) = g2(i) + ((1.0/6.0)*(k1_g2(i) + 2.0*k2_g2(i) + 2.0*k3_g2(i) + k4_g2(i)))
	enddo
	return
end subroutine Runge_Kutta	




