! --------------------------------------------------------------
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
! CODE BY:
! --------   |         |   ---------        /\        |\      |
! |          |         |  |                /  \       | \     |
! |          |         |  |               /    \      |  \    |
! --------   |         |  |   ------|    /------\     |   \   |
!         |  |         |  |         |   /        \    |    \  |
!         |  |         |  |         |  /          \   |     \ |
! ---------   ----------  ----------  /            \  |      \|
! --------------------------------------------------------------

! #########################
! MODULE: system_auxilaries
! LAST MODIFIED: 2 JUNE 2020
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! AUXILARY FUNCTIONS FOR 3D NSE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_auxilaries
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This module has subroutines, that are minor functions used everywhere in the NSE code.
! These is a completely independent module
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_constants
	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT  NONE

  CONTAINS

	SUBROUTINE find_CFL_timestep(dt_max1,dt_max2,delta_t)
	! CALL this to find a nice rounded of time step based on limits
		IMPLICIT  NONE
		DOUBLE PRECISION,INTENT(IN)::dt_max1,dt_max2
		DOUBLE PRECISION,INTENT(OUT)::delta_t

    delta_t = 10.0D0 ** DBLE( FLOOR( DLOG10( dt_max1 ) ) )

    DO WHILE( delta_t .GE. dt_max2 )
      delta_t = 10.0D0 ** DBLE( FLOOR( DLOG10( hf *delta_t ) ) )
    END DO

    IF( delta_t * 5.0D0 .LT. dt_max2 ) THEN
      delta_t = delta_t * 5.0D0
    END IF

	END

	SUBROUTINE find_diss_rate_ref( en, k_int, ds )
	! CALL this to get a rough estimate of dissipation rate
		IMPLICIT  NONE
		DOUBLE PRECISION,INTENT(IN)::en
		DOUBLE PRECISION,INTENT(OUT)::ds
		INTEGER(KIND=4),INTENT(IN)::k_int
		DOUBLE PRECISION::adjust_factor

		adjust_factor = 3.0D0
		! Refined adjusting factor

		ds = twothird * en / C_kolmo
		ds = ds ** 1.5D0
		ds = ds / DBLE(k_int)
		ds = adjust_factor * ds
		! Estimated dissipation rate

	END

	SUBROUTINE step_to_time_convert(step,time,delta_t)
	! CALL this to convert time step into actual time of simulation
		IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(IN)::step
		DOUBLE PRECISION,INTENT(IN)::delta_t
		DOUBLE PRECISION,INTENT(OUT)::time
		time=DBLE( step ) * delta_t
	END

	SUBROUTINE time_to_step_convert(time,step,delta_t)
	! CALL this to convert time step into actual time of simulation
		IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(OUT)::step
		DOUBLE PRECISION,INTENT(IN)::delta_t
		DOUBLE PRECISION,INTENT(IN)::time
		step  = CEILING( time / delta_t )
	END

	SUBROUTINE init_random_seed
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! randomize the seed, to generate random numbers. If not called, everytime it will
	! provide same set of random numbers
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! -----------------------------------------------------------------
    ! THIS PRODUCES NEW RANDOM NUMBERS EVERYTIME CALLED (uniform dist.)
    ! -----------------------------------------------------------------
    INTEGER, DIMENSION(8) :: seed_values
    INTEGER(KIND=4)::seed_size
    ! Declare an assumed shape, dynamic array
    INTEGER, DIMENSION(:),ALLOCATABLE :: seed
    ! gfortran SUBROUTINE to return date and time INFORMATion
    ! from the real time system clock. Works DOwn to milliseconds
    ! and stores the eight return values IN array values.
    CALL DATE_AND_TIME(VALUES=seed_values)
    ! restart the state of the pseuDOranDOm number generator
    ! k = mINimum size of seed (12 on my system)
    seed_size=20
    CALL RANDOM_SEED(size=seed_size)
    ! ALLOCATE memory to seed
    ALLOCATE(seed(seed_size))
    ! assign INFORMATion IN values to seed
    seed(:) = seed_values(:)
    ! seed the ranDOm number generator
    CALL RANDOM_SEED(put=seed)
    ! -----------------------------------------------------------------------

	END

	SUBROUTINE normalization_exponential_spectrum(s_exp,k_I,A_norm)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! gives normalization constant, for Exponentially decayin spectrum. Refer manual.pdf
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
		INTEGER (KIND=4),INTENT(IN)::s_exp,k_I
		DOUBLE PRECISION,INTENT(OUT)::A_norm
		DOUBLE PRECISION::gaussian_integral

		gaussian_integral	=	double_fac( s_exp + 1 ) * DSQRT( hf * two_pi )
		gaussian_integral = gaussian_integral / ( two ** ( hf * DBLE( s_exp ) + one ) )
		gaussian_integral = gaussian_integral / ( ( hf * s_exp ) ** ( hf * s_exp + hf ))

		A_norm = DSQRT( one / ( two_pi * ( k_I ** thr ) * gaussian_integral ) )

	END

	SUBROUTINE normalization_kolmogorov_spectrum(k_I,k_D,A_norm)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! gives normalization constant, for Kolmogorov spectrum. Refer manual.pdf
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
		INTEGER (KIND=4),INTENT(IN)::k_I,k_D
		DOUBLE PRECISION,INTENT(OUT)::A_norm

		A_norm  =  ( ( k_D * k_I ) ** twothird ) / ( k_D ** twothird  - k_I ** twothird )
		A_norm  =  DSQRT( A_norm / ( two_pi * thr ) )

	END

	SUBROUTINE kolmogorov_spectrum_integralscale_subpart(k_ratio,s_exp,factor)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! gives normalization constant, for Kolmogorov spectrum. Refer manual.pdf
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
		INTEGER (KIND=4),INTENT(IN)::s_exp
		DOUBLE PRECISION,INTENT(IN)::k_ratio
		DOUBLE PRECISION,INTENT(OUT)::factor
		DOUBLE PRECISION::c_L

		c_L = 6.78D0

		factor = ( k_ratio / DSQRT( k_ratio ** two + c_L) ) ** ( hf * ( fivthird + s_exp ) )

	END

	SUBROUTINE kolmogorov_spectrum_dissipationscale_subpart(k_ratio,factor)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! gives normalization constant, for Kolmogorov spectrum. Refer manual.pdf
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
		DOUBLE PRECISION,INTENT(IN)::k_ratio
		DOUBLE PRECISION,INTENT(OUT)::factor
		DOUBLE PRECISION::beta,c_eta

		c_eta = 0.4D0
		beta  = 5.2D0

		factor = DEXP( - hf * beta * ( ( k_ratio ** 4.0D0 + c_eta ** 4.0D0 ) ** qtr - c_eta ) )

	END

	SUBROUTINE get_simulation_name( sim_char )
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! This provides a string  of format
	!     'run_d121120_t104022' for run dated 12/11/2020 timed 10:40:22
	!
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CHARACTER(LEN=8)::year_char,month_char,date_char
		CHARACTER(LEN=8)::hour_char,min_char,sec_char
		INTEGER,DIMENSION(8)::values
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CHARACTER(LEN=*),INTENT(OUT)::sim_char

		CALL DATE_AND_TIME(VALUES=values)
		! Gets the date and time as integer array

		values(1)   =   MOD(values(1),2000)
		! writing only last two digits of year

		WRITE(year_char,'(I2)')			 values(1)
		WRITE(month_char,'(I2.2)')	 values(2)
		WRITE(date_char,'(I2.2)')		 values(3)
		WRITE(hour_char,'(I2.2)')		 values(5)
		WRITE(min_char,'(I2.2)')	   values(6)
		WRITE(sec_char,'(I2.2)')		 values(7)
		! Self-explained

		sim_char    =  'run_dt_'//TRIM(ADJUSTL(date_char))//TRIM(ADJUSTL(month_char))//&
		TRIM(ADJUSTL(year_char))//'_t_'//TRIM(ADJUSTL(hour_char))//&
		TRIM(ADJUSTL(min_char))//TRIM(ADJUSTL(sec_char))

	END

	FUNCTION double_fac(num)
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS FUNCTION TO:
	! return the double factorial of an even integer
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		INTEGER(KIND=4)::double_fac
		INTEGER(KIND=4)::i_num , num

		i_num      = 1
		double_fac = 1

		DO WHILE ( i_num .LE. num + 1 )
			double_fac = double_fac * i_num
			i_num      = i_num + 2
		END DO

	END FUNCTION double_fac

	DOUBLE PRECISION FUNCTION dot(vec_1,vec_2)
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS FUNCTION TO:
	! return the dot product of two vectors
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! Calculates the dot product.
		DOUBLE PRECISION,INTENT(IN),DIMENSION(3)::vec_1,vec_2

		dot = SUM( vec_1 * vec_2 )

	END FUNCTION dot

	FUNCTION normal_dist( Z, avg, std ) RESULT (dist)
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS FUNCTION TO:
	! Return a Z vector filled with normal distribution of mean 'avg'
	! and standard deviation 'std'
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		DOUBLE PRECISION,INTENT(IN)::avg, std
		INTEGER(KIND = 4),INTENT(IN):: Z
		INTEGER(KIND = 4):: ind
		DOUBLE PRECISION,DIMENSION( Z )::dist
		DOUBLE PRECISION,DIMENSION( Z )::u_dist_1, u_dist_2

		CALL RANDOM_NUMBER( u_dist_1 )
		CALL RANDOM_NUMBER( u_dist_2 )

		DO ind = 1, Z

			dist( ind ) = DSQRT( - two * DLOG( u_dist_1( ind ) ) ) * DCOS( two_pi * u_dist_2( ind ) )
			dist( ind ) = dist( ind ) * std + avg  

		END DO

	END FUNCTION normal_dist

	SUBROUTINE compute_cross_product(vec_a,vec_b,cross_product)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! gives normalization constant, for Exponentially decayin spectrum. Refer manual.pdf
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
		! Calculates the cross product.
		DOUBLE PRECISION,INTENT(IN),DIMENSION(3)::vec_a,vec_b
		DOUBLE PRECISION,INTENT(OUT),DIMENSION(3)::cross_product

		cross_product( 1 ) = vec_a( 2 ) * vec_b( 3 ) - vec_a( 3 ) * vec_b( 2 )
		cross_product( 2 ) = vec_a( 3 ) * vec_b( 1 ) - vec_a( 1 ) * vec_b( 3 )
		cross_product( 3 ) = vec_a( 1 ) * vec_b( 2 ) - vec_a( 2 ) * vec_b( 1 )

	END

END MODULE system_auxilaries
