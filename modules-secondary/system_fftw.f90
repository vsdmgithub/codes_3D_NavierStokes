module system_fftw
	! HEADER FILES/MODULES INCLUSION
    ! ----------------------
    use,intrinsic::iso_c_binding ! Standard module which defines the equivalent of C types in fortran
    implicit none
    include 'fftw3.f03'  ! Fortran interface files for all of the C routines for FFTW operation
    ! C VARIABLES DECLARATION
    ! -----------------------
    integer(C_INT)::N0
    type(C_PTR)::plan_r2c,plan_c2r,cdata_r2c_in,cdata_r2c_out,cdata_c2r_in,cdata_c2r_out ! all fftw plans are of this datatype in FORTRAN
    complex(C_DOUBLE_COMPLEX),pointer::data_r2c_out(:,:,:)
    real(C_DOUBLE),pointer::data_r2c_in(:,:,:)
    complex(C_DOUBLE_COMPLEX),pointer::data_c2r_in(:,:,:)
    real(C_DOUBLE),pointer::data_c2r_out(:,:,:)
    ! -----------------------------------------------
    ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    ! INFO: We have two subroutines, for vector fourier
    ! and inverse fourier transforms (DFT)
    contains
    subroutine fft_r2c(in_x,in_y,in_z,M,Mh,out_x,out_y,out_z)
    ! Call this with real input array 'in' and get spectral output array 'out'
        implicit none
        integer(kind=4),intent(in)::M,Mh
        double precision,dimension(0:M-1,0:M-1,0:M-1),intent(in)::in_x,in_y,in_z
        double complex,dimension(0:Mh,-Mh:Mh-1,-Mh:Mh-1),intent(out)::out_x,out_y,out_z
        integer(kind=4)::j_x,j_y,j_z
        double precision::M3
        N0=M
        M3=DBLE(M*M*M)! Normalization factor in the FFT
        ! ALLOCATE ARRAYS - DYNAMIC
        ! NOTE:- The array dimensions are in reverse order for FORTRAN
        ! ---------------------------------------
        cdata_r2c_in=fftw_alloc_real(int(N0*N0*N0,C_SIZE_T))
        call c_f_pointer(cdata_r2c_in,data_r2c_in,[N0,N0,N0])
        cdata_r2c_out=fftw_alloc_complex(int((N0/2+1)*N0*N0,C_SIZE_T))
        call c_f_pointer(cdata_r2c_out,data_r2c_out,[(N0/2+1),N0,N0])
        ! PLAN FOR OUT-PLACE FORWARD DFT R2C
        ! -----------------------------------
        plan_r2c=fftw_plan_dft_r2c_3d(N0,N0,N0,data_r2c_in,data_r2c_out,FFTW_ESTIMATE)
        ! INITIALIZE INPUT DATA
        ! ---------------------
        data_r2c_in=in_x
         ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_r2c(plan_r2c,data_r2c_in,data_r2c_out)
        ! WRITE OUTPUT (in format of first Brillouin zone format)
        ! -----------
        do j_x=0,Mh
        do j_y=-Mh,-1
        do j_z=-Mh,-1
        out_x(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+M+1)/M3
        end do
        end do
        do j_y=0,Mh-1
        do j_z=0,Mh-1
        out_x(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+1)/M3
        end do
        end do
        do j_y=0,Mh-1
        do j_z=-Mh,-1
        out_x(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+M+1)/M3
        end do
        end do
        do j_y=-Mh,-1
        do j_z=0,Mh-1
        out_x(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+1)/M3
        end do
        end do
        end do
        ! INITIALIZE INPUT DATA
        ! ---------------------
        data_r2c_in=in_y
        ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_r2c(plan_r2c,data_r2c_in,data_r2c_out)
        ! WRITE OUTPUT  (in format of first Brillouin zone format)
        ! -----------
        do j_x=0,Mh
        do j_y=-Mh,-1
        do j_z=-Mh,-1
        out_y(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+M+1)/M3
        end do
        end do
        do j_y=0,Mh-1
        do j_z=0,Mh-1
        out_y(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+1)/M3
        end do
        end do
        do j_y=0,Mh-1
        do j_z=-Mh,-1
        out_y(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+M+1)/M3
        end do
        end do
        do j_y=-Mh,-1
        do j_z=0,Mh-1
        out_y(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+1)/M3
        end do
        end do
        end do
        ! INITIALIZE INPUT DATA
        ! ---------------------
        data_r2c_in=in_z
        ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_r2c(plan_r2c,data_r2c_in,data_r2c_out)
        ! WRITE OUTPUT  (in format of first Brillouin zone format)
        ! -----------
        do j_x=0,Mh
        do j_y=-Mh,-1
        do j_z=-Mh,-1
        out_z(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+M+1)/M3
        end do
        end do
        do j_y=0,Mh-1
        do j_z=0,Mh-1
        out_z(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+1)/M3
        end do
        end do
        do j_y=0,Mh-1
        do j_z=-Mh,-1
        out_z(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+M+1)/M3
        end do
        end do
        do j_y=-Mh,-1
        do j_z=0,Mh-1
        out_z(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+1)/M3
        end do
        end do
        end do
        ! DESTROY PLANS
        ! -------------
        call fftw_destroy_plan(plan_r2c)
        call fftw_free(cdata_r2c_in)
        call fftw_free(cdata_r2c_out)
    end
    subroutine fft_c2r(in_x,in_y,in_z,M,Mh,out_x,out_y,out_z)
    ! Call this with complex input array 'in' and get real output array 'out'
        implicit none
        integer(kind=4),intent(in)::M,Mh
        double complex,dimension(0:Mh,-Mh:Mh-1,-Mh:Mh-1),intent(in)::in_x,in_y,in_z
        double precision,dimension(M,M,M),intent(out)::out_x,out_y,out_z
        integer(kind=4)::j_x,j_y,j_z
        double precision::M3
        N0=M
        M3=1.0D0 ! Normalization factor in the FFT
        ! ALLOCATE ARRAYS - DYNAMIC
        ! NOTE:- The array dimensions are in reverse order for FORTRAN
        ! ---------------------------------------
        cdata_c2r_in=fftw_alloc_complex(int((N0/2+1)*N0*N0,C_SIZE_T))
        call c_f_pointer(cdata_c2r_in,data_c2r_in,[(N0/2+1),N0,N0])
        cdata_c2r_out=fftw_alloc_real(int(N0*N0*N0,C_SIZE_T))
        call c_f_pointer(cdata_c2r_out,data_c2r_out,[N0,N0,N0])
        ! PLAN FOR OUT-PLACE FORWARD DFT R2C
        ! -----------------------------------
        plan_c2r=fftw_plan_dft_c2r_3d(N0,N0,N0,data_c2r_in,data_c2r_out,FFTW_MEASURE)
        ! INITIALIZE INPUT DATA  (in format of FFTW from first Brillouin zone format)
        ! ---------------------
        do j_x=0,Mh
        do j_y=-Mh,-1
        do j_z=-Mh,-1
        data_c2r_in(j_x+1,j_y+M+1,j_z+M+1)=in_x(j_x,j_y,j_z)
        end do
        end do
        do j_y=0,Mh-1
        do j_z=0,Mh-1
        data_c2r_in(j_x+1,j_y+1,j_z+1)=in_x(j_x,j_y,j_z)
        end do
        end do
        do j_y=0,Mh-1
        do j_z=-Mh,-1
        data_c2r_in(j_x+1,j_y+1,j_z+M+1)=in_x(j_x,j_y,j_z)
        end do
        end do
        do j_y=-Mh,-1
        do j_z=0,Mh-1
        data_c2r_in(j_x+1,j_y+M+1,j_z+1)=in_x(j_x,j_y,j_z)
        end do
        end do
        end do
        ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_c2r(plan_c2r,data_c2r_in,data_c2r_out)
        ! WRITE OUTPUT
        ! ------------
        out_x=data_c2r_out/M3
        ! INITIALIZE INPUT DATA (in format of FFTW from first Brillouin zone format)
        ! ---------------------
        do j_x=0,Mh
        do j_y=-Mh,-1
        do j_z=-Mh,-1
        data_c2r_in(j_x+1,j_y+M+1,j_z+M+1)=in_y(j_x,j_y,j_z)
        end do
        end do
        do j_y=0,Mh-1
        do j_z=0,Mh-1
        data_c2r_in(j_x+1,j_y+1,j_z+1)=in_y(j_x,j_y,j_z)
        end do
        end do
        do j_y=0,Mh-1
        do j_z=-Mh,-1
        data_c2r_in(j_x+1,j_y+1,j_z+M+1)=in_y(j_x,j_y,j_z)
        end do
        end do
        do j_y=-Mh,-1
        do j_z=0,Mh-1
        data_c2r_in(j_x+1,j_y+M+1,j_z+1)=in_y(j_x,j_y,j_z)
        end do
        end do
        end do
        ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_c2r(plan_c2r,data_c2r_in,data_c2r_out)
        ! WRITE OUTPUT
        ! -----------
        out_y=data_c2r_out/M3
        ! INITIALIZE INPUT DATA (in format of FFTW from first Brillouin zone format)
        ! ---------------------
        do j_x=0,Mh
        do j_y=-Mh,-1
        do j_z=-Mh,-1
        data_c2r_in(j_x+1,j_y+M+1,j_z+M+1)=in_z(j_x,j_y,j_z)
        end do
        end do
        do j_y=0,Mh-1
        do j_z=0,Mh-1
        data_c2r_in(j_x+1,j_y+1,j_z+1)=in_z(j_x,j_y,j_z)
        end do
        end do
        do j_y=0,Mh-1
        do j_z=-Mh,-1
        data_c2r_in(j_x+1,j_y+1,j_z+M+1)=in_z(j_x,j_y,j_z)
        end do
        end do
        do j_y=-Mh,-1
        do j_z=0,Mh-1
        data_c2r_in(j_x+1,j_y+M+1,j_z+1)=in_z(j_x,j_y,j_z)
        end do
        end do
        end do
        ! EXECUTE DFT
        ! -----------
        call fftw_execute_dft_c2r(plan_c2r,data_c2r_in,data_c2r_out)
        ! WRITE OUTPUT
        ! -----------
        out_z=data_c2r_out/M3
        ! DESTROY PLANS
        ! -------------
        call fftw_destroy_plan(plan_c2r)
        call fftw_free(cdata_c2r_in)
        call fftw_free(cdata_c2r_out)
    end
  subroutine fft_r2c_scalar(in_sc,M,Mh,out_sc)
  ! Call this with real input array 'in' and get spectral output array 'out'
      implicit none
      integer(kind=4),intent(in)::M,Mh
      double precision,dimension(0:M-1,0:M-1,0:M-1),intent(in)::in_sc
      double complex,dimension(0:Mh,-Mh:Mh-1,-Mh:Mh-1),intent(out)::out_sc
      integer(kind=4)::j_x,j_y,j_z
      double precision::M3
      N0=M
      M3=DBLE(M*M*M)! Normalization factor in the FFT
      ! ALLOCATE ARRAYS - DYNAMIC
      ! NOTE:- The array dimensions are in reverse order for FORTRAN
      ! ---------------------------------------
      cdata_r2c_in=fftw_alloc_real(int(N0*N0*N0,C_SIZE_T))
      call c_f_pointer(cdata_r2c_in,data_r2c_in,[N0,N0,N0])
      cdata_r2c_out=fftw_alloc_complex(int((N0/2+1)*N0*N0,C_SIZE_T))
      call c_f_pointer(cdata_r2c_out,data_r2c_out,[(N0/2+1),N0,N0])
      ! PLAN FOR OUT-PLACE FORWARD DFT R2C
      ! -----------------------------------
      plan_r2c=fftw_plan_dft_r2c_3d(N0,N0,N0,data_r2c_in,data_r2c_out,FFTW_ESTIMATE)
      ! INITIALIZE INPUT DATA
      ! ---------------------
      data_r2c_in=in_sc
       ! EXECUTE DFT
      ! -----------
      call fftw_execute_dft_r2c(plan_r2c,data_r2c_in,data_r2c_out)
      ! WRITE OUTPUT (in format of first Brillouin zone format)
      ! -----------
      do j_x=0,Mh
      do j_y=-Mh,-1
      do j_z=-Mh,-1
      out_sc(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+M+1)/M3
      end do
      end do
      do j_y=0,Mh-1
      do j_z=0,Mh-1
      out_sc(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+1)/M3
      end do
      end do
      do j_y=0,Mh-1
      do j_z=-Mh,-1
      out_sc(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+1,j_z+M+1)/M3
      end do
      end do
      do j_y=-Mh,-1
      do j_z=0,Mh-1
      out_sc(j_x,j_y,j_z)=data_r2c_out(j_x+1,j_y+M+1,j_z+1)/M3
      end do
      end do
      end do
      ! DESTROY PLANS
      ! -------------
      call fftw_destroy_plan(plan_r2c)
      call fftw_free(cdata_r2c_in)
      call fftw_free(cdata_r2c_out)
  end

end module system_fftw
