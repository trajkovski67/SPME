!============================================================
! MODULE: PME_main
!============================================================
! Implementation of Particle Mesh Ewald (PME) and
! Smooth Particle Mesh Ewald (SPME) algorithms for periodic
! electrostatic energy and force calculations.
!
! References:
!   [1] P. P. Ewald, Ann. Phys. 64, 253–287 (1921).
!   [2] T. Darden, D. York, and L. Pedersen,
!       "Particle mesh Ewald: An N·log(N) method for Ewald sums
!       in large systems," J. Chem. Phys. 98, 10089 (1993).
!   [3] U. Essmann et al.,
!       "A smooth particle mesh Ewald method,"
!       J. Chem. Phys. 103, 8577 (1995).
!   [4] OpenMP Architecture Review Board,
!       "OpenMP Application Programming Interface Version 4.5",
!       https://www.openmp.org/specifications/
!
! Notes:
!   - The code uses OpenMP for shared-memory parallelization.
!   - FFTs are computed via FFTW (Fastest Fourier Transform in the West).
!   - Written for double precision (wp = selected_real_kind(15)).
!
! Author: Aleksandar Trajkovski
!============================================================

module PME_main
    use iso_fortran_env, only : output_unit, error_unit
    use print_matrix,   only : write_vector, write_matrix
    use, intrinsic :: iso_c_binding, only: c_loc, c_double_complex, c_int, c_ptr, c_double
    use omp_lib
    implicit none
    private
    public :: PME_prog

    integer, parameter :: wp = selected_real_kind(15)

contains

!=======================================================================
subroutine PME_prog(input)
    implicit none
    integer, intent(in) :: input
    real(wp) :: a,b,c,alpha,beta,gamma,volume,ewald_coeff,pi
    real(wp) :: rms_global, rms_madelung, max_atom_error, max_atom_madelung_error
    real(wp) :: direct_sum,energy,benchmark
    real(wp) :: a1(3),rx(3),ry(3),rz(3),t1,t2
    integer :: io,i,j,no_ion,k1,k2,k3,Rcut,Gcut,p
    real(wp), allocatable :: ion_frac_coord(:,:),ion_coord(:,:)
    real(wp), allocatable :: coord_rec(:,:),forces(:,:),forces_spme(:,:),forces_conv(:,:)
    real(wp), allocatable :: atom_madelung_rec(:),atom_madelung_real(:),atom_madelung_conv_rec(:)
    real(wp) :: coord(3,3)
    real(wp), allocatable :: atom_diff(:), atom_full_norm(:), atom_madelung_diff(:), atom_madelung_full(:)
    character(len=100) :: line,title,key
    logical :: force_switch
    !============================================
    ! Header
    write(*,'(a)') "============================================================"
    write(*,'(a)') "        PARTICLE MESH EWALD (PME) SIMULATION PROGRAM"
    write(*,'(a)') "============================================================"
    write(*,*)

    !============================================
    ! Input parsing
    read(input,'(A)') line; read(line,*,iostat=io) title
    write(*,'(a,a)') ">> Simulation Title: ", trim(title)
    read(input,'(A)') line; read(line,*,iostat=io) a,b,c
    read(input,'(A)') line; read(line,*,iostat=io) alpha,beta,gamma

    pi=3.141592653589793_wp
    alpha=alpha*pi/180.0_wp; beta=beta*pi/180.0_wp; gamma=gamma*pi/180.0_wp

    write(*,'(a)') "------------------------------------------------------------"
    write(*,'(a)') " Lattice Parameters (Å, degrees):"
    write(*,'(a,3f12.6)') "   a,b,c              = ", a,b,c
    write(*,'(a,3f12.6)') "   alpha,beta,gamma   = ", alpha*180/pi,beta*180/pi,gamma*180/pi
    write(*,'(a)') "------------------------------------------------------------"

    read(input,'(A)') line
    read(line,*,iostat=io) no_ion
    allocate(ion_frac_coord(no_ion,4))
    do i=1,no_ion
       read(input,'(A)') line
       read(line,*,iostat=io) ion_frac_coord(i,1:4)
    enddo

    write(*,'(a,i6)') " Number of ions: ", no_ion
    write(*,*)

    coord=0.0_wp
    coord(1,:)=(/a,0.0_wp,0.0_wp/)
    coord(2,:)=(/b*cos(gamma),b*sin(gamma),0.0_wp/)
    coord(3,:)=(/c*cos(beta), c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma), &
                c*sqrt(1.0_wp-cos(beta)**2-((cos(alpha)-cos(beta)*cos(gamma))/sin(gamma))**2)/)

    read(input,'(A)') line; read(line,*,iostat=io) k1,k2,k3
    read(input,'(A)') line; read(line,*,iostat=io) ewald_coeff,p,Rcut,Gcut,key
    force_switch = .false.
    if (io == 0) then
        if (trim(adjustl(key)) == 'force' .or. trim(adjustl(key)) == 'FORCE') then
            force_switch = .true.
        end if
    end if

    write(*,'(a)') " Grid / PME parameters:"
    write(*,'(a,3i6)') "   Grid size (k1,k2,k3): ", k1,k2,k3
    write(*,'(a,f12.3)') "   Ewald coefficient    : ", ewald_coeff
    write(*,'(a,i3)') "   B-spline order (p)   : ", p
    write(*,'(a,i3)') "   Real-space cutoff    : ", 5/ewald_coeff, "Bohr"
    write(*,'(a,i3)') "   Reciprocal cutoff    : ", Gcut
    write(*,'(a)') "------------------------------------------------------------"
    write(*,*)

    allocate(forces(no_ion,3),forces_spme(no_ion,3),forces_conv(no_ion,3))
    allocate(atom_madelung_rec(no_ion),atom_madelung_real(no_ion),atom_madelung_conv_rec(no_ion))
    atom_madelung_real=0.0_wp; atom_madelung_conv_rec=0.0_wp

    call compute_cross_product(coord(:,1),coord(:,2),a1)
    volume=a*b*c*sqrt(1.0_wp-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+ &
             2.0_wp*cos(alpha)*cos(beta)*cos(gamma))
    write(*,'(a,f15.6)') " Unit cell volume: ", volume
    write(*,*)

    allocate(coord_rec(3,3)); coord_rec=0.0_wp
    call compute_cross_product(coord(:,1),coord(:,2),coord_rec(:,3))
    call compute_cross_product(coord(:,1),coord(:,3),coord_rec(:,2))
    call compute_cross_product(coord(:,2),coord(:,3),coord_rec(:,1))
    coord_rec=abs(coord_rec)/volume
    rx=coord(:,1); ry=coord(:,2); rz=coord(:,3)

    allocate(ion_coord(no_ion,3))
    ion_coord=0.0_wp
    do i=1,no_ion
       do j=1,3
          ion_coord(i,:)=ion_coord(i,:)+ion_frac_coord(i,j)*coord(:,j)
       enddo
    enddo

    allocate(atom_diff(no_ion), atom_full_norm(no_ion), atom_madelung_diff(no_ion), atom_madelung_full(no_ion))

    write(*,'(a)') "============================================================"
    write(*,'(a)') "                 BEGINNING CALCULATIONS"
    write(*,'(a)') "============================================================"
    call cpu_time(t1)

    direct_sum=0.0_wp
    call direct_space_sum(no_ion,real(Rcut,wp),ion_coord,ion_frac_coord,coord,ewald_coeff, direct_sum,atom_madelung_real,forces,force_switch)

    do i=1,no_ion
       direct_sum=direct_sum-ewald_coeff*(ion_frac_coord(i,4))**2/sqrt(pi)
       atom_madelung_real(i)=atom_madelung_real(i)-2.0_wp*ewald_coeff*ion_frac_coord(i,4)/sqrt(pi)
    enddo
    call cpu_time(t2)
    write(*,'(a,f10.4,a)') " Direct-space Ewald completed in ", t2-t1, " s"
    write(*,*)

    call  reciprocal_ewald(no_ion,ion_frac_coord,ion_coord,coord_rec,Gcut,ewald_coeff,volume,benchmark,atom_madelung_conv_rec,forces_conv,force_switch)
    call  spme_energy(ewald_coeff,[k1,k2,k3],p,rx,ry,rz,ion_coord,ion_frac_coord(:,4),no_ion,energy,atom_madelung_rec,forces_spme,force_switch)

    forces_spme=forces_spme+forces
    forces_conv=forces_conv+forces

    !=========================================================
    ! RMS force & Madelung error analysis
    !=========================================================
    do i=1,no_ion
        atom_diff(i)         = sqrt(sum((forces_spme(i,:) - forces_conv(i,:))**2))
        atom_full_norm(i)    = sqrt(sum(forces_conv(i,:)**2))
        atom_madelung_full(i)= atom_madelung_conv_rec(i)
        atom_madelung_diff(i)= atom_madelung_rec(i) - atom_madelung_conv_rec(i)
    end do

    rms_global   = sqrt(sum(atom_diff(:)**2) / sum(atom_full_norm(:)**2))
    rms_madelung = sqrt(sum(atom_madelung_diff(:)**2) / sum(atom_madelung_full(:)**2))
    max_atom_error          = maxval(atom_diff(:))
    max_atom_madelung_error = maxval(abs(atom_madelung_diff(:)))

    write(*,'(a)') "------------------------------------------------------------"
    write(*,'(a)') " ERROR ANALYSIS (Deserno & Holm metrics)"
    write(*,'(a)') "------------------------------------------------------------"
    write(*,'(a,4i6)') " Grid (k1,k2,k3,p): ", k1,k2,k3,p
    write(*,'(a,f14.8)') " RMS force error (relative):   ", rms_global
    write(*,'(a,f14.8)') " Max atom force error:         ", max_atom_error
    write(*,'(a,f14.8)') " RMS Madelung error (relative):", rms_madelung
    write(*,'(a,f14.8)') " Max atom Madelung error:      ", max_atom_madelung_error
    write(*,'(a)') "------------------------------------------------------------"
    write(*,*)

    write(*,'(a)') " ENERGY SUMMARY (Hartree units)"
    write(*,'(a)') "------------------------------------------------------------"
    write(*,'(a,f18.8)') " Direct-space energy          :", direct_sum
    write(*,'(a,f18.8)') " Reciprocal-space (PME)       :", energy
    write(*,'(a,f18.8)') " Reciprocal-space (Ewald)     :", benchmark
    write(*,'(a,f18.8)') " Total PME Electrostatic      :", direct_sum + energy
    write(*,'(a,f18.8)') " Total Ewald Electrostatic    :", direct_sum + benchmark
    write(*,'(a,f18.8)') " PME - Ewald difference       :", energy - benchmark
    write(*,'(a)') "------------------------------------------------------------"
    write(*,*)
    write(*,'(a)') " Calculation completed successfully"
    write(*,'(a)') "============================================================"
   if (force_switch) then
    write(*,'(a)') "------------------------------------------------------------"
    write(*,'(a)') " PER-ATOM FORCES (SPME vs Ewald, in Hartree/Bohr)"
    write(*,'(a)') "------------------------------------------------------------"
    write(*,'(a)') "   Atom        Fconv(1:3)                Fspme(1:3)              ΔF magnitude"
    write(*,'(a)') "------------------------------------------------------------"
    do i = 1, no_ion
        write(*,'(i6,2x,3f10.6,3x,3f10.6,3x,f10.6)') i, &
             forces_conv(i,1:3), forces_spme(i,1:3), &
             sqrt(sum((forces_spme(i,:) - forces_conv(i,:))**2))
    end do
    write(*,'(a)') "------------------------------------------------------------"
    write(*,'(a,3f14.8)') " Net force (SPME):  ", sum(forces_spme(:,1)), sum(forces_spme(:,2)), sum(forces_spme(:,3))
    write(*,'(a,3f14.8)') " Net force (Ewald): ", sum(forces_conv(:,1)), sum(forces_conv(:,2)), sum(forces_conv(:,3))
    write(*,'(a)') "------------------------------------------------------------"
    write(*,*) 
else
    write(*,'(a)') " Force computation disabled (no 'force' keyword in input)."
    write(*,*)
end if
write(*,'(a)') "------------------------------------------------------------"
write(*,'(a)') " REFERENCES:"
write(*,'(a)') "   Ewald (1921), Ann. Phys. 64, 253–287"
write(*,'(a)') "   Darden et al. (1993), J. Chem. Phys. 98, 10089 – PME"
write(*,'(a)') "   Essmann et al. (1995), J. Chem. Phys. 103, 8577 – SPME"
write(*,'(a)') "   OpenMP ARB (2015), OpenMP API 4.5 Specification"
write(*,'(a)') "------------------------------------------------------------"
write(*,*)
 
end subroutine PME_prog

        !============================================================
subroutine compute_cross_product(a,b,c)
 real(wp),intent(in)::a(3),b(3); real(wp),intent(out)::c(3)
 c(1)=a(2)*b(3)-a(3)*b(2)
 c(2)=a(3)*b(1)-a(1)*b(3)
 c(3)=a(1)*b(2)-a(2)*b(1)
end subroutine

pure function cross_product(u,v) result(w)
 real(wp),intent(in)::u(3),v(3); real(wp)::w(3)
 w(1)=u(2)*v(3)-u(3)*v(2)
 w(2)=u(3)*v(1)-u(1)*v(3)
 w(3)=u(1)*v(2)-u(2)*v(1)
end function

subroutine invert_3x3_matrix(A,Ainv)
 real(wp),intent(in)::A(3,3); real(wp),intent(out)::Ainv(3,3); real(wp)::det
 Ainv(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
 Ainv(1,2)=A(1,3)*A(3,2)-A(1,2)*A(3,3)
 Ainv(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
 Ainv(2,1)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
 Ainv(2,2)=A(1,1)*A(3,3)-A(1,3)*A(3,1)
 Ainv(2,3)=A(1,3)*A(2,1)-A(1,1)*A(2,3)
 Ainv(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
 Ainv(3,2)=A(1,2)*A(3,1)-A(1,1)*A(3,2)
 Ainv(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
 det=dot_product(A(1,:),Ainv(1,:)); Ainv=Ainv/det
end subroutine

!============================================================
! Direct sum (OpenMP)
subroutine direct_space_sum(no_ion, Rcut, ion_coord, ion_frac_coord, coord, ewald_coeff, direct_sum, atom_m, forces,force_switch)
  use omp_lib
  implicit none
  !=====================
  ! INPUT
  !=====================
  integer, intent(in)              :: no_ion
  real(wp), intent(in)             :: Rcut
  real(wp), intent(in)             :: ion_coord(no_ion,3)        ! cartesian coords
  real(wp), intent(in)             :: ion_frac_coord(no_ion,4)   ! frac coords + charges
  real(wp), intent(in)             :: coord(3,3)                 ! lattice vectors
  real(wp), intent(in)             :: ewald_coeff
  !=====================
  ! OUTPUT
  !=====================
  real(wp), intent(inout)          :: atom_m(no_ion)             ! Madelung accumulators
  real(wp), intent(out)            :: direct_sum                 ! direct space energy
  real(wp), intent(inout)          :: forces(no_ion,3)           ! forces array
  !=====================
  ! LOCALS
  !=====================
  integer :: ntot, kk, i, j, k, n1, n2, n3, tid, nthreads
  integer :: ncell_max
  real(wp) :: Rcut_bohr, rr, qq, fscalar
  real(wp), allocatable :: trans(:,:)
  real(wp), allocatable :: mad_acc(:,:)
  real(wp) :: dist(3), tvec(3), vij(3), pi
  logical :: force_switch

  ! ---- constants ----
  pi = acos(-1.0_wp)

  ! ---- number of real-space cells ----
  ncell_max = int(Rcut)
  Rcut_bohr = 5.0_wp/ewald_coeff   ! typical definition, adjust as needed

  ! ---- allocate translation vectors ----
  allocate(trans(3,(2*ncell_max+1)**3 - 1))
  kk = 0
  do n1 = -ncell_max, ncell_max
    do n2 = -ncell_max, ncell_max
      do n3 = -ncell_max, ncell_max
        if(n1==0 .and. n2==0 .and. n3==0) cycle
        tvec = n1*coord(:,1) + n2*coord(:,2) + n3*coord(:,3)
        rr   = sqrt(sum(tvec**2))
        kk   = kk + 1
        trans(:,kk) = tvec
      end do
    end do
  end do
  ntot = kk

  ! ---- initialize accumulators ----
  direct_sum = 0.0_wp
  forces     = 0.0_wp
  nthreads   = omp_get_max_threads()
  allocate(mad_acc(no_ion,nthreads))
  mad_acc = 0.0_wp

  !$omp parallel default(none) &
  !$omp shared(no_ion,ion_coord,ion_frac_coord,trans,ntot,ewald_coeff,mad_acc,forces,pi) &
  !$omp private(i,j,k,dist,tvec,rr,qq,fscalar,vij,tid) &
  !$omp reduction(+:direct_sum)

    tid = omp_get_thread_num() + 1

    ! === pair interactions ===
    !$omp do schedule(dynamic)
    do i = 1, no_ion
      do j = i+1, no_ion
        dist = ion_coord(i,:) - ion_coord(j,:)
        qq   = ion_frac_coord(i,4) * ion_frac_coord(j,4)

        ! --- home cell ---
        rr = sqrt(sum(dist**2))
        if( rr <= Rcut_bohr) then
          direct_sum = direct_sum + qq * erfc(ewald_coeff*rr)/rr
          mad_acc(i,tid) = mad_acc(i,tid) + ion_frac_coord(j,4)*erfc(ewald_coeff*rr)/rr
          mad_acc(j,tid) = mad_acc(j,tid) + ion_frac_coord(i,4)*erfc(ewald_coeff*rr)/rr

          fscalar = qq * ( erfc(ewald_coeff*rr)/rr**2 + &
                     2.0_wp*ewald_coeff/sqrt(pi) * exp(-(ewald_coeff*rr)**2)/rr )
          vij = fscalar * dist/rr
          !$omp atomic
          forces(i,1) = forces(i,1) + vij(1)
          !$omp atomic
          forces(i,2) = forces(i,2) + vij(2)
          !$omp atomic
          forces(i,3) = forces(i,3) + vij(3)
          !$omp atomic
          forces(j,1) = forces(j,1) - vij(1)
          !$omp atomic
          forces(j,2) = forces(j,2) - vij(2)
          !$omp atomic
          forces(j,3) = forces(j,3) - vij(3)
        end if

        ! --- periodic images ---
        do k = 1, ntot
          tvec = trans(:,k)
          rr   = sqrt(sum((dist + tvec)**2))
          if( rr <= Rcut_bohr) then
            direct_sum = direct_sum + qq*erfc(ewald_coeff*rr)/rr
            mad_acc(i,tid) = mad_acc(i,tid) + ion_frac_coord(j,4)*erfc(ewald_coeff*rr)/rr
            mad_acc(j,tid) = mad_acc(j,tid) + ion_frac_coord(i,4)*erfc(ewald_coeff*rr)/rr

            fscalar = qq * ( erfc(ewald_coeff*rr)/rr**2 + &
                       2.0_wp*ewald_coeff/sqrt(pi) * exp(-(ewald_coeff*rr)**2)/rr )
            vij = fscalar * (dist+tvec)/rr
            !$omp atomic
            forces(i,1) = forces(i,1) + vij(1)
            !$omp atomic
            forces(i,2) = forces(i,2) + vij(2)
            !$omp atomic
            forces(i,3) = forces(i,3) + vij(3)
            !$omp atomic
            forces(j,1) = forces(j,1) - vij(1)
            !$omp atomic
            forces(j,2) = forces(j,2) - vij(2)
            !$omp atomic
            forces(j,3) = forces(j,3) - vij(3)
          end if
        end do
      end do
    end do
    !$omp end do

    ! === self-image interactions (no forces, only Madelung correction) ===
    !$omp do schedule(dynamic)
    do i = 1, no_ion
      qq = ion_frac_coord(i,4)**2
      do k = 1, ntot
        tvec = trans(:,k)
        rr   = sqrt(sum(tvec**2))
        if( rr <= Rcut_bohr) then
          direct_sum = direct_sum + 0.5_wp * qq * erfc(ewald_coeff*rr)/rr
          mad_acc(i,tid) = mad_acc(i,tid) + ion_frac_coord(i,4) * erfc(ewald_coeff*rr)/rr
        end if
      end do
    end do
    !$omp end do

  !$omp end parallel

  ! collect Madelung accumulators across threads
  do i = 1, no_ion
    atom_m(i) = atom_m(i) + sum(mad_acc(i,:))
  end do

  deallocate(trans, mad_acc)
end subroutine direct_space_sum
!============================================================
! Reciprocal Ewald
    subroutine  reciprocal_ewald(no_ion,ion_frac_coord,ion_coord,coord_rec,Gcut,ewald_coeff,volume,resultat,atom_madelung_conv_rec,forces,force_switch)
        implicit none
        real(wp), intent(in) :: ion_coord(:,:), ion_frac_coord(:,:), coord_rec(:,:), ewald_coeff, volume
        integer, intent(in) :: no_ion, Gcut
        real(wp), intent(out) :: resultat,atom_madelung_conv_rec(no_ion),forces(no_ion,3)
        integer :: m1, m2, m3, i
        real(wp) :: g1(3), g2(3), g3(3), rec_vector(3), g_squared,pi = 2.0_wp*acos(0.0_wp),t1,t2
        complex(wp), parameter :: imaginary = (0.0_wp, 1.0_wp)
        complex(wp) :: S, resultat_complex,S_atom
        logical :: force_switch

        forces=0.0_wp
        g1 = coord_rec(1:3,1)
        g2 = coord_rec(1:3,2)
        g3 = coord_rec(1:3,3)
        resultat_complex = (0.0_wp, 0.0_wp)
        call CPU_TIME(t1)
        do m1 = -Gcut, Gcut
            do m2 = -Gcut, Gcut
                do m3 = -Gcut, Gcut
                    if (m1 /= 0 .or. m2 /= 0 .or. m3 /= 0) then
                        rec_vector = m1*g1 + m2*g2 + m3*g3
                        g_squared = dot_product(rec_vector, rec_vector)
                        S = (0.0_wp, 0.0_wp)
                        do i = 1, no_ion
                            S = S + ion_frac_coord(i,4) * exp(2.0_wp * pi * imaginary * dot_product(rec_vector, ion_coord(i,1:3)))
                        end do
                        resultat_complex = resultat_complex + 0.5_wp* exp(-pi**2 * g_squared / ewald_coeff**2) * (S * conjg(S)) / (volume * pi * g_squared)
                        do i = 1, no_ion
                            forces(i,1:3) = forces(i,1:3) + (1.0_wp * ion_frac_coord(i,4) / (volume*pi)) * &
                            exp(-pi**2 * g_squared / ewald_coeff**2) / g_squared *( 2.0_wp*pi * rec_vector ) * &
                            aimag( exp(2.0_wp*pi*imaginary*dot_product(rec_vector,ion_coord(i,1:3))) * conjg(S) )

                            atom_madelung_conv_rec(i) = atom_madelung_conv_rec(i) + &
                            exp(-pi**2 * g_squared / ewald_coeff**2) / (volume * pi * g_squared) * &
                            real( exp(-2.0_wp * pi * imaginary * dot_product(rec_vector, ion_coord(i,1:3))) * S )
                        end do
                    end if
                end do
            end do
        end do
        call cpu_time(t2)
        print *, 'Elapsed CPU time (normal ewald): ', t2-t1, ' seconds'
        resultat = real(resultat_complex, wp)
    end subroutine reciprocal_ewald

!============================================================
subroutine FFT_INVERSE(nx,ny,nz,C,psi)
 use, intrinsic :: iso_c_binding
 integer(c_int),intent(in)::nx,ny,nz
 complex(c_double_complex),target,intent(in)::C(0:nx-1,0:ny-1,0:nz-1)
 complex(c_double_complex),target,intent(out)::psi(0:nx-1,0:ny-1,0:nz-1)
 type(c_ptr)::plan_inv
 integer(c_int),parameter::FFTW_BACKWARD=1,FFTW_ESTIMATE=64
 interface
   function fftw_plan_dft_3d(n1,n2,n3,in,out,sign,flags) bind(C,name="fftw_plan_dft_3d")
     use iso_c_binding
     integer(c_int),value::n1,n2,n3,sign,flags
     type(c_ptr),value::in,out; type(c_ptr)::fftw_plan_dft_3d
   end function
   subroutine fftw_execute_dft(plan,in,out) bind(C,name="fftw_execute_dft")
     use iso_c_binding; type(c_ptr),value::plan,in,out
   end subroutine
   subroutine fftw_destroy_plan(plan) bind(C,name="fftw_destroy_plan")
     use iso_c_binding; type(c_ptr),value::plan
   end subroutine
 end interface
 plan_inv=fftw_plan_dft_3d(nx,ny,nz,c_loc(C(0,0,0)),c_loc(psi(0,0,0)),FFTW_BACKWARD,FFTW_ESTIMATE)
 call fftw_execute_dft(plan_inv,c_loc(C(0,0,0)),c_loc(psi(0,0,0)))
 call fftw_destroy_plan(plan_inv)
end subroutine

!============================================================

subroutine FFT_forward_3d(nx, ny, nz, C, psi)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: nx, ny, nz
    complex(c_double_complex), target, intent(in)  :: C(0:nx-1,0:ny-1,0:nz-1)
    complex(c_double_complex), target, intent(out) :: psi(0:nx-1,0:ny-1,0:nz-1)
    type(c_ptr) :: plan

    interface
        function fftw_plan_dft_3d(n1,n2,n3,in,out,sign,flags) bind(C,name="fftw_plan_dft_3d")
            use iso_c_binding
            integer(c_int), value :: n1,n2,n3,sign,flags
            type(c_ptr), value :: in,out
            type(c_ptr) :: fftw_plan_dft_3d
        end function

        subroutine fftw_execute_dft(plan, in, out) bind(C,name="fftw_execute_dft")
            use iso_c_binding
            type(c_ptr), value :: plan, in, out
        end subroutine

        subroutine fftw_destroy_plan(plan) bind(C,name="fftw_destroy_plan")
            use iso_c_binding
            type(c_ptr), value :: plan
        end subroutine
    end interface

    ! Plan with FFTW_MEASURE to avoid segfaults
    plan = fftw_plan_dft_3d(nx, ny, nz, c_loc(C(0,0,0)), c_loc(psi(0,0,0)), -1, 0)

    ! Execute
    call fftw_execute_dft(plan, c_loc(C(0,0,0)), c_loc(psi(0,0,0)))

    ! Free plan
    call fftw_destroy_plan(plan)
end subroutine

!============================================================
subroutine fft_forward(input_complex,output_complex,nfft)
 use, intrinsic :: iso_c_binding
 integer(c_int),intent(in)::nfft
 complex(c_double_complex),target,intent(in)::input_complex(0:nfft-1)
 complex(c_double_complex),target,intent(out)::output_complex(0:nfft-1)
 interface
   function fftw_plan_dft_1d(n,in,out,sign,flags) bind(C,name="fftw_plan_dft_1d")
     use iso_c_binding; integer(c_int),value::n,sign,flags
     type(c_ptr),value::in,out; type(c_ptr)::fftw_plan_dft_1d
   end function
   subroutine fftw_execute_dft(plan,in,out) bind(C,name="fftw_execute_dft")
     use iso_c_binding; type(c_ptr),value::plan,in,out
   end subroutine
   subroutine fftw_destroy_plan(plan) bind(C,name="fftw_destroy_plan")
     use iso_c_binding; type(c_ptr),value::plan
   end subroutine
 end interface
 type(c_ptr)::plan
 integer(c_int),parameter::FFTW_FORWARD=-1,FFTW_ESTIMATE=64
 plan=fftw_plan_dft_1d(nfft,c_loc(input_complex(0)),c_loc(output_complex(0)),FFTW_FORWARD,FFTW_ESTIMATE)
 call fftw_execute_dft(plan,c_loc(input_complex(0)),c_loc(output_complex(0)))
 call fftw_destroy_plan(plan)
end subroutine

!============================================================
function one_pass_bspline(array,val,n) result(outarray)
 real(wp),intent(in)::array(:),val; integer,intent(in)::n
 real(wp)::outarray(0:size(array)-1)
 real(wp)::div; integer::j
 outarray=array; div=1.0_wp/real(n-1,kind=wp)
 outarray(n-1)=div*val*outarray(n-2)
 do j=1,n-2
   outarray(n-j-1)=div*((val+real(j,kind=wp))*outarray(n-j-2)+ &
                 (real(n-j,kind=wp)-val)*outarray(n-j-1))
 enddo
 outarray(0)=outarray(0)*div*(1.0_wp-val)
end function

function build_spline(order,u) result(spline)
 integer,intent(in)::order; real(wp),intent(in)::u
 real(wp)::spline(0:order-1),array(0:order-1)
 integer::m
 array=0.0_wp; array(0)=1.0_wp-u; array(1)=u
 do m=1,order-2; array=one_pass_bspline(array,u,m+2); enddo
 spline=array
end function

function build_spline_derivative(order,u) result(dspline)
 integer,intent(in)::order; real(wp),intent(in)::u
 real(wp)::dspline(0:order-1),array(0:order-1),prev(0:order-1)
 integer::m
 array=0.0_wp; array(0)=1.0_wp-u; array(1)=u
 do m=1,order-3; array=one_pass_bspline(array,u,m+2); enddo
 prev=array; array=one_pass_bspline(array,u,order)
 dspline=0.0_wp; dspline(0)=-prev(0)
 do m=1,order-2; dspline(m)=prev(m-1)-prev(m); enddo
 dspline(order-1)=prev(order-2)
end function

function dftmod(order,n) result(modulus)
 integer,intent(in)::order,n
 real(wp)::modulus(0:n-1)
 complex(wp)::transform(0:n-1),spline_complex(0:n-1)
 real(wp)::spline(0:n-1); integer::i
 spline=0.0_wp; spline(0:order-1)=build_spline(order,0.0_wp)
 spline_complex=cmplx(spline,0.0_wp,kind=wp)
 call fft_forward(spline_complex,transform,n)
 do i=0,n-1; modulus(i)=real(transform(i),kind=wp)**2+aimag(transform(i))**2; enddo
end function


subroutine spme_energy(ewald_coeff,grid_size,splineorder,a,b,c,coords,charges,natoms,energy,atom_m,forces,force_switch)
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none

  ! Arguments
  real(wp),intent(in) :: ewald_coeff,a(3),b(3),c(3)
  integer,intent(in)  :: grid_size(3),splineorder,natoms
  real(wp),intent(in) :: coords(natoms,3),charges(natoms)
  real(wp),intent(out):: energy,atom_m(natoms)
  real(wp),intent(out),allocatable::forces(:,:)
  logical :: force_switch
  ! Locals
  integer :: nx,ny,nz,i,j,k,atom,ix,iy,iz,shift_x,shift_y,shift_z
  real(wp) :: volume,pi2_a2,q,phi_val
  real(wp) :: recip_a(3),recip_b(3),recip_c(3),kvec(3),pi
  real(wp) :: cell(3,3),inv_cell(3,3),t11,t22
  real(wp),allocatable :: red(:,:),spl(:,:,:),dspl(:,:,:),Qgrid(:,:,:),phi(:,:,:)
  real(wp),allocatable :: splmod_x(:),splmod_y(:),splmod_z(:)
  integer,allocatable :: gridstart(:,:)
  real(wp),allocatable :: gF(:,:,:),m2(:,:,:)
  complex(wp),allocatable :: Sm(:,:,:),conv(:,:,:),phi_c(:,:,:),Qgrid_complex(:,:,:)
  real(wp) :: du(3),fcart(3)

  pi = 3.141592653589793_wp
  nx = grid_size(1); ny = grid_size(2); nz = grid_size(3)

  ! Cell setup
  volume = dot_product(a,cross_product(b,c))
  recip_a = cross_product(b,c)/volume
  recip_b = cross_product(c,a)/volume
  recip_c = cross_product(a,b)/volume
  pi2_a2 = merge(pi**2/ewald_coeff**2,1.0e50_wp,ewald_coeff/=0.0_wp)

  call cpu_time(t11)

  allocate(red(natoms,3),spl(natoms,3,0:splineorder-1),dspl(natoms,3,0:splineorder-1))
  allocate(gridstart(natoms,3),Qgrid(0:nz-1,0:ny-1,0:nx-1)); Qgrid=0.0_wp

  cell(:,1)=a; cell(:,2)=b; cell(:,3)=c
  call invert_3x3_matrix(cell,inv_cell)

  ! Fractional coords + splines
  do atom=1,natoms
     red(atom,:)=matmul(inv_cell,coords(atom,:))
     gridstart(atom,1)=floor(nx*red(atom,1),kind=wp)
     gridstart(atom,2)=floor(ny*red(atom,2),kind=wp)
     gridstart(atom,3)=floor(nz*red(atom,3),kind=wp)
     spl(atom,1,:)=build_spline(splineorder,nx*red(atom,1)-gridstart(atom,1))
     spl(atom,2,:)=build_spline(splineorder,ny*red(atom,2)-gridstart(atom,2))
     spl(atom,3,:)=build_spline(splineorder,nz*red(atom,3)-gridstart(atom,3))
     dspl(atom,1,:)=build_spline_derivative(splineorder,nx*red(atom,1)-gridstart(atom,1))
     dspl(atom,2,:)=build_spline_derivative(splineorder,ny*red(atom,2)-gridstart(atom,2))
     dspl(atom,3,:)=build_spline_derivative(splineorder,nz*red(atom,3)-gridstart(atom,3))
  enddo

  ! Charge spreading
  do atom=1,natoms
     q=charges(atom)
     do i=0,splineorder-1
         do j=0,splineorder-1
                do k=0,splineorder-1
                Qgrid(mod(gridstart(atom,3)+k,nz),mod(gridstart(atom,2)+j,ny),mod(gridstart(atom,1)+i,nx)) = &
                Qgrid(mod(gridstart(atom,3)+k,nz),mod(gridstart(atom,2)+j,ny),mod(gridstart(atom,1)+i,nx)) + &
                q*spl(atom,1,i)*spl(atom,2,j)*spl(atom,3,k)
              enddo
          enddo
      enddo
  enddo

  ! FFT forward
  allocate(Sm(0:nz-1,0:ny-1,0:nx-1))
  allocate(Qgrid_complex(0:nx-1,0:ny-1,0:nz-1))
  Qgrid_complex = cmplx(Qgrid,0.0_wp,wp)
  call FFT_forward_3d(nx,ny,nz,Qgrid_complex,Sm)

  ! Structure factors
  allocate(splmod_x(0:nx-1)); splmod_x = 1.0_wp/dftmod(splineorder,nx)
  allocate(splmod_y(0:ny-1)); splmod_y = 1.0_wp/dftmod(splineorder,ny)
  allocate(splmod_z(0:nz-1)); splmod_z = 1.0_wp/dftmod(splineorder,nz)
  allocate(m2(0:nz-1,0:ny-1,0:nx-1),gF(0:nz-1,0:ny-1,0:nx-1))

  do i=0,nx-1
     shift_x = i; if(i>=nx/2) shift_x=i-nx
     do j=0,ny-1
        shift_y = j; if(j>=ny/2) shift_y=j-ny
        do k=0,nz-1
           shift_z = k; if(k>=nz/2) shift_z=k-nz
           kvec = shift_x*recip_a + shift_y*recip_b + shift_z*recip_c
           m2(k,j,i) = dot_product(kvec,kvec)
           if(i==0.and.j==0.and.k==0) m2(k,j,i)=huge(1.0_wp)
           gF(k,j,i) = splmod_x(i)*splmod_y(j)*splmod_z(k)*exp(-pi2_a2*m2(k,j,i))/m2(k,j,i)
        enddo
     enddo
  enddo

  allocate(conv(0:nz-1,0:ny-1,0:nx-1)); conv=Sm*cmplx(gF,0.0_wp,wp)
  allocate(phi_c(0:nz-1,0:ny-1,0:nx-1))
  call FFT_inverse(nx,ny,nz,conv,phi_c)
  allocate(phi(0:nz-1,0:ny-1,0:nx-1)); phi=real(phi_c,wp)

  ! Energy + forces
  energy=0.0_wp; atom_m=0.0_wp
  allocate(forces(natoms,3)); forces=0.0_wp

  do atom=1,natoms
     q=charges(atom)
     du=0.0_wp
     do i=0,splineorder-1
         do j=0,splineorder-1;
             do k=0,splineorder-1
                ix=mod(gridstart(atom,1)+i,nx)
                iy=mod(gridstart(atom,2)+j,ny)
                iz=mod(gridstart(atom,3)+k,nz)
                phi_val=phi(iz,iy,ix)
        ! fractional spline weights
                if (force_switch) then
                du(1) = du(1) + dspl(atom,1,i)*nx * spl(atom,2,j)*spl(atom,3,k)*phi_val
                du(2) = du(2) + dspl(atom,2,j)*ny * spl(atom,1,i)*spl(atom,3,k)*phi_val
                du(3) = du(3) + dspl(atom,3,k)*nz * spl(atom,1,i)*spl(atom,2,j)*phi_val
                end if
                energy = energy + q*spl(atom,1,i)*spl(atom,2,j)*spl(atom,3,k)*phi_val/(2.0_wp*pi*volume)
                atom_m(atom) = atom_m(atom) + spl(atom,1,i)*spl(atom,2,j)*spl(atom,3,k)*phi_val/(pi*volume)
             enddo
         enddo
     enddo
     ! Transform fractional gradient to Cartesian force
     if (force_switch) then
     fcart = matmul(transpose(inv_cell),du)
     forces(atom,:) = forces(atom,:) - q*fcart/(pi*volume)
     end if
  enddo
  call cpu_time(t22); print*,"SPME time=",t22-t11
end subroutine spme_energy

!============================================================
subroutine sort_vectors_by_distance_merge(trans, dist, n)
  implicit none
  integer, intent(in) :: n
  real(wp), intent(inout) :: trans(3,n), dist(n)
  integer :: i, j
  real(wp) :: tmpd, tmpv(3)

  do i = 2, n
    tmpd = dist(i)
    tmpv = trans(:,i)
    j = i - 1
    do while (j >= 1 .and. dist(j) > tmpd)
      dist(j+1) = dist(j)
      trans(:,j+1) = trans(:,j)
      j = j - 1
    end do
    dist(j+1) = tmpd
    trans(:,j+1) = tmpv
  end do
end subroutine sort_vectors_by_distance_merge

end module PME_main
