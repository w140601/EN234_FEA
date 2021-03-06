subroutine user_print(n_steps)
!  use Types
!  use ParamIO
!  use Globals, only : TIME, DTIME
!!  use Mesh
!  use Printparameters, only : n_user_print_files                  ! No. files specified by the user
!  use Printparameters, only : n_user_print_parameters             ! No. user supplied print parameters
!  use Printparameters, only : user_print_units                    ! Unit numbers
!  use Printparameters, only : user_print_parameters
!  !use Printparameters, only : J_integral_value               ! List of user supplied parameters
!  use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element
!  implicit none
!
!  integer, intent(in) :: n_steps                                 ! Current step number
!
!  integer ::  lmn
!  integer ::  status
!  integer ::  n_state_vars_per_intpt                                         ! No. state variables per integration point
!  real (prec) ::   vol_averaged_strain(6)                                    ! Volume averaged strain in an element
!  real (prec), allocatable ::   vol_averaged_state_variables(:)              ! Volume averaged state variables in an element
!  real (prec),allocatable :: J_integral_value
!
!
!!
!!  Use this file to process or print time histories of the solution, or to print a non-standard mesh.
!!
!!  As an example, this subroutine computes the volume averaged infinitesimal strain and the volume average
!!  element state variables (if they exist) in an element.   The element is specified by the user.
!!  The first six state variables (which are usually the stresses) are printed along with the strains.
!!
!!
!
!   allocate(J_integral_value, stat=status)
!
!   if (status/=0) then
!      write(IOW,*) ' Error in subroutine user_print'
!      write(IOW,*) ' Unable to allocate memory for state variables '
!      stop
!   endif
!
!  !lmn = int(user_print_parameters(1))     ! The element number
!
!   !call compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_variables,length_state_variable_array, &
!                                          !             n_state_vars_per_intpt)
!   call compute_J_integral(J_integral_value)
!   write(user_print_units(1),'(A)')'J_integral_value'
!   write(user_print_units(1),'(1(d12.5))')J_integral_value
!   write(*,*) J_integral_value
!
!!    if (TIME<1.d-12) then
!!      if (n_state_vars_per_intpt<6) then
!!        write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23'
!!      else
!!         write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23'
!!      endif
!!    endif
!!
!!   if (n_state_vars_per_intpt<6) then
!!      write(user_print_units(1),'(7(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6)
!!   else
!!      vol_averaged_state_variables(1:3) = vol_averaged_state_variables(1:3) + vol_averaged_state_variables(7)
!!      write(user_print_units(1),'(13(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6),vol_averaged_state_variables(1:6)
!!   endif
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
!  use Mesh
  use Printparameters, only : n_user_print_files                  ! No. files specified by the user
  use Printparameters, only : n_user_print_parameters             ! No. user supplied print parameters
  use Printparameters, only : user_print_units                    ! Unit numbers
  use Printparameters, only : user_print_parameters               ! List of user supplied parameters
  use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element
  implicit none
  
  integer, intent(in) :: n_steps                                 ! Current step number
  
  integer ::  lmn
  integer ::  status
  integer ::  n_state_vars_per_intpt                                         ! No. state variables per integration point
  real (prec) ::   vol_averaged_strain(6)                                    ! Volume averaged strain in an element
  real (prec) ::   vol_averaged_stress(6),results(2)                                    ! Volume averaged stress in an element
  real (prec), allocatable ::   vol_averaged_state_variables(:)              ! Volume averaged state variables in an element



!
!  Use this file to process or print time histories of the solution, or to print a non-standard mesh.
!
!  As an example, this subroutine computes the volume averaged infinitesimal strain and the volume average
!  element state variables (if they exist) in an element.   The element is specified by the user.
!  The first six state variables (which are usually the stresses) are printed along with the strains.
!
!

  allocate(vol_averaged_state_variables(length_state_variable_array), stat=status)

   if (status/=0) then
      write(IOW,*) ' Error in subroutine user_print'
      write(IOW,*) ' Unable to allocate memory for state variables '
      stop
   endif

   lmn = int(user_print_parameters(1))     ! The element number

   call compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_variables,length_state_variable_array, &
                                                       n_state_vars_per_intpt,vol_averaged_stress)
   results=0.d0
  results(1)=vol_averaged_strain(1)
  results(2)=vol_averaged_stress(1)

    if (TIME<1.d-12) then
      if (n_state_vars_per_intpt<6) then
        write(user_print_units(1),'(A)') 'VARIABLES = e11, s11'
      else
         write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23'
      endif
    endif

   if (n_state_vars_per_intpt<6) then
      write(user_print_units(1),'(2(1x,D12.5))')vol_averaged_strain(1),vol_averaged_stress(1)

   else
      vol_averaged_state_variables(1:3) = vol_averaged_state_variables(1:3) + vol_averaged_state_variables(7)
      write(user_print_units(1),'(13(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6),vol_averaged_state_variables(1:6)
   endif


end subroutine user_print

subroutine compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_vars,length_output_array, &
    n_state_vars_per_intpt,vol_averaged_stress)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent ( in )      :: length_output_array

    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent (out)   ::  vol_averaged_stress(6)
    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)

    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_deltaifier                              ! Flag deltaifying node type
    integer    :: element_deltaifier                           ! Flag deltaifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i,j,flag
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol,volume
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  strain(6), straintotal(6),ee(6),e(3,3),strainmatrix(3,3),zero(6)                       ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  volumestrain                      ! volumetric strain
    real (prec)  ::  Es,Et,eij,stressE,strainE,e0,s0,nn,k,dstressE
    real (prec)  :: e_dyadic_e(6,6),matrix1(6,6),matrix2(6,6)
    real (prec)  :: stress(6),stressmatrix(3,3),D(6,6)
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)

    if (status/=0) then
        write(IOW,*) ' Error in subroutine compute_volume_average_3D'
        write(IOW,*) ' Unable to allocate memory for element variables '
        stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_deltaifier,n_nodes,node_list,n_properties,element_properties, &
        n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_deltaifier,n_coords,x(1:3,i),n_dof, &
            dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    vol_averaged_strain = 0.d0
    vol_averaged_state_vars = 0.d0
    el_vol = 0.d0
    n_state_vars_per_intpt = n_state_variables/n_points
    vol_averaged_stress=0.d0
    straintotal=0.d0
    s0 = element_properties(1)
    e0 = element_properties(2)
    nn=element_properties(3)
    k=element_properties(4)

    if (n_state_vars_per_intpt>size(vol_averaged_state_vars)) then
        write(IOW,*) ' Error detected in subroutine compute_element_volume_average_3d '
        write(IOW,*) ' The element contains ',n_state_vars_per_intpt
        write(IOW,*) ' but the array storing averaged state variables has length ',size(vol_averaged_state_vars)
        stop
    endif
    !     --  Loop over integration points
    volume=0.d0
    dNbardx=0.d0
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        do i=1,n_nodes
            do j=1,3
                dNbardx(i,j)=dNbardx(i,j)+dNdx(i,j)*w(kint)*determinant
            end do
        end do
        volume=volume+w(kint)*determinant
    end do
    dNbardx=dNbardx/volume


    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)

        iof = n_state_vars_per_intpt*(kint-1)+1
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)+1.d0/3.d0*(dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))
        B(1,2:3*n_nodes-1:3)=1.d0/3.d0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        B(1,3:3*n_nodes:3)=1.d0/3.d0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        B(2,1:3*n_nodes-2:3)=1.d0/3.d0*(dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)+1.d0/3.d0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        B(2,3:3*n_nodes:3)=1.d0/3.d0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        B(3,1:3*n_nodes-2:3)=1.d0/3.d0*(dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))
        B(3,2:3*n_nodes-1:3)=1.d0/3.d0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)+1.d0/3.d0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B(1:6,1:3*n_nodes),dof_total(1:3*n_nodes))
        dstrain = matmul(B(1:6,1:3*n_nodes),dof_increment(1:3*n_nodes))

        vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + (strain(1:6)+dstrain(1:6))*w(kint)*determinant

        if (n_state_vars_per_intpt>0) then
            vol_averaged_state_vars(1:n_state_vars_per_intpt) = vol_averaged_state_vars(1:n_state_vars_per_intpt) &
                + updated_state_variables(iof:iof+n_state_vars_per_intpt-1)*w(kint)*determinant
        endif

        el_vol = el_vol + w(kint)*determinant

        !calculate stress
        straintotal=strain+dstrain
        strainmatrix(1,1)=straintotal(1)
        strainmatrix(2,2)=straintotal(2)
        strainmatrix(3,3)=straintotal(3)
        strainmatrix(1,2)=straintotal(4)/2.d0
        strainmatrix(2,1)=straintotal(4)/2.d0
        strainmatrix(1,3)=straintotal(5)/2.d0
        strainmatrix(3,1)=straintotal(5)/2.d0
        strainmatrix(2,3)=straintotal(6)/2.d0
        strainmatrix(3,2)=straintotal(6)/2.d0
        !write(*,*)straintotal

        ! calculate D matrix
        volumestrain=0.d0
        volumestrain=straintotal(1)+straintotal(2)+straintotal(3)
        e=0.d0
        ! calculate devitoric strain
        do i=1,3
            do j=1,3
                if(i==j)then
                    e(i,j)=strainmatrix(i,j)-1.d0/3.d0*volumestrain
                else
                    e(i,j)=strainmatrix(i,j)
                end if
            end do
        end do

        eij=0.d0
        do i=1,3
            do j=1,3
                eij=eij+e(i,j)*e(i,j)
            end do
        end do
        ! calculate stressE
        strainE=dsqrt((2.d0/3.d0)*eij)

        if(strainE<=e0)then
            stressE=(dsqrt((1.d0+nn**2.d0)/(nn-1.d0)**2.d0-(nn/(nn-1.d0)-strainE/e0)**2.d0)-1.d0/(nn-1.d0))*s0
        else
            stressE=((strainE/e0)**(1.d0/nn))*s0
        end if
        !calculate stress vector
        stressmatrix=0.d0
        stress=0.d0
        if (strainE>1.d-010)then
        do i = 1,3
            do j=1,3
                if(i==j)then
                    stressmatrix(i,j) = 2.d0/3.d0*stressE*e(i,j)/strainE+K*(straintotal(1)+straintotal(2)+straintotal(3))
                else
                    stressmatrix(i,j)= 2.d0/3.d0*stressE*e(i,j)/strainE
                end if
            end do
        end do
        else
         do i=1,3
        do j=1,3
        if(i==j)then
       stressmatrix(i,j)=K*(straintotal(1)+straintotal(2)+straintotal(3))
       end if
       end do
       end do
        end if
        stress(1)=stressmatrix(1,1)
        stress(2)=stressmatrix(2,2)
        stress(3)=stressmatrix(3,3)
        stress(4)=stressmatrix(1,2)
        stress(5)=stressmatrix(1,3)
        stress(6)=stressmatrix(2,3)
        ! calculate Et and Es
        Es=stressE/strainE
        if(strainE<=e0) then
            Et=(2.d0*dsqrt(nn/(nn-1.d0)-strainE/e0)/e0)*s0
        else
            Et=(1.d0/nn*((strainE/e0)**(1.d0/nn-1))/e0)*s0
        end if
        ! calculate total deviatoric strain
        ee(1)=e(1,1)
        ee(2)=e(2,2)
        ee(3)=e(3,3)
        ee(4)=e(1,2)
        ee(5)=e(1,3)
        ee(6)=e(2,3)

        e_dyadic_e=0.d0
        e_dyadic_e = spread(ee,dim=2,ncopies=6)*spread(ee,dim=1,ncopies=6)
        matrix1=0.d0
        matrix2=0.d0
        do i=1,3
            matrix1(i,i)=2.d0
            matrix1(i+3,i+3)=1.d0
        end do
        do i=1,3
            do j=1,3
                matrix2(i,j)=1.d0
            end do
        end do
        D=0.d0
        zero=0.d0
        flag=1
        do i=1,6
            if(straintotal(i)/=zero(i))then
                flag=0
            end if
        end do
        ! write(*,*)flag
        if(flag==1)then
            D(1:6,1:6)=Et/3.d0*matrix1(1:6,1:6)+(K-2.d0*Et/9.d0)*matrix2(1:6,1:6)
        else
            D(1:6,1:6)=4.d0/(9.d0*strainE**2.d0)*(Et-Es)*e_dyadic_e(1:6,1:6)+Es/3.d0*matrix1(1:6,1:6)&
                +(K-2.d0*Es/9.d0)*matrix2(1:6,1:6)
        end if
        !stress = matmul(D,strain+dstrain)
        vol_averaged_stress(1:6) = vol_averaged_stress(1:6) + stress(1:6)*w(kint)*determinant



    end do

    vol_averaged_strain(1:6) = vol_averaged_strain(1:6)/el_vol
    vol_averaged_state_vars = vol_averaged_state_vars/el_vol
    vol_averaged_stress(1:6) = vol_averaged_stress(1:6)/volume


    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)

    return




end subroutine compute_element_volume_average_3D

subroutine compute_J_integral(J_integral_value)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use Mesh, only : zone,zone_list
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    real (prec), intent( out )  ::  J_integral_value


    ! Local variables
    integer      ::  n_points,kint
    !real (prec)  ::  B(3 ,length_dof_array)
    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  stressmatrix(2,2)                     ! Stress matrix
    real (prec)  ::  strainmatrix(2,2)                     ! Strain matrix
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  E, xnu, D44, D11, D12             ! Material properties
    real (prec)  ::  delta(2,2)                         ! deltaity
    real (prec)  :: NN(2)
    real (prec)  ::  r0, r
    real (prec)  :: sx(2)
    real (prec)  :: ww                                  ! energy

    integer    :: node_deltaifier                              ! Flag deltaifying node type
    integer    :: element_deltaifier                           ! Flag deltaifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element
    integer    :: n_coords                                     ! No. coords for a node
    integer    :: n_dof                                        ! No. DOFs for a node

    integer      :: status
    integer      :: iof
    integer      :: lmn               ! Element number
    integer      :: lmn_start,lmn_end ! First and last crack tip element
    integer      ::  i, i_count , j_count, n_count                ! Loop counter

!   The arrays below have to be given dimensions large enough to store the data. It doesnt matter if they are too large.

    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)
    !real (prec),  allocatable   :: sx(:)
    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    real (prec), allocatable  ::  B(:,:)                                   ! strain = B*(dof_total+dof_increment)
    !
    !
    !  The variables specifying the sizes of the arrays listed below are determined while reading the input file
    !  They specify the array dimensions required to store the relevant variables for _any_ element or node in the mesh
    !  The actual data will vary depending on the element or node selected
    !
    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(2,length_coord_array/2), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(3,length_dof_array), stat=status)
    !allocate(sx(2),stat=status)
  !  Write your code to calculate the J integral here

  !  You will need to loop over the crack tip elements, and sum the contribution to the J integral from each element.
  !
  !  You can access the first and last crack tip element using

  lmn_start = zone_list(2)%start_element
  lmn_end = zone_list(2)%end_element
write(*,*) lmn_start, lmn_end
     J_integral_value = 0.d0
     delta = 0.d0
     delta(1:2,1:2)=0.d0
     delta(1,1) = 1.d0
     delta(2,2) = 1.d0


  !  The two subroutines below extract data for elements and nodes (see module Mesh.f90 for the source code for these subroutines)
do lmn=lmn_start, lmn_end
    call extract_element_data(lmn,element_deltaifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    if (n_nodes == 6) n_points = 9
    if (n_nodes == 8) n_points = 9
    call initialize_integration_points(n_points, n_nodes, xi, w)

    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1.D0+xnu)
    d11 = (1.D0-xnu)*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
    d12 = xnu*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
    D(1,2) = d12
    D(2,1) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d44


    do i = 1, n_nodes
        iof = 2*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_deltaifier,n_coords,x(1:2,i),n_dof, &
                                                 dof_increment(iof:iof+1),dof_total(iof:iof+1))
    end do
do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)

        B = 0.d0
        B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(3,2:2*n_nodes:2) = dNdx(1:n_nodes,1)


        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        strain = strain + dstrain
        stress = matmul(D,strain)

        stressmatrix(1,1) = stress(1)
        stressmatrix(2,2) = stress(2)
        stressmatrix(1,2) = stress(3)
        stressmatrix(2,1) = stress(3)
        strainmatrix(1,1) = strain(1)
        strainmatrix(2,2) = strain(2)
        strainmatrix(1,2) = strain(3)/2.d0
        strainmatrix(2,1) = strain(3)/2.d0
        sx = 0.d0

        do i = 1, n_nodes
        sx(1) = sx(1)+ N(i)*x(1,i)
        sx(2)=sx(2)+N(i)*x(2,i)

        end do
        r0 = 0.0006d0
        r = dsqrt(dot_product(sx(1:2),sx(1:2)))
        ww=0.d0
       do i_count=1,2
       do j_count=1,2
          ww=ww+stressmatrix(i_count,j_count)*strainmatrix(i_count,j_count)/2.d0
       end do
       end do
       NN=0.d0
       do n_count =1,n_nodes
             NN(1:2)=NN(1:2)+dNdx(n_count,2)*(dof_total(2*n_count-1:2*n_count)+dof_increment(2*n_count-1:2*n_count))
       end do
       do i_count = 1,2
         do j_count = 1,2

!             if(r.gt.r0) then
!             J_integral_value=J_integral_value+0
!             else
             J_integral_value = J_integral_value +(stressmatrix(i_count,j_count)*NN(i_count)-0.5*ww*delta(j_count,2))&
                                *(-sx(j_count)/(r0*r))*w(kint)*determinant
              !end if
           end do
         end do
       end do


end do

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)


    return




end subroutine compute_J_integral

