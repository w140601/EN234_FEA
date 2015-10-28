!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_hypoelast_3dbasic ==============================
subroutine el_hypoelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D
    implicit none
    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,i,j,flag

    real (prec)  ::  strain(6), dstrain(6),straintotal(6),e(3,3),strainmatrix(3,3),zero(6) ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6) ,stressmatrix(3,3),ee(6)                        ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  volume    ! Material properties
    real (prec)  ::  Es,Et,eij,stressE,strainE,e0,s0,nn,k,dstressE
    real (prec)  ::  volumestrain
    real (prec)  ::  e_dyadic_e(6,6), matrix1(6,6), matrix2(6,6)
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    s0 = element_properties(1)
    e0 = element_properties(2)
    nn=element_properties(3)
    k=element_properties(4)
  !write(*,*) s0
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
       do i=1,n_nodes
          do j=1,3
        dNbardx(i,j)=dNbardx(i,j)/volume
       end do
       end do
    ! dNbardx=dNbardx(1:n_nodes,1:3)/volume
    !write(*,*) volume
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
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
       ! calculate strain
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        straintotal=0.d0
        straintotal(1:6)=strain(1:6)+dstrain(1:6)
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

      if(strainE.le.e0)then
         stressE=(dsqrt((1.d0+nn**2.d0)/(nn-1.d0)**2.d0-(nn/(nn-1.d0)-strainE/e0)**2.d0)-1.d0/(nn-1.d0))*s0
      else
         stressE=((strainE/e0)**(1.d0/nn))*s0
      end if
      !calculate stress vector
      stressmatrix=0.d0
      stress=0.d0
      if(strainE>1.d-010)then
        do i = 1,3
          do j=1,3
          if(i==j)then
           stressmatrix(i,j) = 2.d0/3.d0*stressE*e(i,j)/strainE&
                               +K*(straintotal(1)+straintotal(2)+straintotal(3))
          else
           stressmatrix(i,j)= 2.d0/3.d0*stressE*e(i,j)/strainE
           end if
        end do
        end do
      else
      do i=1,3
        do j=1,3
        if(i==j)then
       stressmatrix(i,j)=0
       end if
       end do
       end do
      end if
write(*,*)'strainE=',strainE
write(*,*)'stress11=',stressmatrix(1,1)
        stress(1)=stressmatrix(1,1)
        stress(2)=stressmatrix(2,2)
        stress(3)=stressmatrix(3,3)
        stress(4)=stressmatrix(1,2)
        stress(5)=stressmatrix(1,3)
        stress(6)=stressmatrix(2,3)
      ! calculate Et and Es
         Es=stressE/strainE
      if(strainE.le.e0) then
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


       element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do
    ! write (*,*) B
    return
end subroutine el_hypoelast_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_hypoelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6),straintotal(6)! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
	
    D = 0.d0
   E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
       ! calculate strain
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      ! calculate D matrix
      straintotal=strain+dstrain
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_hypoelast_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hypoelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx=>vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,i,j,kk,flag

    real (prec)  ::  strain(6), dstrain(6),straintotal(6),e(3,3),strainmatrix(3,3) ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6) ,zero(6),ee(6),stressmatrix(3,3)                        ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  volume    ! Material properties
    real (prec)  ::  Es,Et,eij,stressE,strainE,e0,s0,nn,k,dstressE
    real (prec)  ::  volumestrain
    real (prec)  :: e_dyadic_e(6,6),matrix1(6,6),matrix2(6,6)
    real (prec)  :: p, smises                         ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))
write(*,*)1
    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
    D = 0.d0
    s0 = element_properties(1)
    e0 = element_properties(2)
    nn=element_properties(3)
    k=element_properties(4)
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
       do i=1,n_nodes
          do j=1,3
        dNbardx(i,j)=dNbardx(i,j)/volume
       end do
       end do

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
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
       ! calculate strain
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
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
!write(*,*)'strainE=',strainE
!write(*,*)'stress11=',stress(1)
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
    !    stress = matmul(D,strain+dstrain)
      !calculate stress vector
!      do i = 1,3
!           stress(i) = 2.d0/3.d0*stressE*e(i)/strainE+K*(straintotal(1)+straintotal(2)+straintotal(3))
!        end do
!      do i = 4,6
!           stress(i) = 2.d0/3.d0*stressE*e(i)/strainE
!      end do
      ! calculate Et and Es


        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do kk = 1,n_field_variables
            if (strcmp(field_variable_names(kk),'S11',3) ) then
                nodal_fieldvariables(kk,1:n_nodes) = nodal_fieldvariables(kk,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(kk),'S22',3) ) then
                nodal_fieldvariables(kk,1:n_nodes) = nodal_fieldvariables(kk,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(kk),'S33',3) ) then
                nodal_fieldvariables(kk,1:n_nodes) = nodal_fieldvariables(kk,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(kk),'S12',3) ) then
                nodal_fieldvariables(kk,1:n_nodes) = nodal_fieldvariables(kk,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(kk),'S13',3) ) then
                nodal_fieldvariables(kk,1:n_nodes) = nodal_fieldvariables(kk,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(kk),'S23',3) ) then
                nodal_fieldvariables(kk,1:n_nodes) = nodal_fieldvariables(kk,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(kk),'SMISES',6) ) then
                nodal_fieldvariables(kk,1:n_nodes) = nodal_fieldvariables(kk,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_hypoelast_3dbasic

