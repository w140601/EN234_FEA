!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_hyperelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
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
    integer      :: n_points,kint,i,j,k,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6),stressmatrix(3,3)       ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant,detf           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12 ,volume      ! Material properties
    real (prec)  :: F(3,3),finv(3,3),doff(length_dof_array/3,3),II(6),BB(6),BBmatrix(3,3)
    real (prec)  :: matrix1(6,6),nv,bulk,detbb,G(6,9),Bstar(9,3*n_nodes)
    real (prec)  :: bbmatrixinv(3,3),detb,ii_dyadic_bbinv(6,6),bbinv(6),ii_dyadic_ii(6,6),bb_dyadic_bbinv(6,6)
    real (prec)  :: s(3,length_dof_array/3),Pvec(3*n_nodes),pmat(3*n_nodes,3*n_nodes),svec(3*n_nodes)
    real (prec)  :: smat(3*n_nodes,3*n_nodes),Sigma(3*n_nodes,3*n_nodes)
!    real (prec)  :: dNdy(n_nodes,3)
    !
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
    nv = element_properties(1)
    bulk = element_properties(2)
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        doff=transpose(reshape((dof_total+dof_increment),(/3,length_dof_array/3/)))
  ! write(*,*)doff
     F=0.d0
     F(1,1)=1.d0
     F(2,2)=1.d0
     F(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F(i,j)=F(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+dof_increment(3*(a-1)+i))
                end do
            end do
        end do
        write (*,*) 'f(1,1)=',f(1,1)
! calculate inverse F and J
        call invert_small(F,finv,detf)

!write(*,*) dndx
!        do a=1,n_nodes
!           do i=1,3
!              do k=1,3
!              dNdy(a,i)=dNdy(a,i)+finv(k,i)*dNdx(a,k)
!              end do
!           end do
!        end do
      dNdy=matmul(dNdx,finv)

!write(*,*)dNdy

        B=0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        II=[1.d0,1.d0,1.d0,0.d0,0.d0,0.d0]

       BBmatrix=matmul(F,transpose(F))
!        bbmatrix=0.d0
!        do i=1,3
!            do j=1,3
!                do a=1,3
!                    BBmatrix(i,j)=BBmatrix(i,j)+F(i,a)*F(j,a)
!                end do
!            end do
!        end do
         bb=0.d0
        BB(1)=BBmatrix(1,1)
        BB(2)=BBmatrix(2,2)
        BB(3)=BBmatrix(3,3)
        BB(4)=BBmatrix(1,2)
        BB(5)=BBmatrix(1,3)
        BB(6)=BBmatrix(2,3)
! calculate inverse left cauchy green tensor
      call invert_small(BBmatrix,BBmatrixinv,detBB)
!      do i=1,6
!          bbinv(i)=1.d0/bb(i)
!          end do
        bbinv=0.d0
        BBinv(1)=bbmatrixinv(1,1)
        BBinv(2)=BBmatrixinv(2,2)
        BBinv(3)=BBmatrixinv(3,3)
        BBinv(4)=BBmatrixinv(1,2)
        BBinv(5)=BBmatrixinv(1,3)
        BBinv(6)=BBmatrixinv(2,3)
write (*,*) 'bbinv = ' , BBmatrixinv(1,1)
!
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        ! calculate D
        ! calculate matrix1
        matrix1=0.d0
        do i=1,3
            matrix1(i,i)=1.d0
            matrix1(i+3,i+3)=0.5d0
        end do

        II_dyadic_BBinv=spread(II,dim=2,ncopies=6)*spread(BBinv,dim=1,ncopies=6)
        II_dyadic_II=spread(II,dim=2,ncopies=6)*spread(II,dim=1,ncopies=6)
        BB_dyadic_BBinv=spread(BB,dim=2,ncopies=6)*spread(BBinv,dim=1,ncopies=6)


        D=nv/(detf**(2.d0/3.d0))*matrix1+nv/(3.d0*detf**(2.d0/3.d0))*((BB(1)+BB(2)+BB(3))/3.d0*II_dyadic_BBinv&
           -II_dyadic_II-BB_dyadic_BBinv)+bulk*detf*(detf-0.5d0)*II_dyadic_BBinv
write (*,*) 'detf',detf
write(*,*)'BB: ', BB(1)
write(*,*)'II_dya_bbinv: ', II_dyadic_BBinv
write(*,*)'II_dyadic_II: ', II_dyadic_II
write(*,*)'BB_dyadic_BBinv: ', BB_dyadic_BBinv
        !stress = matmul(D,strain+dstrain)
!write(*,*)bbmatrix
!write(*,*)0000000000000000000000
        G=0.d0
        G(1,1:9)=[2.d0*BBmatrix(1,1),0.d0,0.d0,2.d0*BBmatrix(1,2),0.d0,2.d0*BBmatrix(1,3),0.d0,0.d0,0.d0]
        G(2,1:9)=[0.d0,2.d0*BBmatrix(2,2),0.d0,0.d0,2.d0*BBmatrix(1,2),0.d0,0.d0,2.d0*BBmatrix(2,3),0.d0]
        G(3,1:9)=[0.d0,0.d0,2.d0*BBmatrix(3,3),0.d0,0.d0,0.d0,2.d0*BBmatrix(1,3),0.d0,2.d0*BBmatrix(1,3)]
        G(4,1:9)=[2.d0*BBmatrix(1,2),2.d0*BBmatrix(1,2),0.d0,2.d0*BBmatrix(2,2),2.d0*BBmatrix(1,1),&
            2.d0*bbmatrix(2,3),0.d0,2.d0*BBmatrix(1,3),0.d0]
        G(5,1:9)=[2.d0*BBmatrix(1,3),0.d0,2.d0*BBmatrix(1,3),2.d0*BBmatrix(2,3),0.d0,2.d0*BBmatrix(3,3)&
            ,2.d0*BBmatrix(1,1),0.d0,2.d0*BBmatrix(1,2)]
        G(6,1:9)=[0.d0,2.d0*BBmatrix(2,3),2.d0*BBmatrix(2,3),0.d0,2.d0*BBmatrix(1,3),0.d0,2.d0*BBmatrix(1,2),&
            2.d0*BBmatrix(3,3),2.d0*BBmatrix(2,2)]

write(*,*)g(6,9)
           !stressmatrix=0.d0
            do i=1,3
                do j=1,3
                    if(i==j)then
                        stressmatrix(i,j)=nv/(detf**(5.d0/3.d0))*(BBmatrix(i,j)-1.d0/3.d0*(BBmatrix(1,1)+BBmatrix(2,2)&
                            +BBmatrix(3,3)))+bulk*(detf-1.d0)
                    else
                        stressmatrix(i,j)=nv/(detf**(5.d0/3.d0))*BBmatrix(i,j)
                    end if
                end do
            end do
            stressmatrix=stressmatrix*detf
!write(*,*) stressmatrix
            stress=0.d0
            stress(1)=stressmatrix(1,1)
            stress(2)=stressmatrix(2,2)
            stress(3)=stressmatrix(3,3)
            stress(4)=stressmatrix(1,2)
            stress(5)=stressmatrix(1,3)
            stress(6)=stressmatrix(2,3)
                  bstar=0.d0
                  Bstar(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
                  Bstar(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
                  Bstar(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
                  Bstar(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
                  Bstar(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
                  Bstar(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
                  Bstar(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
                  Bstar(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
                  Bstar(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
write(*,*) 'bstar= ',Bstar(1,1)
write(*,*) 'stressmatrix=',stressmatrix
                  S = reshape(matmul(transpose(B),stress),(/3,length_dof_array/3/))
                  do i = 1,n_nodes
                   Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
                   Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
                   Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
                   Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
                  end do
                  Sigma = Pmat*transpose(Smat)
!write(*,*)sigma(1,1)
!write(*,*)G(1,1)
!write(*,*)D(1,1)
                  element_residual = element_residual-matmul(transpose(B),stress)&
                      *w(kint)*determinant

                  element_stiffness = element_stiffness+ (matmul(transpose(B),matmul(D,matmul(G,&
                      bstar)))-sigma)*w(kint)*determinant

    end do
    ! write (*,*) B
    return
end subroutine el_hyperelast_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_hyperelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
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

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_hyperelast_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hyperelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
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
  
    integer      :: n_points,kint,k,i,j,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6),stressmatrix(3,3)       ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises ,volume                        ! Pressure and Mises stress
    real (prec)  :: F(3,3),finv(3,3),doff(length_dof_array/3,3),II(6),BB(6),BBmatrix(3,3),detf
    real (prec)  :: nv,bulk
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
    D = 0.d0
    nv = element_properties(1)
    bulk = element_properties(2)

    !     --  Loop over integration points
   do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        doff=transpose(reshape((dof_total+dof_increment),(/3,length_dof_array/3/)))

     F=0.d0
     F(1,1)=1.d0
     F(2,2)=1.d0
     F(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F(i,j)=F(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+dof_increment(3*(a-1)+i))
                end do
            end do
        end do

! calculate inverse F and J
        call invert_small(F,finv,detf)
!        do a=1,n_nodes
!           do i=1,3
!              do k=1,3
!              dNdy(a,i)=dNdy(a,i)+finv(k,i)*dNdx(a,k)
!              end do
!           end do
!        end do
    dNdy=matmul(dNdx,finv)

        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        II=[1.d0,1.d0,1.d0,0.d0,0.d0,0.d0]
        BBmatrix=matmul(F,transpose(F))

! calculate inverse left cauchy green tensor


        !
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        ! calculate D
        ! calculate matrix1

            do i=1,3
                do j=1,3
                    if(i==j)then
        stressmatrix(i,j)=nv/(detf**(5.d0/3.d0))*(BBmatrix(i,j)-1.d0/3.d0*(BBmatrix(1,1)+BBmatrix(2,2)&
                         +BBmatrix(3,3)))+bulk*(detf-1.d0)
                    else
                        stressmatrix(i,j)=nv/(detf**(5.d0/3.d0))*BBmatrix(i,j)
                    end if
                end do
            end do
         !stressmatrix=stressmatrix*detf

            stress(1)=stressmatrix(1,1)
            stress(2)=stressmatrix(2,2)
            stress(3)=stressmatrix(3,3)
            stress(4)=stressmatrix(1,2)
            stress(5)=stressmatrix(1,3)
            stress(6)=stressmatrix(2,3)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_hyperelast_3dbasic

