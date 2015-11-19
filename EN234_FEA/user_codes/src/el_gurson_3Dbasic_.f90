!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================

 subroutine gurson(element_properties,n_properties,n_state_variables,initial_state_variables, &
updated_state_variables,dstrain,dRot,stress1)
     use Types
     use ParamIO
     use Globals, only : TIME, DTIME
     use Element_Utilities, only : invert_small

     implicit none

     integer, intent( in ) :: n_properties
     integer, intent( in ) :: n_state_variables

     real (prec), intent( in ) :: element_properties(n_properties)
     real (prec), intent( in ) :: initial_state_variables(n_state_variables)
     real (prec), intent( in ) :: dstrain(3,3)
     real (prec), intent( in ) :: dRot(3,3)
     real (prec), intent( out ) :: stress1(6)
     real (prec), intent( out ) :: updated_state_variables(n_state_variables)
    integer :: i,j,k,l
    real (prec)  :: stress0(6), ematrix, Vf,phisq,eps,pp,sigma,yy
    real (prec)  :: taomat(3,3), taomatdiv(3,3), dstraindiv(3,3), pstar, Sstar(3,3), phi, sigmastar, fstar
    real (prec)  :: ffbar, tol, taomat1(3,3),e,xnu,e0dot,mm,q1,q2,q3,fn,een,ssn,fc,ff,dee,dev
    real (prec)  :: ddeevv(2),dphidse,dphidp,dphiddee,dphiddev,dtempddee,dtempddev,se,matrix(2,2)
    real (prec)  :: matrix_inv(2,2),df1ddee,df1ddev,df2ddee,df2ddev,dpddev,dseddee,dt,f1,f2,line(2)
    real (prec)  :: matrix_det,p,temp,error,vf1,ematrix1,dematrix

    E = element_properties(1)
    xnu = element_properties(2)
    YY=element_properties(3)
    e0dot=element_properties(4)
    mm=element_properties(5)
    q1=element_properties(6)
    q2=element_properties(7)
    q3=element_properties(8)
    fn=element_properties(9)
    een=element_properties(10)
    ssn=element_properties(11)
    fc=element_properties(12)
    ff=element_properties(13)

    stress0 = initial_state_variables(1:6) ! Stress at start of increment
    ematrix = initial_state_variables(7)
    Vf = initial_state_variables(8)

! define 3*3 matrix taomat
    taomat = 0.d0
    taomat(1,1) = stress0(1)
    taomat(2,2) = stress0(2)
    taomat(3,3) = stress0(3)
    taomat(1,2) = stress0(4)
    taomat(2,1) = stress0(4)
    taomat(1,3) = stress0(5)
    taomat(3,1) = stress0(5)
    taomat(2,3) = stress0(6)
    taomat(3,2) = stress0(6)

! compute taomatdiv, dstraindiv, pstar and Sij star
    taomatdiv = 0.d0
    do i = 1,3
        do j = 1,3
        if(i==j)then
        taomatdiv(i,j) = taomat(i,j)-(taomat(1,1)+taomat(2,2)+taomat(3,3))/3.d0
        else
        taomatdiv(i,j)=taomat(i,j)
        end if
        end do
    end do

    dstraindiv= 0.d0
    do i = 1,3
        do j = 1,3
            if (i == j) then
            dstraindiv(i,j) = dstrain(i,j)-(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.d0
            else
            dstraindiv(i,j) = dstrain(i,j)
            end if
        end do
    end do

    pstar = (taomat(1,1)+taomat(2,2)+taomat(3,3))/3.d0+E/(3.d0*(1.d0-2.d0*xnu))* &
                     (dstrain(1,1)+ dstrain(2,2)+ dstrain(3,3))

   Sstar = E/(1.d0+xnu)*dstraindiv + matmul(dRot,matmul(taomatdiv,transpose(dRot)))
   !write(*,*)'drot',drot
  !  write(*,*)'sstar',sstar

! define fstar
    ffbar = (q1+dsqrt(q1**2.d0-q3))/q3

    if (Vf < fc) then
    fstar = Vf
    else if (fc < Vf .AND. Vf < ff) then
    fstar = fc + (ffbar-fc)/(ff-fc)*(Vf-fc)
    else if (Vf > ff) then
    fstar = ffbar
    end if

! define sigmastar and phi

    sigmastar = dsqrt(3.d0/2.d0*(Sstar(1,1)**2.d0+Sstar(2,2)**2.d0+Sstar(3,3)**2.d0+ &
                       2.d0*Sstar(1,2)**2.d0+2.d0*Sstar(2,3)**2.d0+2.d0*Sstar(1,3)**2.d0))
  !write(*,*) 's1',sigmastar
    Phisq = sigmastar**2.d0/YY**2.d0+2.d0*q1*fstar*cosh(3.d0/2.d0*q2*pstar/YY)-(1.d0+q3*fstar**2.d0)

! Determine if Phi is greater or smaller than tolerance
   eps = 10.d0**(-8.d0)
   taomat1=0.d0
    if (Phisq < eps) then
        do i = 1,3
            do j = 1,3
            if (i == j) then
            taomat1(i,j) = Sstar(i,j)+pstar
            else
            taomat1(i,j) = Sstar(i,j)
            end if
            end do
         end do
       dematrix=0.d0
       ematrix1=ematrix
    else

!write(*,*) ematrix1

        dee = 0.d0
        dev = 0.d0
        ddeevv=1.d0
        error=1.d0
        do while(error>eps)
        dt=dtime
        !write(*,*)dt
        se = sigmastar - 1.5d0*E/(1.d0+xnu)*dee
   !     write(*,*)'s2',sigmastar
        p = pstar - E/3.d0/(1.d0-2.d0*xnu)*dev
        phi = se**2.d0/YY**2.d0 + 2.d0*q1*fstar*cosh(1.5d0*q2*p/Yy) &
             - (1.d0+q3*fstar**2.d0)
        if (phi<0.d0) then
                dee=dee/10.d0
                dev=dev/10.d0
                cycle
            endif
        phi=dsqrt(phi)
        dphidse = se / phi / Yy**2.d0
        dphidp = 1.5d0 * q1 * q2 * fstar / phi / Yy * sinh(1.5d0 * q2 * p / Yy)
        dseddee = -1.5d0*e/(1+xnu)
        dpddev= -e/(3.d0*(1-2.d0*xnu))
        dphiddee = dphidse * dseddee
        dphiddev = dphidp * dpddev
        temp = dsqrt(dphidse**2.d0 + 2.d0/9.d0*dphidp**2.d0)
        dtempddee = (2.d0*se*dseddee/phi/phi/Yy**4.d0 - se**2.d0/Yy**4.d0*2.d0/phi**3.d0*dphiddee &
                    + 2.d0/9.d0*dphidp**2.d0*(-2.d0)/phi*dphiddee)/2.d0/temp
        dtempddev = (q1**2.d0*q2**2.d0*fstar**2.d0/phi**2.d0/Yy**2.d0*sinh(1.5*q2/Yy*p)*cosh(1.5*q2/Yy*p) &
                    * 1.5*q2/Yy*dpddev + 2.d0/9.d0*dphidp**2.d0*(-2.d0)/phi*dphiddev-se**2.d0/yy**4.d0*2.d0&
                    /phi**3.d0*dphiddev)/2.d0/temp
        F1 = temp * dee /dt/e0dot - dphidse * phi**mm
        F2 = temp * dev /dt/e0dot - dphidp * phi**mm

        dF1ddee = temp/dt/e0dot - dseddee / Yy**2.d0 * phi**(mm-1.d0) - &
                    se/Yy**2.d0*(mm-1.d0)*phi**(mm-2.d0)*dphiddee + dee/dt/e0dot*dtempddee
        dF1ddev = - se/Yy**2.d0*(mm-1.d0)*phi**(mm-2.d0)*dphiddev + dee/dt/e0dot*dtempddev
        dF2ddee = dev/dt/e0dot * dtempddee - dphidp*(phi**mm)*(mm-1.d0)/(phi)*dphiddee
        dF2ddev = temp/dt/e0dot - dphidp*(phi**mm)*(mm-1.d0)/(phi)*dphiddev - 1.5*q1*q2*fstar/Yy*phi**(mm-1.d0)* &
                    cosh(1.5d0*q2*p/Yy) * dpddev*1.5d0*q2/Yy + dev/dt/e0dot * dtempddev

            matrix(1,1:2)=[dF1ddee,df1ddev]
            matrix(2,1:2)=[df2ddee,df2ddev]
            call invert_small(matrix,matrix_inv,matrix_det)
            line=[-f1,-f2]
            ddeevv=matmul(matrix_inv,line)
           ! write(*,*) 'm', matrix
           ! write(*,*) ddeevv
            error=dsqrt(ddeevv(1)**2.d0+ddeevv(2)**2.d0)
            dee=dee+ddeevv(1)
            dev=dev+ddeevv(2)


        end do

     taomat1=0.d0
        do i=1,3
            do j=1,3
                if(i==j) then
                    if(sigmastar/=0.d0) then
                    taomat1(i,j)=sstar(i,j)-dee*e/(1.d0+xnu)*1.5d0*sstar(i,j)/sigmastar+(pstar-e/(3.d0*(1.d0-2.d0*xnu))*dev)
                    else
                     taomat1(i,j)=sstar(i,j)+(pstar-e/(3.d0*(1.d0-2.d0*xnu))*dev)
                     end if
                else
                     if(sigmastar/=0.d0) then
                    taomat1(i,j)=sstar(i,j)-dee*e/(1.d0+xnu)*1.5d0*sstar(i,j)/sigmastar
                    else
                     taomat1(i,j)=sstar(i,j)
                     end if
                end if
            end do
        end do

        dematrix=e0dot*dtime/(1.d0-vf)*phi**mm*(dphidse**2.d0+2.d0/9.d0*dphidp**2.d0)**(-1.d0/2.d0)&
        *(dphidse*se+1.d0/3.d0*dphidp*p)

        ematrix1=ematrix+dematrix

end if

vf1=1.d0+(vf-1)*exp(-dev)+fn*dematrix/(ssn*dsqrt(2.d0*3.1415d0))*exp(-0.5d0*((ematrix-een)/ssn)**2.d0)

stress1(1)=taomat1(1,1)
stress1(2)=taomat1(2,2)
stress1(3)=taomat1(3,3)
stress1(4)=taomat1(1,2)
stress1(5)=taomat1(1,3)
stress1(6)=taomat1(2,3)
 !stressDn=stress0
updated_state_variables(1:6)=stress1(1:6)
updated_state_variables(7)=ematrix1
updated_state_variables(8)=vf1

 end subroutine gurson


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_gurson_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardy => vol_avg_shape_function_derivatives_3D
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
    integer      :: n_points,kint,i,j,k,a,l

    real (prec)  ::  strain(6), dstrain(3,3)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: fff_increment(3,3), fff_mid(3,3),fff_mid_inv(3,3),llbar(3,3)
    real (prec)  :: jj,ww(3,3),ll(3,3),llsum,y,volume,strain_increment(3,3),i_invert(3,3)
    real (prec)  :: yita,yita_increment,e0dot,mm,q1,q2,q3,fn,een,ssn,fc,ff,identity(3,3)
    real (prec)  :: iii,rr_increment(3,3),summ,temp(3,3)
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
	

    E = element_properties(1)
    xnu = element_properties(2)
    Y=element_properties(3)
    e0dot=element_properties(4)
    mm=element_properties(5)
    q1=element_properties(6)
    q2=element_properties(7)
    q3=element_properties(8)
    fn=element_properties(9)
    een=element_properties(10)
    ssn=element_properties(11)
    fc=element_properties(12)
    ff=element_properties(13)

    volume=0.d0
    yita=0.d0
    yita_increment=0.d0
    fff_mid=0.d0
    dNbardy=0.d0
    fff_increment=0.d0
    do kint=1,n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        fff_increment=0.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    fff_increment(i,j)=fff_increment(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)

                end do
            end do
        end do
        fff_mid=0.d0
        fff_mid(1,1)=1.d0
        fff_mid(2,2)=1.d0
        fff_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
              fff_mid(i,j)=fff_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do

        call invert_small(fff_mid,fff_mid_inv,jj)
!write (*,*)1,jj
        yita=yita+jj*w(kint)*determinant
        volume=volume+w(kint)*determinant

        ll=matmul(fff_increment,fff_mid_inv)
!        do i=1,3
!            do j=1,3
!                do k=1,3
!                    ll(i,j)=ll(i,j)+fff_increment(i,k)*fff_mid_inv(k,j)
!                end do
!            end do
!        end do
        llsum=ll(1,1)+ll(2,2)+ll(3,3)
        yita_increment=yita_increment+jj*w(kint)*determinant*llsum
       dNdy=matmul(dNdx,fff_mid_inv)
       dNbardy=dNbardy+jj*dNdy*w(kint)*determinant
!        do a=1,n_nodes
!          do i=1,3
!        dNbardy(a,i)=dNbardy(a,i)+jj*dNdy(a,i)*w(kint)*determinant
!           end do
!           end do
    end do
    yita=yita/volume
    yita_increment=yita_increment/volume/yita
    dNbardy=dNbardy/yita/volume

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        fff_increment=0.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    fff_increment(i,j)=fff_increment(i,j)+dNdx(a,j)*(dof_increment(3*(a-1)+i))

                end do
            end do
        end do
        fff_mid=0.d0
        fff_mid(1,1)=1.d0
        fff_mid(2,2)=1.d0
        fff_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    fff_mid(i,j)=fff_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do
        call invert_small(fff_mid,fff_mid_inv,jj)

           llbar=0.d0
           summ=0.d0
       llbar=matmul(fff_increment,fff_mid_inv)
       summ=llbar(1,1)+llbar(2,2)+llbar(3,3)
       do i=1,3
       do j=1,3
         if(i==j)then
          llbar(i,j)=llbar(i,j)+(yita_increment-summ)/3.d0
          end if
          end do
          end do

!        do i=1,3
!            do j=1,3
!                if (i==j) then
!                    do k=1,3
!                        llbar(i,j)=llbar(i,j)+fff_increment(i,k)*fff_mid_inv(k,j)
!                     end do
!                        summ=0.d0
!                        do k=1,3
!                        do l=1,3
!                            summ=summ+fff_increment(l,k)*fff_mid_inv(k,l)
!                        end do
!                        end do
!                      llbar(i,j)=llbar(i,j)+(yita_increment-summ)/3.d0
!                else
!                    do k=1,3
!
!                        llbar(i,j)=llbar(i,j)+fff_increment(i,k)*fff_mid_inv(k,j)
!
!                    end do
!                end if
!            end do
!        end do
  !write(*,*)'llbar', llbar
   strain_increment=(llbar+transpose(llbar))/2.d0
   ww=(llbar-transpose(llbar))/2.d0
!        do i=1,3
!            do j=1,3
!                strain_increment(i,j)=(llbar(i,j)+llbar(j,i))/2.d0
!                ww(i,j)=(llbar(i,j)-llbar(j,i))/2.d0
!            end do
!        end do
        identity=0.d0
        identity(1,1)=1.d0
        identity(2,2)=1.d0
        identity(3,3)=1.d0
        temp=identity-ww/2.d0
        call invert_small(temp,i_invert,iii)

        rr_increment=matmul(i_invert,(identity+ww/2.d0))
    ! write(*,*)'ww',ww
    !write(*,*)'i_invert' ,i_invert
     !write(*,*) 'rr_increment',rr_increment
      !write(*,*) 'strain_increment',strain_increment
        call gurson(element_properties,n_properties,8,initial_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),updated_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),strain_increment,rr_increment,stress)
 !  write(*,*)'stress',stress
         if(updated_state_variables((kint-1)*8+7)>ff)then
            element_deleted= .true.
         else
            element_deleted= .false.
         end if

         B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)+1.d0/3.d0*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(1,2:3*n_nodes-1:3)=1.d0/3.d0*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(1,3:3*n_nodes:3)=1.d0/3.d0*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(2,1:3*n_nodes-2:3)=1.d0/3.d0*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)+1.d0/3.d0*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(2,3:3*n_nodes:3)=1.d0/3.d0*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(3,1:3*n_nodes-2:3)=1.d0/3.d0*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(3,2:3*n_nodes-1:3)=1.d0/3.d0*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)+1.d0/3.d0*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

    element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

end do
  
    return
end subroutine el_gurson_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_gurson_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNbardyi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardy => vol_avg_shape_function_derivatives_3D
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
    real( prec ), intent( inout )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,i,j,a,l

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises ,volume                        ! Pressure and Mises stress
    real (prec)  :: fff_increment(3,3), fff_mid(3,3),fff_mid_inv(3,3),llbar(3,3)
    real (prec)  :: jj,ww(3,3),ll(3,3),llsum,y,strain_increment(3,3),i_invert(3,3)
    real (prec)  :: yita,yita_increment,e0dot,mm,q1,q2,q3,fn,een,ssn,fc,ff,identity(3,3)
    real (prec)  :: iii,rr_increment(3,3),summ,temp(3,3),vvf
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
	
E = element_properties(1)
    xnu = element_properties(2)
    Y=element_properties(3)
    e0dot=element_properties(4)
    mm=element_properties(5)
    q1=element_properties(6)
    q2=element_properties(7)
    q3=element_properties(8)
    fn=element_properties(9)
    een=element_properties(10)
    ssn=element_properties(11)
    fc=element_properties(12)
    ff=element_properties(13)

    volume=0.d0
    yita=0.d0
    yita_increment=0.d0
    fff_mid=0.d0
    dNbardy=0.d0
    fff_increment=0.d0
    do kint=1,n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        fff_increment=0.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    fff_increment(i,j)=fff_increment(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)

                end do
            end do
        end do
        fff_mid=0.d0
        fff_mid(1,1)=1.d0
        fff_mid(2,2)=1.d0
        fff_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
              fff_mid(i,j)=fff_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do

        call invert_small(fff_mid,fff_mid_inv,jj)
!write (*,*)1,jj
        yita=yita+jj*w(kint)*determinant
        volume=volume+w(kint)*determinant

        ll=matmul(fff_increment,fff_mid_inv)
!        do i=1,3
!            do j=1,3
!                do k=1,3
!                    ll(i,j)=ll(i,j)+fff_increment(i,k)*fff_mid_inv(k,j)
!                end do
!            end do
!        end do
        llsum=ll(1,1)+ll(2,2)+ll(3,3)
        yita_increment=yita_increment+jj*w(kint)*determinant*llsum
       dNdy=matmul(dNdx,fff_mid_inv)
       dNbardy=dNbardy+jj*dNdy*w(kint)*determinant
!        do a=1,n_nodes
!          do i=1,3
!        dNbardy(a,i)=dNbardy(a,i)+jj*dNdy(a,i)*w(kint)*determinant
!           end do
!           end do
    end do
    yita=yita/volume
    yita_increment=yita_increment/volume/yita
    dNbardy=dNbardy/yita/volume

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        fff_increment=0.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    fff_increment(i,j)=fff_increment(i,j)+dNdx(a,j)*(dof_increment(3*(a-1)+i))

                end do
            end do
        end do
        fff_mid=0.d0
        fff_mid(1,1)=1.d0
        fff_mid(2,2)=1.d0
        fff_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    fff_mid(i,j)=fff_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do
        call invert_small(fff_mid,fff_mid_inv,jj)

           llbar=0.d0
           summ=0.d0
       llbar=matmul(fff_increment,fff_mid_inv)
       summ=llbar(1,1)+llbar(2,2)+llbar(3,3)
       do i=1,3
       do j=1,3
         if(i==j)then
          llbar(i,j)=llbar(i,j)+(yita_increment-summ)/3.d0
          end if
          end do
          end do

!        do i=1,3
!            do j=1,3
!                if (i==j) then
!                    do k=1,3
!                        llbar(i,j)=llbar(i,j)+fff_increment(i,k)*fff_mid_inv(k,j)
!                     end do
!                        summ=0.d0
!                        do k=1,3
!                        do l=1,3
!                            summ=summ+fff_increment(l,k)*fff_mid_inv(k,l)
!                        end do
!                        end do
!                      llbar(i,j)=llbar(i,j)+(yita_increment-summ)/3.d0
!                else
!                    do k=1,3
!
!                        llbar(i,j)=llbar(i,j)+fff_increment(i,k)*fff_mid_inv(k,j)
!
!                    end do
!                end if
!            end do
!        end do
  !write(*,*)'llbar', llbar
   strain_increment=(llbar+transpose(llbar))/2.d0
   ww=(llbar-transpose(llbar))/2.d0
!        do i=1,3
!            do j=1,3
!                strain_increment(i,j)=(llbar(i,j)+llbar(j,i))/2.d0
!                ww(i,j)=(llbar(i,j)-llbar(j,i))/2.d0
!            end do
!        end do
        identity=0.d0
        identity(1,1)=1.d0
        identity(2,2)=1.d0
        identity(3,3)=1.d0
        temp=identity-ww/2.d0
        call invert_small(temp,i_invert,iii)

        rr_increment=matmul(i_invert,(identity+ww/2.d0))
    ! write(*,*)'ww',ww
    !write(*,*)'i_invert' ,i_invert
     !write(*,*) 'rr_increment',rr_increment
      !write(*,*) 'strain_increment',strain_increment
        call gurson(element_properties,n_properties,8,initial_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),updated_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),strain_increment,rr_increment,stress)
        vvf=updated_state_variables((kint-1)*8+8)
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
            else if (strcmp(field_variable_names(k),'Vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + vvf*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_gurson_3dbasic

