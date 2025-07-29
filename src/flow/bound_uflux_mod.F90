MODULE bound_uflux_mod
    USE core_mod
    USE fields_mod
    USE flowcore_mod

    IMPLICIT NONE(type, external)
    PRIVATE

   
    TYPE, EXTENDS(bound_t) :: bound_uflux_t
    CONTAINS
        PROCEDURE, NOPASS :: front  => bfront
        PROCEDURE, NOPASS :: back   => bfront
        PROCEDURE, NOPASS :: right  => bright
        PROCEDURE, NOPASS :: left   => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top    => bbottom
    END TYPE bound_uflux_t

    TYPE(bound_uflux_t) :: bound_uflux

    PUBLIC :: bound_uflux

CONTAINS

    ! --------------------------------------------------------------------
    ! bfront
    !   Handles flux prolongation at the 'front' (or 'back') boundary 
    !   for x-direction faces. 
    !   Distributes the incoming buffer flux onto the corresponding 
    !   fine-grid faces using area-based weights.
    ! --------------------------------------------------------------------
    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i2, i3, i4, istag1, istag2, dir
        REAL(realk)   :: area1, area2, area3, area4, area5, area6
        REAL(realk)   :: arecvtot, qtot, w1, w2, w3, w4, w5, w6
        REAL(realk), POINTER, CONTIGUOUS :: quo_x(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: quo_x_buf(:, :, :)

        ! Get pointers to field arrays
        CALL f1%get_ptr(quo_x, igrid)
        CALL f1%buffers%get_buffer(quo_x_buf, igrid, iface)
        CALL get_fieldptr(bp, "BP", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        
        SELECT CASE (iface)
        CASE (1)  ! front
            i2 = 2 
            i3 = 3 
            i4 = 4 
            istag1 = 1 
            istag2 = 2 
        CASE (2)  ! back
            i2 = ii - 1
            i3 = ii - 2
            i4 = ii - 3
            istag1 = ii - 1
            istag2 = ii - 2
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Loop over the receiver cells at the interface
        DO j = 3, jj-4, 2
            DO k = 3, kk-3, 2
                ! Read coarse-grid total flux from buffer, with orientation
                qtot = quo_x_buf(k, j, 2)

                area1 = 0.5_realk * bp(k,   j,   i3)
                area2 = 0.5_realk * bp(k+1, j,   i3)
                area3 =             bp(k,   j+1, i3)
                area4 =             bp(k+1, j+1, i3)
                area5 = 0.5_realk * bp(k,   j+2, i3)
                area6 = 0.5_realk * bp(k+1, j+2, i3)

                arecvtot = area1 + area2 + area3 + area4 + area5 + area6
                IF (arecvtot == 0.0_realk) CYCLE

                ! Compute weights for distributing flux proportionally
                w1 = divide0(area1, arecvtot); w2 = divide0(area2, arecvtot)
                w3 = divide0(area3, arecvtot); w4 = divide0(area4, arecvtot)
                w5 = divide0(area5, arecvtot); w6 = divide0(area6, arecvtot)

                ! Deposit distributed flux to fine-grid faces
                quo_x(k  , j  , istag2) = quo_x(k  , j  , istag2) + w1*qtot
                quo_x(k+1, j  , istag2) = quo_x(k+1, j  , istag2) + w2*qtot
                quo_x(k  , j+1, istag2) = quo_x(k  , j+1, istag2) + w3*qtot
                quo_x(k+1, j+1, istag2) = quo_x(k+1, j+1, istag2) + w4*qtot
                quo_x(k  , j+2, istag2) = quo_x(k  , j+2, istag2) + w5*qtot
                quo_x(k+1, j+2, istag2) = quo_x(k+1, j+2, istag2) + w6*qtot
            END DO
        END DO

    END SUBROUTINE bfront


    ! --------------------------------------------------------------------
    ! bright
    !   Handles flux prolongation at the 'right' (or 'left') boundary 
    !   for y-direction faces.
    ! --------------------------------------------------------------------
    SUBROUTINE bright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, i, j2, j3, j4, jstag1, jstag2, dir
        REAL(realk)   :: area1, area2, area3, area4, area5, area6
        REAL(realk)   :: arecvtot, qtot, w1, w2, w3, w4, w5, w6

        REAL(realk), POINTER, CONTIGUOUS :: quo_y(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: quo_y_buf(:, :, :)

        CALL f2%get_ptr(quo_y, igrid)
        CALL f2%buffers%get_buffer(quo_y_buf, igrid, iface)

        CALL get_fieldptr(bp, "BP", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (3)  ! right
            j2 = 2 
            j3 = 3 
            j4 = 4 
            jstag1 = 1 
            jstag2 = 2 
           
            
        CASE (4)  ! left
            j2 = jj - 1
            j3     = jj - 2
            j4 = jj - 3
            jstag2 = jj - 2
            
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT


        DO i = 3, ii-4, 2
            DO k = 3, kk-3, 2
                qtot = quo_y_buf(k, i, 2)

                area1 = 0.5_realk * bp(k,   j3, i  )
                area2 = 0.5_realk * bp(k+1, j3, i  )
                area3 =            bp(k,   j3, i+1)
                area4 =            bp(k+1, j3, i+1)
                area5 = 0.5_realk * bp(k,   j3, i+2)
                area6 = 0.5_realk * bp(k+1, j3, i+2)

                arecvtot = area1 + area2 + area3 + area4 + area5 + area6
                IF (arecvtot == 0.0_realk) CYCLE

                w1 = divide0(area1, arecvtot); w2 = divide0(area2, arecvtot)
                w3 = divide0(area3, arecvtot); w4 = divide0(area4, arecvtot)
                w5 = divide0(area5, arecvtot); w6 = divide0(area6, arecvtot)

                quo_y(k  , jstag2, i  ) = quo_y(k  , jstag2, i  ) + w1*qtot
                quo_y(k+1, jstag2, i  ) = quo_y(k+1, jstag2, i  ) + w2*qtot
                quo_y(k  , jstag2, i+1) = quo_y(k  , jstag2, i+1) + w3*qtot
                quo_y(k+1, jstag2, i+1) = quo_y(k+1, jstag2, i+1) + w4*qtot
                quo_y(k  , jstag2, i+2) = quo_y(k  , jstag2, i+2) + w5*qtot
                quo_y(k+1, jstag2, i+2) = quo_y(k+1, jstag2, i+2) + w6*qtot
            END DO
        END DO

    END SUBROUTINE bright


    ! --------------------------------------------------------------------
    ! bbottom
    !   Handles flux prolongation at the 'bottom' (or 'top') boundary 
    !   for z-direction faces.
    !   Similar logic as bfront/bright but for z-faces.
    ! --------------------------------------------------------------------

    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph


        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: j, i, k2, k3, k4, kstag1, kstag2, dir
        REAL(realk)   :: area1, area2, area3, area4, area5, area6
        REAL(realk)   :: arecvtot, qtot, w1, w2, w3, w4, w5, w6
        REAL(realk), POINTER, CONTIGUOUS :: quo_z(:, :, :), bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: quo_z_buf(:, :, :)

        CALL f3%get_ptr(quo_z, igrid)
        CALL f3%buffers%get_buffer(quo_z_buf, igrid, iface)

        CALL get_fieldptr(bp, "BP", igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (5)
            k2 = 2
            k3 = 3
            k4 = 4
            kstag1 = 1
            kstag2 = 2
          
        CASE (6)
            k2 = kk - 1
            k3 = kk - 2
            k4 = kk - 3
            kstag2 = kk - 2
            
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT


        DO i = 3, ii-4, 2
            DO j = 3, jj-3, 2

                qtot = quo_z_buf(j, i, 2)

                area1 = 0.5_realk * bp(k3, j,   i  )
                area2 = 0.5_realk * bp(k3, j+1, i  )
                area3 =            bp(k3, j,   i+1)
                area4 =            bp(k3, j+1, i+1)
                area5 = 0.5_realk * bp(k3, j,   i+2)
                area6 = 0.5_realk * bp(k3, j+1, i+2)

                arecvtot = area1 + area2 + area3 + area4 + area5 + area6
                IF (arecvtot == 0.0_realk) CYCLE

                w1 = divide0(area1, arecvtot); w2 = divide0(area2, arecvtot)
                w3 = divide0(area3, arecvtot); w4 = divide0(area4, arecvtot)
                w5 = divide0(area5, arecvtot); w6 = divide0(area6, arecvtot)

                quo_z(kstag2, j  , i  ) = quo_z(kstag2, j  , i  ) + w1*qtot
                quo_z(kstag2, j+1, i  ) = quo_z(kstag2, j+1, i  ) + w2*qtot
                quo_z(kstag2, j  , i+1) = quo_z(kstag2, j  , i+1) + w3*qtot
                quo_z(kstag2, j+1, i+1) = quo_z(kstag2, j+1, i+1) + w4*qtot
                quo_z(kstag2, j  , i+2) = quo_z(kstag2, j  , i+2) + w5*qtot
                quo_z(kstag2, j+1, i+2) = quo_z(kstag2, j+1, i+2) + w6*qtot
            END DO
        END DO

    END SUBROUTINE bbottom

END MODULE bound_uflux_mod
