
MODULE mem_sfcanl

  TYPE sfcanl_fields
    REAL, POINTER, DIMENSION(:, :)    ::  u, v, pr, t, td, vv, rh, pm, tad &
                                         , th, the, hi, ps, vor, qm, qcon, div, thad &
                                         , qad, spd, css, vis, fwx, tgd
  END TYPE

  TYPE(sfcanl_fields) :: sfcanl

! Pointers for renaming arrays in lapsvanl
  REAL, POINTER, DIMENSION(:, :) :: &
    u_a, v_a, p_a, t, td, vv, rh, hi, mslp, tadv, theta, thetae, psfc &
    , vort, q, qcon, div, thadv, qadv, spd, cssi, vis, fire, tgd_k

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE alloc_sfcanl_arrays(nxl, nyl)

    IMPLICIT NONE
    INTEGER :: nxl, nyl

    ALLOCATE (sfcanl%u(nxl, nyl))
    ALLOCATE (sfcanl%v(nxl, nyl))
    ALLOCATE (sfcanl%pr(nxl, nyl))
    ALLOCATE (sfcanl%t(nxl, nyl))
    ALLOCATE (sfcanl%td(nxl, nyl))
    ALLOCATE (sfcanl%vv(nxl, nyl))
    ALLOCATE (sfcanl%rh(nxl, nyl))
    ALLOCATE (sfcanl%pm(nxl, nyl))
    ALLOCATE (sfcanl%tad(nxl, nyl))
    ALLOCATE (sfcanl%th(nxl, nyl))
    ALLOCATE (sfcanl%the(nxl, nyl))
    ALLOCATE (sfcanl%hi(nxl, nyl))
    ALLOCATE (sfcanl%ps(nxl, nyl))
    ALLOCATE (sfcanl%vor(nxl, nyl))
    ALLOCATE (sfcanl%qm(nxl, nyl))
    ALLOCATE (sfcanl%qcon(nxl, nyl))
    ALLOCATE (sfcanl%div(nxl, nyl))
    ALLOCATE (sfcanl%thad(nxl, nyl))
    ALLOCATE (sfcanl%qad(nxl, nyl))
    ALLOCATE (sfcanl%spd(nxl, nyl))
    ALLOCATE (sfcanl%css(nxl, nyl))
    ALLOCATE (sfcanl%vis(nxl, nyl))
    ALLOCATE (sfcanl%fwx(nxl, nyl))
    ALLOCATE (sfcanl%tgd(nxl, nyl))

    RETURN
  END SUBROUTINE alloc_sfcanl_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE point_sfcanl_arrays()

    IMPLICIT NONE

    u_a => sfcanl%u
    v_a => sfcanl%v
    p_a => sfcanl%pr
    t => sfcanl%t
    td => sfcanl%td
    vv => sfcanl%vv
    rh => sfcanl%rh
    hi => sfcanl%pm
    mslp => sfcanl%tad
    tadv => sfcanl%th
    theta => sfcanl%the
    thetae => sfcanl%hi
    psfc => sfcanl%ps
    vort => sfcanl%vor
    q => sfcanl%qm
    qcon => sfcanl%qcon
    div => sfcanl%div
    thadv => sfcanl%thad
    qadv => sfcanl%qad
    spd => sfcanl%spd
    cssi => sfcanl%css
    vis => sfcanl%vis
    fire => sfcanl%fwx
    tgd_k => sfcanl%tgd

    RETURN
  END SUBROUTINE point_sfcanl_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE deallocate_sfcanl_arrays()

    IMPLICIT NONE

    IF (ASSOCIATED(sfcanl%u)) DEALLOCATE (sfcanl%u)
    IF (ASSOCIATED(sfcanl%v)) DEALLOCATE (sfcanl%v)
    IF (ASSOCIATED(sfcanl%pr)) DEALLOCATE (sfcanl%pr)
    IF (ASSOCIATED(sfcanl%t)) DEALLOCATE (sfcanl%t)
    IF (ASSOCIATED(sfcanl%td)) DEALLOCATE (sfcanl%td)
    IF (ASSOCIATED(sfcanl%vv)) DEALLOCATE (sfcanl%vv)
    IF (ASSOCIATED(sfcanl%rh)) DEALLOCATE (sfcanl%rh)
    IF (ASSOCIATED(sfcanl%pm)) DEALLOCATE (sfcanl%pm)
    IF (ASSOCIATED(sfcanl%tad)) DEALLOCATE (sfcanl%tad)
    IF (ASSOCIATED(sfcanl%th)) DEALLOCATE (sfcanl%th)
    IF (ASSOCIATED(sfcanl%the)) DEALLOCATE (sfcanl%the)
    IF (ASSOCIATED(sfcanl%hi)) DEALLOCATE (sfcanl%hi)
    IF (ASSOCIATED(sfcanl%ps)) DEALLOCATE (sfcanl%ps)
    IF (ASSOCIATED(sfcanl%vor)) DEALLOCATE (sfcanl%vor)
    IF (ASSOCIATED(sfcanl%qm)) DEALLOCATE (sfcanl%qm)
    IF (ASSOCIATED(sfcanl%qcon)) DEALLOCATE (sfcanl%qcon)
    IF (ASSOCIATED(sfcanl%div)) DEALLOCATE (sfcanl%div)
    IF (ASSOCIATED(sfcanl%thad)) DEALLOCATE (sfcanl%thad)
    IF (ASSOCIATED(sfcanl%qad)) DEALLOCATE (sfcanl%qad)
    IF (ASSOCIATED(sfcanl%spd)) DEALLOCATE (sfcanl%spd)
    IF (ASSOCIATED(sfcanl%css)) DEALLOCATE (sfcanl%css)
    IF (ASSOCIATED(sfcanl%vis)) DEALLOCATE (sfcanl%vis)
    IF (ASSOCIATED(sfcanl%fwx)) DEALLOCATE (sfcanl%fwx)
    IF (ASSOCIATED(sfcanl%tgd)) DEALLOCATE (sfcanl%tgd)

    RETURN
  END SUBROUTINE deallocate_sfcanl_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE nullify_sfcanl_arrays()

    IMPLICIT NONE

    IF (ASSOCIATED(sfcanl%u)) NULLIFY (sfcanl%u)
    IF (ASSOCIATED(sfcanl%v)) NULLIFY (sfcanl%v)
    IF (ASSOCIATED(sfcanl%pr)) NULLIFY (sfcanl%pr)
    IF (ASSOCIATED(sfcanl%t)) NULLIFY (sfcanl%t)
    IF (ASSOCIATED(sfcanl%td)) NULLIFY (sfcanl%td)
    IF (ASSOCIATED(sfcanl%vv)) NULLIFY (sfcanl%vv)
    IF (ASSOCIATED(sfcanl%rh)) NULLIFY (sfcanl%rh)
    IF (ASSOCIATED(sfcanl%pm)) NULLIFY (sfcanl%pm)
    IF (ASSOCIATED(sfcanl%tad)) NULLIFY (sfcanl%tad)
    IF (ASSOCIATED(sfcanl%th)) NULLIFY (sfcanl%th)
    IF (ASSOCIATED(sfcanl%the)) NULLIFY (sfcanl%the)
    IF (ASSOCIATED(sfcanl%hi)) NULLIFY (sfcanl%hi)
    IF (ASSOCIATED(sfcanl%ps)) NULLIFY (sfcanl%ps)
    IF (ASSOCIATED(sfcanl%vor)) NULLIFY (sfcanl%vor)
    IF (ASSOCIATED(sfcanl%qm)) NULLIFY (sfcanl%qm)
    IF (ASSOCIATED(sfcanl%qcon)) NULLIFY (sfcanl%qcon)
    IF (ASSOCIATED(sfcanl%div)) NULLIFY (sfcanl%div)
    IF (ASSOCIATED(sfcanl%thad)) NULLIFY (sfcanl%thad)
    IF (ASSOCIATED(sfcanl%qad)) NULLIFY (sfcanl%qad)
    IF (ASSOCIATED(sfcanl%spd)) NULLIFY (sfcanl%spd)
    IF (ASSOCIATED(sfcanl%css)) NULLIFY (sfcanl%css)
    IF (ASSOCIATED(sfcanl%vis)) NULLIFY (sfcanl%vis)
    IF (ASSOCIATED(sfcanl%fwx)) NULLIFY (sfcanl%fwx)
    IF (ASSOCIATED(sfcanl%tgd)) NULLIFY (sfcanl%tgd)

    RETURN
  END SUBROUTINE nullify_sfcanl_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE
