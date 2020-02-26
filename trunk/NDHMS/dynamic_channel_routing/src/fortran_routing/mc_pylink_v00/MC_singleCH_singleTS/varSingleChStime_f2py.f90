module var
    implicit none
    save

    real :: dt, dx, qup, quc, qdp, qdc, ql, Bw, Tw, TwCC, nCC, Cs, So, n, z, vel, depth
    real :: bfd, WPC, AREAC, C1, C2, C3, C4

    integer :: ntim
    integer :: ncomp0, ncomp, uslinkID, linkID !,nlinks

    real,allocatable,dimension(:) :: qlat !real,allocatable,dimension(:,:) :: qlat
    real,allocatable,dimension(:,:) :: vela, deptha
    real,allocatable,dimension(:,:,:) :: Qd

end module var


