
!=======================================================
! HYDROTHERMAL ALTERATION OF THERMAL DIFFUSIVITY
!=======================================================
subroutine ReadHydro()
include 'precision.inc'
include 'params.inc'

call AdvanceToNextInputLine( 4 ) ! add to the last line of the input file and param.inc !Liu 7/2017                                                                               
read(4,*) if_hydro, nusselt, xmaxdepth,xmaxt

return
end
!common /hydroth/ xenhc,xwidth,xmaxstr,xmaxdepth,xmaxt,ixh1t
! See Lavier and Buck 2002 JGR for Fig.7 detial explaination !Liu
! hydrothermal cooling is efficient in fault zone
!open( 9, file='hydrother.inp',status='old',err=2001 )
!if_hydro = 1
!read (9,*) ixh1t,xmaxdepth,xenhc,xwidth,xmaxstr,xmaxt
!close( 9 )

!return

!2001 if_hydro = 0
!return

!end

function HydroDiff( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
!common /hydroth/ xenhc,xwidth,xmaxstr,xmaxdepth,xmaxt,ixh1t

    iph = iphase(j,i)
    diff = conduct(iph)/den(iph)/cp(iph) 
!   if(i.gt.ixh1t+1.and.i.lt.nx-ixh1t-1) then
    do  i= 1,nx-1
       do j= 1,nz-1
    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    yc = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
!    if(aps(j,k).gt.xmaxstr) write(*,*) tmpr,yc,aps(j,k)
if( tmpr.le.xmaxt .and.yc.le.xmaxdepth) then
!    xx0 = 0.25*(cord(j,i,1)+cord(j+1,i,1)+cord(j,i+1,1)+cord(j+1,i+1,1))
!    xx = 0.25*(cord(j,k,1)+cord(j+1,k,1)+cord(j,k+1,1)+cord(j+1,k+1,1))
!        xx=abs(xx0-xx)
!        xbets=xenhc*exp(-xx/xwidth)
        HydroDiff = diff * nusselt
        write(*,*) HydroDiff,nusselt
else 
        HydroDiff = diff 
    endif

      end do
   end do
!  continue

!    do 20 k= i+ixh1t,i,-1
!    tmpr = 0.25*(temp(j,k)+temp(j+1,k)+temp(j,k+1)+temp(j+1,k+1))
!    yc = 0.25*(cord(j,k,2)+cord(j+1,k,2)+cord(j,k+1,2)+cord(j+1,k+1,2))
!if( tmpr.le.xmaxt .and. aps(j,k).gt.xmaxstr.and.yc.ge.xmaxdepth) then
!    xx0 = 0.25*(cord(j,i,1)+cord(j+1,i,1)+cord(j,i+1,1)+cord(j+1,i+1,1))
!    xx = 0.25*(cord(j,k,1)+cord(j+1,k,1)+cord(j,k+1,1)+cord(j+1,k+1,1))
!        xx=abs(xx0-xx)
!        xbets=xenhc*exp(-xx/xwidth)
!        HydroDiff = diff * (1.+xbets)
!        write(*,*) HydroDiff,xbets
!       return
!    endif
!20  continue
!    endif
!
!        HydroDiff = diff
!      write(*,*) ixh1t,xmaxdepth,xenhc,xwidth,xmaxstr,xmaxt,HydroDiff,xbets
return
end                      

!=======================================================
! REMESHING STRESS FOR MODE 3 - put hydrostatic pressure to new elements at the bottom
!=======================================================
subroutine rem_stress_alt()
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

common /remeshing/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt

do i = 1,nx-1
    do j = 1,nz-1
        if( cnew(j,i,2) .lt. cold(j,i,2) ) then
            yc = cnew(j,i,2)
            press_norm = pisos-(den(iphsub)+drosub)*g*(yc-rzbo)
            do ii = 1,4 
                stress0(j,i,1,ii) = - press_norm 
                stress0(j,i,2,ii) = - press_norm 
                stress0(j,i,3,ii) = 0.
                stress0(j,i,4,ii) = - press_norm 
            end do
        endif
    end do
end do

return
end
