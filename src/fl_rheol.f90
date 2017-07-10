
!  Rheology (Update stresses depending on rheology)
!  Calculate total finite strain and plastic strain  
    
subroutine fl_rheol
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension depl(4)
dimension s11p(4),s22p(4),s12p(4),s33p(4),s11v(4),s22v(4),s12v(4),s33v(4)
logical rh_sel
real :: rate_inject_later
!for Roger_stress
!real, dimension(1:nelem_inject) :: stress_roger, pressure_roger
!real :: dz_elem = 1000, density_dike = 2800

!if( mod(nloop,10).eq.0 .OR. ireset.eq.1 ) then
!    rh_sel = .true.
!else
!    rh_sel = .false.
!endif
rh_sel = .true.

!XXX: irh==11, or irh>=11?
irh=irheol(mphase)
if(irh.eq.11) call init_visc

!if(iynts.eq.1) call init_temp

! Initial stress boundary condition
! Accretional Stresses
if (ny_inject.gt.0) then
         sarc1 = 0.
         sarc2 = 0. 
         if (ny_inject.eq.1) iinj = 1
         if (ny_inject.eq.2) iinj = nx/2 
         !write (*,*) iinj
         !nelem_inject = nz-1 (Tian: Comment this hardwired dike depth)
         !average dx for injection:
         dxinj = 0.
         do jinj = 1,nelem_inject
            !Tian Add temperature to the dike as 1300 C
            !if (mod(int(time*0.5*3.171d-8*1.d-6*10),2).eq.0) then
            !   temp(jinj,iinj) = 1100 
            !endif
            !temp(jinj,iinj + 1) = 1200
            if (mod(int(time*0.5*3.171d-8*1.d-6),2).eq.0) then  !(Tian: added for change color dike in different period)
               iphase(jinj,iinj) = 6 !(Tian:add a line of hardwired dike phase) 
!               iphase(jinj,iinj+1) = 8 !(Tian:add a line of hardwired dike phase) 
!               iphase(jinj,iinj+2) = 6 !(Tian:add a line of hardwired dike phase)  For HD1
!               iphase(jinj,iinj+3) = 6 !(Tian:add a line of hardwired dike phase)  For HD1
!               iphase(jinj,iinj-2) = 6 !Tian, wider_dike                           For HD1
!               iphase(jinj,iinj-1) = 6 !Tian, wider_dike
            elseif (mod(int(time*0.5*3.171d-8*1.d-6),2).eq.1) then !(Tian: added for change color dike in different period 2Myr
               iphase(jinj,iinj) = 7 !(Tian:add a line of hardwired dike phase)
!               iphase(jinj,iinj+1) = 8 !(Tian:add a line of hardwired dike phase)
!               iphase(jinj,iinj+2) = 7 !(Tian:add a line of hardwired dike phase) For HD1
!               iphase(jinj,iinj+3) = 7 !(Tian:add a line of hardwired dike phase) For HD1
!               iphase(jinj,iinj-1) = 7 !Tian, wider_dike
!               iphase(jinj,iinj-2) = 7 !Tian, wider_dike                          For HD1
            endif
            ! iphase(jinj,iinj) = 6 !(Tian:add a line of hardwired dike phase)
            iph=iphase(jinj,iinj)
!            dxinj=dxinj+cord(jinj,iinj+1,1)-cord(jinj,iinj,1)
!            dxinj=dxinj+cord(jinj,iinj+1,1)-cord(jinj,iinj-1,1) !Tian, wider_dike 
            dxinj=dxinj+cord(jinj,iinj+1,1)-cord(jinj,iinj,1) !Tian, wider_dike 
         enddo
         dxinj = dxinj/nelem_inject 
!         dxinj = dxinj/nelem_inject * 3 !Tian, wider_dike
!         dxinj = dxinj/nelem_inject * 6 !Tian, wider_dike For HD1
         ! Constants Elastic:
         poiss = 0.5*rl(iph)/(rl(iph)+rm(iph))
         young = rm(iph)*2.*(1.+poiss)   
         ! Additional Stress:
         !rate_inject_later = 0.5 * rate_inject !Tian
         !if (time*3.171d-8*1.d-6 > 3 .and. time*3.171d-8*1.d-6 < 6) then !(Tian: stop infill after 5Myr, and change M )
         !   rate_inject = rate_inject_later !(Tian: change M)                                   
         !else
         !   rate_inject = rate_inject_later * 2
         !endif !(Tian)
         sarc1 = -young/(1.-poiss*poiss)*rate_inject/dxinj*dt
         sarc2 = sarc1*poiss/(1.-poiss)
         !write(*,*) sarc1,sarc2
endif

!Tian2017 added for kinematic bc phase of diking (begin)
!if (ny_inject.eq.0) then
!   iinj = 1
!   do jinj = 1,nelem_inject
!      if (mod(int(time*0.5*3.171d-8*1.d-6),2).eq.0) then  !(Tian: added for change color dike in different period)
!         iphase(jinj,iinj) = 6 !(Tian:add a line of hardwired dike phase) 
!      elseif (mod(int(time*0.5*3.171d-8*1.d-6),2).eq.1) then !(Tian: added for change color dike in different period 2Myr
!         iphase(jinj,iinj) = 7 !(Tian:add a line of hardwired dike phase)
!      endif
!      iph=iphase(jinj,iinj)
!   enddo
!   write(*,*) "iinj=",iinj,"iph=",iph
!endif
!Tian2017 added for kinematic bc phase of diking (end)


irh_mark = 0

! max. deviatoric strain and area change of current time step
curr_devmax = devmax
curr_dvmax = dvmax

!$OMP Parallel Private(i,j,k,iph,irh,bulkm,rmu,coh,phi,psi, &
!$OMP                  stherm,hardn,vis, &
!$OMP                  de11,de22,de12,de33,dv, &
!$OMP                  s11p,s22p,s12p,s33p, &
!$OMP                  s11v,s22v,s12v,s33v, &
!$OMP                  depl,ipls,diss, &
!$OMP                  sII_plas,sII_visc, &
!$OMP                  quad_area,s0a,s0b,s0) &
!$OMP firstprivate(irh_mark)
!$OMP do schedule(guided) reduction(max: curr_devmax, curr_dvmax)
do 3 i = 1,nx-1
    do 3 j = 1,nz-1
        ! iphase (j,i) is number of a phase NOT a rheology
        iph = iphase(j,i)
        irh = irheol(iph)

!        if(ny_inject.gt.0.and.j.le.nelem_inject) then
!        if(i.eq.iinj.or.i.eq.iinj-1) irh_mark = 1
!        if(i.eq.iinj) irh = 3 
!        endif

        ! Elastic modules & viscosity & plastic properties
        bulkm = rl(iph) + 2.*rm(iph)/3.
        rmu   = rm(iph)

        ! Thermal stresses (alfa_v = 3.e-5 1/K)
        stherm = 0.
        if (istress_therm.gt.0) stherm = -alfa(iph)*bulkm*(temp(j,i)-temp0(j,i))


        ! Preparation of plastic properties
        if (irh.eq.6 .or. irh.ge.11) call pre_plast(i,j,coh,phi,psi,hardn)
              
        ! Re-evaluate viscosity
        if (irh.eq.3 .or. irh.eq.12) then 
            if( mod(nloop,ifreq_visc).eq.0 .OR. ireset.eq.1 ) visn(j,i) = Eff_visc(j,i)
!            if (ny_inject.gt.0.and.i.eq.iinj) visn(j,i) = v_min  
!            if (ny_inject.gt.0.and.i.eq.iinj.and.j.le.nelem_inject) visn(j,i) = 1e18  !Tian1607, continuous plate-->broken plate Turcotte 2010
            if (ny_inject.gt.0.and.i.eq.iinj) visn(j,i) = 1e17  !Tian1607, continuous plate-->broken plate Turcotte 2010
!            if (ny_inject.gt.0.and.i.eq.iinj) visn(j,i) = 1e15  !Tian2017 to solve the grid size sigma_xx high at dike problem, continuous plate-->broken plate Turcotte 2010

!            if (ny_inject.gt.0.and.i.eq.iinj+2.and.j.le.nelem_inject) visn(j,i) = 1e18  !Tian1607, continuous plate-->broken plate Turcotte 2010 For HD1
!            if (ny_inject.gt.0.and.i.eq.iinj+3.and.j.le.nelem_inject) visn(j,i) = 1e18  !Tian1607, continuous plate-->broken plate Turcotte 2010 For HD1
!            if (ny_inject.gt.0.and.i.eq.iinj-2.and.j.le.nelem_inject) visn(j,i) = 1e18  !Tian1607, continuous plate-->broken plate Turcotte 2010 For HD1
!            if (ny_inject.gt.0.and.i.eq.iinj+1.and.j.le.nelem_inject) visn(j,i) = 1e18  !Tian1607, continuous plate-->broken plate Turcotte 2010  wider_dike
!            if (ny_inject.gt.0.and.i.eq.iinj-1.and.j.le.nelem_inject) visn(j,i) = 1e18  !Tian1607, continuous plate-->broken plate Turcotte 2010 wider_dike
!            if(visn(j,i).le.1e20) visn(j,i) = 1e15  !Tian2017
        endif
        vis = visn(j,i)

        ! Cycle by triangles
        do k = 1,4

            ! Incremental strains
            de11 = strainr(1,k,j,i)*dt
            de22 = strainr(2,k,j,i)*dt
            de12 = strainr(3,k,j,i)*dt
            de33 = 0.
            dv = dvol(j,i,k)
            s11p(k) = stress0(j,i,1,k) + stherm 
            s22p(k) = stress0(j,i,2,k) + stherm 
            if(ny_inject.gt.0.and.j.le.nelem_inject) then
                !XXX: iinj is un-init'd if ny_inject is not 1 or 2.
                if(i.eq.iinj) then
!               if(i.eq.iinj.or.i.eq.iinj-2.or.i.eq.iinj-1.or.i.eq.iinj+1.or.i.eq.iinj+2.or.i.eq.iinj+3) then !Tian wider_dike
!               if(i.eq.iinj.or.i.eq.iinj-1.or.i.eq.iinj+1) then !Tian wider_dike
                    s11p(k) = stress0(j,i,1,k) + stherm +sarc1
                    s22p(k) = stress0(j,i,2,k) + stherm +sarc2
                    !!            irh = 1
                endif
            endif
! Roger_stress begins (Tian Stress_Roger 20160810)
!            if(ny_inject.gt.0.and.j.le.nelem_inject) then
!                !XXX: iinj is un-init'd if ny_inject is not 1 or 2.
!               if(i.eq.iinj-1) then
!                  stress_roger(j) = 0.5 * (0.25 * (stress0(j,i,1,1) + stress0(j,i,1,2) +stress0(j,i,1,3)+ &
!                       stress0(j,i,1,4)) + 0.25*( stress0(j,i+2,1,1) + stress0(j,i+2,1,2) +stress0(j,i+2,1,3) +&
!                       stress0(j,i+2,1,4)))
!                  stress_roger(j) = - (density_dike * 10 * (j-0.5) * dz_elem) - stress_roger(j)
!                  stress_roger(j) = stress_roger(j) * 0.5 ! damping
!                  if(stress_roger(j).le.-1.e6) stress_roger(j) = -1.e6
!                  if(stress_roger(j).gt.1.e6) stress_roger(j) = 1.e6
                 
!                  write(*,*) 'j=',j, 'pressure=',- (2800 * 10 * (j-0.5) * dz_elem), 'stress_roger=', stress_roger(j)
!               endif
!               if(i.eq.iinj) then
!!               if(i.eq.iinj.or.i.eq.iinj-1.or.i.eq.iinj+1) then !Tian wider_dike
!                    s11p(k) = stress0(j,i,1,k) + stherm +sarc1+ stress_roger(j)
!                    s22p(k) = stress0(j,i,2,k) + stherm +sarc2+ stress_roger(j) * poiss/(1.-poiss)
!                    !!            irh = 1
!               endif
!            endif
! Roger_stress ends
            s12p(k) = stress0(j,i,3,k) 
            s33p(k) = stress0(j,i,4,k) + stherm
            s11v(k) = s11p(k)
            s22v(k) = s22p(k)
            s12v(k) = s12p(k)
            s33v(k) = s33p(k)
!!            if(abs(sarc11).gt.0.) write(*,*) i,j,sarc11,sarc22
            if (irh.eq.1) then
                ! elastic
                call elastic(bulkm,rmu,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de12)
                irheol_fl(j,i) = 0  
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.eq.3) then
                ! viscous
                call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),de11,de22,de33,de12,dv,&
                     ndim,dt,curr_devmax,curr_dvmax)
                irheol_fl(j,i) = -1  
                stress0(j,i,1,k) = s11v(k)
                stress0(j,i,2,k) = s22v(k)
                stress0(j,i,3,k) = s12v(k)
                stress0(j,i,4,k) = s33v(k)

            elseif (irh.eq.6) then
                ! plastic
                call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                     ten_off,ndim,irh_mark)
                irheol_fl(j,i) = 1
                stress0(j,i,1,k) = s11p(k)
                stress0(j,i,2,k) = s22p(k)
                stress0(j,i,3,k) = s12p(k)
                stress0(j,i,4,k) = s33p(k)

            elseif (irh.ge.11) then 
                ! Mixed rheology (Maxwell or plastic)
                if( rh_sel ) then
                    call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                        s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                        ten_off,ndim,irh_mark)
                    call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                        de11,de22,de33,de12,dv,&
                        ndim,dt,curr_devmax,curr_dvmax)
                else ! use previously defined rheology
                    if( irheol_fl(j,i) .eq. 1 ) then
                        call plastic(bulkm,rmu,coh,phi,psi,depl(k),ipls,diss,hardn,&
                            s11p(k),s22p(k),s33p(k),s12p(k),de11,de22,de33,de12,&
                            ten_off,ndim,irh_mark)
                        stress0(j,i,1,k) = s11p(k)
                        stress0(j,i,2,k) = s22p(k)
                        stress0(j,i,3,k) = s12p(k)
                        stress0(j,i,4,k) = s33p(k)
                    else  ! irheol_fl(j,i) = -1
                        call maxwell(bulkm,rmu,vis,s11v(k),s22v(k),s33v(k),s12v(k),&
                            de11,de22,de33,de12,dv,&
                            ndim,dt,curr_devmax,curr_dvmax)
                        stress0(j,i,1,k) = s11v(k)
                        stress0(j,i,2,k) = s22v(k)
                        stress0(j,i,3,k) = s12v(k)
                        stress0(j,i,4,k) = s33v(k)
                    endif
                endif
            endif
        enddo

        if( irh.ge.11 .AND. rh_sel ) then
            ! deside - elasto-plastic or viscous deformation
            sII_plas = (s11p(1)+s11p(2)+s11p(3)+s11p(4)-s22p(1)-s22p(2)-s22p(3)-s22p(4))**2 &
                     + 4*(s12p(1)+s12p(2)+s12p(3)+s12p(4))**2

            sII_visc = (s11v(1)+s11v(2)+s11v(3)+s11v(4)-s22v(1)-s22v(2)-s22v(3)-s22v(4))**2 &
                     + 4*(s12v(1)+s12v(2)+s12v(3)+s12v(4))**2

            if (sII_plas .lt. sII_visc) then
                do k = 1, 4
                    stress0(j,i,1,k) = s11p(k)
                    stress0(j,i,2,k) = s22p(k)
                    stress0(j,i,3,k) = s12p(k)
                    stress0(j,i,4,k) = s33p(k)
                end do
                irheol_fl (j,i) = 1
            else 
                do k = 1, 4
                    stress0(j,i,1,k) = s11v(k)
                    stress0(j,i,2,k) = s22v(k)
                    stress0(j,i,3,k) = s12v(k)
                    stress0(j,i,4,k) = s33v(k)
                end do
                irheol_fl (j,i) = -1
            endif
        endif


        ! Averaging of isotropic stresses for pair of elements
        if (mix_stress .eq. 1 ) then
        
            ! For A and B couple:
            ! area(n,it) is INVERSE of "real" DOUBLE area (=1./det)
            quad_area = 1./(area(j,i,1)+area(j,i,2))
            s0a=0.5*(stress0(j,i,1,1)+stress0(j,i,2,1))
            s0b=0.5*(stress0(j,i,1,2)+stress0(j,i,2,2))
            s0=(s0a*area(j,i,2)+s0b*area(j,i,1))*quad_area
            stress0(j,i,1,1) = stress0(j,i,1,1) - s0a + s0
            stress0(j,i,2,1) = stress0(j,i,2,1) - s0a + s0
            stress0(j,i,1,2) = stress0(j,i,1,2) - s0b + s0
            stress0(j,i,2,2) = stress0(j,i,2,2) - s0b + s0

            ! For C and D couple:
            quad_area = 1./(area(j,i,3)+area(j,i,4))
            s0a=0.5*(stress0(j,i,1,3)+stress0(j,i,2,3))
            s0b=0.5*(stress0(j,i,1,4)+stress0(j,i,2,4))
            s0=(s0a*area(j,i,4)+s0b*area(j,i,3))*quad_area
            stress0(j,i,1,3) = stress0(j,i,1,3) - s0a + s0
            stress0(j,i,2,3) = stress0(j,i,2,3) - s0a + s0
            stress0(j,i,1,4) = stress0(j,i,1,4) - s0b + s0
            stress0(j,i,2,4) = stress0(j,i,2,4) - s0b + s0
        endif

        if (irh.eq.6 .or. irh.ge.11) then
            !  ACCUMULATED PLASTIC STRAIN
            ! Average the strain for pair of the triangles
            ! Note that area (n,it) is inverse of double area !!!!!
            aps(j,i) = aps(j,i) &
                 + 0.5*( depl(1)*area(j,i,2)+depl(2)*area(j,i,1) ) / (area(j,i,1)+area(j,i,2)) &
                 + 0.5*( depl(3)*area(j,i,4)+depl(4)*area(j,i,3) ) / (area(j,i,3)+area(j,i,4))
            if( aps(j,i) .lt. 0. ) aps(j,i) = 0.

!            if( aps(j,i) .ne. 0. ) aps(j,i) = 0.            ! Tian 2016 Jan29 make aps = 0 all the time

            !	write(*,*) depl(1),depl(2),depl(3),depl(4),area(j,i,1),area(j,i,2),area(j,i,3),area(j,i,4)

            ! LINEAR HEALING OF THE PLASTIC STRAIN
            if (tau_heal .ne. 0.) &
                 aps (j,i) = aps (j,i)/(1.+dt/tau_heal)
!            if (ny_inject.gt.0.and.i.eq.iinj) aps (j,i) = 0.
            if (ny_inject.gt.0.and.i.lt.iinj+5) aps (j,i) = 0. !Tian wider zero plastic strain zero zone
!            if (ny_inject.gt.0.and.i.eq.iinj.or.i.eq.iinj-1) aps (j,i) = 0. !Tian wider_dike
!            if (ny_inject.gt.0.and.i.eq.iinj.or.i.eq.iinj+1) aps (j,i) = 0. !Tian wider_dike
!            if (ny_inject.gt.0.and.i.eq.iinj.or.i.eq.iinj+2) aps (j,i) = 0. !Tian wider_dike For HD1
!            if (ny_inject.gt.0.and.i.eq.iinj.or.i.eq.iinj+3) aps (j,i) = 0. !Tian wider_dike For HD1
!            if (ny_inject.gt.0.and.i.eq.iinj.or.i.eq.iinj-2) aps (j,i) = 0. !Tian wider_dike For HD1
        end if

        ! TOTAL FINITE STRAIN
        strain(j,i,1) = strain(j,i,1) + 0.25*dt*(strainr(1,1,j,i)+strainr(1,2,j,i)+strainr(1,3,j,i)+strainr(1,4,j,i))
        strain(j,i,2) = strain(j,i,2) + 0.25*dt*(strainr(2,1,j,i)+strainr(2,2,j,i)+strainr(2,3,j,i)+strainr(2,4,j,i))
        strain(j,i,3) = strain(j,i,3) + 0.25*dt*(strainr(3,1,j,i)+strainr(3,2,j,i)+strainr(3,3,j,i)+strainr(3,4,j,i))

3 continue
!$OMP end do
!$OMP end parallel

devmax = max(devmax, curr_devmax)
dvmax = max(dvmax, curr_dvmax)

return
end
