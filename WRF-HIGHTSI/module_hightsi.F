MODULE module_hightsi
  USE module_model_constants, only : CP, R_D, XLF, XLV, RHOWATER, STBOLT
  use module_sf_noahlsm, only : RD, SIGMA, CPH2O, CPICE, LSUBF, EMISSI_S

  PUBLIC  HIGHTSI
  PRIVATE CSNOW
  PRIVATE PENMAN
  PRIVATE SNOPAC
  PRIVATE SNOWPACK
  PRIVATE SNOWZ0
  PRIVATE SNOW_NEW

  INTEGER, PRIVATE :: ILOC
  INTEGER, PRIVATE :: JLOC
!$omp threadprivate(iloc, jloc)

  REAL, PARAMETER, PRIVATE :: TFREEZ = 273.15

CONTAINS


  SUBROUTINE emb(t,et,de)
! Calculate the saturation pressure of water vapour et(mb) at temperature t
  IMPLICIT NONE
  REAL, INTENT(IN)  :: t ! temperature (K)
  REAL, INTENT(OUT) :: et  ! water vapour (mb)
  REAL, INTENT(OUT) :: de  ! derivative term of et

  IF(t.GE.273.15) THEN
    et = EXP((-6763.6/t)-4.9283*log(t)+54.23)
    de = (exp((-6763.6/t)-4.9283*log(t)+54.23))*(6763.6/t/t-4.9283/t)
  ELSEIF(t.LT.273.15) THEN
    et = exp((-6141.0/t)+24.3)
    de=(exp((-6141.0/t)+24.3))*(6141.0/t/t)
  ELSE
    et=6.25
    de=0.0
  ENDIF
  END SUBROUTINE emb


  SUBROUTINE wvap(t,q,es,qs,e)
!Calculate water vapour pressure: e; specific humidities: qs, q
  IMPLICIT NONE
  REAL, INTENT(IN) :: t   !skin temperature (K)
  REAL, INTENT(IN) :: q   !specific humidity
  REAL, INTENT(OUT) :: qs   !specific humidity
  REAL, INTENT(OUT) :: es  !saturation pressure at temperature of a surface
  REAL, INTENT(OUT) :: e   !water vapour pressure: e
  REAL :: de,pe

  pe = 1013.0
! 1. Calculate saturation pressure at temperature of a surface: es=emb(tsfc)
  CALL emb(t,es,de)
!2. Calculate water vapour pressure in air: e
  e=q*pe/(0.622+0.378*q)

!3. Specific humidities qs and q calculated from es and e
  qs=0.622*es/(pe-0.378*es)

  END SUBROUTINE wvap

  SUBROUTINE coe(tsfc,tair,uz,z,qs,q,ch,ce)
!Calculate the drag coefficient cd, Stanton number ce
! Dalton number ch; ce=ch
!friction velocity uf  and  z/L
  IMPLICIT NONE
  REAL, INTENT(IN) :: tsfc ! surface temperature [K]
  REAL, INTENT(IN) :: tair ! surface air temperature [K]
  REAL, INTENT(IN) :: uz   ! surface wind speed [m/s]
  REAL, INTENT(IN) :: z    ! height of the atmospheric parameters [m]
! REAL, INTENT(IN) :: zo   ! aerodynamic roughness length zo  [m]
! REAL, INTENT(IN) :: rat  ! roughness ratio zo/zt [-]
  REAL, INTENT(IN) :: qs   ! specific humidity
  REAL, INTENT(IN) :: q    ! specific humidity
! REAL, INTENT(OUT) :: vr  ! reference level wind speed [m/s]
  REAL, INTENT(OUT) :: ch  ! turbulent exchange coefficient
  REAL, INTENT(OUT) :: ce  ! ##
  REAL :: ct  ! z/L  [ ]
  REAL :: z0  ! roughness [m]
  REAL :: uf  ! friction velocity [m/s]
  REAL :: zt  ! roughness [m]
  REAL :: zq  ! roughness [m]
  REAL :: cd  ! drag coefficient

  REAL :: qt, zz
  REAL :: t0,cdn,cen,chn,lnz0,lnzt,lnzq,rdm,rdt,rdq,rlnz0,rz0,re,rb0,rb1,&
          rb2,b,rz,ufm,ufh,ufe,tem1,tem2,tem3,unfm,unfh,unfe,r1,r2
  REAL :: kc  ! obs. surface roughness [cm]
  REAL :: rk  ! Von karman constant [-]
  REAL :: g   ! gravity constant
  REAL :: pi  ! universal constant

!1.  Define constant coefficients
  kc = 3.0
  rk = 0.405
  g  = 9.81
  pi = 3.1415926
  r1 = 16.0
  r2 = 16.0 ! original Businger-Dyer constant
  t0 = 0.5*(tsfc+tair)
!2. Calculate the roughness z0, zt, zq

            cdn=(1.10+0.072*kc)*10.0**(-3)
            rlnz0=log(10.0)-rk*cdn**(-0.5)
            rz0=exp(rlnz0)
            re=rz0*cdn**(0.5)*uz/((0.9065*tair-112.7)/10.0**7)

               if(re.gt.0.0.and.re.lt.0.135) then
                  rb0=1.43
                  rb1=0.0
                  rb2=0.0
               elseif(re.gt.0.135.and.re.lt.2.5) then
                  rb0=0.25
                  rb1=-0.589
                  rb2=0.0
               elseif(re.gt.2.5.and.re.lt.1000) then
                  rb0=0.356
                  rb1=-0.538
                  rb2=-0.181
               endif

            lnz0=log(10.0)-rk*cdn**(-0.5)
            lnzt=lnz0+rb0+rb1*log(re)+rb2*(log(re))**2.0
            lnzq=lnzt

            rdm=log(z)-lnz0
            rdt=log(z)-lnzt
            rdq=log(z)-lnzq

! Roughness  z0,zt,zq

            z0=exp(lnz0)
            zt=exp(lnzt)
            zq=zt

!3. Calculate the Bulk-Richardson number and z/L
!       b=1+0.61*t0*(q-qs)/(tair-tsfc)

! Cal. taking moisture effect "b" into account

!       rz=(z*g*(tair-tsfc))/(t0*uz**2.0)*b
        rz=(z*g*(tair-tsfc))/(t0*uz**2.0)

!3.1  Unstable condation  rz < 0.0, i.e. z/L < 0.0
       if (rz.lt.0.0) then

         if(r1.eq.r2) then
             ct=(rdm*rdm/rdt-0.55)*rz   !:(Launiainen1995, eq(11))         
         else
             ct=(rdm*rdm/rdt-0.9)*rz    !:(Launiainen1995, eq(11))         
         endif

         ufm=(1-r1*ct)**(-1.0/4.0)
         ufh=(1-r2*ct)**(-1.0/2.0)
         ufe=ufh

         tem1=2*log(0.5*(1+ufm**(-1.0)))
         tem2=log(0.5*(1+ufm**(-2.0)))
         tem3=2*atan(ufm**(-1.0))

! Universal functions (B-D) as follows

         unfm=tem1+tem2-tem3+pi/2.0
         unfh=2*log(0.5*(1+ufh**(-1.0)))
         unfe=unfh

! 3.2  Stable condition  rz > 0.0, i.e. z/L > 0.0

       elseif(rz.gt.0.0) then

         ct=(1.89*rdm+44.2)*rz*rz+(1.18*rdm-1.5*log(z0/zt)-1.37)*rz

! above ct: (Launiainen1995, eq(9))

! Universal functions (Holtslag) as follows

         unfm=-0.7*ct-0.75*(ct-14.2857)*exp(-0.35*ct)-10.7143
         unfh=unfm
         unfe=unfh

       endif

! 4. Finnal coefficients cd,ch,ce are calculated as follow

       cd=rk*rk*(rdm-unfm)**(-2.0)
       ch=rk*rk/(rdm-unfm)/(rdt-unfh)
       ce=rk*rk/(rdm-unfm)/(rdq-unfe)
       uf=uz*rk/(rdm-unfm)     ! friction velocity

  RETURN
  END SUBROUTINE coe

  REAL FUNCTION rl1(t)
!Calculate the enthalpy of vaporiztion/fusion [J]
  IMPLICIT NONE
  REAL :: t  !temperature [K]
  REAL :: abst

  abst = 273.15

        if(t.ge.abst) then

        rl1=(2500.-2.375*(t-abst))*1000.0

        elseif(t.lt.abst) then

        rl1=(2500.-2.375*(t-abst))*1000.+335.*10.**3

        endif

  return
  END FUNCTION rl1

  SUBROUTINE inice(cn,dih,nvi,totalq,md,qin,rsi)
!The solar radiation penetrated into the ice layer
  IMPLICIT NONE
  REAL, INTENT(IN) :: cn      !cloudiness factor
  REAL, INTENT(IN) :: dih     !dih=hi/nvi
  REAL, INTENT(IN) :: totalq  !total available Qs
  INTEGER, INTENT(IN) :: nvi  !number of layers in ice
  INTEGER, INTENT(IN) :: md   !method of solar radiation penetration in ice
!solar radiation flux distributed within the ice
  REAL, DIMENSION(100), INTENT(OUT) ::  qin
  REAL, DIMENSION(100), INTENT(OUT) ::  rsi  !  derivative of qin

  REAL :: ag,bg,cg,dg,ek,vzi
     !  ag,             ! parameters    
     !  bg,             ! parameters    
     !  cg,             ! parameters
     !  dg,             ! parameters
     !  ek,             ! extinction coefficient of ice
     !  vzi             ! increasing segment of depth in ice

  INTEGER :: i

  DO i = 1,nvi
  vzi = dih*i

          if(md.eq.1)then                    ! exp + exp(Launiainen and Cheng, 1998)
           if(vzi.le.0.1)then
           ek = 17*(1-cn)+10.5*cn
           qin(i)= totalq*exp(-ek*vzi)
           rsi(i) = ek*totalq*exp(-ek*vzi)
           elseif(vzi.gt.0.1)then
           ek = 1.5
           qin(i) = totalq*(0.1826*(1-cn)+0.3499*cn)&
                    *exp(-ek*(vzi-0.1))
           rsi(i)  = ek*totalq*(0.1826*(1-cn)+0.3499*cn)&
                *exp(-ek*(vzi-0.1))    
           endif
        
          elseif(md.eq.2)then                   ! linear + exp (Sahlberg 1988)
           if(vzi.le.0.1)then
           qin(i) = totalq*(1-(8.174*(1-cn)+6.501*cn)*vzi)
           rsi(i)  = totalq*(8.174*(1-cn)+6.501*cn)
           elseif(vzi.ge.0.1)then
           ek = 1.5
           qin(i) = totalq*(0.1826*(1-cn)+0.3499*cn)&
                    *exp(-ek*(vzi-0.1))
           rsi(i)  = ek*totalq*(0.1826*(1-cn)+0.3499*cn)&
                *exp(-ek*(vzi-0.1))    
           endif
          elseif(md.eq.3)then           ! exp (Grenfell and Maykut, 1979)
           ag  = 0.3894*(1-cn)+0.3456*cn
           bg  = 12.39*(1-cn)+10.30*cn
           cg  = 0.350*(1-cn)+0.590*cn
           dg  = 1.578*(1-cn)+1.315*cn
        
           ek  = ag*exp(-bg*vzi)+cg*exp(-dg*vzi)
           qin(i) = totalq*ek
           rsi(i)  = totalq*(ag*bg*exp(-bg*vzi)+cg*dg*exp(-dg*vzi))
        
          elseif(md.eq.4)then           ! exp(Ebery and Curry,1993)
           ek  = 1.5
           qin(i) = totalq*(0.1826*(1-cn)+0.3499*cn)*exp(-ek*vzi)
           rsi(i)  = ek*totalq*(0.1826*(1-cn)+0.3499*cn)&
                *exp(-ek*vzi)
          endif

  END DO
  RETURN
  END SUBROUTINE inice

  SUBROUTINE insnow1(extink,hs,totalq,qin,method)
!Calculate the distribution of solar radiation for surface heat balance
!Single thin snow layer as a whole (heat conduction equation is not applied)
!Snow thickness  Hs:  0< Hs <=  0.01m = 1cm
!This code only consider snow thickness less than 1.0 cm
  IMPLICIT NONE
  REAL, INTENT(IN) :: extink,hs,totalq
  REAL, DIMENSION(100), INTENT(OUT) :: qin
  CHARACTER(3), INTENT(IN) :: method
!extink  extinction coefficient of snow
!hs      Snow thickness
!totalq  downward solar radiation at the surface
!qin     solar radaition penetrated into snow/ice
!        only calculate qin(1) for thin snow layer(hs)(<0.01m)
!method  represent the way of penetration (linear/exponential)

  IF (method.EQ.'exp') THEN
! exponential attenuate approximation
  qin(1) = totalq*EXP(-extink*hs)
  ELSEIF (method.EQ.'lin') THEN
! linear attenuate approximation
  qin(1) = totalq*(1-9*hs)
  ENDIF

  RETURN
  END SUBROUTINE insnow1

  SUBROUTINE insnow2(extink,dsh,nvs,totalq,md,qb,qin,rss)
!Snow layer is larger than 0.01m (heat conduction equation is considered)
!The solar radiation penetrated into the snow layer
!parametered according to different approach
  IMPLICIT NONE
  REAL, INTENT(IN) :: dsh,totalq,extink
  INTEGER, INTENT(IN) :: md,nvs
  REAL, DIMENSION(100), INTENT(OUT) :: qin,rss
  REAL, INTENT(OUT) :: qb
  INTEGER :: i
  REAL :: vzs
!qin     solar flux distributed within snow layer
!rss     This is the deriative of "qin"
!dsh     dsh = hs/nvs
!hs      snow thickness
!nvs     umber of layers in snow
!totalq  downward solar radiation at the surface
!qb      solar radaitive flux attenuaded at the snow/ice interface
!extink  extinction coefficient of snow
!md      method of solar radiation penetration(It is more meaningful for ice)
!nvs     number of layers in snow
!vzs     increasing segment of depth in snow, refers to the grid numbers

  DO i = 1,nvs

  vzs = dsh*i

  IF (md.eq.1) THEN  !exp + exp (Launiainen and Cheng 1998)
  qin(i) = totalq*EXP(-extink*vzs)
  rss(i) = extink*totalq*EXP(-extink*vzs)

  ELSEIF (md.eq.2) THEN  ! linear + exp (Sahlberg 1988)
    IF (vzs.LE.0.1) THEN ! linear variation within first 10 cm
    qin(i)  = totalq*(1-9*vzs)
    rss(i)   = totalq*9
    ELSEIF (vzs.GT.0.1) THEN !exponential variation below 10 cm
    qin(i)  = totalq*0.1*exp(-extink*(vzs-0.1))
    rss(i)   = extink*totalq*0.1*exp(-extink*(vzs-0.1))
    ENDIF

  ELSEIF (md.eq.3) THEN  ! exp (Grenfell and Maykut 1979)
  qin(i)  = totalq*EXP(-extink*vzs)
  rss(i)  = extink*totalq*exp(-extink*vzs)

  ELSEIF (md.eq.4) THEN  ! exp (Ebery and Curry, 1993)
  qin(i)  = totalq*EXP(-extink*vzs)
  rss(i)  = extink*totalq*exp(-extink*vzs)

  ENDIF

  END DO

  qb = qin(nvs)

  RETURN
  END SUBROUTINE insnow2

  SUBROUTINE ts1(tsur,tair,q,uz,z,dp,qsum,ql,&
              sss,kkkk,dih,dsh,ns,tice,tsnow,rki,rks)
  IMPLICIT NONE
!Calculate the surface temperature with iterative method
  REAL, INTENT(INOUT) :: tsur
  REAL, INTENT(IN) :: tair,q,uz,z,dp,sss,dih,dsh,rki,rks
  REAL, INTENT(IN) :: qsum,ql
  INTEGER, INTENT(IN) :: kkkk,ns
  REAL, DIMENSION(100), INTENT(IN) :: tice,tsnow ! [C]
  REAL :: tp,tp1,ti,es,qs,e,ch,ce,des,de
  REAL :: deh,dele,deql,cgm,detc,dfdts,sumd,dte,cp,h,rle,abst
  INTEGER :: i
  CHARACTER(255) ::  message

  cp  = 1004.0
  cgm = 5.67e-8
  abst = 273.15


!Do loop for surface and snow/ice interface temperature calculation

  tp = tsur

  DO i = 1,15

!Calculate snow/ice interface temp. (ti)
! tint(tsfc,sss,kkkk,dih,dsh,ns,tice,tsnow,rki,rks,ti)
  CALL tint(tp,sss,kkkk,dih,dsh,ns,tice,tsnow,rki,rks,ti)

!Calculate the water vapour pressure es, e; qs; q;
  CALL wvap(tp,q,es,qs,e)

!Calculate drag coefficients cd, ch, and ce=ch
!coe(tsfc,tair,uz,z,qs,q,ch,ce)
  CALL coe(tp,tair,uz,z,qs,q,ch,ce)

!Calculate turbulent heat fluxes
  h   = -dp*cp*ch*(tp-tair)*uz          ! sensible heat
  rle = -dp*ce*(qs-q)*uz*rl1(tp)         ! latent heat

!Calculate long-wave radiation
! CALL qlw1(qlin,tp)

!Calculate derivative term of TS (TP) (Norton-iterative method)
  deh = -dp*cp*ch*uz       ! derivative term of sensible heat

!derivative term of latent heat
  CALL emb(tp,des,de)
  dele = -dp*ce*uz*rl1(tp)*0.622*1013.0/&
         ((1013.0-0.378*des)**2)*de

!derivative term of long wave
  deql = -0.97*4*cgm*tp**3

!derivative term of thermal conductivity
  if(sss.le.0.005)then      !: without snow
  detc = -rki/dih
  elseif(sss.gt.0.005.and.sss.lt.0.01)then  !: with snow
  detc = rks**2*dih/(sss*(rks*dih+rki*sss))-rks/sss
  elseif(sss.ge.0.01)then
    if(kkkk.le.15)then
    detc = rks**2*dih/(sss*(rks*dih+rki*sss))-rks/sss
    else
    detc = -rks/dsh
    endif
  endif

!sum of derivative terms
   dfdts = deh+dele+deql+detc

  if(sss.le.0.005)then
  sumd = qsum+h+rle+ql+rki/dih*((tice(1)+abst)-tp)

  elseif(sss.gt.0.005.and.sss.lt.0.01)then
  sumd = qsum+h+rle+ql+rks/sss*(ti-tp)

  elseif(sss.ge.0.01)then
    if(kkkk.le.15)then
    sumd = qsum+h+rle+ql+rks/sss*(ti-tp)
    else
    sumd = qsum+h+rle+ql+rks/dsh*((tsnow(1)+abst)-tp)
    endif
  endif

  tp1 = tp-sumd/dfdts  ! surface temperature at next step
  dte = abs(tp1-tp)
  IF (ISNAN(tp1)) THEN
!     WRITE(message, *),'qsum,h,rle,ql,sss,',qsum,h,rle,ql,sss
!     WRITE(message,*),'dp,ch,ce,uz,rl1,sumd,dfdts,',dp,ch,ce,uz,rl1(tp),sumd,dfdts
      WRITE(message, *),'sss,dih,tice(1),ti,qsum',sss,dih,tice(1),ti,qsum
      CALL wrf_message(message)
      WRITE(message, *),'ql,qs,q,tp,uz',ql,qs,q,tp,uz
      CALL wrf_message(message)
      WRITE(message, *),'tair,ce,sumd,dfdts,',tair,ce,sumd,dfdts
      CALL wrf_error_fatal(message)
  ENDIF

  IF (dte.LT.0.01) GOTO 333
  IF (tp1 .LT. 100) THEN
      WRITE(message, *),'tp1 too low at iter',i,tp1,'use appr'
      CALL wrf_message(message)
      CALL ts1_appr(tsur,tp1,tair,q,uz,z,qsum,ql,rki,rks,dih,dsh,&
                    dp,sss,tice(1),tsnow(1),kkkk)
      WRITE(message, *),'in appr: tsur,tp1,tair',tsur,tp1,tair
      CALL wrf_message(message)
      WRITE(message, *),'in appr: sss,dih',sss,dih
      CALL wrf_message(message)
      WRITE(message, *),'in appr: uz,dp,',uz,dp
      CALL wrf_message(message)
      WRITE(message, *),'in appr: tsnow,tice',tsnow(1),tice(1)
      CALL wrf_message(message)
      GOTO 333
!     CALL wrf_error_fatal(message)
  ENDIF
  tp = tp1
  END DO

!Set surface temperature equal the air temperature
!if iterative method is divergence in do loop 100 (15 steps)

  tsur = tair
  GOTO 444

333 CONTINUE
  tsur = tp1

444 CONTINUE

  RETURN
  END SUBROUTINE ts1

  SUBROUTINE ts1_appr(tp,tp1,tair,q,uz,z,qsum,ql,rki,rks,dih,dsh,&
                      dp,sss,tice1,tsnow1,kkkk)
  IMPLICIT NONE
  REAL, INTENT(IN)  :: tp,tair,q,uz,z,qsum,ql,tice1,tsnow1
  REAL, INTENT(IN)  :: rks,rki,dih,dsh,dp,sss
  INTEGER, INTENT(IN) :: kkkk
  REAL, INTENT(OUT) :: tp1
  REAL :: es,qs,e,ch,ce,ddtt,h,rle
  REAL, PARAMETER :: cgm = 5.67e-8
  REAL, PARAMETER :: abst = 273.15

! cp  = 1004.0
! cgm = 5.67e-8
!Calculate the water vapour pressure es, e; qs; q;
  CALL wvap(tp,q,es,qs,e)

!Calculate drag coefficients cd, ch, and ce=ch
!coe(tsfc,tair,uz,z,qs,q,ch,ce)
  CALL coe(tp,tair,uz,z,qs,q,ch,ce)

!Calculate turbulent heat fluxes
  h   = -dp*cp*ch*(tp-tair)*uz          ! sensible heat
  rle = -dp*ce*(qs-q)*uz*rl1(tp)         ! latent heat

  IF (sss.LT.0.005) THEN
     ddtt = (h+rle+qsum+ql-0.97*cgm*tp**4+rki/dih*(tice1+abst-tp)) &
            /(rki/dih+4*0.97*cgm**3)
     tp1  = tp+ddtt
  ELSEIF (sss.GT.0.005.AND.sss.LE.0.01) THEN
     ddtt = (h+rle+qsum+ql-0.97*cgm*tp**4+rks*rki*(tice1+abst-tp) &
            /(rks*dih+sss*rki))/(4*0.97*cgm*tp**3+rks*rki/(rks*dih+sss*rki))
     tp1  = tp+ddtt
  ELSEIF (sss.GT.0.01) THEN
     IF (kkkk.LT.15) THEN
        ddtt = (h+rle+qsum+ql-0.97*cgm*tp**4+rks*rki*(tice1+abst-tp) &
               /(rks*dih+sss*rki))/(4*0.97*cgm*tp**3+rks*rki/(rks*dih+sss*rki))
        tp1  = tp+ddtt
     ELSE
        ddtt = (h+rle+qsum+ql-0.97*cgm*tp**4+rks/dsh*(tsnow1+abst-tp)) &
               /(rks/dsh+4*0.97*cgm*tp**3)
        tp1  = tp+ddtt
     ENDIF
  ENDIF

  RETURN
  END SUBROUTINE ts1_appr

  SUBROUTINE qsall1(toqs,qsum,extink,cn,dih,dsh,nvi,nvs,md,sss,kkkk,&
                    qini,qins,rsi,rss)
!Calculate the distribution of solar radiation for surface heat balance
  IMPLICIT NONE
  REAL, INTENT(IN) :: extink
  REAL, INTENT(IN) :: toqs
  REAL, INTENT(IN) :: cn,dih,dsh,sss
  INTEGER, INTENT(IN) :: nvi,nvs,md,kkkk
  REAL, INTENT(OUT) :: qsum
!extink    extinction coefficient for snow
!toqsum   net short-wave radiation at surface [W/m2]
!qsum     energy for the surface heat balance

  REAL :: qb,toqsum
  REAL, DIMENSION(100), INTENT(OUT) :: qini,rsi,qins,rss
  toqsum = toqs

  IF(sss.LE.0.005) THEN
  CALL inice(cn,dih,nvi,toqsum,md,qini,rsi)
!qini(1): solar energy available at the interface of surface layer after being
!attenuated in this layer according to the Bouguer-Lambert law

  qsum = toqsum -qini(1)

  ELSEIF(sss.GT.0.005.AND.sss.LT.0.01) THEN
  CALL insnow1(extink,sss,toqsum,qins,'exp')
  qsum   = toqsum - qins(1)
  toqsum = qins(1)
  CALL inice(cn,dih,nvi,toqsum,md,qini,rsi)
        
  ELSEIF(sss.GE.0.01) THEN
    IF(kkkk.LE.15) THEN
    CALL insnow1(extink,sss,toqsum,qins,'exp')
    qsum   = toqsum - qins(1)
    toqsum = qins(1)
    CALL inice(cn,dih,nvi,toqsum,md,qini,rsi)
        
    ELSEIF(kkkk.GT.15) THEN
    CALL insnow2(extink,dsh,nvs,toqsum,md,qb,qins,rss)
    qsum   = toqsum - qins(1)
    toqsum = qb
    CALL inice(cn,dih,nvi,toqsum,md,qini,rsi)
    ENDIF
  ENDIF

  RETURN
  END SUBROUTINE qsall1

  SUBROUTINE tint(tsfc,sss,kkkk,dih,dsh,ns,tice,tsnow,rki,rks,ti)
!Calculate snow/ice interface temperature
  REAL, INTENT(IN) :: tsfc  ! surface temperature [K]
  REAL, INTENT(IN) :: sss,dih,dsh,rks,rki
  INTEGER, INTENT(IN) :: kkkk,ns
  REAL, DIMENSION(100), INTENT(IN) :: tice, tsnow
  REAL, INTENT(OUT) :: ti   ! snow/ice interface temp. [K]
  REAL :: abst
! rks    thermal conductivity of snow [W/mK]
! rki    thermal conductivity of ice [W/mK]
! kkkk   loop account of snow part
! sss    snow thickness [m]
! tice   ice temp.  [deg.]
! tsnow  snow temp. [deg.]

    abst = 273.15


        if(sss.le.0.005)then
        ti=tsfc

        elseif(sss.gt.0.005.and.sss.lt.0.01)then
        ti=(rks*dih*tsfc+rki*sss*(tice(1)+abst))/(rks*dih+rki*sss) ! new(2)

        elseif(sss.ge.0.01)then
           if(kkkk.le.15)then
           ti=(rks*dih*tsfc+rki*sss*(tice(1)+abst))/(rks*dih+rki*sss)
           else
           ti=(rks*dih*(tsnow(ns)+abst)+rki*dsh*(tice(1)+abst))&    ! (3)
               /(rks*dih+rki*dsh)
           endif
        endif

        IF ( ISNAN(ti) ) THEN
        PRINT *, 'NAN ti' ,sss,tsfc
        STOP
        ENDIF

  RETURN
  END SUBROUTINE tint


  SUBROUTINE dhs(tsur,ti,tair,q,uz,z,dt,dp,ql,qsum,sss,hhh,rks,rki,dsh,dih,&
                 tsnow,tice,nvi,nvs,kkkk,sdth0)
!Calculate the snow or ice melting at the surface
!This is considered in a classical way, i.e. if Tsfc >0, then Tsfc=0, and
!melting is calculated for the surface layer(thickness between real surface
!and first inner grid point in snow or ice) only.
  IMPLICIT NONE
  REAL, INTENT(IN) :: tsur,q,tair,uz,z,ql,qsum,dp,rki,rks,dt,ti
  REAL, DIMENSION(100), INTENT(IN) :: tice,tsnow
  INTEGER, INTENT(IN) :: nvi,kkkk,nvs
  REAL, INTENT(INOUT) :: sss,hhh,dih,dsh
  REAL, INTENT(OUT) :: sdth0
  REAL :: fcontsi,sumd
  REAL :: es,qs,e,ch,ce,h,rle
  REAL :: abst,qii,qsi,cp

  abst = 273.15
  qii  = 302e6       ! heat capacity of ice [J/m3]
  qsi  = 110e6       ! heat capacity of snow [J/m3]
  cp   = 1004.0
!fcontsi  heat conductive flux
!rl       rl(T) enthalpy of vaporization [J]
!sumd     sum of the fluxes
!sdth0    mass change at the surface [m]

!Calculate the ice melt on the surface

!Calculate surface heat fluxes

  CALL wvap(tsur,q,es,qs,e)

  CALL coe(tsur,tair,uz,z,qs,q,ch,ce)

! CALL qlw1(qlin,tsur)

  h=-dp*cp*ch*(tsur-tair)*uz       ! Sensible heat
  rle=-dp*ce*(qs-q)*uz*rl1(tsur)    ! Latent heat

  IF(sss.LE.0.005)THEN

  fcontsi = rki/dih*((tice(1)+abst)-tsur) ! new version for ice
  sumd = qsum+h+rle+ql- fcontsi            ! new version for ice
  sdth0 = -dt/qii*sumd                   ! surface melt value(ice)
  hhh   = hhh+sdth0                     ! surface melt value(ice)
  dih=hhh/nvi
  PRINT *,sdth0,qsum,h,rle,ql,fcontsi

  ELSEIF(sss.GT.0.005.and.sss.LT.0.01) THEN

  fcontsi = -1.0*rki/dih*(tice(1)+abst-ti)
  sumd  =  qsum+h+rle+ql+ fcontsi
  sdth0 = -dt/qii*sumd                   ! surface melt value(ice)
  hhh   = hhh+sdth0                     ! surface melt value(ice)
  dih=hhh/nvi
! PRINT *,sdth0,qsum,h,rle,ql,fcontsi
! PRINT *,tice(1)+abst,ti,tsur

  fcontsi = rks/sss*(ti-tsur)
  sumd  =  qsum+h+rle+ql- fcontsi

  sdth0 = -dt/qsi*sumd    ! surface melt value(snow)
  sss   = sss+sdth0                     ! surface melt value(snow)
  dsh=sss/nvs

  ELSEIF(sss.GE.0.01) THEN
    IF(kkkk.LE.15) THEN
    fcontsi = rks/sss*(ti-tsur)
    sumd = qsum+h+rle+ql- fcontsi
    ELSE 
    fcontsi = rks/dsh*((tsnow(1)+abst)-tsur)
    sumd = qsum+h+rle+ql- fcontsi
    ENDIF

  sdth0 = -dt/qsi*sumd    !: surface melt value(snow)
  sss   = sss+sdth0                     ! surface melt value(snow)
  dsh=sss/nvs

  ENDIF

  RETURN
  END SUBROUTINE dhs

  SUBROUTINE nume1(nl,hsi,tsfc,tbn,tsnow,tice,y,vrs,chr,rk0,dt,srk,srki0,srkin)
!Calculate the ice temperature
!Numerical method as a whole
  CHARACTER(1), INTENT(IN) :: chr
  INTEGER, INTENT(IN) :: nl
  REAL, INTENT(IN) :: hsi,tsfc,tbn,rk0,dt
  REAL, DIMENSION(100), INTENT(IN) :: vrs,tsnow,tice
  REAL, DIMENSION(100), INTENT(OUT) ::  y,srk
  REAL, INTENT(OUT) :: srki0,srkin
!chr    represent snow or ice case
!nl     the number of vertical layer
!nh
!hsi    thickness of the medium(snow or ice
!tsfc   surface temp. [K]
!tbn    bottom temp.  [K]
!vrs    heat resource item
!y      result(Tsnow or Tice)
!srk
  INTEGER :: n,i,j,k,i1,j1,k1,if1,nh
  REAL, DIMENSION(100) :: a,reta,c,p,f1,f2,df,x,si,hd
  REAL, DIMENSION(100,100) :: al,br,l,r
  REAL :: eta0,etanv,tx,dh,s,w1,w2,w3,z1,s0,sn,ci,cs,di,st

  s0 = 3.0
  sn = 6.2
  ci = 2093.0
  cs = 2093.0
  di = 915.0
  st = 0.5

  n = nl-1

  IF(chr.EQ.'s') THEN
  eta0  = rks
  etanv = rks

  ELSEIF(chr.EQ.'i') THEN

! eta0  = rk0 + bta*s0 / tsfc
! Those are for sea-ice
! etanv = rk0 + bta*sn / tbn
! Those are for sea-ice
! if(eta0.lt.1.43)  eta0  =1.43  ! a threshold value
! if(etanv.lt.1.43) etanv =1.43  ! a threshold value


  eta0  = rk0
  etanv = rk0

  ENDIF

  dh = hsi / nl
  tx = dt / (dh*dh)     ! ratio of time step to vertical ice length

!     Parameterization of salinity effect
!----   si=14.24-19.39*hi, hi < 0.57 m
!----   si=3.2             hi > 0.57 m

  DO  k =1,n
  hd(k)=dh*k
  IF(hd(k).LT.0.57) THEN
  si(k)=14.24-19.39*hd(k)
  ELSEIF(hd(k).ge.0.57) THEN
  si(k)=3.2
  ENDIF
  END DO

  DO k=1,n
  s=(s0+sn)*0.5

  IF(chr.EQ.'s') THEN

  reta(k) = rks
  c(k)    = ds*cs

  ELSEIF(chr.EQ.'i') THEN

  reta(k) = rk0           ! no 'si','Tice' effect
  c(k)    = di*ci

  ENDIF
  ENDDO

  a(1)  = (eta0+reta(1))/2.0
  a(nl) = (etanv+reta(n))/2.0

  IF(chr.EQ.'s') THEN
  srk(1)=a(1)     ! rks
  srk(nl)=a(nl)   ! rks

  ELSEIF(chr.EQ.'i') THEN
!----   This is one way to define Ki, i.e. (Ki1+Ki2)/2
  srk(1)=a(1)
  srk(nl)=a(nl)
  ENDIF

  DO i=2,n
  a(i)=(reta(i-1)+reta(i))/2.0    

  IF(chr.EQ.'s') THEN
  srk(i) = rks                            
                
  ELSEIF(chr.EQ.'i') THEN

  srk(i)=a(i)

  IF(srk(i).LT.1.43) srk(i)=1.43

  ENDIF
  ENDDO

  IF(chr.EQ.'s') THEN
  srki0 =rks
  srkin =rks

  ELSEIF(chr.EQ.'i') THEN
  srki0=a(1)
  srkin=a(nl)

  IF(srki0.LT.1.43) srki0=1.43
  IF(srkin.LT.1.43) srkin=1.43

  ENDIF

! Create the matrixes in each side of the heat conduction equation
!               al*temp(k+1)=br*temp(k)+b
!   where b=f1+f2 see below.

  DO i=1,n
  DO j=1,n
    IF(i.EQ.j) THEN
    al(i,j)=c(i)+tx*st*(a(j)+a(j+1))
    br(i,j)=c(i)-tx*(1-st)*(a(j)+a(j+1))
    ELSEIF((i-j).EQ.1) THEN   ! lower part
    al(i,j)=-tx*st*a(i)
    br(i,j)=tx*(1-st)*a(i)
    ELSEIF((i-j).EQ.-1) THEN  ! upper part
    al(i,j)=-tx*st*a(j)
    br(i,j)=tx*(1-st)*a(j)
    ELSEIF(abs(i-j).EQ.1) THEN
    al(i,j)=0.0
    br(i,j)=0.0
    ENDIF
   ENDDO
   ENDDO

! Separate matrix al to l and r where l is lower triangular matrix
! and r is upper triangular matrix.
! al is the matrix on the left side, e.g.  al*temp(k+1)

  DO i=1,n
  r(i,i)=1
    DO k=1,n
    w1=0.0
      DO j=1,i-1
      w1=w1+r(i,j)*l(j,k)
      ENDDO
    l(i,k)=al(i,k)-w1
    ENDDO

    DO if1=i+1,n
    z1=0.0
      DO nh=1,i-1
      z1=z1+r(if1,nh)*l(nh,i)
      ENDDO
    r(if1,i)=(al(if1,i)-z1)/l(i,i)
    ENDDO
  ENDDO

! Calculate first item of the right side e.g.  br*temp(k)=p

  DO j=1,n
  p(j)=0.0
    DO k=1,n
    IF(chr.EQ.'s') THEN
    p(j)=p(j)+br(j,k)*tsnow(k)
    ELSEIF(chr.EQ.'i') THEN
    p(j)=p(j)+br(j,k)*tice(k)
    ENDIF
    ENDDO
  ENDDO

! Calculate the right side of the equation as a whole
!       p  is the item equal p=br*temp(k)  ! first item
!       f1 is the boundary item            ! second item        
!       f2 is the radiation item           ! second item

!       df is the total right item as a whole

  f1(1)=tsfc*tx*a(1)
  f1(n)=tbn*tx*a(nl)

  DO i=2,n-1
  f1(i)=0.0
  ENDDO

  DO i=1,n
  f2(i)=dt*vrs(i)
  ENDDO

  DO k=1,n
  df(k)=p(k)+f1(k)+f2(k)
  ENDDO


! Calculate the x ------>  l*x=df

  x(1)=df(1)/l(1,1)
  DO i=2,n
  w2=0.0
    DO k=1,i-1
    w2=w2+l(k,i)*x(k)
    ENDDO
  x(i)=(df(i)-w2)/l(i,i)
  ENDDO

! Calculate the y ------>  r*y=x

  y(n)=x(n)/r(n,n)
  DO k=n-1,1,-1
  w3=0.0
    DO j=k+1,n
    w3=w3+r(j,k)*y(j)
    ENDDO
  y(k)=x(k)-w3
  ENDDO

  RETURN
  END SUBROUTINE nume1

  SUBROUTINE nume2(nl,hsi,tsfc,tbn,tsnow,tice,y,vrs,chr,rk0,dt,srk,srki0,srkin)
!----   Calculate in-snow temperature with consideration of snowlayering
!----   effect: fit to CHINARE2003 expedition:
!----   Initially, there was a 22-24cm old snow with an average density of
!----   about 350(assumed). According general knowledge there should be laye
! ring inside
!----   snow, but it is difficult to summarize in the code below,
!----   therefore, we only define 2 layers of snow: new snowfall and old 22cm
!----   initial snow thickness, accordingly:
!       new snow fall  : external forcing [ds=225kgm-3, so rks=    ]
!       old snow       : 22cm [ds=450kgm-3, so rks= 0.347]
!----   This subroutine is only for snow part, i.e. chr.eq.'s'
  IMPLICIT NONE
  CHARACTER(1), INTENT(IN) :: chr
  INTEGER, INTENT(IN) :: nl !,nh
  REAL, INTENT(IN)  ::  hsi,tsfc,tbn,rk0,dt
  REAL, DIMENSION(100), INTENT(IN) :: vrs,tsnow,tice
  REAL, DIMENSION(100), INTENT(OUT) :: y,srk
  REAL, INTENT(OUT) :: srki0,srkin
!chr              represent snow or ice case
!nl               the number of vertical layer
!nh
!hsi              thickness of the medium(snow or ice)
!tsfc             surface temp. [K]
!tbn              bottom temp.  [K]
!vrs              heat resource item
!y                result(Tsnow or Tice)
!srk
  INTEGER :: n,i,j,k,i1,j1,k1,if1,nh
  REAL, DIMENSION(100) :: a,reta,c,p,f1,f2,df,x,si,hd
  REAL, DIMENSION(100,100) :: al,br,l,r
  REAL :: eta0,etanv,tx,dh,s,w1,w2,w3,z1,bta,s0,sn,ci,cs,di,st

  bta = 0.12
  s0  = 3.0
  sn  = 6.2
  ci  = 2093.0
  cs  = 2093.0
  di  = 915.0
  st  = 0.5

  n = nl-1

  IF(chr.EQ.'s') THEN

  eta0  = 10**(0.002650*380-1.652)      ! ds= 380
  etanv = 10**(0.002650*550-1.652)      ! ds= 550(snow bottom, assumed)

!----   This means that if snow thickness is less than 0.22(initial)
!       so the ds=350 and eta0==rks= 0.35
!       etanv always indicate the snow/ice interface(snow bottom)
    IF(hsi.LE.0.22) THEN  ! 0.22: snow thickness before snow fall
    eta0  = 10**(0.002650*450-1.652)      ! ds= 450 intermediate layer
    etanv = 10**(0.002650*550-1.652)      ! ds= 550 snow bottom (assumed)
    ENDIF

!----   This means that if snow is further melting away, until it belows
!       15cm assumed to be very old snow, so ds=550-600 and eta0==rks =

    IF(hsi.LE.0.15) THEN
    eta0  = 10**(0.002650*550-1.652)      ! ds= 550 snow bottom (assumed)
    etanv = 10**(0.002650*600-1.652)      ! ds= 600 snow bottom (assumed)
    ENDIF

  ELSEIF(chr.EQ.'i') THEN
  eta0  = rk0 + bta*s0 / tsfc
  etanv = rk0 + bta*sn / tbn
  IF(eta0.LT.1.43)  eta0  =1.43
  IF(etanv.LT.1.43) etanv =1.43

  ENDIF

  dh=hsi / nl
  tx=dt / (dh*dh)     ! ratio of time step to vertical ice length

!----   si=14.24-19.39*hi   hi <0.57m
!----   si=3.2              hi >0.57m
  DO k=1,n
  hd(k)=dh*k
  IF(hd(k).LT.0.57) THEN
  si(k)=14.24-19.39*hd(k)
  ELSEIF(hd(k).ge.0.57) THEN
  si(k)=3.2
  ENDIF
  END DO

  DO k=1,n
  s=(s0+sn)/2.0
  hd(k)=dh*k

  IF(chr.EQ.'s') THEN
        
    IF(hd(k).LE.(hsi-0.22)) THEN
    reta(k) =10**(0.002650*380-1.652)  ! ds=380
    c(k)    = 380*cs

    ELSEIF(hd(k).GE.(hsi-0.22).AND.hd(k).LE.(hsi-0.15)) THEN
    reta(k) = 10**(0.002650*450-1.652) ! ds=450
    c(k)    = 450*cs

    ELSEIF(hd(k).gt.(hsi-0.15)) THEN
    reta(k) = 10**(0.002650*550-1.652) ! ds=550
    c(k)    = 550*cs
    ENDIF

!----   For the constant heat conductivity and density of snow

  ELSEIF(chr.EQ.'i') THEN

  reta(k) = rk0
  c(k)    = di*ci

  ENDIF

  ENDDO

  a(1)  = (eta0+reta(1))/2.0
  a(nl) = (etanv+reta(n))/2.0


  IF(chr.EQ.'s') THEN
  srk(1)=a(1)     ! rks
  srk(nl)=a(nl)   ! rks

  ELSEIF(chr.EQ.'i') THEN

!----   This is one way to define Ki, i.e. (Ki1+Ki2)/2 
  srk(1)=a(1)
  srk(nl)=a(nl)

  ENDIF

  do i=2,n
  a(i)=(reta(i-1)+reta(i))/2.0    

  if(chr.eq.'s')then
  srk(i) = a(i)                 
                
  elseif(chr.eq.'i') then

  srk(i)=a(i)

  if(srk(i).lt.1.43) srk(i)=1.43
  endif

  enddo

  if(chr.eq.'s')then
  srki0 =a(1)     ! rks
  srkin =a(nl)    ! rks

  elseif(chr.eq.'i'.or.chr.eq.'I')then

  srki0 =a(1)     !
  srkin =a(nl)    !

  if(srki0.lt.1.43) srki0=1.43
  if(srkin.lt.1.43) srkin=1.43

  endif

! Create the matrixes in each side of the heat conduction equation
!               al*temp(k+1)=br*temp(k)+b
! where b=f1+f2 see below.

  do i=1,n
  do j=1,n
  if(i.eq.j) then
  al(i,j)=c(i)+tx*st*(a(j)+a(j+1))
  br(i,j)=c(i)-tx*(1-st)*(a(j)+a(j+1))
  elseif((i-j).eq.1) then   ! lower part
  al(i,j)=-tx*st*a(i)
  br(i,j)=tx*(1-st)*a(i)
  elseif((i-j).eq.-1) then  ! upper part
  al(i,j)=-tx*st*a(j)
  br(i,j)=tx*(1-st)*a(j)
  elseif(abs(i-j).gt.1) then
  al(i,j)=0.0
  br(i,j)=0.0
  endif
  enddo
  enddo

! Separate matrix al to l and r where l is lower triangular matrix
! and r is upper triangular matrix.
! al is the matrix on the left side, e.g.  al*temp(k+1)

  do i=1,n
  r(i,i)=1
    do k=1,n
    w1=0.0
      do j=1,i-1
      w1=w1+r(i,j)*l(j,k)
      enddo
    l(i,k)=al(i,k)-w1
    enddo

    do if1=i+1,n
    z1=0.0
      do nh=1,i-1
      z1=z1+r(if1,nh)*l(nh,i)
      enddo
      r(if1,i)=(al(if1,i)-z1)/l(i,i)
    enddo
  enddo

! Calculate first item of the right side e.g.  br*temp(k)=p

  do j=1,n
  p(j)=0.0
    do k=1,n
    if(chr.eq.'s')then
    p(j)=p(j)+br(j,k)*tsnow(k)
    elseif(chr.eq.'i')then
    p(j)=p(j)+br(j,k)*tice(k)
    endif
  enddo
  enddo

! Calculate the right side of the equation as a whole
!       p  is the item equal p=br*temp(k)  ! first item
!       f1 is the boundary item            ! second item        
!       f2 is the radiation item           ! second item

!       df is the total right item as a whole

  f1(1)=tsfc*tx*a(1)
  f1(n)=tbn*tx*a(nl)

  do i=2,n-1
  f1(i)=0.0
  enddo

  do i=1,n
  f2(i)=dt*vrs(i)
  enddo

  do k=1,n
  df(k)=p(k)+f1(k)+f2(k)
  enddo

! Calculate the x ------>  l*x=df

  x(1)=df(1)/l(1,1)
  do i=2,n
  w2=0.0
    do k=1,i-1
    w2=w2+l(k,i)*x(k)
    enddo
  x(i)=(df(i)-w2)/l(i,i)
  enddo

! Calculate the y ------>  r*y=x

  y(n)=x(n)/r(n,n)
  do k=n-1,1,-1
  w3=0.0
    do j=k+1,n
    w3=w3+r(j,k)*y(j)
    enddo
  y(k)=x(k)-w3
  enddo

  return
  END SUBROUTINE nume2




  SUBROUTINE HIGHTSI (IILOC, JJLOC, SEAICE_ALBEDO_OPT, SEAICE_ALBEDO_DEFAULT, &    !C
       &                  SEAICE_SNOWDEPTH_OPT, SEAICE_SNOWDEPTH_MAX,      &    !C
       &                  SEAICE_SNOWDEPTH_MIN,FFROZP,DT,ZLVL,               &    !C
       &                  SITHICK,SFCSPD,rki,tsnow,tice,T1,ti,         &
       &                  nvs,nvi,ETA,SHEAT,         &
       &                  LWDN,SWDN,SFCPRS,PRCP,SFCTMP,Q2,               &    !F
       &                  TH2,Q2SAT,DQSDT2,Z0BRD,Z0,CH,RIBB,               &    !I
       &                  TBOT, SNOWH, SNEQV, ALBEDO,rks,dsnow,kkkk, &
       &                  ETA_KINEMATIC,FDOWN,ESNOW,DEW,ETP,SSOIL,FLX1, &
       &                  FLX2,FLX3,SNOMLT,SNCOVR,RUNOFF1,Q1     )

! ----------------------------------------------------------------------
! SUBROUTINE SFLX_SEAICE
! ----------------------------------------------------------------------
! SUB-DRIVER FOR "Noah LSM" FAMILY OF PHYSICS SUBROUTINES FOR A SEA-ICE
! LAND-SURFACE MODEL TO UPDATE ICE TEMPERATURE, SKIN TEMPERATURE,
! SNOWPACK WATER CONTENT, SNOWDEPTH, AND ALL TERMS OF THE SURFACE ENERGY
! BALANCE (EXCLUDING INPUT ATMOSPHERIC FORCINGS OF DOWNWARD RADIATION
! AND PRECIP)
! ----------------------------------------------------------------------
! SFLX_SEAICE ARGUMENT LIST KEY:
! ----------------------------------------------------------------------
!  C  CONFIGURATION INFORMATION
!  F  FORCING DATA
!  I  OTHER (INPUT) FORCING DATA
!  S  SURFACE CHARACTERISTICS
!  H  HISTORY (STATE) VARIABLES
!  O  OUTPUT VARIABLES
!  D  DIAGNOSTIC OUTPUT
! ----------------------------------------------------------------------
! 1. CONFIGURATION INFORMATION (C):
! ----------------------------------------------------------------------
!   DT         TIMESTEP (SEC) (DT SHOULD NOT EXCEED 3600 SECS, RECOMMEND
!                1800 SECS OR LESS)
!   ZLVL       HEIGHT (M) ABOVE GROUND OF ATMOSPHERIC FORCING VARIABLES
!   NSOIL      NUMBER OF SOIL LAYERS (AT LEAST 2, AND NOT GREATER THAN
!                PARAMETER NSOLD SET BELOW)
! ----------------------------------------------------------------------
! 3. FORCING DATA (F):
! ----------------------------------------------------------------------
!   LWDN       LW DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET LONGWAVE)
!   SOLNET     NET DOWNWARD SOLAR RADIATION ((W M-2; POSITIVE)
!   SFCPRS     PRESSURE AT HEIGHT ZLVL ABOVE GROUND (PASCALS)
!   PRCP       PRECIP RATE (KG M-2 S-1) (NOTE, THIS IS A RATE)
!   SFCTMP     AIR TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND
!   TH2        AIR POTENTIAL TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND
!   Q2         MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1)
!   FFROZP     FRACTION OF FROZEN PRECIPITATION
! ----------------------------------------------------------------------
! 4. OTHER FORCING (INPUT) DATA (I):
! ----------------------------------------------------------------------
!   Q2SAT      SAT SPECIFIC HUMIDITY AT HEIGHT ZLVL ABOVE GROUND (KG KG-1)
!   DQSDT2     SLOPE OF SAT SPECIFIC HUMIDITY CURVE AT T=SFCTMP
!                (KG KG-1 K-1)
! ----------------------------------------------------------------------
! 5. CANOPY/SOIL CHARACTERISTICS (S):
! ----------------------------------------------------------------------
!   SNOALB     UPPER BOUND ON MAXIMUM ALBEDO OVER DEEP SNOW (E.G. FROM
!                ROBINSON AND KUKLA, 1985, J. CLIM. & APPL. METEOR.)
!   TBOT       BOTTOM SOIL TEMPERATURE (LOCAL YEARLY-MEAN SFC AIR
!                TEMPERATURE)
!   Z0BRD      Background fixed roughness length (M)
!   Z0         Time varying roughness length (M) as function of snow depth
!
!   EMISSI     Surface emissivity (between 0 and 1)
! ----------------------------------------------------------------------
! 6. HISTORY (STATE) VARIABLES (H):
! ----------------------------------------------------------------------
!  T1          GROUND/CANOPY/SNOWPACK) EFFECTIVE SKIN TEMPERATURE (K)
!  STC(NSOIL)  SOIL TEMP (K)
!  SNOWH       ACTUAL SNOW DEPTH (M)
!  SNEQV       LIQUID WATER-EQUIVALENT SNOW DEPTH (M)
!                NOTE: SNOW DENSITY = SNEQV/SNOWH
!  ALBEDO      SURFACE ALBEDO
!  CH          SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE
        !                (M S-1); NOTE: CH IS TECHNICALLY A CONDUCTANCE SINCE
!                IT HAS BEEN MULTIPLIED BY WIND SPEED.
! ----------------------------------------------------------------------
! 7. OUTPUT (O):
! ----------------------------------------------------------------------
! OUTPUT VARIABLES NECESSARY FOR A COUPLED NWP MODEL.  FOR THIS APPLICATION,
! THE REMAINING OUTPUT/DIAGNOSTIC/PARAMETER BLOCKS BELOW ARE NOT
! NECESSARY.  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES.
!   ETA        ACTUAL LATENT HEAT FLUX (W m-2: NEGATIVE, IF UP FROM
                !              SURFACE)
!  ETA_KINEMATIC actual latent heat flux in Kg m-2 s-1
!   SHEAT      SENSIBLE HEAT FLUX (W M-2: NEGATIVE, IF UPWARD FROM
!              SURFACE)
!   FDOWN      Radiation forcing at the surface (W m-2) = SOLDN*(1-alb)+LWDN
! ----------------------------------------------------------------------
!   ESNOW      SUBLIMATION FROM (OR DEPOSITION TO IF <0) SNOWPACK (W m-2)
!   DEW        DEWFALL (OR FROSTFALL FOR T<273.15) (M)
! ----------------------------------------------------------------------
!   ETP        POTENTIAL EVAPORATION (W m-2)
!   SSOIL      SOIL HEAT FLUX (W M-2: NEGATIVE IF DOWNWARD FROM SURFACE)
! ----------------------------------------------------------------------
!   FLX1       PRECIP-SNOW SFC (W M-2)
!   FLX2       FREEZING RAIN LATENT HEAT FLUX (W M-2)
!   FLX3       PHASE-CHANGE HEAT FLUX FROM SNOWMELT (W M-2)
! ----------------------------------------------------------------------
!   SNOMLT     SNOW MELT (M) (WATER EQUIVALENT)
!   SNCOVR     FRACTIONAL SNOW COVER (UNITLESS FRACTION, 0-1)
! ----------------------------------------------------------------------
!   RUNOFF1    SURFACE RUNOFF (M S-1), NOT INFILTRATING THE SURFACE
! ----------------------------------------------------------------------
! 8. DIAGNOSTIC OUTPUT (D):
! ----------------------------------------------------------------------
!   Q1         Effective mixing ratio at surface (kg kg-1), used for
!              diagnosing the mixing ratio at 2 meter for coupled model
!  Documentation SNOABL2 ?????
!  What categories of arguments do these variables fall into ????
!  Documentation for RIBB ?????
!  What category of argument does RIBB fall into ?????
! ----------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------------------
      integer, intent(in) :: iiloc, jjloc
      INTEGER, INTENT(IN) :: SEAICE_ALBEDO_OPT
      REAL,    INTENT(IN) :: SEAICE_ALBEDO_DEFAULT
      INTEGER, INTENT(IN) :: SEAICE_SNOWDEPTH_OPT
      REAL,    INTENT(IN) :: SEAICE_SNOWDEPTH_MAX
      REAL,    INTENT(IN) :: SEAICE_SNOWDEPTH_MIN

      LOGICAL            ::  FRZGRA, SNOWNG


      REAL, INTENT(IN)   :: DT,DQSDT2,LWDN,PRCP,            &
                            Q2,Q2SAT,SFCPRS,SFCTMP,SFCSPD,   &
                            SWDN,TBOT,TH2,ZLVL,&
                            FFROZP
      REAL, INTENT(INOUT)  :: ALBEDO,RIBB
      REAL, INTENT(INOUT):: CH,SNEQV,SNCOVR,SNOWH,T1,Z0BRD
      REAL               :: EMISSI
      REAL, INTENT(INOUT)   :: SITHICK

      REAL,INTENT(OUT)   :: ETA_KINEMATIC,DEW,ESNOW,ETA,                    &
                            ETP,FLX1,FLX2,FLX3,SHEAT,RUNOFF1,               &
                            SSOIL,                                          &
                            SNOMLT,                                         &
                            FDOWN,Q1,Z0

      REAL :: SNDENS,SNCOND,SN_NEW,RCH,RR,T2V,T24,RHO,SOLNET

!yyao 2014/11/18 ccsm scheme
      REAL :: TT,AHMAX,FHTAN,FH,ALBO,alb_ice_v,alb_ice_i,&
              FT,alb_snow_v,alb_snow_i,asnow,alb_v,alb_i
!yyao 2014/12 hsi
      REAL, INTENT(INOUT) :: rki,rks,dsnow
      INTEGER, INTENT(INOUT) :: kkkk
      INTEGER, INTENT(IN) :: nvs,nvi
      REAL, INTENT(OUT) :: ti
      REAL :: dp,tmelt,dih,dsh,qsum,extink,ql,sdth,es,qs,e,&
              ch1,ce
      REAL :: fcontsi,sumn
      REAL :: srki0,srkin
      REAL, DIMENSION(100), INTENT(INOUT) :: tice,tsnow
      INTEGER :: ns,ni,md
      CHARACTER(1) :: sur
      REAL :: qii,qsi
      REAL :: ts0,tin1,hhhs,hhhi,qsumsubi
      INTEGER :: jsnow,ks,js
      REAL, DIMENSION(100) :: ys,yi, rss,rsi, srks,srki,qini,qins
      INTEGER,DIMENSION(100) :: iarrs,iarri
      REAL :: sumqins,qsumsubs,ffup,ffdn
      REAL ::tb,rk0,ds,dwater,ws
      REAL :: tin2,tmeltice,tu,td,sumqini
      INTEGER :: jice,ki,ji
      REAL :: dsnew,di,dsi,dswater,c1,c2
      REAL :: dth
      REAL :: fw   !ocean heat flux

! ----------------------------------------------------------------------
! DECLARATIONS - PARAMETERS
! ----------------------------------------------------------------------

      REAL, PARAMETER :: LVH2O = 2.501E+6
      REAL, PARAMETER :: LSUBS = 2.83E+6
      REAL, PARAMETER :: R = 287.04

      iloc = iiloc
      jloc = jjloc
! ----------------------------------------------------------------------
!   INITIALIZATION
! ----------------------------------------------------------------------

      RUNOFF1 = 0.0
      SNOMLT = 0.0
      EMISSI = 0.98
      Z0BRD = 0.001 
      SNCOVR = 1.0

! ----------------------------------------------------------------------
!  INITIALIZE PRECIPITATION LOGICALS.
! ----------------------------------------------------------------------

      SNOWNG = .FALSE.
      FRZGRA = .FALSE.

! ----------------------------------------------------------------------
! OVER SEA-ICE, IF S.W.E. (SNEQV) BELOW THRESHOLD LOWER
! BOUND (0.01 M FOR SEA-ICE, 0.10 M FOR GLACIAL-ICE), THEN SET AT LOWER
! BOUND
! ----------------------------------------------------------------------
! FOR SEA-ICE CASE, ASSIGN DEFAULT WATER-EQUIV SNOW ON TOP
! ----------------------------------------------------------------------



      SNDENS = SNEQV / SNOWH
      SNOWH = MAX ( 0.001 , MIN ( SNOWH , 10.0 ) )
      SNEQV = SNOWH * SNDENS 


! ----------------------------------------------------------------------
! IF INPUT SNOWPACK IS NONZERO, THEN COMPUTE SNOW DENSITY "SNDENS" AND
! SNOW THERMAL CONDUCTIVITY "SNCOND"
! ----------------------------------------------------------------------

      IF(SNDENS > 1.0) THEN
         CALL wrf_error_fatal ( 'Physical snow depth is less than snow water equiv.' )
      ENDIF
      CALL CSNOW (SNCOND,SNDENS)

! ----------------------------------------------------------------------
! DETERMINE IF IT'S PRECIPITATING AND WHAT KIND OF PRECIP IT IS.
! IF IT'S PRCPING AND THE AIR TEMP IS COLDER THAN 0 C, IT'S SNOWING!
! IF IT'S PRCPING AND THE AIR TEMP IS WARMER THAN 0 C, BUT THE GRND
! TEMP IS COLDER THAN 0 C, FREEZING RAIN IS PRESUMED TO BE FALLING.
! ----------------------------------------------------------------------

      IF (PRCP > 0.0) THEN
! snow defined when fraction of frozen precip (FFROZP) > 0.5,
! passed in from model microphysics.
         IF (FFROZP .GT. 0.5) THEN
            SNOWNG = .TRUE.
         ELSE
            IF (T1 <= TFREEZ) FRZGRA = .TRUE.
         END IF
      END IF


! ----------------------------------------------------------------------
! IF EITHER PRCP FLAG IS SET, DETERMINE NEW SNOWFALL (CONVERTING PRCP
! RATE FROM KG M-2 S-1 TO A LIQUID EQUIV SNOW DEPTH IN METERS) AND ADD
! IT TO THE EXISTING SNOWPACK.
! ----------------------------------------------------------------------

      IF ( SNOWNG .OR. FRZGRA ) THEN
         SN_NEW = PRCP * DT * 0.001
         SNEQV = SNEQV + SN_NEW

! ----------------------------------------------------------------------
! UPDATE SNOW DENSITY BASED ON NEW SNOWFALL, USING OLD AND NEW SNOW.
! UPDATE SNOW THERMAL CONDUCTIVITY
! ----------------------------------------------------------------------

         CALL SNOW_NEW ( SFCTMP , SN_NEW , SNOWH , SNDENS )
         SNEQV = SNOWH*SNDENS

         CALL CSNOW (SNCOND,SNDENS)
      ENDIF


! ----------------------------------------------------------------------
! ALBEDO OF SEA ICE
! ----------------------------------------------------------------------


!thickness dependence
      AHMAX = 0.5
      FHTAN = ATAN(AHMAX*4.0)
      FH = MIN(ATAN(SITHICK*4.0)/FHTAN,1.0)
      ALBO = 0.1*(1.0-FH)

      alb_ice_v = 0.78*FH+ALBO
      alb_ice_i = 0.36*FH+ALBO

!temperature dependence
      TT = MIN(T1,273.15)
      FT = MIN((273.15-TT)-1.0,0.0)
      alb_ice_v = alb_ice_v+0.075*FT   
      alb_ice_i = alb_ice_i+0.075*FT
      IF (alb_ice_v.GT.1) THEN
         CALL wrf_error_fatal ( ' alb_ice_v > 1 ' )
      END IF
      IF (alb_ice_i.GT.1) THEN
         CALL wrf_error_fatal ( ' alb_ice_i > 1 ' )
      END IF

!avoid negative albedos
      alb_ice_v = MAX(alb_ice_v,0.1)
      alb_ice_i = MAX(alb_ice_i,0.1)

      IF (SNOWH.GT.1.0e-11) THEN
         alb_snow_v = 0.98+0.1*FT
         alb_snow_i = 0.7+0.15*FT
         IF (alb_snow_v.GT.1) THEN
            CALL wrf_error_fatal ( ' alb_snow_v > 1 ' )
         END IF
         IF (alb_snow_i.GT.1) THEN
            CALL wrf_error_fatal ( ' alb_snow_i > 1 ' )
         END IF
         IF (alb_snow_v.LT.0) THEN
            CALL wrf_error_fatal ( ' alb_snow_v < 0 ' )
         END IF
         IF (alb_snow_i.LT.0) THEN
            CALL wrf_error_fatal ( ' alb_snow_i < 0 ' )
         END IF

         asnow = SNOWH/(SNOWH+0.02)

         alb_v = alb_ice_v*(1.0-asnow)+ &
                 alb_snow_v*asnow
         alb_i = alb_ice_i*(1.0-asnow)+ &
                 alb_snow_i*asnow

      ELSE
         alb_v  = alb_ice_v
         alb_i  = alb_ice_i

      END IF

      ALBEDO    = alb_v*0.636+alb_i*0.364
      IF (ALBEDO.GT.1) THEN
         CALL wrf_error_fatal ( ' ALBEDO > 1 ' )
      END IF
      IF (ALBEDO.LT.0) THEN
         CALL wrf_error_fatal ( ' ALBEDO < 0 ' )
      END IF
      SOLNET    = SWDN*(1.0-ALBEDO)

! ----------------------------------------------------------------------
! CALCULATE TOTAL DOWNWARD RADIATION (SOLAR PLUS LONGWAVE) NEEDED IN
! PENMAN EP SUBROUTINE THAT FOLLOWS
! ----------------------------------------------------------------------
      FDOWN =  SOLNET + LWDN
! ----------------------------------------------------------------------
! CALC VIRTUAL TEMPS AND VIRTUAL POTENTIAL TEMPS NEEDED BY SUBROUTINES
! PENMAN.
! ----------------------------------------------------------------------
      T2V = SFCTMP * (1.0+ 0.61 * Q2 )
      T24 = SFCTMP * SFCTMP * SFCTMP * SFCTMP
      RHO = SFCPRS / ( RD * T2V )
      ! RCH = RHO * CP * CH
      RCH = RHO * 1004.6 * CH  ! CP is defined different in subroutine PENMAN.
                               ! Pulling this computation out of PENMAN changed
                               ! the results.  So I'm hard-coding the PENMAN
                               ! value here, but perhaps this should go back
                               ! into PENMAN for now.

!-----------------
! HIGHTSI
!-----------------

! surface heat balance
      qii  = 302e6       ! heat capacity of ice [J/m3]
      qsi  = 110e6       ! heat capacity of snow [J/m3]
      dp   = SFCPRS/(287.05*SFCTMP)
      md   = 3
      extink  = 25.0
      rk0  = 1.83
      dwater    = 1000.0
      dsnew = 320.0
      di    = 915.0
      dsi   = 820.0
      dswater = 1028.0
      c1  = 5.0
      c2  = 21.0
      tmeltice = -0.1355
      fw  = 0.50
      iarrs = 0
      iarri = 0
      dih = SITHICK/nvi
      dsh = SNOWH/nvs
      ns  = nvs-1
      ni  = nvi-1

      IF ( SNOWH .LT. 0.005) THEN
         kkkk = 1
      ELSEIF ( SNOWH .GE. 0.01) THEN
         kkkk = kkkk+1
      ENDIF

      kkkk = MIN(1000,kkkk)

      ds = SNDENS*1000.0
!     ds = dsnow*1000.0
      rks= SNCOND
!     rks=10**(0.002650*ds-1.652)+0.2
!     rks=20**(0.002650*ds-1.652)+0.2

      IF ( SNOWH .LE. 0.005) THEN
         sur = 'i'
         tmelt = 273.15 - 0.1355
         SSOIL = rki/dih*(T1-273.15-tice(1))
      ELSEIF ( SNOWH .GT. 0.005) THEN
         sur = 's'
         tmelt = 273.15 - 0.0
         IF ( SNOWH .LT. 0.01 ) THEN
            SSOIL = rks/SNOWH*(T1-273.15-ti)
         ELSE
            SSOIL = rks/dsh*(T1-273.15-tsnow(1))
         ENDIF
      ENDIF

! ----------------------------------------------------------------------
! CALL PENMAN SUBROUTINE TO CALCULATE POTENTIAL EVAPORATION (ETP), AND
! OTHER PARTIAL PRODUCTS AND SUMS FOR LATER CALCULATIONS.
! ----------------------------------------------------------------------

      CALL PENMAN (SFCTMP,SFCPRS,CH,TH2,PRCP,FDOWN,T24,SSOIL,     &
           Q2,Q2SAT,ETP,RCH,RR,SNOWNG,FRZGRA,                     &
           DQSDT2,FLX2,EMISSI,T1)

      DEW   = 0.0
      ESNOW = 0.0
      FLX1  = 0.0
      FLX3  = 0.0

      IF (SNOWNG) THEN
         FLX1 = CPICE * PRCP * (T1- SFCTMP)
      ELSE
         IF (PRCP >  0.0) FLX1 = CPH2O * PRCP * (T1- SFCTMP)
      END IF

      IF (ETP <= 0.0) THEN
         IF ( ( RIBB >= 0.1 ) .AND. ( FDOWN > 150.0 ) ) THEN
            ETP=(MIN(ETP*(1.0-RIBB),0.)*SNCOVR/0.980 + ETP*(0.980-SNCOVR))/0.980
         ENDIF
         DEW = -ETP * 0.001
      ELSE
         ESNOW  = ETP*SNCOVR
      END IF

      ETA_KINEMATIC =  ESNOW


!Calculate solar radiation at surface and penetrating in snow and ice
      CALL qsall1(SOLNET,qsum,extink,0.0,dih,dsh,nvi,nvs,md,SNOWH,kkkk,&
                    qini,qins,rsi,rss)

      ql = LWDN - EMISSI*5.67e-8*T1**4.0
!Calculate surface temperature
      CALL ts1(T1,SFCTMP,Q2,SFCSPD,ZLVL,dp,qsum,ql,SNOWH,kkkk,&
               dih,dsh,ns,tice,tsnow,rki,rks)
      IF (T1.LT.100.0) THEN
         CALL wrf_error_fatal('ts1 too small')
!     T1 = SFCTMP
      END IF

!If the surface temperature is higher than the tmelt,
!then calculate melting at surface
      IF ( T1 .GT. tmelt) THEN
         T1 = tmelt
         IF (sur.eq.'i') rki=1.73
!rki: heat conductivity of sea ice at surface layer ONLY

         CALL tint(T1,SNOWH,kkkk,dih,dsh,ns,tice,tsnow,rki,rks,ti)
         ql = LWDN - EMISSI*5.67e-8*T1**4.0
!        CALL dhs(T1,ti,SFCTMP,Q2,SFCSPD,ZLVL,DT,dp,ql,qsum,SNOWH,SITHICK,&
!                rks,rki,dsh,dih,tsnow,tice,nvi,nvs,kkkk,sdth0)
         CALL wvap(T1,Q2,es,qs,e)
         CALL coe(T1,SFCTMP,SFCSPD,ZLVL,qs,Q2,ch1,ce)
         SHEAT  = -dp*CP*ch1*(T1-SFCTMP)*SFCSPD
         ETA    = -dp*ce*(qs-Q2)*SFCSPD*rl1(T1)

         IF (sur.EQ.'i') THEN
            IF (tice(1).LT.-0.1355) THEN
               fcontsi = -1*rki/dih*((tice(1)+273.15)-T1)
               sumn    = qsum + SHEAT + ETA +ql+fcontsi
               sdth    = -DT/qii*sumn
               SITHICK = SITHICK+sdth
               IF (ABS(sdth).gt.0.1) &
               CALL wrf_error_fatal('Ice melt too much at surface')
               SITHICK  = MAX(SITHICK,0.1)
               dih      = SITHICK/nvi
            ENDIF

         ELSEIF (sur.EQ.'s') THEN
            IF (SNOWH.GT.0.005.AND.SNOWH.LT.0.01) THEN
               fcontsi = rks/SNOWH*(T1-ti)
               sumn    = qsum + SHEAT + ETA +ql+fcontsi
               sdth    = -DT/qsi*sumn
               sdth    = MIN(0.0,sdth)
               SNOWH   = SNOWH+sdth
               SNOMLT  = SNOMLT-sdth
               SNOWH   = MAX(SNOWH,0.001)
               SNEQV   = SNOWH*SNDENS
            ELSEIF (SNOWH.GE.0.01) THEN
               IF (kkkk.LE.15) THEN
                  fcontsi = rks/SNOWH*(T1-ti)
                  sumn    = qsum + SHEAT + ETA +ql+fcontsi
                  sdth    = -DT/qsi*sumn
                  sdth    = MIN(0.0,sdth)
                  SNOMLT  = SNOMLT-sdth
                  SNOWH   = SNOWH+sdth
                  SNOWH   = MAX(SNOWH,0.001)
                  SNEQV   = SNOWH*SNDENS
               ELSEIF (kkkk.GT.15) THEN
                  IF (tsnow(1).LT.0.0) THEN
                     fcontsi = rks/dsh*(T1-(tsnow(1)+273.15))
                     sumn    = qsum+SHEAT+ETA+ql+fcontsi
                     sdth    = -DT/qsi*sumn
                     sdth    = MIN(0.0,sdth)
                     SNOMLT  = SNOMLT-sdth
                     SNOWH   = SNOWH+sdth
                     SNOWH   = MAX(SNOWH,0.001)
                     SNEQV   = SNOWH*SNDENS
                  ENDIF
               ENDIF
            ENDIF !loop SNOWH > 0.01
            dsh = SNOWH/nvs
         ENDIF    !loop sur = 's'
      ENDIF       !loop T1 > tmelt

!Calculate and output surface heat fluxes

      CALL wvap(T1,Q2,es,qs,e)
      CALL coe(T1,SFCTMP,SFCSPD,ZLVL,qs,Q2,ch1,ce)
      SHEAT  = -dp*CP*ch1*(T1-SFCTMP)*SFCSPD
      ETA    = -dp*ce*(qs-Q2)*SFCSPD*rl1(T1)

!Calculate snow/ice interface temperature
      CALL tint(T1,SNOWH,kkkk,dih,dsh,ns,tice,tsnow,rki,rks,ti)

!Heat conductivity in snow layer(vertical layer=ns)
!This is only considered when sss >= 0.01m and kkkk > 15
!The internal melting of snow is also calculated in this section
      ts0=T1-273.15          !  convert unit from [K] to [Deg.]
      tin1=ti-273.15
      IF (SNOWH.GE.0.01.AND.kkkk.gt.15) THEN
         tsnow(nvs)=tin1
         CALL nume2(nvs,SNOWH,ts0,tin1,tsnow,tice,ys,rss,'s',rk0,DT,srks,srki0,srkin)
         jsnow=0
         DO ks = 1,ns
            iarrs(ks) = 0
            IF ( ys(ks).GT.0.0) THEN
!----   Find how many grid point that Tsnow(grid) > zero ----> jsnow
!       Find the grid location that Tsnow(grid) > zero ----> iarrs
               jsnow        = jsnow+1
               iarrs(jsnow) = ks
!----   Restore Tsnow(grid) = tmeltsnow if it's larger than tmeltsnow = 0.0
               ys(ks)       = 0.0
            ENDIF
            tsnow(ks)      = ys(ks)
         ENDDO
         IF (tsnow(1).LT.-100.0) &
         CALL wrf_error_fatal('tsnow too low') 

!----   From the point of view of absorbed energy and heat conductive flux
!----   Below we calculate the internal melting in snow
         IF (jsnow.GE.1) THEN
            IF (jsnow.EQ.iarrs(jsnow).AND.T1.EQ.tmelt) THEN
               hhhs=dsh*jsnow
               qsumsubs = SOLNET-qins(jsnow)

!-       Heat conductive flux between surface and jsnow's layer
               fcontsi  = rks/hhhs*(ts0-tsnow(jsnow))

               sumn     = qsumsubs+SHEAT+ETA+ql+fcontsi
               sdth     = -DT/qsi*sumn  !: surface melt value of snow (a)
               sdth     = MIN(0.0,sdth)
               SNOMLT   = SNOMLT-sdth
               SNOWH    = SNOWH+sdth
               SNOWH    = MAX(SNOWH,0.001)
               SNEQV    = SNOWH*SNDENS

            ELSEIF (jsnow.EQ.iarrs(jsnow).AND.T1.NE.tmelt) THEN
!----    This case means that the surface temperature is below tmelt, but the
!----    lower (jsnow)'th layers have reached the tmelt 

!----    Calculate the accumulation of absorbed solar radiation in sub-surface
!----    layer with Temperature = tmelt.

               sumqins  = 0.0
               DO js     = 1,jsnow
                  sumqins   = sumqins+(qins(js)-qins(js+1))
               ENDDO
               qsumsubs  = sumqins

!!----    Calculate conductive heat flux adjoined with upper and lower grids
               ffup      = rks/dsh*(ts0-(tsnow(iarrs(1))))
               ffdn      = rks/dsh*(tsnow(iarrs(jsnow)+1)-tsnow(iarrs(jsnow)))

               sumn      = qsumsubs+ffup+ffdn

               sdth      = -DT/qsi*sumn
               sdth      = MIN(0.0,sdth)
               SNOMLT    = SNOMLT-sdth
               SNOWH     = SNOWH+sdth
               SNOWH     = MAX(SNOWH,0.001)
               SNEQV     = SNOWH*SNDENS

            ELSEIF (jsnow.LT.iarrs(jsnow).AND.T1.NE.tmelt) THEN
!----    This case means that the internal snow temperature reaches the 
!----    melting point below the first inner grid point.
               sumqins   = 0.0
               DO js     = 1,jsnow
                  sumqins  = sumqins+(qins(js)-qins(js+1))
               ENDDO
               qsumsubs  = sumqins

               IF (iarrs(1).EQ.1) THEN
                  CALL wrf_message('unsuccessivie melting')
                  ffup      = rks/dsh*(ts0-tsnow(iarrs(1)))
               ELSE
                  ffup      = rks/dsh*(tsnow(iarrs(1)-1)-tsnow(iarrs(1)))
               ENDIF

               ffdn      = rks/dsh*(tsnow(iarrs(jsnow)+1)-tsnow(iarrs(jsnow)))

               sumn      = qsumsubs+ffup+ffdn

               sdth      = -DT/qsi*sumn
               sdth      = MIN(0.0,sdth)
               SNOMLT    = SNOMLT-sdth
               SNOWH     = SNOWH+sdth
               SNOWH     = MAX(SNOWH,0.001)
               SNEQV     = SNOWH*SNDENS

            ENDIF   !loop melt type selection
         ENDIF      !loop internal melt exist

         dsh = SNOWH/nvs
         ti  = (rks*dih*(tsnow(ns)+273.15)+rki*dsh*(tice(1)+273.15))&
               /(rks*dih+rki*dsh)
      ENDIF


!----   Below we calculate the change in snow density
      IF (SNOWH.LE.0.005) THEN
         SNDENS  = dsnew/1000.0  
      ELSE
         CALL SNOWPACK (SNEQV,DT,SNOWH,SNDENS,T1,ti)
!        dsnow = dsnow+DT/3600.0*dsnow*ws*c1*EXP(-c2*dsnow)&
!                *EXP(-0.08*(273.15-SFCTMP))
      ENDIF
      SNDENS = MIN(450.0,SNDENS)
      SNEQV  = SNDENS*SNOWH
      dsnow  = SNDENS
      ds     = dsnow*1000.0
      dsh    = SNOWH/nvs


!3.  Heat conductivity in ice layer


      tin2   =  ti-273.15
      tb     =  TBOT-273.15

!Solve the heat conservation equation
      IF (tice(1).LT.-100.0) THEN
         PRINT *,tice(1),tice(11)
         CALL wrf_error_fatal('tice too low before nume')
      END IF

      IF (SNOWH.LE.0.005) THEN
         CALL nume1(nvi,SITHICK,ts0,tb,tsnow,tice,yi,rsi,'i',rk0,DT,srki,srki0,srkin)
      ELSEIF (SNOWH.GT.0.005) THEN
         CALL nume1(nvi,SITHICK,tin2,tb,tsnow,tice,yi,rsi,'i',rk0,DT,srki,srki0,srkin)
      ENDIF

      rki  = srki0
      tu   = tice(ni)   ! temp. at N layer(previous step[i])
      td   = yi(ni)     ! temp. at N layer(new step[i+1])

!----  Set the new value of ice internal temperature
      jice = 0
      DO ki = 1,ni
         IF (yi(ki).gt.tmeltice) THEN

!----    Find how many grid point that Tice(grid) > tmeltice ----> jice
!        Find the grid location that Tice(grid) > tmeltice   ----> iarri
            jice = jice + 1
            iarri(jice) = ki
!----    Restore Tice(grid) = tmeltice if it is larger than tmeltice
            yi(ki) = tmeltice
         ENDIF
         tice(ki) = yi(ki)
      ENDDO
         IF (tice(1).LT.-100.0) CALL wrf_error_fatal('tice too low')

!----   From the point of view of absorbed energy and heat conductive flux
!----   Below we calculate the internal melting in ice
      IF (jice.GE.1) THEN
         IF (jice.EQ.iarri(jice).AND.T1.EQ.tmelt) THEN
            hhhi  = dih*jice
            qsumsubi = SOLNET-qini(jice)

!-       Heat conductive flux between surface and jice's layer
            fcontsi  = rki/hhhi*(ts0-tice(jice))
            sumn     = qsumsubi+SHEAT+ETA+ql+fcontsi

            sdth     = -DT/qii*sumn
            SITHICK  = SITHICK+sdth
            IF (ABS(sdth).gt.0.1) CALL wrf_error_fatal('surface ice melt too fast')
         ELSEIF (jice.EQ.iarri(jice).AND.T1.NE.tmelt) THEN
!----    This case means that surface temperature is below tmeltice, but the
!----    lower (jice)'th layers have reached the tmeltice
!----    Calculate the accumulation of absorbed solar radiation in
!sub-surface layer with Temperature = tmeltice.

            sumqini    = 0.0
            DO ji = 1,jice
               sumqini    = sumqini+(qini(ji)-qini(ji+1))
            ENDDO
            qsumsubi   = sumqini

!----    Calculate conductive heat flux adjoined with upper and lower grids
            ffup = rki/dih*(ts0-(tice(iarri(1))))
            ffdn = rki/dih*(tice(iarri(jice)+1)-tice(iarri(jice)))
            ffup = MAX(ffup,0.0)
            ffdn = MAX(ffdn,0.0)

            sumn = qsumsubi+ffup+ffdn

            sdth = -DT/qii*sumn
            SITHICK   = SITHICK+sdth
            IF (ABS(sdth).gt.0.1) CALL wrf_error_fatal('ice internal melt too fast')

         ELSEIF (jice.LT.iarri(jice).AND.T1.NE.tmelt) THEN
!----    This case means that the internal ice temperature reaches the 
!----    melting point below the first inner ice grid point.
            sumqini   = 0.0
            DO ji = 1,jice
               sumqini   = sumqini+(qini(ji)-qini(ji+1))
            ENDDO
            qsumsubi  = sumqini

            IF (iarri(1).EQ.1) THEN
               ffup = rki/dih*(tin2-tice(iarri(1)))
            ELSE
               ffup = rki/dih*(tice(iarri(1)-1)-tice(iarri(1)))
            ENDIF
            ffdn = rki/dih*(tice(iarri(jice)+1)-tice(iarri(jice)))
            ffup = MAX(ffup,0.0)
            ffdn = MAX(ffdn,0.0)

            sumn = qsumsubi+ffup+ffdn

            sdth    = -DT/qii*sumn
            SITHICK = SITHICK+sdth
            IF (ABS(sdth).gt.0.1) &
            CALL wrf_error_fatal('lower ice melt too fast')

         ENDIF !loop ice melt type
         SITHICK  = MAX(SITHICK,0.1)
      ENDIF    !loop ice melt exist

!--4.  Ice bottom heat balance

      dih   = SITHICK/nvi
!     dth   = DT/qii*(srki(nvi)*(tb-td)/dih-fw)
      dth   = DT/qii*(2.3*(tb-td)/dih-fw)
      IF (ABS(dth).gt.0.1) CALL wrf_error_fatal('bottom ice change too fast')


      SITHICK  = SITHICK + dth
      SITHICK  = MAX(SITHICK,0.1)

! ----------------------------------------------------------------------
! DETERMINE SURFACE ROUGHNESS OVER SNOWPACK USING SNOW CONDITION FROM
! THE PREVIOUS TIMESTEP.
! ----------------------------------------------------------------------

      CALL SNOWZ0 (asnow,Z0,Z0BRD,SNOWH)

!     Calculate effective mixing ratio at ground level (skin)
      Q1=Q2+ETA_KINEMATIC*CP/RCH

! ----------------------------------------------------------------------
! DETERMINE SENSIBLE HEAT (H) IN ENERGY UNITS (W M-2)
! ----------------------------------------------------------------------

      SHEAT = - (CH * CP * SFCPRS)/ (R * T2V) * ( TH2- T1 )

! ----------------------------------------------------------------------
! CONVERT EVAP TERMS FROM KINEMATIC (KG M-2 S-1) TO ENERGY UNITS (W M-2)
! ----------------------------------------------------------------------

      ESNOW = ESNOW * LSUBS
      ETP = ETP*((1.-SNCOVR)*LVH2O + SNCOVR*LSUBS)
      IF (ETP .GT. 0.) THEN
         ETA = ESNOW
      ELSE
         ETA = ETP
      ENDIF

! ----------------------------------------------------------------------
! CONVERT THE SIGN OF SOIL HEAT FLUX SO THAT:
!   SSOIL>0: WARM THE SURFACE  (NIGHT TIME)
!   SSOIL<0: COOL THE SURFACE  (DAY TIME)
! ----------------------------------------------------------------------

      SSOIL = -1.0* SSOIL

! ----------------------------------------------------------------------
! FOR THE CASE OF SEA-ICE, ADD ANY
! SNOWMELT DIRECTLY TO SURFACE RUNOFF (RUNOFF1) SINCE THERE IS NO
! SOIL MEDIUM, AND THUS NO CALL TO SUBROUTINE SMFLX (FOR SOIL MOISTURE
! TENDENCY).
! ----------------------------------------------------------------------
      RUNOFF1 = SNOMLT/DT

! ----------------------------------------------------------------------
    END SUBROUTINE HIGHTSI
! ----------------------------------------------------------------------

      SUBROUTINE CSNOW (SNCOND,DSNOW)

! ----------------------------------------------------------------------
! SUBROUTINE CSNOW
! FUNCTION CSNOW
! ----------------------------------------------------------------------
! CALCULATE SNOW TERMAL CONDUCTIVITY
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL, INTENT(IN) :: DSNOW
      REAL, INTENT(OUT):: SNCOND
      REAL             :: C
      REAL, PARAMETER  :: UNIT = 0.11631

! ----------------------------------------------------------------------
! SNCOND IN UNITS OF CAL/(CM*HR*C), RETURNED IN W/(M*C)
! CSNOW IN UNITS OF CAL/(CM*HR*C), RETURNED IN W/(M*C)
! BASIC VERSION IS DYACHKOVA EQUATION (1960), FOR RANGE 0.1-0.4
! ----------------------------------------------------------------------
      C = 0.328*10** (2.25* DSNOW)
!      CSNOW=UNIT*C

! ----------------------------------------------------------------------
! DE VAUX EQUATION (1933), IN RANGE 0.1-0.6
! ----------------------------------------------------------------------
!      SNCOND=0.0293*(1.+100.*DSNOW**2)
!      CSNOW=0.0293*(1.+100.*DSNOW**2)

! ----------------------------------------------------------------------
! E. ANDERSEN FROM FLERCHINGER
! ----------------------------------------------------------------------
!      SNCOND=0.021+2.51*DSNOW**2
!      CSNOW=0.021+2.51*DSNOW**2

!      SNCOND = UNIT * C
! double snow thermal conductivity
      SNCOND = 2.0 * UNIT * C
!     SNCOND = 2.01*DSNOW**1.88

! ----------------------------------------------------------------------
  END SUBROUTINE CSNOW
! ----------------------------------------------------------------------

  SUBROUTINE PENMAN (SFCTMP,SFCPRS,CH,TH2,PRCP,FDOWN,T24,SSOIL, &
       &             Q2,Q2SAT,ETP,RCH,RR,SNOWNG,FRZGRA,       &
       &             DQSDT2,FLX2,EMISSI,T1)

! ----------------------------------------------------------------------
! CALCULATE POTENTIAL EVAPORATION FOR THE CURRENT POINT.  VARIOUS
! PARTIAL SUMS/PRODUCTS ARE ALSO CALCULATED AND PASSED BACK TO THE
! CALLING ROUTINE FOR LATER USE.
! ----------------------------------------------------------------------

    IMPLICIT NONE
    LOGICAL, INTENT(IN)     :: SNOWNG, FRZGRA
    REAL, INTENT(IN)        :: CH, DQSDT2, FDOWN, PRCP,           &
         &                     Q2, Q2SAT, SSOIL, SFCPRS, SFCTMP,  &
         &                     TH2,EMISSI
    REAL, INTENT(IN)        :: T1, T24, RCH
    REAL, INTENT(OUT)       :: ETP,FLX2,RR
    REAL                    :: ELCP1, LVS, EPSCA, A, DELTA, FNET, RAD

    REAL, PARAMETER      :: ELCP = 2.4888E+3, LSUBC = 2.501000E+6,CP = 1004.6
    REAL, PARAMETER      :: LSUBS = 2.83E+6

! ----------------------------------------------------------------------
! PREPARE PARTIAL QUANTITIES FOR PENMAN EQUATION.
! ----------------------------------------------------------------------

    IF ( T1 > 273.15 ) THEN
       ELCP1=ELCP
       LVS=LSUBC
    ELSE
       ELCP1  = ELCP*LSUBS/LSUBC
       LVS    = LSUBS
    ENDIF

    FLX2 = 0.0
    DELTA = ELCP1 * DQSDT2
    RR = EMISSI * T24 * 6.48E-8 / (SFCPRS * CH) + 1.0

! ----------------------------------------------------------------------
! ADJUST THE PARTIAL SUMS / PRODUCTS WITH THE LATENT HEAT
! EFFECTS CAUSED BY FALLING PRECIPITATION.
! ----------------------------------------------------------------------

    IF ( PRCP > 0.0 ) THEN
       IF (.NOT. SNOWNG) THEN
          RR = RR + CPH2O * PRCP / RCH
       ELSE
          RR = RR + CPICE * PRCP / RCH
       ENDIF
    ENDIF

! ----------------------------------------------------------------------
! INCLUDE THE LATENT HEAT EFFECTS OF FREEZING RAIN CONVERTING TO ICE ON
! IMPACT IN THE CALCULATION OF FLX2 AND FNET.
! ----------------------------------------------------------------------

    FNET = FDOWN - EMISSI * SIGMA * T24 - SSOIL
    IF (FRZGRA) THEN
       FLX2 = - LSUBF * PRCP
       FNET = FNET - FLX2
    END IF

! ----------------------------------------------------------------------
! FINISH PENMAN EQUATION CALCULATIONS.
! ----------------------------------------------------------------------

    RAD = FNET / RCH + TH2 - SFCTMP
    A = ELCP1 * (Q2SAT - Q2)
    EPSCA = (A * RR + RAD * DELTA) / (DELTA + RR)
    ETP = EPSCA * RCH / LVS

! ----------------------------------------------------------------------
  END SUBROUTINE PENMAN
! ----------------------------------------------------------------------

      SUBROUTINE SNOWPACK (ESD,DTSEC,SNOWH,SNDENS,TSNOW,TSOIL)

! ----------------------------------------------------------------------
! SUBROUTINE SNOWPACK
! ----------------------------------------------------------------------
! CALCULATE COMPACTION OF SNOWPACK UNDER CONDITIONS OF INCREASING SNOW
! DENSITY, AS OBTAINED FROM AN APPROXIMATE SOLUTION OF E. ANDERSON'S
! DIFFERENTIAL EQUATION (3.29), NOAA TECHNICAL REPORT NWS 19, BY VICTOR
! KOREN, 03/25/95.
! ----------------------------------------------------------------------
! ESD     WATER EQUIVALENT OF SNOW (M)
! DTSEC   TIME STEP (SEC)
! SNOWH   SNOW DEPTH (M)
! SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
! TSNOW   SNOW SURFACE TEMPERATURE (K)
! TSOIL   SOIL SURFACE TEMPERATURE (K)

! SUBROUTINE WILL RETURN NEW VALUES OF SNOWH AND SNDENS
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER                :: IPOL, J
      REAL, INTENT(IN)       :: ESD, DTSEC,TSNOW,TSOIL
      REAL, INTENT(INOUT)    :: SNOWH, SNDENS
      REAL                   :: BFAC,DSX,DTHR,DW,SNOWHC,PEXP,           &
                                TAVGC,TSNOWC,TSOILC,ESDC,ESDCX
      REAL, PARAMETER        :: C1 = 0.01, C2 = 21.0, G = 9.81,         &
                                KN = 4000.0
! ----------------------------------------------------------------------
! CONVERSION INTO SIMULATION UNITS
! ----------------------------------------------------------------------
      SNOWHC = SNOWH *100.
      ESDC = ESD *100.
      DTHR = DTSEC /3600.
      TSNOWC = TSNOW -273.15
      TSOILC = TSOIL -273.15

! ----------------------------------------------------------------------
! CALCULATING OF AVERAGE TEMPERATURE OF SNOW PACK
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! CALCULATING OF SNOW DEPTH AND DENSITY AS A RESULT OF COMPACTION
!  SNDENS=DS0*(EXP(BFAC*ESD)-1.)/(BFAC*ESD)
!  BFAC=DTHR*C1*EXP(0.08*TAVGC-C2*DS0)
! NOTE: BFAC*ESD IN SNDENS EQN ABOVE HAS TO BE CAREFULLY TREATED
! NUMERICALLY BELOW:
!   C1 IS THE FRACTIONAL INCREASE IN DENSITY (1/(CM*HR))
!   C2 IS A CONSTANT (CM3/G) KOJIMA ESTIMATED AS 21 CMS/G
! ----------------------------------------------------------------------
      TAVGC = 0.5* (TSNOWC + TSOILC)
      IF (ESDC >  1.E-2) THEN
         ESDCX = ESDC
      ELSE
         ESDCX = 1.E-2
      END IF

!      DSX = SNDENS*((DEXP(BFAC*ESDC)-1.)/(BFAC*ESDC))
! ----------------------------------------------------------------------
! THE FUNCTION OF THE FORM (e**x-1)/x EMBEDDED IN ABOVE EXPRESSION
! FOR DSX WAS CAUSING NUMERICAL DIFFICULTIES WHEN THE DENOMINATOR "x"
! (I.E. BFAC*ESDC) BECAME ZERO OR APPROACHED ZERO (DESPITE THE FACT THAT
! THE ANALYTICAL FUNCTION (e**x-1)/x HAS A WELL DEFINED LIMIT AS
! "x" APPROACHES ZERO), HENCE BELOW WE REPLACE THE (e**x-1)/x
! EXPRESSION WITH AN EQUIVALENT, NUMERICALLY WELL-BEHAVED
! POLYNOMIAL EXPANSION.

! NUMBER OF TERMS OF POLYNOMIAL EXPANSION, AND HENCE ITS ACCURACY,
! IS GOVERNED BY ITERATION LIMIT "IPOL".
!      IPOL GREATER THAN 9 ONLY MAKES A DIFFERENCE ON DOUBLE
!            PRECISION (RELATIVE ERRORS GIVEN IN PERCENT %).
!       IPOL=9, FOR REL.ERROR <~ 1.6 E-6 % (8 SIGNIFICANT DIGITS)
!       IPOL=8, FOR REL.ERROR <~ 1.8 E-5 % (7 SIGNIFICANT DIGITS)
!       IPOL=7, FOR REL.ERROR <~ 1.8 E-4 % ...
! ----------------------------------------------------------------------
      BFAC = DTHR * C1* EXP (0.08* TAVGC - C2* SNDENS)
      IPOL = 4
      PEXP = 0.
!        PEXP = (1. + PEXP)*BFAC*ESDC/REAL(J+1)
      DO J = IPOL,1, -1
         PEXP = (1. + PEXP)* BFAC * ESDCX / REAL (J +1)
      END DO

      PEXP = PEXP + 1.
! ----------------------------------------------------------------------
! ABOVE LINE ENDS POLYNOMIAL SUBSTITUTION
! ----------------------------------------------------------------------
!     END OF KOREAN FORMULATION

!     BASE FORMULATION (COGLEY ET AL., 1990)
!     CONVERT DENSITY FROM G/CM3 TO KG/M3
!       DSM=SNDENS*1000.0

!       DSX=DSM+DTSEC*0.5*DSM*G*ESD/
!    &      (1E7*EXP(-0.02*DSM+KN/(TAVGC+273.16)-14.643))

!  &   CONVERT DENSITY FROM KG/M3 TO G/CM3
!       DSX=DSX/1000.0

!     END OF COGLEY ET AL. FORMULATION

! ----------------------------------------------------------------------
! SET UPPER/LOWER LIMIT ON SNOW DENSITY
! ----------------------------------------------------------------------
      DSX = SNDENS * (PEXP)
      IF (DSX > 0.40) DSX = 0.40
      IF (DSX < 0.05) DSX = 0.05
! ----------------------------------------------------------------------
! UPDATE OF SNOW DEPTH AND DENSITY DEPENDING ON LIQUID WATER DURING
! SNOWMELT.  ASSUMED THAT 13% OF LIQUID WATER CAN BE STORED IN SNOW PER
! DAY DURING SNOWMELT TILL SNOW DENSITY 0.40.
! ----------------------------------------------------------------------
      SNDENS = DSX
      IF (TSNOWC >=  0.) THEN
         DW = 0.13* DTHR /24.
         SNDENS = SNDENS * (1. - DW) + DW
         IF (SNDENS >=  0.40) SNDENS = 0.40
! ----------------------------------------------------------------------
! CALCULATE SNOW DEPTH (CM) FROM SNOW WATER EQUIVALENT AND SNOW DENSITY.
! CHANGE SNOW DEPTH UNITS TO METERS
! ----------------------------------------------------------------------
      END IF
      SNOWHC = ESDC / SNDENS
      SNOWH = SNOWHC *0.01

! ----------------------------------------------------------------------
  END SUBROUTINE SNOWPACK
! ----------------------------------------------------------------------

      SUBROUTINE SNOWZ0 (SNCOVR,Z0, Z0BRD, SNOWH)

! ----------------------------------------------------------------------
! SUBROUTINE SNOWZ0
! ----------------------------------------------------------------------
! CALCULATE TOTAL ROUGHNESS LENGTH OVER SNOW
! SNCOVR  FRACTIONAL SNOW COVER
! Z0      ROUGHNESS LENGTH (m)
! Z0S     SNOW ROUGHNESS LENGTH:=0.001 (m)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL, INTENT(IN)        :: SNCOVR, Z0BRD
      REAL, INTENT(OUT)       :: Z0
      REAL, PARAMETER         :: Z0S=0.001
      REAL, INTENT(IN)        :: SNOWH
      REAL                    :: BURIAL
      REAL                    :: Z0EFF

!m      Z0 = (1.- SNCOVR)* Z0BRD + SNCOVR * Z0S
      BURIAL = 7.0*Z0BRD - SNOWH
      IF(BURIAL.LE.0.0007) THEN
        Z0EFF = Z0S
      ELSE      
        Z0EFF = BURIAL/7.0
      ENDIF
      
      Z0 = (1.- SNCOVR)* Z0BRD + SNCOVR * Z0EFF

! ----------------------------------------------------------------------
  END SUBROUTINE SNOWZ0
! ----------------------------------------------------------------------


      SUBROUTINE SNOW_NEW (TEMP,NEWSN,SNOWH,SNDENS)

! ----------------------------------------------------------------------
! SUBROUTINE SNOW_NEW
! ----------------------------------------------------------------------
! CALCULATE SNOW DEPTH AND DENSITY TO ACCOUNT FOR THE NEW SNOWFALL.
! NEW VALUES OF SNOW DEPTH & DENSITY RETURNED.

! TEMP    AIR TEMPERATURE (K)
! NEWSN   NEW SNOWFALL (M)
! SNOWH   SNOW DEPTH (M)
! SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL, INTENT(IN)        :: NEWSN, TEMP
      REAL, INTENT(INOUT)     :: SNDENS, SNOWH
      REAL                    :: DSNEW, HNEWC, SNOWHC,NEWSNC,TEMPC

! ----------------------------------------------------------------------
! CONVERSION INTO SIMULATION UNITS
! ----------------------------------------------------------------------
      SNOWHC = SNOWH *100.
      NEWSNC = NEWSN *100.

! ----------------------------------------------------------------------
! CALCULATING NEW SNOWFALL DENSITY DEPENDING ON TEMPERATURE
! EQUATION FROM GOTTLIB L. 'A GENERAL RUNOFF MODEL FOR SNOWCOVERED
! AND GLACIERIZED BASIN', 6TH NORDIC HYDROLOGICAL CONFERENCE,
! VEMADOLEN, SWEDEN, 1980, 172-177PP.
!-----------------------------------------------------------------------
      TEMPC = TEMP -273.15
      IF (TEMPC <=  -15.) THEN
         DSNEW = 0.05
      ELSE
         DSNEW = 0.05+0.0017* (TEMPC +15.)**1.5
      END IF
! ----------------------------------------------------------------------
! ADJUSTMENT OF SNOW DENSITY DEPENDING ON NEW SNOWFALL
! ----------------------------------------------------------------------
      HNEWC = NEWSNC / DSNEW
      IF (SNOWHC + HNEWC .LT. 1.0E-3) THEN
         SNDENS = MAX(DSNEW,SNDENS)
      ELSE
         SNDENS = (SNOWHC * SNDENS + HNEWC * DSNEW)/ (SNOWHC + HNEWC)
      ENDIF
      SNOWHC = SNOWHC + HNEWC
      SNOWH = SNOWHC *0.01

! ----------------------------------------------------------------------
  END SUBROUTINE SNOW_NEW
! ----------------------------------------------------------------------

END MODULE module_hightsi
