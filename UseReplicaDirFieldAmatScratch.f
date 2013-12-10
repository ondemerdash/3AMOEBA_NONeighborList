c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
c         do i = 1, npole
         do l1 = 1,npole3b
            i=pnum(l1)
            i2 = 3*(l1-1)
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = i, npole
               dscale_dir(ipole(j)) = 1.0d0
               pscale(ipole(j)) = 1.0d0
              dscale_mut(ipole(j)) = 1.0d0
            end do
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale_dir(ip11(j,ii)) = d1scale
               dscale_mut(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale_dir(ip12(j,ii)) = d2scale
               dscale_mut(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale_dir(ip13(j,ii)) = d3scale
               dscale_mut(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale_dir(ip14(j,ii)) = d4scale
               dscale_mut(ip14(j,ii)) = u4scale
            end do
c            do k = i, npole
            do l3 = l1,npole3b
               k = pnum(l3)
               k2 = 3*(l3-1)
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     scale3_dir = 1.0d0
                     scale5_dir = 1.0d0
                     scale7_dir = 1.0d0

                     scale3_mut = dscale_mut(kk)
                     scale5_mut = dscale_mut(kk)

                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3_dir = 1.0d0 - expdamp
                           scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7_dir = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        scale3_mut = scale3_mut * (1.0d0-expdamp)
                        scale5_mut = scale5_mut * (1.0d0-expdamp
     &                                        *(1.0d0-damp))


                        end if
                     end if
                     rr3_dir = scale3_dir / (r*r2)
                     rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
                     rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)

                     rr3 = scale3_mut / (r*r2)
                     rr5 = 3.0d0 * scale5_mut / (r*r2*r2)

                     Txx = -(-rr3 + xr*xr*rr5)
                     Txy = -(xr*yr*rr5)
                     Txz = -(xr*zr*rr5)
                     Tyx = Txy
                     Tyy = -(-rr3 + yr*yr*rr5)
                     Tyz = -(yr*zr*rr5)
                     Tzx = Txz
                     Tzy = Tyz
                     Tzz = -(-rr3 + zr*zr*rr5)

                     M_tot(i2+1,k2+1) = Txx
                     M_tot(i2+1,k2+2) = Txy
                     M_tot(i2+1,k2+3) = Txz
                     M_tot(i2+2,k2+1) = Tyx
                     M_tot(i2+2,k2+2) = Tyy
                     M_tot(i2+2,k2+3) = Tyz
                     M_tot(i2+3,k2+1) = Tzx
                     M_tot(i2+3,k2+2) = Tzy
                     M_tot(i2+3,k2+3) = Tzz

                     M_tot(k2+1,i2+1) = Txx
                     M_tot(k2+1,i2+2) = Txy
                     M_tot(k2+1,i2+3) = Txz
                     M_tot(k2+2,i2+1) = Tyx
                     M_tot(k2+2,i2+2) = Tyy
                     M_tot(k2+2,i2+3) = Tyz
                     M_tot(k2+3,i2+1) = Tzx
                     M_tot(k2+3,i2+2) = Tzy
                     M_tot(k2+3,i2+3) = Tzz


                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
                     fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
                     fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
                     fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
                     fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
                     fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
                     do j = 1, 3
                        fip(j) = fid(j)
                        fkp(j) = fkd(j)
                     end do
                     if (use_polymer .and. r2 .le. polycut2) then
                        do j = 1, 3
                           fid(j) = fid(j) * dscale_dir(kk)
                           fip(j) = fip(j) * pscale(kk)
                           fkd(j) = fkd(j) * dscale_dir(kk)
                           fkp(j) = fkp(j) * pscale(kk)
                        end do
                     end if
                     do j = 1, 3
c                        field(j,i) = field(j,i) + fid(j)
c                        fieldp(j,i) = fieldp(j,i) + fip(j)
                        field(j,l1) = field(j,l1) + fid(j)
                        fieldp(j,l1) = fieldp(j,l1) + fip(j)

                        if (ii .ne. kk) then
c                           field(j,k) = field(j,k) + fkd(j)
c                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                           field(j,l3) = field(j,l3) + fkd(j)
                           fieldp(j,l3) = fieldp(j,l3) + fkp(j)
                        end if
                     end do
c                  end if
               end do
            end do
         end do
      end if

