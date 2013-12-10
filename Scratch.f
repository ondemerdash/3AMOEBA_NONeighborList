      switchmode = 'MPOLE'
      call switch (switchmode)
         call mpi_barrier(mpi_comm_world,ierr)

                ep3bt_tot2=0.0d0
                ep3bt_tot3=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt_tot2(j,i) = 0.0d0
                      dep3bt_tot3(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt_tot2(i,j)=0.0d0
                      virep3bt_tot3(i,j)=0.0d0
                   end do
                end do
              

             offset=int((nmol-1)/numtasks)
             remainder=mod((nmol-1),numtasks)

             if(taskid.le.remainder-1) then
               start=taskid*offset+1
               do moli1 =start,start+offset-1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                 virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot2=ep3bt_tot2+ep3bt                 
               end do 

                moli1=numtasks*offset+taskid+1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot2=ep3bt_tot2+ep3bt

                  call mpi_reduce(ep3bt_tot2,ep3b2,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot2,dep3b2,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot2,virep3b2,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

             else
               start=taskid*offset+1
               do moli1 =start,start+offset-1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)

                 do i=1,3
                   do j=1,3
                    virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                    dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do

                ep3bt_tot2=ep3bt_tot2+ep3bt
               end do 
                  call mpi_reduce(ep3bt_tot2,ep3b2,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot2,dep3b2,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot2,virep3b2,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if

             offset=int((nmol-2)/numtasks)
             remainder=mod((nmol-2),numtasks)


             if(taskid.le.remainder-1) then
               start=taskid*offset+1
               do moli1 =start,start+offset-1
c               do k1 =start,start+offset-1
c                 moli1=mol3new(k1)
                  call Innerloop3(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                 virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot3=ep3bt_tot3+ep3bt
               end do

                moli1=numtasks*offset+taskid+1
                  call Innerloop3(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot3=ep3bt_tot3+ep3bt

                  call mpi_reduce(ep3bt_tot3,ep3b3,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot3,dep3b3,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot3,virep3b3,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

             else
               start=taskid*offset+1
               do moli1 =start,start+offset-1

                  call Innerloop3(moli1,ep3bt,virep3bt,dep3bt)

                 do i=1,3
                   do j=1,3
                    virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                    dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do

                ep3bt_tot3=ep3bt_tot3+ep3bt
               end do
                  call mpi_reduce(ep3bt_tot3,ep3b3,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot3,dep3b3,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot3,virep3b3,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if


           call mpi_barrier(mpi_comm_world,ierr)          
