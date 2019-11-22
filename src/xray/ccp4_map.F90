#include "../include/assert.fh"
   subroutine write_map_ccp4(unit, extent, grid_size, rho, unit_cell, &
        spacegroup_number, spacegroup_name, num_symmops, symmop, &
        mask, scale)

      implicit none

      real :: r_dummy
      integer, parameter :: real_kind = kind(r_dummy)
      integer, parameter :: rk_ = real_kind

      integer, intent(in) :: unit, extent(2,3)
      integer, intent(in) :: grid_size(3)
      real(real_kind), intent(in) :: rho(0:grid_size(1), &
                             0:grid_size(2), &
                             0:grid_size(3))
      real(real_kind), intent(in) :: unit_cell(6)
      integer, intent(in) :: spacegroup_number, num_symmops
      character(len=11), intent(in) :: spacegroup_name
      real(real_kind), intent(in) :: symmop(3,4,16)
      integer, intent(in), optional :: mask( &
                             0:grid_size(1), &
                             0:grid_size(2), &
                             0:grid_size(3))
      logical, intent(in), optional :: scale

      ! local
      integer, parameter :: stdout=6
      integer, parameter :: MAX_SYMMOPS = 16 ! actually 96
      logical :: qscale
      real(real_kind) :: map_min,map_max,map_avg,map_sigma
      integer :: a, b, c, nr, isym, u, v, w
      integer :: i, j
      integer :: abc(3), as, bs, cs
      real(real_kind) :: sect_avg, sect_sigma, ed
      equivalence(as,abc(1))
      equivalence(bs,abc(2))
      equivalence(cs,abc(3))
      integer :: order_uvw(3)
      integer :: alloc_stat
      real(real_kind), allocatable :: section(:,:)
      real(real_kind), parameter :: r_epsilon = 1e-20
      if (present(scale)) then
         qscale = scale
      else
         qscale = .false.
      end if
      !-----------------------------------------------------
      ! get statistics for the whole map:
      map_min=0
      map_max=0
      map_avg=0
      map_sigma=0
      nr=0
      do c=lbound(rho,3),ubound(rho,3)
         do b=lbound(rho,2),ubound(rho,2)
            do a=lbound(rho,1),ubound(rho,1)
               if (present(mask)) then
                  if (mask(a,b,c) == 0) cycle
               end if
               ed=rho(a,b,c)
               map_avg=map_avg+ed
               map_sigma=map_sigma+ed**2
               if (nr==0) then
                  map_min=ed
                  map_max=ed
                  nr=1
               else
                  map_min=min(map_min,ed)
                  map_max=max(map_max,ed)
                  nr=nr+1
               end if
            end do
         end do
      end do

      ! print overall statistics
      map_avg=map_avg/nr
      map_sigma=sqrt(max(0.0_rk_,map_sigma/nr-map_avg**2))
      write(stdout,'(A,2(F14.5,A))') &
            ' MAP: avg. density in unit cell=', map_avg, &
            ' sigma=',map_sigma,' e/A^3'
      write(stdout,'(2(A,F14.5))') &
            ' MAP: minimum =', map_min, '  maximum =',map_max

      if (qscale.and.map_sigma < r_epsilon) then
         write(stdout,'(A)') &
               ' XMAPX: sigma very small for map.  No scaling performed.'
         map_sigma=1
      end if

      ! This also sets order_uvw()
      call write_header(unit)
      allocate(section(extent(1,order_uvw(1)):extent(2,order_uvw(1)), &
                       extent(1,order_uvw(2)):extent(2,order_uvw(2))), &
                       stat=alloc_stat)
      REQUIRE(alloc_stat==0)
      isym=1
      SLOW_AXIS: do w=extent(1,order_uvw(3)),extent(2,order_uvw(3))
         abc(order_uvw(3))=w
         sect_avg=0
         sect_sigma=0
         nr=0
         section=0
         MEDIUM_AXIS: do v=extent(1,order_uvw(2)),extent(2,order_uvw(2))
            abc(order_uvw(2))=v
            FAST_AXIS: do u=extent(1,order_uvw(1)),extent(2,order_uvw(1))
               abc(order_uvw(1))=u
               ! Iterate isym only if the current one does not match,
               ! because the same symmop will be used in the majority of cases.
               DO_SYMMOP: do i = 1,num_symmops
                  do j=1,3
                     if (modulo(nint(sum(symmop(1:3,j,isym)*abc(:)) &
                           + symmop(j,4,isym)*grid_size(j)),grid_size(j)) &
                           >= grid_size(j)) then
                        isym=modulo(isym,num_symmops)+1
                        cycle DO_SYMMOP
                     end if
                  end do
                  exit DO_SYMMOP
               end do DO_SYMMOP
               if (present(mask)) then
                  if (mask(as,bs,cs) == 0) cycle FAST_AXIS
               end if
               ed=rho(as,bs,cs)
               section(u,v)=ed
               sect_avg=sect_avg+ed
               sect_sigma=sect_sigma+ed**2
               nr=nr+1
            end do FAST_AXIS
         end do MEDIUM_AXIS

         ! print info about the section
         sect_avg=sect_avg/nr
         sect_sigma=sqrt(max(0.0_rk_,sect_sigma/nr-sect_avg**2))
         write(stdout,'(A,I4,A,F14.5,A,F14.5,A)') &
               ' MAP: section #',w,' avg. dens.=',sect_avg, &
               ' sigma=',sect_sigma,' e/A^3'

         ! scale the section
         if (qscale) section=section-map_avg/map_sigma

         write(unit) section

      end do SLOW_AXIS

   contains

      ! NOTE: (U,V,W) === (Fast,Medium,Slow)
      ! UNIT MUST be opened for RAW BINARY I/O !
      subroutine write_header(unit)
         integer, intent(in) :: unit
         !-------------------------------------------------------------------------
         ! MAP header structure:
         integer(4) :: map_header(256)

         integer(4) :: h_map_uvw_extent(3) ! Map Extent
         integer(4) :: h_map_mode
         integer(4) :: h_map_uvw_origin(3) ! Map Origin
         integer(4) :: h_map_xyz_size(3) ! Grid size
         real(4)    :: h_map_unit_cell(6)
         integer(4) :: h_map_order_uvw(3)
         real(4)    :: h_map_min,h_map_max,h_map_mean
         integer(4) :: h_map_spacegroup_number
         integer(4) :: h_map_symmop_num_chars
         integer(4) :: h_map_skew_flag
         real(4)    :: h_map_skew_matrix(3,3), h_map_skew_translation(3)
         integer(4) :: h_map_pad(15)
         integer(4) :: h_map_tag, h_map_arch_stamp
         real(4)    :: h_map_rms
         integer(4) :: h_map_num_labels
         integer(4) :: h_map_int_labels(20,10) ! == 10 * char(len=80)

         equivalence(map_header( 1),h_map_uvw_extent)
         equivalence(map_header( 4),h_map_mode)
         equivalence(map_header( 5),h_map_uvw_origin)
         equivalence(map_header( 8),h_map_xyz_size)
         equivalence(map_header(11),h_map_unit_cell)
         equivalence(map_header(17),h_map_order_uvw)
         equivalence(map_header(20),h_map_min)
         equivalence(map_header(21),h_map_max)
         equivalence(map_header(22),h_map_mean)
         equivalence(map_header(23),h_map_spacegroup_number)
         equivalence(map_header(24),h_map_symmop_num_chars)
         equivalence(map_header(25),h_map_skew_flag)
         equivalence(map_header(26),h_map_skew_matrix)
         equivalence(map_header(35),h_map_skew_translation)
         equivalence(map_header(36),h_map_pad)
         equivalence(map_header(53),h_map_tag)
         equivalence(map_header(54),h_map_arch_stamp)
         equivalence(map_header(55),h_map_rms)
         equivalence(map_header(56),h_map_num_labels)
         equivalence(map_header(57),h_map_int_labels)

         !-------------------------------------------------------------------------
         !  spacegroups axis ordering: 1=YXZ, 2=ZXY
         integer, parameter :: axis_order_index(23) = &
               (/2,2,2,2,1,1,1,1,1,2,1,1,1,1,1,2,2,2,1,2,2,1,2/)
         integer :: axis_order_uvw(3,2) = reshape((/  2,1,3,  3,1,2  /),(/3,2/))
         character(len=80) :: label
         integer(4) :: native_IT, native_FT

         ! integer :: j, l, ld, lu, lt
         integer :: i, sg
         character(len=60) :: filename
         ! character(len=12) :: username
         ! character(len=11) :: date
         ! character(len=8) :: time
         !-------------------------------------------------------------------------
         ! Allow for CCP4-style alternate origin numbers:
         sg = mod(spacegroup_number,1000)

         if (sg>0 .and. sg<size(axis_order_index)) then
            order_uvw = axis_order_uvw(:,axis_order_index(sg))
         else
            order_uvw = axis_order_uvw(:,1)
         end if

         map_header = 0
         h_map_int_labels=ICHAR(' ')
         h_map_order_uvw = order_uvw
         h_map_uvw_origin(:) = extent(1,order_uvw(:))
         h_map_uvw_extent(:) = extent(2,order_uvw(:)) - h_map_uvw_origin(:) + 1
         h_map_mode=2
         h_map_xyz_size(:) = grid_size(:)
         h_map_unit_cell(:) = unit_cell
         h_map_spacegroup_number = spacegroup_number
         h_map_max = map_max
         h_map_min = map_min
         h_map_mean = map_avg
         h_map_rms = map_sigma

         !----------------------------------------------------------
         ! Generate the endianness flag:
         native_IT = IAND(15,transfer(char(4)//char(3)//char(2)//char(1),0_4))
         native_FT = native_IT
         h_map_arch_stamp = transfer( CHAR(native_FT*16+native_FT) &
               // CHAR(native_IT*16+1) &
               // CHAR(0) // CHAR(0), &
               h_map_arch_stamp)
         h_map_tag = transfer('MAP ',0_4)

         !----------------------------------------------------------
         inquire(unit=unit,name=filename)
         label='FILENAME="'//trim(filename)//'"'
         h_map_int_labels(:,1) = transfer(label,h_map_int_labels(:,1))

         !label='DATE:'//date(1:ld)//'  '// &
               !      time(1:lt)//'       created by user: '//trim(user_name)
         !h_map_int_labels(:,2) = transfer(label,h_map_int_labels(:,2))

         !h_map_int_labels(:,3) = transfer(label,h_map_int_labels(:,3))

         ! Add up to 7 lines of remark titles:
         !i=3
         !do while(i-3<ntitle .and. i<10)
         !   i=i+1
         !   label = title(i-3)
         !   h_map_int_labels(:,i) = transfer(label,h_map_int_labels(:,i))
         !end do
         i=1
         h_map_num_labels = i

         label=" "
         i=i+1
         do while(i<10)
            i=i+1
            h_map_int_labels(:,i) = transfer(label,h_map_int_labels(:,i))
         end do

         write(unit) map_header

      end subroutine write_header
   end subroutine write_map_ccp4

