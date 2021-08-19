! <compile=optimized>
#include "../include/assert.fh"

! This is a module to encapsulate all of the AMBER common blocks and
! other global junk into a module namespace.
module xray_common_module
use file_io_dat
end module xray_common_module

module xray_interface_module

   ! Calls from main SANDER code:
   !  dynlib.f:   if (xray_active) call xray_write_md_state(stdout=6)
   !  printe.f:   if (xray_active) call xray_write_min_state(stdout=6)
   !  force.f:    if (xray_active) call xray_get_derivative(x,f,ener)
   !  mdread.f:   call xray_read_mdin(mdin_lun=5) !<== Also does the global init
   !  sander.f:   call xray_read_parm(prmtop_lun=8,stdout=6)
   !  sander.f:   call xray_init()
   !  sander.f:   call xray_fini()

   use xray_globals_module
   use bulk_solvent_mod, only: k_sol, b_sol
   implicit none
   private

   namelist /xray/ &
         pdb_infile, pdb_outfile, &
         fave_outfile, fmtz_outfile, &
         pdb_read_coordinates, &
         pdb_use_segid, &
         pdb_wrap_names, &
         spacegroup_name, &
         reflection_infile, &
         resolution_low, &
         resolution_high, &
         xray_weight, xray_offset, &
         target, &
         solvent_mask_probe_radius, &
         solvent_mask_expand, &
         solvent_mask_outfile, &
         solvent_mask_reflection_outfile, &
         user_fmask, &
         fft_method, &
         fft_grid_size, &
         fft_grid_spacing, &
         fft_bfactor_sharpen, &
         fft_density_tolerance, &
         fft_reflection_tolerance, &
         fft_radius_min, fft_radius_max, &
         bfactor_min, bfactor_max, &
         bfactor_refinement_interval, &
         atom_selection_mask, solute_selection_mask, &
         k_sol, b_sol, k_tot, b_tot, inputscale,  &
         mask_update_frequency, scale_update_frequency, &
         ml_update_frequency, xray_nstep, bulk_solvent_model

   ! Common public entities all have an xray_ prefix.
   ! Others assume localization by the module scope.
   public :: xray_init_globals, xray_init, xray_fini
   public :: xray_get_derivative, xray_write_md_state
   public :: xray_read_mdin, xray_read_parm, xray_read_pdb
   public :: xray_write_options, xray_write_min_state

   !-------------------------------------------------------------------
contains

   subroutine xray_read_mdin(mdin_lun)
      implicit none
      integer, intent(in) :: mdin_lun
      integer :: stat, inerr
      rewind(mdin_lun)
      read(unit=mdin_lun,nml=xray,iostat=stat)
      if (stat /= 0) then
         write(stdout,'(A)') 'Error reading namelist &xray.'
         call mexit(stdout,1)
      end if

      ! some basic input checks:

      inerr = 0
      if( target /= 'ls' .and. target /= 'ml' .and. target /= 'vls' &
           .and. target /= 'wls' ) then
         write( 6, '(a,a)' ) 'Bad value for target: ', target
         inerr = 1
      end if
      if( mask_update_frequency < 1 ) then
         write( 6, '(a)' ) 'mask_update_frequency must be > 0'
         inerr = 1
      end if
      if( scale_update_frequency < 1 ) then
         write( 6, '(a)' ) 'scale_update_frequency must be > 0'
         inerr = 1
      end if
      if( ml_update_frequency < 1 ) then
         write( 6, '(a)' ) 'ml_update_frequency must be > 0'
         inerr = 1
      end if

      if( inerr > 0 ) call mexit(6,1)

      return
   end subroutine xray_read_mdin

   subroutine xray_write_options()
      use bulk_solvent_mod, only: k_sol, b_sol
      implicit none

      write(stdout,'(/,A)') 'X-ray Refinement Parameters:'
      write(stdout,'(5X,2A)') 'PDB InFile: ',trim(pdb_infile)
      if( pdb_outfile .ne. '' ) &
         write(stdout,'(5X,2A)') 'PDB OutFile:',trim(pdb_outfile)
      if( fave_outfile .ne. '' ) &
         write(stdout,'(5X,2A)') 'FCALC_AVE OutFile:',trim(fave_outfile)
      if( fmtz_outfile .ne. '' ) &
         write(stdout,'(5X,2A)') 'FMTZ OutFile:',trim(fmtz_outfile)
      write(stdout,'(5X,A,L1)') 'PDB Read Coordinates: ',pdb_read_coordinates
      write(stdout,'(5X,A,L1)') 'PDB Use SegID: ',pdb_use_segid
      write(stdout,'(5X,A,L1)') 'PDB Wrap Names: ',pdb_wrap_names
      write(stdout,'(5X,2A)') 'Spacegroup: ',trim(spacegroup_name)
      write(stdout,'(5X,2A)') 'Reflection InFile: ',trim(reflection_infile)
      write(stdout,'(5X,2(A,F8.3))') 'Resolution Range: ',resolution_low,',',resolution_high
      write(stdout,'(5X,A,E10.3)') 'X-ray weight: ',xray_weight
      write(stdout,'(5X,A,A4)') 'Use target: ',target
      write(stdout,'(5X,A,I5)') 'Scale update Interval: ',scale_update_frequency
      ! write(stdout,'(5X,A,F8.3)') 'Solvent mask probe radius: ',solvent_mask_probe_radius
      ! write(stdout,'(5X,A,F8.3)') 'Solvent mask expand: ',solvent_mask_expand
      ! write(stdout,'(5X,2A)') 'Solvent Mask OutFile:',trim(solvent_mask_outfile)
      ! write(stdout,'(5X,2A)') 'Solvent Mask Reflection OutFile:',trim(solvent_mask_reflection_outfile)
      write(stdout,'(5X,A,I5)') 'Solvent Mask Update Interval: ',mask_update_frequency
      write(stdout,'(5X,2(A,F8.3))') 'Solvent scale:',k_sol,', B-factor:', b_sol
      write(stdout,'(5X,A,I2)')   'FFT method: ',fft_method
      if( fft_method > 0 ) then
         write(stdout,'(5X,A,3(5X,I5))') 'FFT Grid Size: ',fft_grid_size
         write(stdout,'(5X,A,F9.5)') 'FFT Grid Spacing: ',fft_grid_spacing
         write(stdout,'(5X,A,F8.3)') 'FFT B-factor Sharpen: ',fft_bfactor_sharpen
         write(stdout,'(5X,A,E10.3)') 'FFT Densty Toleranec: ',fft_density_tolerance
         write(stdout,'(5X,A,E10.3)') 'FFT Reflection Tolerance: ',fft_reflection_tolerance
         write(stdout,'(5X,2(A,F8.3))') 'FFT Radius Min:',fft_radius_min,', Max: ',fft_radius_max
      endif
      ! write(stdout,'(5X,2(A,F8.3))') 'B-Factor Min:',bfactor_min,', Max: ',bfactor_max
      ! write(stdout,'(5X,A,I4)') 'B-factor Refinement Interval: ',bfactor_refinement_interval
      write(stdout,'(5X,2A)') 'Atom Selection Mask:   ',trim(atom_selection_mask)
      write(stdout,'(5X,2A)') 'Solute Selection Mask: ',trim(solute_selection_mask)
      return
   end subroutine xray_write_options

   ! Read X-ray data from the PRMTOP file, and also save pointers to global
   ! PRMTOP data.
   subroutine xray_read_parm(prmtop_lun,out_lun)
      use memory_module, only: natom, nres
      use xray_reciprocal_space_module, only: SYMM_TRICLINIC
      implicit none
      integer, intent(in) :: prmtop_lun, out_lun
      ! local
      character(len=32) :: fmt
      integer :: alloc_status, ierr
      logical :: master
#ifdef MPI
#     include "parallel.h"
#else
      integer :: mytaskid = 0
#endif
      master = (mytaskid == 0)

      ! if (pdb_outfile /= '') then
         num_atoms = natom
         num_residues = nres

         allocate(atom_bfactor(natom), atom_occupancy(natom), &
            atom_selection(natom), residue_chainid(nres), residue_icode(nres), &
            atom_element(natom), atom_altloc(natom), residue_number(nres), &
            solute_selection(natom), stat=alloc_status)
         REQUIRE(alloc_status==0)

         call nxtsec_reset()
         call nxtsec(prmtop_lun,STDOUT,0,'(20I4)','RESIDUE_NUMBER',fmt,ierr)
         read(prmtop_lun,fmt) residue_number
         call nxtsec(prmtop_lun,STDOUT,0,'(20A4)','RESIDUE_CHAINID',fmt,ierr)
         read(prmtop_lun,fmt) residue_chainid

         call nxtsec(prmtop_lun,STDOUT,1,'*','RESIDUE_ICODE',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) residue_icode
         else
            residue_icode=' '
         end if

         call nxtsec(prmtop_lun,STDOUT,1,'*','ATOM_ALTLOC',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) atom_altloc
         else
            atom_altloc=' '
         end if

         call nxtsec(prmtop_lun,STDOUT,0,'(20A4)','ATOM_ELEMENT',fmt,ierr)
         read(prmtop_lun,fmt) atom_element
      ! end if

      if (reflection_infile == '') return

      call nxtsec(prmtop_lun,out_lun,0,'(I4)','XRAY_NUM_SCATTER_TYPES',fmt,ierr)
      if (fmt=='*') then
         write(stdout,'(A)') &
            'ERROR: XRAY_NUM_SCATTER_TYPES not found in PRMTOP file.'
         call mexit(stdout,1)
      end if
      read(prmtop_lun,fmt) num_scatter_types

      allocate(atom_scatter_type(natom),  &
            scatter_coefficients(2,scatter_ncoeffs,num_scatter_types), &
            stat=alloc_status)
      REQUIRE(alloc_status==0)

      call nxtsec(prmtop_lun,out_lun,0,'(20I4)','XRAY_ATOM_SCATTER_TYPE_INDEX',fmt,ierr)
      read(prmtop_lun,fmt) atom_scatter_type
      call nxtsec(prmtop_lun,out_lun,0,'(F12.6)','XRAY_SCATTER_COEFFICIENTS',fmt,ierr)
      read(prmtop_lun,fmt) scatter_coefficients
      call nxtsec(prmtop_lun,out_lun,1,'*','XRAY_SYMMETRY_TYPE',fmt,ierr)
      if (ierr==-2) then
         if( master ) write(STDOUT,*) &
               'XRAY_SYMMETRY_TYPE not found in PRMTOP file; assuming P1'
         num_symmops = 1
         spacegroup_number = 1
         spacegroup_name = 'P 1'
         au_type = 1
         system = SYMM_TRICLINIC
         symmop(:,:,1) = reshape((/1,0,0, 0,1,0, 0,0,1, 0,0,0/),(/3,4/))
         symmop_inv(:,:,1) = symmop(:,:,1)
      else
         stop 'ONLY P1 SUPPORTED FOR NOW'
         read(prmtop_lun,fmt) num_symmops, spacegroup_number, au_type, system
         call nxtsec(prmtop_lun,out_lun,1,'*', &
               'XRAY_SYMMETRY_OPERATORS',fmt,ierr)
         ! ...
      end if
   end subroutine xray_read_parm

   subroutine xray_read_pdb(filename)
      use memory_module, only: residue_label,atom_name,coordinate
      use xray_utils_module, only: allocate_lun
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      character(len=4) :: name,resName,segID,element,altLoc,chainID,iCode
      integer :: resSeq
      real(real_kind) :: xyz(3),occupancy,tempFactor
      character(len=80) :: line
      integer :: unit, iostat, iatom, ires, i, j, ndup, nmiss
      real(real_kind), parameter :: MISSING = -999.0_rk_
      logical :: master
#ifdef MPI
#     include "parallel.h"
#else
      integer :: mytaskid = 0
#endif
      master = (mytaskid == 0)
      ! begin
      atom_occupancy(:)=MISSING
      call amopen(allocate_lun(unit),filename,'O','F','R')
      ndup=0
      iatom=1
      ires=1
      do
         read(unit,'(A)',iostat=iostat) line
         if (iostat/=0) exit
         if (line(1:6)=='END   ') exit
         if (line(1:6)=='ATOM  ' .or. line(1:6)=='HETATM') then
            read(line,'(12X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)') &
                  name,altLoc,resName,chainID,resSeq,iCode, &
                  xyz,occupancy,tempFactor,segID,element
            i = find_atom(name,resName,chainID,resSeq,iCode)
            if (i<0) then
               write(stdout,'(A)') 'Atom not found:'
               write(stdout,'(A)') trim(line)
               stop
            end if
            if (atom_occupancy(i) >= 0) then
               ndup=ndup+1
               if (ndup<10) then
                  if( master ) write(stdout,'(3(A,1X),A,I4,A)') 'PDB: Duplicate ATOM:', &
                        name,resName,chainID(1:1),resSeq,iCode(1:1)
               end if
            end if
            if (pdb_read_coordinates) coordinate(1:3,i) = xyz
            atom_bfactor(i) = tempFactor
            atom_occupancy(i) = occupancy
         end if
      end do
      nmiss = count(atom_occupancy==MISSING)
      if (nmiss>0) then
         if( master ) write(stdout,'(A,I4,A)') 'PDB: missing data for ',nmiss,' atoms.'
         j=0
         do i=1,num_atoms
            if (atom_occupancy(i)==MISSING) then
               atom_occupancy(i)=0
               j=j+1
               if (j<=10) then
                  if( master ) write(stdout,'(3(A,1X),A,I4,A)') 'PDB: Missing ATOM:', &
                        atom_name(i),residue_label(i),residue_chainID(i)(1:1),&
                        residue_number(i),residue_iCode(i)(1:1)
               end if
            end if
         end do
      end if
      if (nmiss==0 .and. ndup==0) then
         if( master ) write(stdout,'(A)') 'PDB: All atoms read successfully.'
      end if
      close(unit)
      return
   end subroutine xray_read_pdb

   function find_atom(name,resName,chainID,resSeq,iCode) result(atom_serial)
      use memory_module, only: residue_pointer,residue_label,atom_name
      implicit none
      integer :: atom_serial
      character(len=4), intent(in) :: name, resName, chainID, iCode
      integer, intent(in) :: resSeq
      ! locals
      character(len=4) :: lname
      integer, save :: ires = 1
      integer :: i,j
      lname = adjustl(name)
      ! first find the matching residue: (ignore resName!)
      do i=1,num_residues
         if (resSeq==residue_number(ires) &
               .and. chainID==residue_chainid(ires) &
               .and. iCode==residue_icode(ires)) then
            ! then find the matching atom name:
            do j = residue_pointer(ires),residue_pointer(ires+1)-1
               if (lname==atom_name(j)) then
                  atom_serial = j
                  return
               end if
            end do
            ! Continue searching, just in case there is a residue
            ! that has been split into two parts.
         end if
         ires = ires + 1
      end do
      atom_serial = -1
      return
   end function find_atom

   subroutine xray_write_pdb(filename)
      use xray_common_module, only: owrite, title, title1
      use xray_utils_module, only: allocate_lun
      use memory_module, only: &
            residue_pointer,residue_label,atom_name,coordinate
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      integer :: unit, iatom, ires, ierr
      integer :: first1, last1, first2, last2, ibond
      integer :: iatom_p, ires_p
      character(len=4) :: name
      character(len=8) :: date
      character(len=10) :: time
      ! character(len=1) :: altloc
      ! character(len=4) :: segid
      character(len=3) :: resName
      logical          :: isStandardRes
      character(len=3), parameter :: standard_pdb_residues(28) = (/ &
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE", &
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL", &
            " DG"," DA"," DT"," DC","  G","  A","  U","  C" /)
      ! Amber modres types: "CYX","HID","HIE","HIP",
      character(len=*), parameter :: pdbfmt_MODRES = &
            '("MODRES",1X,A4,1X,A3,1X,A1,1X,I4,A1,1X,A3,2X,A41)'
      ! GMS: Fix for pgf90 compiler
      character(len=4) :: this_residue_chainid

      call amopen(allocate_lun(unit),filename,owrite,'F','R')
      call date_and_time(date,time)
      if (title/='') write(unit,'(2A)') 'REMARK  ', title
      if (title1/='') write(unit,'(2A)') 'REMARK  ', title1
      write(unit,'(12A)') 'REMARK  Written by Amber 20, SANDER, ', &
            date(1:4),'.',date(5:6),'.',date(7:8),'  ', &
            time(1:2),':',time(3:4),':',time(5:6)

      ! Actually '(6A,3F9.3A9,3F7.2,1X,A11,I4)', with last value = Z
      write(unit,'(A6,3F9.3,3F7.2,1X,A11)') &
            'CRYST1',unit_cell,spacegroup_name
#if 0
      do ires = 1,num_residues
         if (residue_chainid(ires)=='*') cycle
         if (residue_label(ires)=='HID') then
            write(unit,pdbfmt_MODRES) &
                  '----','HID',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HE2 ATOM REMOVED'
         else if (residue_label(ires)=='HIE') then
            write(unit,pdbfmt_MODRES) &
                  '----','HIE',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HD1 ATOM REMOVED'
         else if (residue_label(ires)=='HIP') then
            write(unit,pdbfmt_MODRES) &
                  '----','HIP',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HD1 AND HE2 ATOMS REMOVED'
         end if
      end do
#endif

      do ires = 1,num_residues
         if (residue_chainid(ires)=='*') cycle
         do iatom = residue_pointer(ires), residue_pointer(ires+1)-1
            ! ***NOTE***
            ! This code only adds a leading space to give element-alignment
            ! where possible. It is impossible to follow the PDB version 3
            ! "remediated" format correctly, because it has no alignment rules.
            ! Instead, it assumes you have a complete database of all known
            ! residues, and any other residue names are a fatal error.
            name = atom_name(iatom)
            if (atom_element(iatom)(1:1)==' ' &
                  .and. name(1:1) == atom_element(iatom)(2:2)) then
               if (len_trim(name) < 4 .or. pdb_wrap_names) then
                  name = name(4:4)//name(1:3)
               end if
            end if
            resName=residue_label(ires)(1:3)
            resName=adjustr(resName)
            ! GMS: Fix for pgf90 compiler
            this_residue_chainid = residue_chainid(ires)
            ! DRR: PGI does not seem to like any() inside merge() intrinsic.
            isStandardRes = any(resName==standard_pdb_residues)
            ! don't overflow atom or residue numbers:
            iatom_p = mod( iatom, 100000 )
            ires_p = mod( residue_number(ires), 10000 )
            write(unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)')&
                  merge('ATOM  ', 'HETATM', isStandardRes), &
                  iatom_p,name,atom_altloc(iatom)(1:1), &
                  resName,residue_chainid(ires)(1:1), &
                  ires_p,residue_icode(ires)(1:1), &
                  coordinate(1:3,iatom), &
                  atom_occupancy(iatom), &
                  atom_bfactor(iatom), &
                  merge(this_residue_chainid,'    ',pdb_use_segid), &
                  atom_element(iatom)
         end do
      end do
      write(unit,'(A)') 'END'
      close(unit)
      return
   end subroutine xray_write_pdb

   subroutine xray_init()
      use xray_common_module, only: inpcrd
      use AmberNetcdf_mod, only: NC_checkRestart
      use binrestart, only: read_nc_restart_box
      use xray_utils_module, only: allocate_lun
      use xray_reciprocal_space_module, only: derive_cell_info
      use xray_fourier_module, only: get_mss4
      use findmask, only: atommask
      use memory_module, only: natom,nres,ih,m02,m04,m06,ix,i02,x,lcrd
      use ml_mod, only: init_ml, init_scales
      use bulk_solvent_mod, only: init_bulk_solvent, f_mask, f_solvent
      implicit none
      ! local
      integer :: hkl_lun, i, ier, alloc_status, nstlim = 1, NAT_for_mask1
      double precision :: resolution, fabs_solvent, phi_solvent
      double precision :: a,b,c,alpha,beta,gamma
      real(real_kind) :: phi
      logical :: master
      double precision :: time0, time1
#ifdef MPI
#     include "parallel.h"
#else
      integer :: mytaskid = 0
#endif
      master = (mytaskid == 0)

      if (pdb_infile /= '') call xray_read_pdb(trim(pdb_infile))

      if (reflection_infile == '') xray_active = .false.

      ! get the values for ucell:
      if ( NC_checkRestart(inpcrd) ) then
        if( master ) &
          write(6,'(a,a)') ' getting box info from netcdf file ',trim(inpcrd)
        call read_nc_restart_box(inpcrd,a,b,c,alpha,beta,gamma)
      else
        if( master ) &
          write(6,'(a,a)') ' getting box info from bottom of ',trim(inpcrd)
        call peek_ewald_inpcrd(inpcrd,a,b,c,alpha,beta,gamma)
      endif

      if( master ) write(stdout,'(A,3F9.3,3F7.2)') &
            'XRAY: UNIT CELL= ',a, b, c, alpha, beta, gamma
      call derive_cell_info(a, b, c, alpha, beta, gamma)
      ! Ewald/X-ray equivalences:
      !     XRAY                EWALD
      ! orth_to_frac    == transpose(recip)
      ! frac_to_orth    == ucell
      ! volume          == volume
      ! unit_cell(1:3)  == dirlng
      ! recip_cell(1:3) == 1.0/(reclng)
      !
      ! NOTE: orth_to_frac and ewald:recip are the same as PDB SCALEn records

      !--------------------------------------------------------------
      ! Read reflection data
      call amopen(allocate_lun(hkl_lun),reflection_infile,'O','F','R')
      read(hkl_lun,*,end=1,err=1) num_hkl

      allocate(hkl_index(3,num_hkl),abs_Fobs(num_hkl),sigFobs(num_hkl), &
            mSS4(num_hkl),test_flag(num_hkl),d_star_sq(num_hkl), &
            Fcalc(num_hkl), k_scale(num_hkl), f_weight(num_hkl), &
            stat=alloc_status)
      REQUIRE(alloc_status==0)

      if ( user_fmask ) then
         allocate(f_solvent(num_hkl), stat=alloc_status)
         REQUIRE(alloc_status==0)
      endif

      if (fave_outfile /= '') then
         allocate(Fcalc_ave(num_hkl), stat=alloc_status)
         REQUIRE(alloc_status==0)
         Fcalc_ave(:) = cmplx(0._rk_, 0._rk_, rk_)
      endif

      !  each line contains h,k,l two reals, and an r-free flag
      !  if target /= "vls"  reals are Fobs, sigFobs (for diffraction)
      !  if target == "vls", reals are Fobs, phiFobs (for cryoEM)

      !  further, if user_fmask is set, each line has two additional
      !  reals, giving fabs_solvent and phi_solvent

      if( user_fmask ) then
         do i = 1,num_hkl
            read(hkl_lun,*,end=1,err=2) &
               hkl_index(1:3,i),abs_Fobs(i),sigFobs(i),test_flag(i), &
               fabs_solvent, phi_solvent
            phi_solvent = phi_solvent * 0.0174532925d0
            f_solvent(i) = cmplx( fabs_solvent*cos(phi_solvent), &
                                  fabs_solvent*sin(phi_solvent), rk_ )
            test_flag(i) = min(test_flag(i),1)
            f_weight(i) = 1._rk_/(2._rk_*sigFobs(i)**2)
         end do
      else
         do i = 1,num_hkl
            read(hkl_lun,*,end=1,err=2) &
               hkl_index(1:3,i),abs_Fobs(i),sigFobs(i),test_flag(i)
            test_flag(i) = min(test_flag(i),1)
            f_weight(i) = 1._rk_/(2._rk_*sigFobs(i)**2)
         end do
      endif

      ! 'ls' is an unweighted least-squares target; use 'wls' for
      !    weighted least-squares
      if( target(1:2) == 'ls' ) f_weight(:) = 1.0_rk_

      ! set up complex Fobs(:), if vector target is requested
      if( target(1:3) == 'vls' ) then
         allocate(Fobs(num_hkl),stat=alloc_status)
         REQUIRE(alloc_status==0)
!$omp parallel do  private(phi)
         do i = 1,num_hkl
            !  sigFobs() here is assumed to be really phi(), in degrees
            phi = sigFobs(i) * 0.0174532925d0
            Fobs(i) = cmplx( abs_Fobs(i)*cos(phi), abs_Fobs(i)*sin(phi), rk_ )
         end do
!$omp end parallel do
         ! f_weight(:) = 1.d0
      endif

      ! if( fft_method > 0 ) call FFT_setup()

      call atommask(natom=natom,nres=nres,prnlev=0, &
            igraph=ih(m04),isymbl=ih(m06),ipres=ix(i02), &
            lbres=ih(m02),crd=x(lcrd), &
            maskstr=atom_selection_mask,mask=atom_selection)

      NAT_for_mask1 = sum(atom_selection)
      if( master ) write(6,'(a,i6,a,a)') 'Found ',NAT_for_mask1, &
           ' atoms in ', trim(atom_selection_mask)
      !  also ignore any atoms with zero occupancy:
      do i=1,natom
         if( atom_occupancy(i) == 0._rk_) atom_selection(i) = 0
      end do
      NAT_for_mask = sum(atom_selection)
      if( master .and. NAT_for_mask1 /= NAT_for_mask ) &
         write(6,'(a,i4,a)') 'Removing ',NAT_for_mask1 - NAT_for_mask, &
           ' additional atoms with zero occupancy'

      call atommask(natom=natom,nres=nres,prnlev=0, &
            igraph=ih(m04),isymbl=ih(m06),ipres=ix(i02), &
            lbres=ih(m02),crd=x(lcrd), &
            maskstr=solute_selection_mask,mask=solute_selection)
      if( master ) write(6,'(a,i6,a,a)') 'Found ',sum(solute_selection), &
           ' atoms in ', trim(solute_selection_mask)

      call init_ml(target, nstlim, d_star_sq, resolution)
      call init_bulk_solvent(resolution)

      call get_mss4(num_hkl, hkl_index, mSS4 )

      return
      1 continue
      write(stdout,'(A)') 'End-of-file reading HKL file.'
      call mexit(stdout,1)
      2 continue
      write(stdout,'(A)') 'Error reading HKL file.'
      call mexit(stdout,1)
   end subroutine xray_init

   subroutine xray_init_globals()
      pdb_infile = ''
      pdb_outfile = ''
      fave_outfile = ''
      fmtz_outfile = ''
      pdb_read_coordinates = .false.
      pdb_use_segid = .false.
      pdb_wrap_names = .false.
      target = 'ls  '
      spacegroup_name = 'P 1'
      reflection_infile = ''
      resolution_low = 50
      resolution_high = 0
      xray_weight = 1.0
      solvent_mask_probe_radius = 1.0
      solvent_mask_expand = 0.8
      solvent_mask_reflection_outfile = ''
      solvent_mask_outfile = ''
      fft_method = 0
      fft_grid_size = (/0,0,0/)
      fft_grid_spacing = 0.33_rk_
      fft_bfactor_sharpen = 20
      fft_density_tolerance = 1e-4_rk_
      fft_reflection_tolerance = 1e-4_rk_
      fft_radius_min = 1.0
      fft_radius_max = 4.0
      bfactor_min = 1.0
      bfactor_max = 999.0
      bfactor_refinement_interval = 0
      atom_selection_mask = '!@H='
      solute_selection_mask = ':*'
      mask_update_frequency = 100
      scale_update_frequency = 100
      ml_update_frequency = 100
      xray_nstep = 1
      bulk_solvent_model = 'none'
      return
   end subroutine xray_init_globals

   ! Write X-ray output files and deallocate.
   subroutine xray_fini()
      use bulk_solvent_mod, only : k_mask, f_mask
      implicit none
#     include "extra.h"
      ! local
      integer :: dealloc_status, i
      double precision :: phicalc, phimask

      if (.not.xray_active) return

      if (master .and. pdb_outfile /= '') then
         call xray_write_pdb(trim(pdb_outfile))
      end if
      if (master .and. fave_outfile /= '') then
         Fcalc_ave(:) = Fcalc_ave(:)/n_fcalc_ave
         open(20,file=trim(fave_outfile),action='write')
         write(20,'(11a)') 'h',achar(9),'k',achar(9),'l',achar(9), &
            'fobs',achar(9),'real_fcalc_ave',achar(9),'imag_fcalc_ave'
         write(20,'(11a)') '4N',achar(9),'4N',achar(9),'4N',achar(9), &
            '15N',achar(9),'15N',achar(9),'15N'
         do i=1,num_hkl
            write(20,'(i4,a,i4,a,i4,a,f12.3,a,f15.5,a,f15.5)') hkl_index(1,i), &
             achar(9),hkl_index(2,i),achar(9),hkl_index(3,i),achar(9), &
             abs_Fobs(i), achar(9), real(Fcalc_ave(i)), achar(9), &
             aimag(Fcalc_ave(i))
         end do
         close(20)
      endif

      if (master .and. fmtz_outfile /= '') then
         open(20,file=trim(fmtz_outfile),action='write')
         if( target(1:3) == 'vls' ) then
            ! rdb header:
            write(20,'(15a)') 'h', achar(9), 'k', achar(9), 'l', achar(9), &
               'd', achar(9), 'Fobsr', achar(9), 'Fcalcr', achar(9), &
               'Fobsi', achar(9), 'Fcalci' 
            write(20,'(15a)') '4N', achar(9), '4N', achar(9), '4N', achar(9), &
               '15N', achar(9), '15N', achar(9), '15N', achar(9), &
               '15N', achar(9), '15N' 
            do i=1,num_hkl
#  if 1
               write(20, &
               '(i4,a,i4,a,i4,a,f12.3,a,f12.3,a,f12.3,a,f12.3,a,f12.3)') &
                hkl_index(1,i), &
                achar(9),hkl_index(2,i),achar(9),hkl_index(3,i),achar(9), &
                1./sqrt(d_star_sq(i)), achar(9), &
                real(Fobs(i)), achar(9), real(Fcalc(i)), achar(9),  &
                aimag(Fobs(i)), achar(9), aimag(Fcalc(i))
 
#  else
               phi = atan2( Fcalc(i)%im, Fcalc(i)%re ) * 57.2957795d0
               write(20,'(i4,a,i4,a,i4,a,f12.3,a,f12.3)') hkl_index(1,i), &
                achar(9),hkl_index(2,i),achar(9),hkl_index(3,i),achar(9), &
                abs_Fcalc(i), achar(9), phi
#  endif
            end do
         else
            write(20,'(19a)') 'h',achar(9),'k',achar(9),'l',achar(9), &
               'd',achar(9),'fobs',achar(9),'sigfobs',achar(9), &
               'fcalc',achar(9),'phicalc', achar(9), 'rfree-flag', achar(9), &
               'k_scale'
            write(20,'(19a)') '4N',achar(9),'4N',achar(9),'4N',achar(9), &
               '15N',achar(9),'15N',achar(9), '15N',achar(9),'15N',&
               achar(9),'15N',achar(9),'3N',achar(9),'15N'
            do i=1,num_hkl
               phicalc = atan2( aimag(Fcalc(i)), real(Fcalc(i)) ) * 57.2957795d0
               phimask = atan2( aimag(f_mask(i)), real(f_mask(i)) ) * 57.2957795d0
               write(20,&
         '(i4,a,i4,a,i4,a,f8.3,a,f12.3,a,f12.3,a,f12.3,a,f12.3,a,i1,a,f12.3)') &
                hkl_index(1,i), &
                achar(9),hkl_index(2,i), achar(9), hkl_index(3,i), achar(9), &
                1./sqrt(d_star_sq(i)), achar(9),abs_Fobs(i), achar(9), &
                sigFobs(i), achar(9), abs(Fcalc(i)), achar(9), phicalc, &
                achar(9), test_flag(i), achar(9), k_scale(i) 
            end do
         endif
         close(20)
      endif

      deallocate(atom_bfactor,atom_occupancy,atom_scatter_type, &
            atom_selection,residue_chainid,residue_icode, &
            atom_element,atom_altloc,residue_number, &
            scatter_coefficients, &
            hkl_index,abs_Fobs,sigFobs,mSS4,test_flag, &
            stat=dealloc_status)
      REQUIRE(dealloc_status==0)
      if( target(1:3) == 'vls' ) then
         deallocate(Fobs,stat=dealloc_status)
         REQUIRE(dealloc_status==0)
      endif

   end subroutine xray_fini

   ! gets xray_energy and derivatives
   subroutine xray_get_derivative(xyz,dxyz,xray_e,dB)
      use xray_fourier_module
      use xray_utils_module, only: pack_index
      implicit none
      real(real_kind), intent(in) :: xyz(3,num_atoms)
      real(real_kind), intent(inout) :: dxyz(3,num_atoms)
      real(real_kind), intent(out) :: xray_e
      real(real_kind), optional, intent(out) :: dB(num_atoms)
      ! local
      integer, allocatable :: sel_index(:)
      real(real_kind), allocatable :: frac_xyz(:,:)
      real(real_kind), allocatable :: xray_dxyz(:,:), xray_dB(:)
      real(real_kind), allocatable, target :: abs_Fcalc(:)
      complex(real_kind), allocatable :: dF(:)
      real(real_kind) :: phi, gradnorm_amber, gradnorm_xray
      integer :: status, alloc_status, num_selected, dealloc_status
      integer :: i,ierr
      logical, save :: first=.true.
      integer :: free_flag(num_hkl)
#include "def_time.h"
#ifdef MPI
      include 'mpif.h'
#endif

      call timer_start(TIME_XRAY)
      allocate(sel_index(num_atoms),stat=alloc_status)
      REQUIRE(alloc_status==0)
      call pack_index(atom_selection(:)==1 .and. atom_scatter_type(:)>0, sel_index, num_selected)
      allocate(frac_xyz(3,num_selected),dF(num_hkl), &
            abs_Fcalc(num_hkl), &
            xray_dxyz(3,num_selected),xray_dB(num_selected), stat=alloc_status)
      REQUIRE(alloc_status==0)

      frac_xyz=modulo(matmul(transpose(orth_to_frac),xyz(:,sel_index(1:num_selected))),1.0_rk_)
      
      call timer_start(TIME_IHKL)
      if( fft_method == 0 ) then
         call fourier_Fcalc(num_hkl,hkl_index,Fcalc,mSS4, &
            num_selected,frac_xyz, &
            atom_bfactor(sel_index(1:num_selected)), &
            atom_scatter_type(sel_index(1:num_selected)) , &
            atom_occupancy(sel_index(1:num_selected)) )
#if 0
      else
         call FFT_Fcalc(num_hkl,Fcalc,test_flag-1, &
            num_selected,frac_xyz, &
            atom_bfactor(sel_index(1:num_selected)), &
            atom_occupancy(sel_index(1:num_selected)), &
            atom_scatter_type(sel_index(1:num_selected)), &
            num_scatter_types,scatter_ncoeffs,scatter_coefficients)
#endif
      endif
#ifdef MPI
      ! need reduction on Fcalc so all nodes can subsequently calculate
      !   the target energy, etc.
      call mpi_allreduce( MPI_IN_PLACE, Fcalc, num_hkl, &
           MPI_DOUBLE_COMPLEX, mpi_sum, commsander, ierr)
#endif
      call timer_stop(TIME_IHKL)

      if( target(1:3) == 'vls' ) then
         call dTargetV_dF(xyz, deriv=dF, residual=r_work, xray_energy=xray_energy)
      else if( target(1:2) == 'ls' .or. target(1:3) == 'wls' ) then
         call dTargetLS_dF(xyz, selected=test_flag, deriv=dF, &
                           xray_energy=xray_energy)
      else if(target(1:2) == 'ml' ) then
         call dTargetML_dF(xyz, deriv=dF, xray_energy=xray_energy)
      else
         write(6,*) 'Bad target: ', target
         call mexit(6,1)
      endif

      ! now we can get r_work/r_free since target functions may have scaled
      !   Fcalc, or added a bulk solvent contribution, etc.:
      if( target(1:3) /= 'vls' ) then
         abs_Fcalc(:) = abs(Fcalc(:))
         call get_residual(num_hkl,abs_Fobs,abs_Fcalc,r_work,selected=test_flag)
         free_flag(:) = test_flag(:) - 1
         call get_residual(num_hkl,abs_Fobs,abs_Fcalc,r_free,selected=free_flag)
      else
         abs_Fcalc(:) = abs(Fcalc(:))
         call get_residual(num_hkl,abs_Fobs,abs_Fcalc,r_free)
      endif

      if (xray_weight == 0._rk_) then  ! skip remaining calculations
         xray_energy = 0._rk_
         xray_e = 0._rk_
      else
         xray_energy = xray_weight * xray_energy - xray_offset
         xray_e = xray_energy  ! energy get sent back both through the argument
                               ! list and via xray_globals, to accommodate
                               ! both pmemd and sander interfaces
         dF = xray_weight * dF

         call timer_start(TIME_DHKL)
         if( fft_method == 0 ) then
            if( present(dB) ) then
               call fourier_dTarget_dXYZBQ(num_hkl,hkl_index,dF,mSS4, &
               num_selected,frac_xyz, &
               atom_bfactor(sel_index(1:num_selected)), &
               atom_scatter_type(sel_index(1:num_selected)), &
               dxyz=xray_dxyz, d_tempFactor=xray_dB )
            else
               call fourier_dTarget_dXYZBQ(num_hkl,hkl_index,dF,mSS4, &
               num_selected,frac_xyz, &
               atom_bfactor(sel_index(1:num_selected)), &
               atom_scatter_type(sel_index(1:num_selected)), &
               dxyz=xray_dxyz )
            endif
#if 0
         else
            call FFT_dXYZBQ_dF(num_hkl,dF,test_flag-1, &
               num_selected,frac_xyz, &
               atom_bfactor(sel_index(1:num_selected)), &
               atom_occupancy(sel_index(1:num_selected)), &
               atom_scatter_type(sel_index(1:num_selected)), &
               num_scatter_types,scatter_ncoeffs,scatter_coefficients,xray_dxyz)
#endif
         endif
         call timer_stop(TIME_DHKL)

         ! Convert xray_dxyz() back to orthogonal coordinates: 
         xray_dxyz(:,:) = matmul(orth_to_frac,xray_dxyz(:,:))

#ifndef CLANG  /* flang does not (yet) support norm2()  */
#ifndef MPI
         ! compute norm of gradient from Amber, and from xray: this
         !   information could be used to estimate xray_weight:
         if ( first ) then
            gradnorm_amber = norm2(dxyz(:,sel_index(1:num_selected)))
            gradnorm_xray  = norm2(xray_dxyz(:,:))
            write(6,'(a,3e12.5)') '| gradient norms, amber/xray: ', &
               gradnorm_amber, gradnorm_xray, gradnorm_amber/gradnorm_xray
            first = .false.
         endif
#endif
#endif

         ! combine xray forces (-xray_dxyz) into the force array:
         dxyz(:,sel_index(1:num_selected)) = dxyz(:,sel_index(1:num_selected)) &
             - xray_dxyz(:,:)
         if( present(dB) ) dB(sel_index(1:num_selected)) = - xray_dB(:)

      end if  ! from skipping target + derivatives if xray_weight is zero

      ! DAC: why not allocate/deallocate just once, or make these
      !    automatic variables?  Or is this not important?
      deallocate(frac_xyz,dF, &
            abs_Fcalc, xray_dxyz, xray_dB, stat=dealloc_status)
      REQUIRE(dealloc_status==0)
      deallocate(sel_index,stat=dealloc_status)
      REQUIRE(dealloc_status==0)

      call timer_stop(TIME_XRAY)

   end subroutine xray_get_derivative

   subroutine xray_write_md_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f14.4))') &
        'Exray  = ', xray_energy, ' Rwork   = ', r_work, ' Rfree      = ', r_free
   end subroutine xray_write_md_state

   subroutine xray_write_min_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f13.4))') &
        'Exray   = ', xray_energy, ' Rwork   = ', r_work, ' Rfree      = ', r_free
   end subroutine xray_write_min_state

end module xray_interface_module
