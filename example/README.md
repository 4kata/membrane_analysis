## about system
DPPC : 100 * 2 = 200  
Chol : 11  * 2 = 22  
water : 11100 * 2 = 22200  

## about trajctory
dt = 0.02 ps  
flame = 50000  
total_time = 1 ns  
temp = 303K  

## about toppar
Place itpfile.  
In the case of charmm-gui, simply copy or link to the file as it is.  
If you use other tools to create the system, put ${RESNAME}.itp.  

## about MDinfo
&system_info  
nojumpfile = "*xxx*.xtc" ! trajctory name, available only xtc file  
wrapfile = "*xxx*.xtc" ! same as nojumpfile  
nstep = *xxx*  ! enter number of flames  
dt = *xxx* ! writing frequency [ps]  
/  

## about MOLinfo
*N* ! total molecular species in the system  
[ *name* ] ! name   
*"order in the system"* *"resname"* *"number of atoms in a molecule"* *"number of molecules"* *"COM or NOT"*  

#COM is used to calculate the center of mass in the membrane

## about IDinfo
[ *name* ] ! The name used in MOLinfo cannot be used.  
*resname*  
*total number of atoms in the group*
*ID1 ID2 ID3*.....  

# usage programs

bindir=../programs/bin # or write PATH(bindir) in .bashrc  

## phb : Fig.S14, S15
arg1= #donor side oxygen name  
arg2= #acceptor oxygen (atom or group) name  
$bindir/phb $arg1 $arg2  

## pmf : Fig.5, 6, S6-S12
arg1= #donor side oxygen name  
arg2= #acceptor oxygen (atom or group) name  
$bindir/pmf $arg1 $arg2  

## pz2 : Fig.3
arg1= #centering_atom_name  
arg2= #group_name to  distrution  
$bindir/pz2 $arg1 $arg2  

## tcf_transition : Fig.5
arg1= #centering_atom_name  
arg2= #group_name to  distrution  
arg3= #z coordinate which splits region between 1 and 2  
arg4= #z coordinate which splits region between 2 and 3  
$bindir/tcf_transition $arg1 $arg2 $arg3 $arg4  

## zp :Fig.1
arg1= #group_name  
$bindir/zp $arg1  

## zp_near :Fig.4(a)(b)
arg1= #centering_atom_name  
arg2= #group_name to  distrution  
$bindir/zp_near $arg1 $arg2  

## zphb_near :Fig.4(c)(d)
arg1= #centering_atom_name  
arg2= #donor side oxygen name  
arg3(...argn)= #acceptor oxygen (atom or group) name  
$bindir/zp_near $arg1 $arg2 $arg3 ($arg4...$argn)  

