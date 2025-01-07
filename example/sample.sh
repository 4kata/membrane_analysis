#!bin/bash
### please read README to custom arguments

bindir=../programs/bin # or write PATH(bindir) in .bashrc
exe=$1  # enter executable file name


export OMP_NUM_THREAD=10



if [ "$exe" = 'phb' ]; then
### phb
  arg1=ow 
  arg2=ow 
  $bindir/phb $arg1 $arg2
elif [ $exe = 'pmf' ]; then
### pmf
  arg1=ow 
  arg2=ow 
  $bindir/pmf $arg1 $arg2
elif [ $exe = 'pz2' ]; then
### pz2
  arg1=phosphate
  arg2=ow
  $bindir/pz2 $arg1 $arg2
elif [ $exe = 'tcf_transition' ]; then
### tcf_transition
  arg1=phosphate
  arg2=ow
  arg3= -0.25
  arg4= 0.80
  $bindir/tcf_transition $arg1 $arg2 $arg3 $arg4
elif [ $exe = 'zp' ]; then
## zp
  arg1=ow
  $bindir/zp $arg1
elif [ $exe = 'zp_near' ]; then
#zp_near
  arg1=phosphate
  arg2=ow
  $bindir/zp_near $arg1 $arg2
elif [ $exe = 'zphb_near' ]; then
#zphb_near
  arg1=phosphate
  arg2=ow
  arg3=olip
  arg4=ochl
  $bindir/zp_near $arg1 $arg2 $arg3 $arg4
fi
