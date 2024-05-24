FC=ifort
HOME=.
XDRLIB=../../xdr/lib
XDRINTERFACE=../mod
MODDIR=../mod
EXEF=tcf_transition
TRGT=tcf_transition
OBJ=${XDRINTERFACE}/xdr.F90 \
    ${MODDIR}/info.f90 \
    ${MODDIR}/calc.f90 \
    ${HOME}/${TRGT}.f90

all   : ${TRGT}
clean :
	rm ${SRC}*.o *.mod ${EXEF}
opt:
	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -qopenmp
#	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -O3
	rm -f *.mod

dbg:
	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -check all -traceback
	rm -f *.mod

${TRGT} : ${OBJ} 
	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -qopenmp -O3
	cp ${EXEF} ../bin/
	rm -f *.mod
