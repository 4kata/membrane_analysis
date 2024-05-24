FC=ifort
HOME=.
XDRLIB=../../xdr/lib
XDRINTERFACE=../mod
MODDIR=../mod
TRGT=phb
EXEF=${TRGT}
OBJ=${XDRINTERFACE}/xdr.F90 \
    ${MODDIR}/info.f90 \
    ${MODDIR}/calc.f90 \
    ${HOME}/${EXEF}.f90

all   : ${TRGT}
clean :
	rm ${SRC}*.o *.mod ${EXEF}
opt:
	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -qopenmp -O3
#	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -O3
	rm -f *.mod

dbg:
	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -check all -traceback
	rm -f *.mod

${TRGT} : ${OBJ} 
	${FC} -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -qopenmp -O3
	cp ${EXEF} ../bin/
	rm -f *.mod
