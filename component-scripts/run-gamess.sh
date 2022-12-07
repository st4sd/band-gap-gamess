#!/usr/bin/env sh

molecule=$1
version=$2
cpus=$3

PATH_RUNGMS_WRAPPER=${PATH_RUNGMS:-/usr/local/bin/rungms}
PATH_GAMESS=${PATH_GAMESS:-/usr/local/bin/gamess}

PATH_MY_GAMESS=${PATH_MY_GAMESS:-/tmp/gamess}
# VV: This Scratch dir should point to a High Performance filesystem
# for example, the ephemeral space of a container or a SSD/NVME.
# You probably do not want to use PVCs or GPFS for these intermediate files.
GAMESS_SCRATCH_DIR=${GAMESS_SCRATCH_DIR:-${PATH_MY_GAMESS}/scratch}

here=`pwd`
mkdir -p "${PATH_MY_GAMESS}"
mkdir -p "${GAMESS_SCRATCH_DIR}"

sed -e "s#set USERSCR=/workspace/restart#set USERSCR=${here}#g" \
    -e "s#set currentdir=\`pwd\`#set currentdir=${PATH_GAMESS}#g" \
    -e "s#set SCR=\`pwd\`/scratch#set SCR=${GAMESS_SCRATCH_DIR}#g" \
    -e "s#TARGET=mpi#TARGET=ga#g" \
    "${PATH_GAMESS}/rungms" >"${PATH_MY_GAMESS}/run-gamess.sh"

cp /usr/local/bin/gamess/install.info "${PATH_GAMESS}/install.info"


# The NVidia Image Features version 00 ONLY and target=GA ONLY

if [[ "${version}" != "00" ]]; then
    echo "Overriding VERSION to 00"
    version=00
fi

chmod +x ${PATH_MY_GAMESS}/run-gamess.sh
# cat ${PATH_MY_GAMESS}/run-gamess.sh
"${PATH_MY_GAMESS}"/run-gamess.sh "${molecule}" "${version}" "${cpus}"