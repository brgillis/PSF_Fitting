# /bin/sh

export TINYTIM=/home/brg/Program_Files/tinytim-7.5

# fits file
# x/y coords
export x=${2}
export y=${3}

PS_SIZE=74
# subsampling factor
export SUBSAMP=${4}
export TMPFILE_BASE=${1}
export FOCUS=${5}

# set of parameters that will be fed into tinytim

# CAM 15 is ACS/WFC
export CAM=${7}
# chip no should be computed from y coord instead?
# CHIP=1|2
export CHIP=${6}
export FILTER=${8}
# for blackbody source
export SOURCE_TEMP_K=5000
export PSF_WIDTH=3.0

# create a psf at the CCD position of the star
# only run tiny1 tiny2 here

${TINYTIM}/tiny1 ${TMPFILE_BASE}.par << EOF
${CAM}
${CHIP}
${x}
${y}
${FILTER}
2
${SOURCE_TEMP_K}
${PSF_WIDTH}
${FOCUS}
${TMPFILE_BASE}
EOF

${TINYTIM}/tiny2 ${TMPFILE_BASE}.par

# run tiny3 too, in order to get charge diffusion kernel in the header
${TINYTIM}/tiny3 ${TMPFILE_BASE}.par SUB=${SUBSAMP}

