#!/usr/bin/env bash
# Test suite for schmutzi. Invoked via "make test" (quick) or "make test-full".
#
#   quick : runs contDeam.pl on both test BAMs and checks that the deamination
#           based contamination estimates land where they are expected
#   full  : additionally runs the whole iterative schmutzi.pl pipeline on the
#           simulated dataset (~20-40 minutes)
#
# Requires the test BAMs, get them with "make testdata".

set -u

BASEDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )
SRC=${BASEDIR}/src
TESTDATA=${BASEDIR}/testdata
OUT=${OUTDIR:-${BASEDIR}/testoutput}
FULL=0
THREADS=${THREADS:-$( (nproc 2>/dev/null || echo 1) )}

if [ "${1:-}" = "--full" ]; then
    FULL=1
fi

FAILED=0
SKIPPED=0

pass(){ echo "  [PASS] $1"; }
fail(){ echo "  [FAIL] $1"; FAILED=$((FAILED+1)); }
skip(){ echo "  [SKIP] $1"; SKIPPED=$((SKIPPED+1)); }

# is $1 within [$2,$3] ? using awk to avoid depending on bc
inrange(){
    awk -v v="$1" -v lo="$2" -v hi="$3" 'BEGIN{ exit !(v+0>=lo && v+0<=hi) }'
}

echo
echo "=== schmutzi test suite ==="
echo "output directory: ${OUT}"
echo

##############################################################################
# 1) are the binaries and scripts there and do they run
##############################################################################
echo "[1] executables"

BINARIES="endoCaller mtCont damage2profile log2freq log2fasta contDeam msa2freq \
bam2prof insertSize logs2pos countRecords jointFreqDeaminated \
jointFreqDeaminatedDouble msa2log log2ConsensusLog logdiff logmask filterlog \
msa2singlefreq subSampleBAM addXACircular splitEndoVsCont/poshap2splitbam"

for b in ${BINARIES}; do
    if [ -x "${SRC}/${b}" ]; then
        pass "${b} built"
    else
        fail "${b} is missing, run make"
    fi
done

for s in contDeam.pl schmutzi.pl; do
    if perl -c "${SRC}/${s}" >/dev/null 2>&1; then
        pass "${s} compiles"
    else
        fail "${s} does not compile"
    fi
done

# endoCaller with no argument prints its usage and exits non zero, we only
# check that it does not segfault
"${SRC}/endoCaller" >/dev/null 2>&1
if [ $? -lt 128 ]; then
    pass "endoCaller runs"
else
    fail "endoCaller crashed"
fi

##############################################################################
# 2) reference data shipped with the package
##############################################################################
echo
echo "[2] shipped data"

for f in share/schmutzi/refs/human_MT_wrapped.fa share/schmutzi/refs/human_MT.fa; do
    if [ -s "${BASEDIR}/${f}" ]; then
        pass "${f}"
    else
        fail "${f} is missing"
    fi
done

FREQDIR=${BASEDIR}/share/schmutzi/alleleFreqMT/197/freqs
NFREQ=$(ls "${FREQDIR}"/*.freq 2>/dev/null | wc -l)
if [ "${NFREQ}" -gt 0 ]; then
    pass "contaminant database (${NFREQ} freq files)"
else
    fail "no *.freq file found in ${FREQDIR}"
fi

##############################################################################
# 3) test BAM files
##############################################################################
echo
echo "[3] test data"

HAVEBAM=1
for f in mezB9687.bam mezB9687.bam.bai simulation.bam simulation.bam.bai; do
    if [ -s "${TESTDATA}/${f}" ]; then
        pass "testdata/${f}"
    else
        fail "testdata/${f} is missing, run \"make testdata\""
        HAVEBAM=0
    fi
done

if [ "${HAVEBAM}" -eq 0 ]; then
    echo
    echo "=== no test data, stopping here ==="
    exit 1
fi

mkdir -p "${OUT}"

if ! command -v Rscript >/dev/null 2>&1; then
    echo
    echo "  NOTE: Rscript was not found, the plots (*.pdf) cannot be produced."
    echo "        The tests below still check the numerical estimates."
    HAVER=0
else
    HAVER=1
fi

##############################################################################
# 4) countRecords / bam2prof on the test BAMs
##############################################################################
echo
echo "[4] reading the test BAM files"

NREC=$("${SRC}/countRecords" "${TESTDATA}/simulation.bam" 2>/dev/null | tr -dc '0-9')
if [ -n "${NREC}" ] && [ "${NREC}" -gt 0 ]; then
    pass "countRecords on simulation.bam (${NREC} records)"
else
    fail "countRecords returned nothing for simulation.bam"
fi

"${SRC}/bam2prof" -length 5 -endo -double \
    -5p "${OUT}/prof.5p" -3p "${OUT}/prof.3p" \
    "${TESTDATA}/simulation.bam" >/dev/null 2>&1

if [ -s "${OUT}/prof.5p" ]; then
    # the C->T rate on the first base of the 5' end, simulated data is heavily
    # deaminated so this has to be clearly above zero
    CT=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="C>T") c=i} NR==2{print $c}' "${OUT}/prof.5p")
    if [ -n "${CT}" ] && inrange "${CT}" 0.01 1; then
        pass "bam2prof 5' C>T on simulation.bam = ${CT}"
    else
        fail "bam2prof 5' C>T on simulation.bam is '${CT}', expected a rate above 0.01"
    fi
else
    fail "bam2prof produced no 5' profile"
fi

##############################################################################
# 5) contDeam.pl, deamination based contamination estimate
##############################################################################
echo
echo "[5] contDeam.pl"

# $1 name, $2 bam, $3 library, $4 expected low, $5 expected high
run_contdeam(){
    local name=$1 bam=$2 lib=$3 lo=$4 hi=$5
    local log=${OUT}/${name}.contDeam.log

    "${SRC}/contDeam.pl" --library "${lib}" --ref "${BASEDIR}/share/schmutzi/refs/human_MT_wrapped.fa" \
        --out "${OUT}/${name}" "${bam}" > "${log}" 2>&1
    local rc=$?

    if [ ! -s "${OUT}/${name}.cont.est" ]; then
        fail "contDeam.pl on ${name} produced no ${name}.cont.est (see ${log})"
        return
    fi

    if [ ${rc} -ne 0 ]; then
        if [ "${HAVER}" -eq 0 ]; then
            skip "contDeam.pl on ${name} exited non zero, plotting needs R"
        else
            fail "contDeam.pl on ${name} exited with ${rc} (see ${log})"
        fi
    else
        pass "contDeam.pl on ${name} ran"
    fi

    local est low high
    est=$(cut -f1 "${OUT}/${name}.cont.est")
    low=$(cut -f2 "${OUT}/${name}.cont.est")
    high=$(cut -f3 "${OUT}/${name}.cont.est")

    if inrange "${est}" "${lo}" "${hi}"; then
        pass "${name} contamination estimate ${est} [${low}-${high}], expected ${lo}-${hi}"
    else
        fail "${name} contamination estimate is ${est} [${low}-${high}], expected ${lo}-${hi}"
    fi

    if [ "${HAVER}" -eq 1 ]; then
        if [ -s "${OUT}/${name}.cont.pdf" ]; then
            pass "${name}.cont.pdf was produced"
        else
            fail "${name}.cont.pdf was not produced"
        fi
    fi
}

# simulated data, 20% contamination, double stranded library
run_contdeam sim "${TESTDATA}/simulation.bam" double 0.15 0.25
# empirical data (Mezmaiskaya), single stranded library, deamination based
# estimate for this sample sits around 0.5
run_contdeam mez "${TESTDATA}/mezB9687.bam" single 0.40 0.60

##############################################################################
# 6) endoCaller, endogenous consensus
##############################################################################
echo
echo "[6] endoCaller"

"${SRC}/endoCaller" -seq "${OUT}/endo.fa" -log "${OUT}/endo.log" \
    -deam5p "${OUT}/sim.endo.5p.prof" -deam3p "${OUT}/sim.endo.3p.prof" \
    "${BASEDIR}/share/schmutzi/refs/human_MT_wrapped.fa" \
    "${TESTDATA}/simulation.bam" > "${OUT}/endoCaller.log" 2>&1

if [ -s "${OUT}/endo.fa" ]; then
    LEN=$(grep -v '^>' "${OUT}/endo.fa" | tr -d '\n-' | wc -c)
    # the consensus is called against the wrapped reference so it has the
    # length of that reference, give or take a few indels
    REFLEN=$(grep -v '^>' "${BASEDIR}/share/schmutzi/refs/human_MT_wrapped.fa" | tr -d '\n' | wc -c)
    if inrange "${LEN}" $((REFLEN-100)) $((REFLEN+100)); then
        pass "endoCaller consensus is ${LEN} bp, reference is ${REFLEN} bp"
    else
        fail "endoCaller consensus is ${LEN} bp, expected about ${REFLEN}"
    fi

    NS=$(grep -v '^>' "${OUT}/endo.fa" | tr -cd 'Nn' | wc -c)
    if [ "${NS}" -lt 1000 ]; then
        pass "endoCaller consensus has ${NS} N"
    else
        fail "endoCaller consensus has ${NS} N, the consensus is mostly empty"
    fi
else
    fail "endoCaller produced no consensus (see ${OUT}/endoCaller.log)"
fi

##############################################################################
# 7) log2fasta / log2freq on the consensus produced above
##############################################################################
echo
echo "[7] log file utilities"

if [ -s "${OUT}/endo.log" ]; then
    if "${SRC}/log2fasta" -indel 50 "${OUT}/endo.log" > "${OUT}/fromlog.fa" 2>/dev/null &&
       [ -s "${OUT}/fromlog.fa" ]; then
        pass "log2fasta"
    else
        fail "log2fasta produced no fasta"
    fi

    if "${SRC}/log2freq" "${OUT}/endo.log" > "${OUT}/fromlog.freq" 2>/dev/null &&
       [ -s "${OUT}/fromlog.freq" ]; then
        pass "log2freq"
    else
        fail "log2freq produced no freq file"
    fi
else
    skip "log2fasta/log2freq, endoCaller produced no log"
fi

##############################################################################
# 8) the R scripts
##############################################################################
echo
echo "[8] R scripts"

if [ "${HAVER}" -eq 0 ]; then
    skip "R scripts, Rscript is not installed"
else
    # $1 script, $2 package it needs (or "-"), $3 expected output file, rest: arguments
    run_rscript(){
        local script=$1 pkg=$2 expected=$3
        shift 3
        local log=${OUT}/${script}.log

        if [ "${pkg}" != "-" ]; then
            local missing=""
            for p in ${pkg}; do
                Rscript -e "quit(status=!requireNamespace('${p}',quietly=TRUE))" >/dev/null 2>&1 ||
                    missing="${missing} ${p}"
            done
            if [ -n "${missing}" ]; then
                skip "${script}, missing R package(s):${missing}"
                return
            fi
        fi

        rm -f "${expected}"
        if ! "${SRC}/${script}" "$@" > "${log}" 2>&1; then
            fail "${script} exited non zero (see ${log})"
            return
        fi

        if [ -s "${expected}" ]; then
            pass "${script}"
        else
            fail "${script} produced no ${expected} (see ${log})"
        fi
    }

    # posterior probability plot, takes the contDeam.pl log
    run_rscript posteriorDeam.R - "${OUT}/rtest.cont.pdf" \
        "${OUT}/sim.cont.deam" "${OUT}/rtest.cont.pdf" "test"

    # log normal fit of the fragment length distribution, needs fitdistrplus.
    # This is what schmutzi.pl --uselength relies on, it has to print the three
    # line table that contDeam.pl parses
    "${SRC}/insertSize" "${TESTDATA}/simulation.bam" 2>/dev/null | gzip > "${OUT}/rtest.size.gz"
    if Rscript -e "quit(status=!requireNamespace('fitdistrplus',quietly=TRUE))" >/dev/null 2>&1; then
        "${SRC}/approxDist.R" "${OUT}/rtest.size.gz" > "${OUT}/rtest.size.param" 2>"${OUT}/approxDist.R.log"
        NLINE=$(grep -c "" "${OUT}/rtest.size.param")
        MEANLOG=$(awk 'NR==2{print $1}' "${OUT}/rtest.size.param")
        if [ "${NLINE}" -eq 3 ] && [ -n "${MEANLOG}" ] && inrange "${MEANLOG}" 1 10; then
            pass "approxDist.R, meanlog ${MEANLOG}"
        else
            fail "approxDist.R gave ${NLINE} line(s), meanlog '${MEANLOG}', expected 3 lines and a meanlog between 1 and 10"
        fi
    else
        skip "approxDist.R, missing R package: fitdistrplus (needed by --uselength)"
    fi

    # damage plot, takes a bam2prof profile
    run_rscript deamProf2pdf.R "gplots ggplot2 reshape" "${OUT}/rtest.deam.pdf" \
        "${OUT}/prof.5p" "${OUT}/rtest.deam.pdf" "test"

    # consensus plots, take the endoCaller log. This one names its own output
    if [ -s "${OUT}/endo.log" ]; then
        cp "${OUT}/endo.log" "${OUT}/rtest.cons.log"
        if "${SRC}/mitoConsPDF.R" "${OUT}/rtest.cons.log" > "${OUT}/mitoConsPDF.R.log" 2>&1; then
            NPDF=0
            for p in qual cov div; do
                [ -s "${OUT}/rtest.cons.log.${p}.pdf" ] && NPDF=$((NPDF+1))
            done
            if [ "${NPDF}" -eq 3 ]; then
                pass "mitoConsPDF.R, 3 plots"
            else
                fail "mitoConsPDF.R produced ${NPDF} of the 3 expected plots"
            fi
        else
            fail "mitoConsPDF.R exited non zero (see ${OUT}/mitoConsPDF.R.log)"
        fi
    else
        skip "mitoConsPDF.R, endoCaller produced no log"
    fi

    # chi squared test on the joint deamination frequencies
    if "${SRC}/jointFreqDeaminated" "${TESTDATA}/simulation.bam" > "${OUT}/rtest.joint" 2>/dev/null &&
       [ -s "${OUT}/rtest.joint" ]; then
        if "${SRC}/jointFreqDeaminated.R" "${OUT}/rtest.joint" > "${OUT}/jointFreqDeaminated.R.log" 2>&1 &&
           grep -q "X-squared" "${OUT}/jointFreqDeaminated.R.log"; then
            pass "jointFreqDeaminated.R"
        else
            fail "jointFreqDeaminated.R gave no chi squared result (see ${OUT}/jointFreqDeaminated.R.log)"
        fi
    else
        fail "jointFreqDeaminated produced no output to feed jointFreqDeaminated.R"
    fi
fi

##############################################################################
# 9) full pipeline, only with --full
##############################################################################
echo
echo "[9] schmutzi.pl"

if [ "${FULL}" -eq 0 ]; then
    skip "iterative pipeline, run \"make test-full\" to include it"
else
    echo "  running the iterative pipeline on the simulated data with ${THREADS} thread(s),"
    echo "  this takes a while ..."
    LOG=${OUT}/sim.schmutzi.log
    "${SRC}/schmutzi.pl" --ref "${BASEDIR}/share/schmutzi/refs/human_MT_wrapped.fa" \
        --t "${THREADS}" "${OUT}/sim" \
        "${BASEDIR}/share/schmutzi/alleleFreqMT/197/freqs/" \
        "${TESTDATA}/simulation.bam" > "${LOG}" 2>&1
    rc=$?

    if [ ${rc} -ne 0 ] && [ "${HAVER}" -eq 0 ]; then
        skip "schmutzi.pl exited non zero, plotting needs R"
    elif [ ${rc} -ne 0 ]; then
        fail "schmutzi.pl exited with ${rc} (see ${LOG})"
    else
        pass "schmutzi.pl ran"
    fi

    for f in _final.cont.est _final_endo.fa _final_cont.fa; do
        if [ -s "${OUT}/sim${f}" ]; then
            pass "sim${f} was produced"
        else
            fail "sim${f} was not produced (see ${LOG})"
        fi
    done

    if [ -s "${OUT}/sim_final.cont.est" ]; then
        FEST=$(cut -f1 "${OUT}/sim_final.cont.est")
        if inrange "${FEST}" 0.15 0.25; then
            pass "final contamination estimate ${FEST}, expected 0.15-0.25"
        else
            fail "final contamination estimate is ${FEST}, expected 0.15-0.25"
        fi
    fi

    if [ -s "${OUT}/sim_final_endo.fa" ]; then
        LEN=$(grep -v '^>' "${OUT}/sim_final_endo.fa" | tr -d '\n-' | wc -c)
        if inrange "${LEN}" 16400 16700; then
            pass "endogenous genome is ${LEN} bp"
        else
            fail "endogenous genome is ${LEN} bp, expected about 16569"
        fi
    fi
fi

##############################################################################
echo
echo "=== ${FAILED} failure(s), ${SKIPPED} skipped ==="
echo
exit $(( FAILED > 0 ))
