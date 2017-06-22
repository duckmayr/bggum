#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP ggumR_acceptanceAlpha(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_acceptanceDelta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_acceptanceTau(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_acceptanceTheta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_d4beta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_d_4beta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_dlst(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_d_lst(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_dtrunclst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_d_trunclst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_dtruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_getPriorAlpha(SEXP);
extern SEXP ggumR_getPriorDelta(SEXP);
extern SEXP ggumR_getPriorTaus(SEXP);
extern SEXP ggumR_getPriorTheta(SEXP);
extern SEXP ggumR_ggumMCMC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_ggumProbability(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_loglikelihoodCol(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_loglikelihoodRow(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_p4beta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_plst(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_p_lst(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_proposerAlpha(SEXP, SEXP);
extern SEXP ggumR_proposerDelta(SEXP, SEXP);
extern SEXP ggumR_proposerTau(SEXP, SEXP);
extern SEXP ggumR_proposerTheta(SEXP, SEXP);
extern SEXP ggumR_ptrunclst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_ptruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_q4beta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_qlst(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_qtrunclst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_qtruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_r4beta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_r_4beta(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_rlst(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_rtrunclst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_r_trunclst(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_rtruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggumR_r_truncnorm(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ggumR_acceptanceAlpha",  (DL_FUNC) &ggumR_acceptanceAlpha,  6},
    {"ggumR_acceptanceDelta",  (DL_FUNC) &ggumR_acceptanceDelta,  6},
    {"ggumR_acceptanceTau",    (DL_FUNC) &ggumR_acceptanceTau,    7},
    {"ggumR_acceptanceTheta",  (DL_FUNC) &ggumR_acceptanceTheta,  6},
    {"ggumR_d4beta",           (DL_FUNC) &ggumR_d4beta,           5},
    {"ggumR_d_4beta",          (DL_FUNC) &ggumR_d_4beta,          5},
    {"ggumR_dlst",             (DL_FUNC) &ggumR_dlst,             4},
    {"ggumR_d_lst",            (DL_FUNC) &ggumR_d_lst,            4},
    {"ggumR_dtrunclst",        (DL_FUNC) &ggumR_dtrunclst,        6},
    {"ggumR_d_trunclst",       (DL_FUNC) &ggumR_d_trunclst,       6},
    {"ggumR_dtruncnorm",       (DL_FUNC) &ggumR_dtruncnorm,       5},
    {"ggumR_getPriorAlpha",    (DL_FUNC) &ggumR_getPriorAlpha,    1},
    {"ggumR_getPriorDelta",    (DL_FUNC) &ggumR_getPriorDelta,    1},
    {"ggumR_getPriorTaus",     (DL_FUNC) &ggumR_getPriorTaus,     1},
    {"ggumR_getPriorTheta",    (DL_FUNC) &ggumR_getPriorTheta,    1},
    {"ggumR_ggumMCMC",         (DL_FUNC) &ggumR_ggumMCMC,         5},
    {"ggumR_ggumProbability",  (DL_FUNC) &ggumR_ggumProbability,  5},
    {"ggumR_loglikelihoodCol", (DL_FUNC) &ggumR_loglikelihoodCol, 5},
    {"ggumR_loglikelihoodRow", (DL_FUNC) &ggumR_loglikelihoodRow, 5},
    {"ggumR_p4beta",           (DL_FUNC) &ggumR_p4beta,           5},
    {"ggumR_plst",             (DL_FUNC) &ggumR_plst,             4},
    {"ggumR_p_lst",            (DL_FUNC) &ggumR_p_lst,            4},
    {"ggumR_proposerAlpha",    (DL_FUNC) &ggumR_proposerAlpha,    2},
    {"ggumR_proposerDelta",    (DL_FUNC) &ggumR_proposerDelta,    2},
    {"ggumR_proposerTau",      (DL_FUNC) &ggumR_proposerTau,      2},
    {"ggumR_proposerTheta",    (DL_FUNC) &ggumR_proposerTheta,    2},
    {"ggumR_ptrunclst",        (DL_FUNC) &ggumR_ptrunclst,        6},
    {"ggumR_ptruncnorm",       (DL_FUNC) &ggumR_ptruncnorm,       5},
    {"ggumR_q4beta",           (DL_FUNC) &ggumR_q4beta,           5},
    {"ggumR_qlst",             (DL_FUNC) &ggumR_qlst,             4},
    {"ggumR_qtrunclst",        (DL_FUNC) &ggumR_qtrunclst,        6},
    {"ggumR_qtruncnorm",       (DL_FUNC) &ggumR_qtruncnorm,       5},
    {"ggumR_r4beta",           (DL_FUNC) &ggumR_r4beta,           5},
    {"ggumR_r_4beta",          (DL_FUNC) &ggumR_r_4beta,          4},
    {"ggumR_rlst",             (DL_FUNC) &ggumR_rlst,             4},
    {"ggumR_rtrunclst",        (DL_FUNC) &ggumR_rtrunclst,        6},
    {"ggumR_r_trunclst",       (DL_FUNC) &ggumR_r_trunclst,       5},
    {"ggumR_rtruncnorm",       (DL_FUNC) &ggumR_rtruncnorm,       5},
    {"ggumR_r_truncnorm",      (DL_FUNC) &ggumR_r_truncnorm,      4},
    {NULL, NULL, 0}
};

void R_init_ggumR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
