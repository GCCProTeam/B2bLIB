#include "./rtklib.h"


/* B2bLIB function */
extern int satSlot2Sat(int SatSlot);
extern int sat2Slot(int sat);
extern int add_eph(nav_t* nav, const eph_t* eph);
extern initB2b(nav_t* navs);
extern void readB2bType1(FILE* fp, b2bsat_t* b2bsat, b2bsat_t* b2bsat_pre, gtime_t obstime);
extern void readB2bType2(FILE* fp, b2bsat_t* b2bsat, gtime_t obstime);
extern void readB2bType3(FILE* fp, b2bsat_t* b2bsat, gtime_t obstime);
extern void readB2bType4(FILE* fp, b2bsat_t* b2bsat, gtime_t obstime);
extern int readRinex4Nav(const char* file, nav_t* nav);
extern int satpos_B2b(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    double* rs, double* dts, double* var, int* svh);
extern int ephpos_B2b(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    int iode, double* rs, double* dts, double* var, int* svh);
extern eph_t* selB2beph(gtime_t time, int sat, int iodn, const nav_t* nav);
extern void eph2pos_CNAV(gtime_t time, const eph_t* eph, double* rs, double* dts,
    double* var);
extern double varUraB2b(double* ura);
extern int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    int iode, double* rs, double* dts, double* var, int* svh);
extern double var_uraeph(int ura);
extern double prange_dualfrequency(const obsd_t* obs, const nav_t* nav, double* var,
    const double* dantr, const double* dants, const prcopt_t* opt);




