#pragma once
//GCC
#define MaskBDS 63       // Number of BDS masks
#define MaskGPS 37       // Number of GPS masks
#define MaskGAL 37       // Number of Galileo masks
#define MaskGLO 37       // Number of GLONASS masks
#define MaskNSAT (MaskBDS+MaskGPS+MaskGAL+MaskGLO) // Total satellite masks in B2b_Type1
#define NMODESINGAL 16   // Number of signal tracking modes
#define NINB2bTYPE 86400 // Number of B2b_Type messages per second (one per second)
#define NSAT 168         // max satellite number (1 to MAXSAT) 


// Incremental parameters for B2b message types
#define MAXORBITCOR 6+1  // Max orbit corrections per Type2 message (+1 to prevent overflow)
#define MAXCLOCKCOR 23   // Max clock corrections per Type4 message

// Validity periods for B2b corrections
#define MAXAGEB2b 96.0          // Max age of B2b orbit/URA corrections (seconds)
#define MAXAGEB2b_CBIAS 86400   // Max age of B2b code bias corrections (seconds)
#define MAXAGEB2b_CLOCK 12.0    // Max age of B2b clock corrections (seconds)

// Reference semi-major axes for BDS satellites
#define ArefMEO_BDS 27906100     // Reference semi-major axis (A) for BDS MEO satellites
#define ArefIGSOGEO_BDS 42162200 // Reference semi-major axis (A) for BDS GEO/IGSO satellites


typedef struct {        /* B2b message type */
    int week, tow;       /* reception time (GPS week, time of week) */
    int prn;            /* B2b satellite PRN number */
    int type;           /* message type */
    int iB2b;           /* index within all messages */
    unsigned char msg[61]; /* B2b message (486 bits) */
} b2bmsg_t;

typedef struct {        /* B2b messages container type */
    int n, nmax;        /* current number of messages / max allocated */
    b2bmsg_t* msgs;     /* array of B2b messages */
} b2b_t;

typedef struct {        /* B2b Type1 - Satellite slot information */
    double TodBDT;      /* BDT time of day (TOD) */
    gtime_t t0;         /* Reference time (GPST) */
    int IodSsr;         /* SSR Issue of Data for Type1 */
    int Iodp;           /* Issue of Data for phase */
    int iB2b;           /* Index in message sequence */
} B2bType1_t;

typedef struct {        /* B2b Type2 - Orbit correction */
    double TodBDT;      /* BDT time of day (TOD) */
    gtime_t t0;         /* Reference time (GPST) */
    int IodSsr;         /* SSR Issue of Data */
    int IODN;           /* Issue of Data for navigation */
    int IodCorr;        /* Issue of Data for correction */
    double OrbCorr[3];  /* Orbit corrections: 0:Radial, 1:Along-track, 2:Cross-track */
    double UraClass;    /* URA classification */
    double UraValue;    /* URA value (m) */
    int iB2b;           /* Index in message sequence */
} B2bType2_t;

typedef struct {        /* B2b Type3 - DCB correction */
    int IodSsr;         /* SSR Issue of Data */
    double TodBDT;      /* BDT time of day (TOD) */
    gtime_t t0;         /* Reference time (GPST) */
    int iB2b;           /* Index in message sequence */
    double SatDCB[NMODESINGAL]; /* Satellite DCB corrections:
                                0: B1I L1C/A G1C/A; 1: B1C(D) L1_P G1_P E1_B;
                                2: B1C(P) G2C/A E1C; 3: null;
                                4: B2a(D) L1C(P) E5aQ; 5: B2a(P) L1C(D+P) E5aI;
                                6: null; 7: B2b_I L2C(L) E5b_Q;
                                8: B2b_Q L2C(M+L) E5b_Q; 9: null;
                                10: null; 11: L5I E6_C;
                                12: B3I L5_Q; 13: L5_I+Q */
} B2bType3_t;

typedef struct {        /* B2b Type4 - Clock correction */
    int Iodp;           /* Issue of Data for phase */
    int IodSsr;         /* SSR Issue of Data */
    double TodBDT;      /* BDT time of day (TOD) */
    gtime_t t0;         /* Reference time (GPST) */
    double C0;          /* Clock correction parameter */
    int IodCorr;        /* Issue of Data for correction */
    int iB2b;           /* Index in message sequence */
} B2bType4_t;

typedef struct {        /* B2b correction for current epoch */
    int sat;            /* Satellite index */
    B2bType1_t b2btype1;/* Type1 correction data */
    B2bType2_t b2btype2;/* Type2 correction data */
    B2bType3_t b2btype3;/* Type3 correction data */
    B2bType4_t b2btype4;/* Type4 correction data */
} b2bsatp_t;

typedef struct {        /* B2b satellite corrections container */
    int nsat;           /* Number of satellites */
    int SatSlot[MaskNSAT]; /* Satellite slot numbers:
                            0-62: BDS; 63-99: GPS;
                            100-136: Galileo; 137-173: GLONASS */
    b2bsatp_t b2bsats[MAXSAT]; /* Array of satellite corrections */
} b2bsat_t;
