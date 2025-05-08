#include "rtklib.h"
#include "B2bLIB.h"

#define SQR(x)      ((x)*(x))
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */
#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */
#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */
#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */


static const double ura_eph[] = {         /* ura values (ref [3] 20.3.3.3.1.1) */
	2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
	3072.0,6144.0,0.0
};
// Convert satellite slot number to satellite number
extern int satSlot2Sat(int SatSlot)
{
    int sat, sys;

    if (SatSlot <= 63 && SatSlot >= 1) { // BDS satellites (1-63)
        sys = SYS_CMP;
        sat = satno(sys, SatSlot);
    }
    else if (SatSlot >= 64 && SatSlot <= 100) { // GPS satellites (64-100)
        sys = SYS_GPS;
        sat = satno(sys, SatSlot - MaskBDS);
    }
    else if (SatSlot >= 101 && SatSlot <= 137) { // Galileo satellites (101-137)
        sys = SYS_GAL;
        sat = satno(sys, SatSlot - MaskBDS - MaskGPS);
    }
    else if (SatSlot >= 138 && SatSlot <= 174) { // GLONASS satellites (138-174)
        sys = SYS_GLO;
        sat = satno(sys, SatSlot - MaskBDS - MaskGPS - MaskGLO);
    }
    else if (SatSlot >= 175 || SatSlot <= 0) { // Invalid slot range
        sat = 0;
    }

    return sat;
}

// Convert satellite number to satellite slot index
extern int sat2Slot(int sat)
{
    int prn = 0, sys = 0, SatSlot = 0;

    sys = satsys(sat, &prn); // Extract system type and PRN

    if (sys == SYS_CMP) {
        SatSlot = prn;
    }
    else if (sys == SYS_GPS) {
        SatSlot = prn + MaskBDS;
    }
    else if (sys == SYS_GAL) {
        SatSlot = prn + MaskBDS + MaskGPS;
    }
    else if (sys == SYS_GLO) {
        SatSlot = prn + MaskBDS + MaskGPS + MaskGLO;
    }
    else {
        SatSlot = 0;
    }

    return SatSlot;
}

// Adjust time to be within ¡À3.5 days of reference epoch
static gtime_t adjweek(gtime_t t, gtime_t t0)
{
    double tt = timediff(t, t0);

    if (tt < -302400.0) return timeadd(t, 604800.0); // Add one GPS week
    if (tt > 302400.0) return timeadd(t, -604800.0); // Subtract one GPS week
    return t;
}

// Convert URA value (meters) to URA index according to standard thresholds
static int uraindex(double value)
{
    int i;
    for (i = 0; i < 15; i++) {
        if (ura_eph[i] >= value) break;
    }
    return i;
}

/* Decode GPS LNAV ephemeris and BeiDou D1 D2 ephemeris -----------------------*/
static int decode_LDeph(int sat, gtime_t toc, const double* data, eph_t* eph)
{
    eph_t eph0 = { 0 };
    int sys;

    sys = satsys(sat, NULL);

    if (!(sys & (SYS_GPS | SYS_GAL | SYS_QZS | SYS_CMP))) {
        //printf("ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    *eph = eph0; // Initialize ephemeris structure to zero

    eph->sat = sat; // Set the satellite number
    eph->toc = toc; // Set the reference time

    eph->f0 = data[0];
    eph->f1 = data[1];
    eph->f2 = data[2];

    eph->A = SQR(data[10]); eph->e = data[8]; eph->i0 = data[15]; eph->OMG0 = data[13];
    eph->omg = data[17]; eph->M0 = data[6]; eph->deln = data[5]; eph->OMGd = data[18];
    eph->idot = data[19]; eph->crc = data[16]; eph->crs = data[4]; eph->cuc = data[7];
    eph->cus = data[9]; eph->cic = data[12]; eph->cis = data[14];

    if (sys == SYS_GPS) { // GPS LNAV message
        eph->iode = (int)data[3]; // IODE (Issue of Data, Ephemeris) converted to integer
        eph->iodc = (int)data[26]; // IODC (Issue of Data, Clock) converted to integer
        eph->toes = data[11]; // TOE (Time of Ephemeris) in seconds within the GPS week
        eph->week = (int)data[21]; //GPS week number
        eph->toe = adjweek(gpst2time(eph->week, data[11]), toc);// Adjust the ephemeris reference time
        eph->ttr = adjweek(gpst2time(eph->week, data[27]), toc);// Information transmission time

        eph->code = (int)data[20]; // GPS: Codes on L2 channel 
        eph->svh = (int)data[24]; // Satellite health status
        eph->sva = uraindex(data[23]); //URA (User Range Accuracy) index
        eph->flag = (int)data[22]; // GPS: L2 P data flag

        eph->tgd[0] = data[25]; // TGD (Group Delay)
        eph->fit = data[28]; // Fit interval, the curve fitting interval for GPS ephemerides
    }
    if (eph->iode < 0 || 1023 < eph->iode) {
        printf("rinex nav invalid: sat=%2d iode=%d\n", sat, eph->iode);
        return 0;
    }
    if (eph->iodc < 0 || 1023 < eph->iodc) {
        printf("rinex nav invalid: sat=%2d iodc=%d\n", sat, eph->iodc);
        return 0;
    }
    return 1;
}
int BDTWeek;    //Storing Compass Week in RINEX 4.0 D1D2
/* decode ephemeris GPS CNAV CNV2 AND BeiDou CNV123---*/
static int decode_CNVeph(int sat, gtime_t toc, const double* data, eph_t* eph)
{
    eph_t eph0 = { 0 };
    int sys;

    sys = satsys(sat, NULL);

    if (!(sys & (SYS_GPS | SYS_GAL | SYS_QZS | SYS_CMP))) {
        //printf("ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    *eph = eph0; // Initialize ephemeris structure to zero
    // CNAV CNV123: Same information storage
    eph->sat = sat; // Set satellite number
    eph->toc = toc; // Set reference time

    eph->f0 = data[0];
    eph->f1 = data[1];
    eph->f2 = data[2];
    eph->toc = bdt2gpst(eph->toc);
    eph->dotA = data[3]; eph->crs = data[4]; eph->deln = data[5]; eph->M0 = data[6]; eph->cuc = data[7];
    eph->e = data[8]; eph->cus = data[9];  eph->A = SQR(data[10]); eph->toes = data[11]; eph->cic = data[12];
    eph->OMG0 = data[13]; eph->cis = data[14]; eph->i0 = data[15]; eph->crc = data[16]; eph->omg = data[17];
    eph->OMGd = data[18]; eph->idot = data[19]; eph->dotn = data[20]; eph->satType = data[21];
    eph->tgd[0] = data[29];      /* TGD1 B1Cp */
    eph->tgd[1] = data[30];      /* TGD2 B2ap */
    eph->tgd[2] = data[27];      //b1cd
    eph->svh = (int)data[32];      /* satH1 */
    eph->iode = (int)data[38];      /* AODE */
    eph->iodc = (int)data[34];      /* AODC */

    eph->week = BDTWeek;      /* bdt week */
    eph->toe = bdt2gpst(bdt2time(eph->week, data[11])); /* bdt -> gpst */
    eph->ttr = bdt2gpst(bdt2time(eph->week, data[35])); /* bdt -> gpst */
    eph->toe = adjweek(eph->toe, toc);
    eph->ttr = adjweek(eph->ttr, toc);
   
    if (eph->iode < 0 || 1023 < eph->iode) {
        printf("rinex nav invalid: sat=%2d iode=%d\n", sat, eph->iode);
        return 0;
    }
    if (eph->iodc < 0 || 1023 < eph->iodc) {
        printf("rinex nav invalid: sat=%2d iodc=%d\n", sat, eph->iodc);
        return 0;
    }
    return 1;
}

/* read rinex 4.0 navigation data, store information consistent with message 71 and 72 --
   Information like STO, EOP, ION are not used at the moment. Positioning has been reserved
   for future storage in respective structures, you can add the necessary structures here later.
   
   INPUT: 
   file: rinex 4.0 navigation file

   OUTPUT:
   struct: nav_t* nav (a pointer to a structure nav_t that stores navigation information extracted from the RINEX 4.0 navigation file))

   Copyright (C) GCC Group

*/
extern int readRinex4Nav(const char* file, nav_t* nav)
{
    nav->eph = NULL; nav->n = nav->nmax = 0;
    eph_t ephg = { 0 }; eph_t ephb = { 0 };    gtime_t toc;
    FILE* fp;
    double ver = 2.10, data[64] = { 0 }; char type1 = ' ';
    int  i, j, sat = 0, mask, stat = 0, sys, LD, CNV;
    char buff[MAXRNXLEN], id[8] = "", * p, * label = buff + 60;

    if (!nav)
        return 0;
    // Zero out the ephemeris structure
    memset(&ephg, 0, sizeof(eph_t));
    memset(&ephb, 0, sizeof(eph_t));
    if ((fp = fopen(file, "r")) == NULL) {
        printf("*** ERROR: open Rinex4.0 nav file failed, please check it!\n");
        return 0;
    }
    // set system mask
    mask = SYS_ALL;
    // Retrieve BeiDou week from D1D2
    while (fgets(buff, MAXRNXLEN, fp)) {
        if ((strstr(buff, "D1") || strstr(buff, "D2")) && (strstr(buff, "> EPH"))) {
            for (i = 0; i < 6; i++)
            {
                fgets(buff, MAXRNXLEN, fp);
            }
            BDTWeek = (int)str2num(buff, 43, 62);
            break;
        }
        continue;
    }
    rewind(fp);
    // Convert epoch to time; time to GPS; or directly use toc
    while (fgets(buff, MAXRNXLEN, fp)) {
        i = 0;     memset(data, 0, sizeof(data));// Zero out the data array
        if (buff[0] == '\0' || buff[0] == '\n') continue;
        if (strstr(label, "RINEX VERSION / TYPE")) {
            ver = str2num(buff, 0, 9);
            type1 = *(buff + 20);
            if (ver < 4.0)
                return 0;
            while (fgets(buff, MAXRNXLEN, fp)) {
                // Header section, only store leap seconds
                if (strstr(label, "LEAP SECONDS")) { 
                    // Optional
                    if (nav) nav->leaps = (int)str2num(buff, 0, 6);
                }
                else if (strstr(label, "END OF HEADER")) {
                    // Read one more line to enter the body of the file
                    fgets(buff, MAXRNXLEN, fp); break;
                }
                continue;
            }
        }
        // Body of the file
        if (ver >= 4.0)
        {
            //if ((strstr(buff, "LNAV"))) ephg->Inmt = 0;
            //else if ((strstr(buff, "CNVX"))) eph->Inmt = 1;
            //else if ((strstr(buff, "D1D2"))) eph->Inmt = 2;
            //else if ((strstr(buff, "SBAS"))) eph->Inmt = 3;
            if ((strstr(buff, "> STO"))) {
                //system time offset
                strncpy(id, buff + 6, 3);
                // Calculate the corresponding satellite number in the navigation message for RTKLIB
                sat = satid2no(id);              
                // Determine the satellite system based on the satellite number
                sys = satsys(sat, NULL);    
                for (i = 0; i < 2; i++)
                {// Read and discard the next two lines
                    fgets(buff, MAXRNXLEN, fp);
                }
                continue;
            }
            if (strstr(buff, "> ION")) {//ionosphere      
                strncpy(id, buff + 6, 3);
                sat = satid2no(id);
                sys = satsys(sat, NULL);
                if (strstr(buff, "IFNV"))
                {
                    for (i = 0; i < 2; i++)
                    {// Read and discard the next two lines
                        fgets(buff, MAXRNXLEN, fp);
                    }
                    continue;
                }
                else
                {
                    for (i = 0; i < 3; i++)
                    {// Read and discard the next three lines
                        fgets(buff, MAXRNXLEN, fp);
                    }
                    continue;
                }
            }
            if ((strstr(buff, "> EOP"))) {// Earth orientation parameters
                strncpy(id, buff + 6, 3);
                // Calculate the corresponding satellite number in the navigation message for RTKLIB
                sat = satid2no(id);              
                // Determine the satellite system based on the satellite number
                sys = satsys(sat, NULL);    
                for (i = 0; i < 3; i++)
                {   // Read and discard the next three lines
                    fgets(buff, MAXRNXLEN, fp);
                }
                continue;
            }
            if ((strstr(buff, "> EPH")) && (strstr(buff, "CNV1") || strstr(buff, "LNAV"))) { // Read LNAV CNV1 ephemerides
                strncpy(id, buff + 6, 3);
                // Calculate the corresponding satellite number in the navigation message for RTKLIB
                sat = satid2no(id);              
                // Determine the satellite system based on the satellite number
                sys = satsys(sat, NULL);    
                if (sys == SYS_GPS || sys == SYS_CMP)
                {
                    //message71
                    if ((strstr(buff, "LNAV"))) ephg.Enmt = 1;
                    //else if ((strstr(buff, "CNAV"))) ephg->Enmt = 2;
                    //else if ((strstr(buff, "D1")))   ephb->Enmt = 3;
                    //else if ((strstr(buff, "D2")))   ephb->Enmt = 4;
                    else if ((strstr(buff, "CNV1"))) ephb.Enmt = 5;//message72
                    //else if ((strstr(buff, "CNV2"))) ephb->Enmt = 5;
                    //else if ((strstr(buff, "CNV3"))) eph->Enmt = 6;
                    //else if ((strstr(buff, "SBAS"))) eph->Enmt = 7;
                    while (fgets(buff, MAXRNXLEN, fp)) { // Read ephemeris parameters
                        if (buff[0] == '\0' || buff[0] == '\n') continue;
                        if (i == 0) {
                            if (str2time(buff + 4, 0, 19, &toc)) {          /* Decode time field */
                                printf("rinex nav toc error: %23.23s\n", buff);
                                return 0;
                            }
                            for (j = 0, p = buff + 4 + 19; j < 3; j++, p += 19) { /* Decode data fields */
                                data[i++] = str2num(p, 0, 19);
                            }
                        }
                        else
                        {
                            for (j = 0, p = buff + 4; j < 4; j++, p += 19) {       /* Decode subsequent data fields */
                                data[i++] = str2num(p, 0, 19);
                            }
                            if ((ephg.Enmt == 1) && i >= 29) { // LNAV D1D2 together contains 31 data points
                                LD = decode_LDeph(sat, toc, data, &ephg);
                                if (LD)
                                {
                                    add_eph(nav, &ephg);
                                }
                                break;
                            }
                            else if ((ephb.Enmt == 5) && i >= 39) { // CNV1 CNV2 together contains 39 data points 
                                CNV = decode_CNVeph(sat, toc, data, &ephb);
                                if (CNV)
                                {
                                    add_eph(nav, &ephb);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(fp); return 1;
}

// Initialize B2b satellite data structure
extern  initB2b(nav_t* navs)
{
    int i, j;

    // Initialize current epoch B2b data
    navs->b2bsat.nsat = 0;
    for (i = 0; i < MaskNSAT; i++) {
        navs->b2bsat.SatSlot[i] = 0;
    }
    for (i = 0; i < MAXSAT; i++) {
        navs->b2bsat.b2bsats[i].sat = -1;

        // B2b Type1 Initialization
        navs->b2bsat.b2bsats[i].b2btype1.Iodp = -1;
        navs->b2bsat.b2bsats[i].b2btype1.iB2b = -1;
        navs->b2bsat.b2bsats[i].b2btype1.IodSsr = -1;
        navs->b2bsat.b2bsats[i].b2btype1.t0.time = 0;
        navs->b2bsat.b2bsats[i].b2btype1.t0.sec = 0;
        navs->b2bsat.b2bsats[i].b2btype1.TodBDT = 0.0;

        // B2b Type2 Initialization
        navs->b2bsat.b2bsats[i].b2btype2.IodSsr = -1;
        navs->b2bsat.b2bsats[i].b2btype2.IODN = -1;
        navs->b2bsat.b2bsats[i].b2btype2.iB2b = -1;
        navs->b2bsat.b2bsats[i].b2btype2.UraClass = -1;
        navs->b2bsat.b2bsats[i].b2btype2.UraValue = -1;
        navs->b2bsat.b2bsats[i].b2btype2.IodCorr = -1;
        navs->b2bsat.b2bsats[i].b2btype2.t0.time = 0;
        navs->b2bsat.b2bsats[i].b2btype2.t0.sec = 0;
        navs->b2bsat.b2bsats[i].b2btype2.TodBDT = 0.0;
        for (j = 0; j < 3; j++) {
            navs->b2bsat.b2bsats[i].b2btype2.OrbCorr[j] = 0.0;
        }

        // B2b Type3 Initialization
        navs->b2bsat.b2bsats[i].b2btype3.TodBDT = 0.0;
        navs->b2bsat.b2bsats[i].b2btype3.IodSsr = -1;
        navs->b2bsat.b2bsats[i].b2btype3.iB2b = -1;
        navs->b2bsat.b2bsats[i].b2btype3.t0.time = 0;
        navs->b2bsat.b2bsats[i].b2btype3.t0.sec = 0;
        for (j = 0; j < NMODESINGAL; j++) {
            navs->b2bsat.b2bsats[i].b2btype3.SatDCB[j] = 0.0;
        }

        // B2b Type4 Initialization
        navs->b2bsat.b2bsats[i].b2btype4.C0 = 0.0;
        navs->b2bsat.b2bsats[i].b2btype4.IodSsr = -1;
        navs->b2bsat.b2bsats[i].b2btype4.Iodp = -1;
        navs->b2bsat.b2bsats[i].b2btype4.IodCorr = -1;
        navs->b2bsat.b2bsats[i].b2btype4.TodBDT = 0.0;
        navs->b2bsat.b2bsats[i].b2btype4.t0.time = 0;
        navs->b2bsat.b2bsats[i].b2btype4.t0.sec = 0;
        navs->b2bsat.b2bsats[i].b2btype4.iB2b = -1;
    }

    // Initialize previous epoch B2b data
    navs->b2bsat_pre.nsat = 0;
    for (i = 0; i < MaskNSAT; i++) {
        navs->b2bsat_pre.SatSlot[i] = 0;
    }
    for (i = 0; i < MAXSAT; i++) {
        navs->b2bsat_pre.b2bsats[i].sat = -1;

        navs->b2bsat_pre.b2bsats[i].b2btype1.Iodp = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype1.iB2b = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype1.IodSsr = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype1.t0.time = 0;
        navs->b2bsat_pre.b2bsats[i].b2btype1.t0.sec = 0;
        navs->b2bsat_pre.b2bsats[i].b2btype1.TodBDT = 0.0;

        navs->b2bsat_pre.b2bsats[i].b2btype2.IodSsr = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype2.IODN = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype2.iB2b = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype2.UraClass = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype2.UraValue = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype2.IodCorr = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype2.t0.time = 0;
        navs->b2bsat_pre.b2bsats[i].b2btype2.t0.sec = 0;
        navs->b2bsat_pre.b2bsats[i].b2btype2.TodBDT = 0.0;
        for (j = 0; j < 3; j++) {
            navs->b2bsat_pre.b2bsats[i].b2btype2.OrbCorr[j] = 0.0;
        }

        navs->b2bsat_pre.b2bsats[i].b2btype3.TodBDT = 0.0;
        navs->b2bsat_pre.b2bsats[i].b2btype3.IodSsr = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype3.iB2b = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype3.t0.time = 0;
        navs->b2bsat_pre.b2bsats[i].b2btype3.t0.sec = 0;
        for (j = 0; j < NMODESINGAL; j++) {
            navs->b2bsat_pre.b2bsats[i].b2btype3.SatDCB[j] = 0.0;
        }

        navs->b2bsat_pre.b2bsats[i].b2btype4.C0 = 0.0;
        navs->b2bsat_pre.b2bsats[i].b2btype4.IodSsr = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype4.Iodp = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype4.IodCorr = -1;
        navs->b2bsat_pre.b2bsats[i].b2btype4.TodBDT = 0.0;
        navs->b2bsat_pre.b2bsats[i].b2btype4.t0.time = 0;
        navs->b2bsat_pre.b2bsats[i].b2btype4.t0.sec = 0;
        navs->b2bsat_pre.b2bsats[i].b2btype4.iB2b = -1;
    }
}

// Get satellite number from subtype index
extern int B2bSubtype2Sat(const int ind, b2bsat_t* b2bsat)
{
    int i, count = 0, sat = -1;

    for (i = 0; i < MaskNSAT; i++) {
        if (b2bsat->SatSlot[i]) {
            count++;
        }
        if (count == ind + 1) {
            sat = satSlot2Sat(i + 1);
            return sat;
        }
    }
    return -1; // Not found
}


extern void readB2bType1(FILE* fp, b2bsat_t* b2bsat, b2bsat_t* b2bsat_pre, gtime_t obstime)
{
    static long last_file_pos = 0; // Static variable to store the position of the last read in the file
    int Week, Sow, Tod, SSRGap, IODSSR, IODP, i = 0, slot = 0, count = 0, t1 = 0;
    char bdsmask[211], gpsmask[211], galmask[211], glomask[211];
    gtime_t time = { 0 };

    // Check if the file pointer is valid
    if (!fp) {
        perror("Failed to open file B2btype1");
        exit(EXIT_FAILURE);
    }

    // If there is a saved position from the last read, seek to that position to resume reading from there
    if (last_file_pos != 0) {
        if (fseek(fp, last_file_pos, SEEK_SET) != 0) {
            perror("Failed to seek to last position");
            exit(EXIT_FAILURE);
        }
    }
    char line[512];

    // Read the file line by line
    while (fgets(line, sizeof(line), fp)) {
        // Record the starting position of the current line at the beginning of each loop
        long current_line_pos = ftell(fp) - strlen(line);


        // Parse the contents of the line
        if (sscanf(line, "%d %d %d %d %d %d %s %s %s %s %s", &Week, &Sow, &Tod, &SSRGap, &IODSSR, &IODP, bdsmask, gpsmask, galmask, glomask) != 10) {
            fprintf(stderr, "An error occurred reading the B2btype1 file\n");
            continue;
        }
        // If the IODP value differs from the previous, update the b2bsat_pre structure
        if (b2bsat->b2bsats->b2btype1.Iodp != IODP)
        {
            memcpy(b2bsat_pre, b2bsat, sizeof(b2bsat_t));
        }
        // Convert the week and seconds of week into BDT time
        time = bdt2time(Week, Sow + 14);

        double ep[6];
        time2epoch(time, ep);
        // If the time difference is greater than a tolerance (DTTOL), save the file position and break the loop
        if (timediff(time, obstime) > DTTOL) // Save the file position of this line for the next read
        {
            last_file_pos = current_line_pos;
            break;
        }
        // Adjust the time based on the seconds of the week and time of day
        t1 = (Sow % 86400) - Tod < 0 ? (Sow % 86400) - Tod + 86400 : (Sow % 86400) - Tod;
        time.time = time.time - t1 + 14;
        // Process the satellite masks for different satellite systems
        // 0~62 BDS mask, 63~99 GPS mask, 100~136 Galileo mask, 137~173 GLONASS mask
        for (i = 0; i < MaskNSAT; ++i) {
            if (i <= 62)
            {
                slot = bdsmask[i] == '1' ? 1 : 0;
            }
            else if (i <= 99)
            {
                slot = gpsmask[i - 63] == '1' ? 1 : 0;
            }
            else if (i <= 136)
            {
                slot = galmask[i - 100] == '1' ? 1 : 0;
            }
            else if (i <= 173)
            {
                slot = glomask[i - 137] == '1' ? 1 : 0;
            }
            else {
                slot = 0; // Invalid index.
            }
            // Set the slot status for the satellite
            b2bsat->SatSlot[i] = slot;
            // Skip if the slot is empty (0)
            if (!slot) continue;
            if (slot < 0 || slot>1) {
                continue;

            }
            // Store the satellite information in the b2bsats array
            b2bsat->b2bsats[count].sat = satSlot2Sat(i + 1);
            b2bsat->b2bsats[count].b2btype1.Iodp = IODP;
            b2bsat->b2bsats[count].b2btype1.IodSsr = IODSSR;
            b2bsat->b2bsats[count].b2btype1.TodBDT = Tod;
            //b2bsat->b2bsats[count].b2btype1.iB2b = msg->iB2b;
            b2bsat->b2bsats[count++].b2btype1.t0 = time;
        }
        // Update the number of satellites
        b2bsat->nsat = count;
    }

    fclose(fp);

}

extern void readB2bType2(FILE* fp, b2bsat_t* b2bsat, gtime_t obstime)
{
    static long last_file_pos = 0; // Static variable to store the position of the last read in the file
    int Week, Sow, Tod, SSRGap, IODSSR, IODN, SatSlot, SatSlot1, IODCorr, URAclass, URAvalue, sat = -1, i = 0, t1, j = 0, ind = -1;
    double Rcor, Tcor, Ncor;
    BOOL isture = 1;
    gtime_t time = { 0 };

    // Check if the file pointer is valid
    if (!fp) {
        perror("Failed to open file B2btype2");
        exit(EXIT_FAILURE);
    }

    // If there is a saved position from the last read, seek to that position to resume reading from there
    if (last_file_pos != 0) {
        if (fseek(fp, last_file_pos, SEEK_SET) != 0) {
            perror("Failed to seek to last position");
            exit(EXIT_FAILURE);
        }
    }
    char line[512];

    // Read the file line by line
    while (fgets(line, sizeof(line), fp)) {
        // Record the starting position of the current line at the beginning of each loop
        long current_line_pos = ftell(fp) - strlen(line);


        // Parse the contents of the line
        if (sscanf(line, "%d %d %d %d %d %d %d %d %d %lf %lf %lf %d %d", &SatSlot1, &Week, &Sow, &Tod, &SSRGap, &IODSSR, &IODN, &SatSlot, &IODCorr, &Rcor, &Tcor, &Ncor, &URAclass, &URAvalue) != 14) {
            fprintf(stderr, "An error occurred reading the B2btype2 file\n");
            continue;
        }
        // Convert the week and seconds of week into BDT time
        time = bdt2time(Week, Sow + 14);

        double ep[6];
        time2epoch(time, ep);

        // If the time difference is greater than a tolerance (DTTOL), save the file position and break the loop
        if (timediff(time, obstime) > DTTOL) 
        {
            last_file_pos = current_line_pos; // Save the file position of this line for the next read
            break;
        }
        // Adjust the time based on the seconds of the week and time of day
        t1 = (Sow % 86400) - Tod < 0 ? (Sow % 86400) - Tod + 86400 : (Sow % 86400) - Tod;
        time.time = time.time - t1 + 14;
        // Check for valid satellite slot
        if (SatSlot < 0 || SatSlot>255) {
            continue;
        }
        sat = satSlot2Sat(SatSlot);
        // Find the index of the satellite in the b2bsat structure
        for (j = 0; j < b2bsat->nsat; j++) {
            if (sat != b2bsat->b2bsats[j].sat) continue;
            ind = j;
            break;
        }
        // If the satellite is not found, skip the current line
        if (ind == -1) {
            continue;
        }
        // Prepare the data structure for B2bType2
        B2bType2_t B2bType2;
        B2bType2.t0 = time;
        B2bType2.IodSsr = IODSSR;
        B2bType2.TodBDT = Tod;
        B2bType2.IODN = IODN;
        B2bType2.IodCorr = IODCorr;
        // Validate the IODCorr value
        if (B2bType2.IodCorr < 0 || B2bType2.IodCorr>7) {
            isture = 0;
            continue;
        }
        // Radial correction value (Rcor)
        B2bType2.OrbCorr[0] = Rcor;
        // Validate the radial correction
        if (fabs(B2bType2.OrbCorr[0]) - 26.2128 > 1E-6) {
            isture = 0;
            continue;
        }
        // Tangential correction value (Tcor) ¨C 13 bits
        B2bType2.OrbCorr[1] = Tcor;
        // Normal correction value (Ncor) ¨C 13 bits
        B2bType2.OrbCorr[2] = Ncor;
        // Validate the tangential and normal corrections
        for (int j = 1; j < 3; j++) {
            if (fabs(B2bType2.OrbCorr[j]) - 26.208 > 1E-6) {
                isture = 0;
                continue;
            }
        }
        // URA class and value validation
        B2bType2.UraClass = URAclass;
        if (B2bType2.UraClass < 0 || B2bType2.UraClass>7) {
            isture = 0;
            continue;
        }

        B2bType2.UraValue = URAvalue;
        if (B2bType2.UraValue < 0 || B2bType2.UraValue>7) {
            isture = 0;
            continue;
        }
        // B2bType2.iB2b = msg->iB2b;// If all values are valid, update the b2bsat structure
        if (isture) {
            b2bsat->b2bsats[ind].b2btype2 = B2bType2;
        }
        // Reset the flags for the next iteration
        else continue;
        ind = -1;
        isture = 1;
    }


    fclose(fp);

}

extern void readB2bType3(FILE* fp, b2bsat_t* b2bsat, gtime_t obstime)
{
    static long last_file_pos = 0; // Static variable to store the last file read position
    int SatSlot, Week, Sow, Tod, SSRGap, IODSSR, SatSlot1, CodeNum, sat = -1, i = 0, t1, j = 0, k = 0, ind = -1, slot = 0;
    double corr[16];
    int model[16];
    gtime_t time = { 0 };

    if (!fp) {
        perror("Failed to open file B2btype3");
        exit(EXIT_FAILURE);
    }

    // If there is a saved file position, start reading from there
    if (last_file_pos != 0) {
        if (fseek(fp, last_file_pos, SEEK_SET) != 0) {
            perror("Failed to seek to last position");
            exit(EXIT_FAILURE);
        }
    }
    char line[512];

    while (fgets(line, sizeof(line), fp)) {
        // Record the start position of the current line at the beginning of each loop
        long current_line_pos = ftell(fp) - strlen(line);

        // Parse the file content
        if (sscanf(line,
            "%d %d %d %d %d %d %d %d "
            "%d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf ",
            &SatSlot, &Week, &Sow, &Tod, &SSRGap, &IODSSR, &SatSlot1, &CodeNum,
            &model[0], &corr[0], &model[1], &corr[1], &model[2], &corr[2],
            &model[3], &corr[3], &model[4], &corr[4], &model[5], &corr[5],
            &model[6], &corr[6], &model[7], &corr[7], &model[8], &corr[8]
        ) != 24) {
            fprintf(stderr, "An error occurred reading the B2btype3 file\n");
            continue;
        }
        time = bdt2time(Week, Sow + 14);

        double ep[6];
        time2epoch(time, ep);

        // If the time difference between the current data and the observation time exceeds the tolerance (DTTOL),
        // save the current file position so that the next time we can resume from where we left off.
        if (timediff(time, obstime) > DTTOL) {
            last_file_pos = current_line_pos;
            break;
        }
        t1 = (Sow % 86400) - Tod < 0 ? (Sow % 86400) - Tod + 86400 : (Sow % 86400) - Tod;
        time.time = time.time - t1 + 14;
        if (SatSlot == 0) continue;
        sat = satSlot2Sat(SatSlot);

        for (j = 0; j < b2bsat->nsat; j++) {
            if (sat != b2bsat->b2bsats[j].sat) continue;
            ind = j;
            break;
        }
        for (i = 0; i < CodeNum; ++i) {
            // Update the B2btype3 structure for the corresponding satellite
            b2bsat->b2bsats[ind].b2btype3.IodSsr = IODSSR;
            b2bsat->b2bsats[ind].b2btype3.t0 = time;
            b2bsat->b2bsats[ind].b2btype3.TodBDT = Tod;

            if (ind == -1) {
                continue;
            }
            b2bsat->b2bsats[ind].b2btype3.SatDCB[model[i]] = corr[i];
        }
        ind = -1;
    }
    fclose(fp);
}

extern void readB2bType4(FILE* fp, b2bsat_t* b2bsat, gtime_t obstime)
{
    static long last_file_pos = 0; // Static variable to store the last file read position
    int SubType, SubType1, Week, Sow, Week1, Sow1, Tod, SSRGap, IODSSR, IODP, SatNum, sat = -1, i = 0, t1, j = 0, k = 0, ind = -1, slot = 0;
    double C0[MAXCLOCKCOR];
    int IODCorr[MAXCLOCKCOR];
    gtime_t time = { 0 };
    b2bsat_t* b2bsat_pre;
    b2bsat_pre = b2bsat;

    if (!fp) {
        perror("Failed to open file B2btype4");
        exit(EXIT_FAILURE);
    }

    // If there is a saved file position, start reading from there
    if (last_file_pos != 0) {
        if (fseek(fp, last_file_pos, SEEK_SET) != 0) {
            perror("Failed to seek to last position");
            exit(EXIT_FAILURE);
        }
    }
    char line[1024];

    while (fgets(line, sizeof(line), fp)) {
        // Record the start position of the current line at the beginning of each loop
        long current_line_pos = ftell(fp) - strlen(line);


        // Parse the file content
        if (sscanf(line,
            "%d %d %d %d %d %d %d %d %d "
            "%d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf "
            "%d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf "
            "%d %lf %d %lf %d %lf "
            "%d %d",
            &SubType, &Week, &Sow, &Tod, &SSRGap, &IODSSR, &IODP, &SubType, &SatNum,
            &IODCorr[0], &C0[0], &IODCorr[1], &C0[1], &IODCorr[2], &C0[2],
            &IODCorr[3], &C0[3], &IODCorr[4], &C0[4], &IODCorr[5], &C0[5],
            &IODCorr[6], &C0[6], &IODCorr[7], &C0[7], &IODCorr[8], &C0[8],
            &IODCorr[9], &C0[9], &IODCorr[10], &C0[10], &IODCorr[11], &C0[11],
            &IODCorr[12], &C0[12], &IODCorr[13], &C0[13], &IODCorr[14], &C0[14],
            &IODCorr[15], &C0[15], &IODCorr[16], &C0[16], &IODCorr[17], &C0[17],
            &IODCorr[18], &C0[18], &IODCorr[19], &C0[19], &IODCorr[20], &C0[20],
            &IODCorr[21], &C0[21], &IODCorr[22], &C0[22],
            &Week1, &Sow1) != 57) {
            fprintf(stderr, "An error occurred reading the B2btype4 file\n");
            continue;
        }
        time = bdt2time(Week, Sow + 14);

        double ep[6];
        time2epoch(time, ep);

        // If the time difference between the current data and the observation time exceeds the tolerance (DTTOL),
        // save the current file position so that the next time we can resume from where we left off.
        if (timediff(time, obstime) > DTTOL) // Save the file position of this line for the next read
        {
            last_file_pos = current_line_pos;
            break;
        }
        t1 = (Sow % 86400) - Tod < 0 ? (Sow % 86400) - Tod + 86400 : (Sow % 86400) - Tod;
        time.time = time.time - t1 + 14;
        if (IODSSR < 0 || IODSSR>3) continue;
        //IODP
        if (IODP < 0 || IODP>15) continue;
        //SubType1
        if (SubType < 0 || SubType>31)continue;
        j = SubType * MAXCLOCKCOR;
        // Clock correction values for 23 satellites, followed by reading the next 68 lines to complete the message
        for (i = 0; i < MAXCLOCKCOR; ++i) {
            // i + j = 0 to 22: First 23 satellites where the mask is set to 1
            // i + j = 23 to 45: Next 23 satellites where the mask is set to 1
            // i + j = 46 to 68: Next 23 satellites where the mask is set to 1
            // i + j = 69 to 91: Next 23 satellites where the mask is set to 1
            sat = B2bSubtype2Sat(i + j, b2bsat);
            for (k = 0; k < b2bsat->nsat; k++) {
                if (sat != b2bsat->b2bsats[k].sat) continue;
                ind = k;
                break;
            }
            if (ind == -1) {
                continue;
            }
            // Update the satellite information for the specific satellite index
            //b2bsat->b2bsats[ind].b2btype4.iB2b = msg->iB2b;
            b2bsat->b2bsats[ind].b2btype4.Iodp = IODP;
            b2bsat->b2bsats[ind].b2btype4.TodBDT = Tod;
            b2bsat->b2bsats[ind].b2btype4.IodSsr = IODSSR;
            b2bsat->b2bsats[ind].b2btype4.t0 = time;
            //IOD Corr
            b2bsat->b2bsats[ind].b2btype4.IodCorr = IODCorr[i];
            if (b2bsat->b2bsats[ind].b2btype4.IodCorr < 0 || b2bsat->b2bsats[ind].b2btype4.IodCorr>7) {
                b2bsat = b2bsat_pre;
                continue;
            }
            //C0(15bit)
            b2bsat->b2bsats[ind].b2btype4.C0 = C0[i];
            if (fabs(b2bsat->b2bsats[ind].b2btype4.C0) - 27 > 1E-6) {
                b2bsat = b2bsat_pre;
                continue;
            }
            ind = -1;

        }

    }
    fclose(fp);

}
// Calculate variance from B2b URA value
extern double varUraB2b(double* ura)
{
    double URA;

    URA = (pow(3, ura[0]) * (1 + 0.25 * ura[1]) - 1) * 1E-3; // URA calculation formula
    return SQR(URA); // Return variance (square of URA)
}

// Compute satellite position and clock bias from CNAV ephemeris
extern void eph2pos_CNAV(gtime_t time, const eph_t* eph, double* rs, double* dts, double* var)
{
    double tk, M, E, Ek, sinE, cosE, u, r, i, O, sin2u, cos2u, x, y, sinO, cosO, cosi, mu, omge;
    double xg, yg, zg, sino, coso;
    int n, sys, prn;
    double Ak, deltna, na, n0;

    if (eph->A <= 0.0) {
        rs[0] = rs[1] = rs[2] = *dts = *var = 0.0;
        return;
    }

    // Time difference from ephemeris reference epoch
    tk = timediff(time, eph->toe);
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;

    Ak = eph->A + eph->dotA * tk; // Corrected semi-major axis

    // Select system constants
    switch ((sys = satsys(eph->sat, &prn))) {
    case SYS_GAL: mu = MU_GAL; omge = OMGE_GAL; break;
    case SYS_CMP: mu = MU_CMP; omge = OMGE_CMP; break;
    default:      mu = MU_GPS; omge = OMGE;     break;
    }

    // Compute mean motion
    n0 = sqrt(mu / (eph->A * eph->A * eph->A));
    deltna = eph->deln + eph->dotn * tk / 2;
    M = eph->M0 + (n0 + deltna) * tk; // Mean anomaly

    // Solve Kepler's equation for eccentric anomaly
    for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER && n < MAX_ITER_KEPLER; n++) {
        Ek = E;
        E -= (E - eph->e * sin(E) - M) / (1.0 - eph->e * cos(E));
    }
    if (n >= MAX_ITER_KEPLER) return; // Kepler iteration failed

    sinE = sin(E);
    cosE = cos(E);

    // Compute satellite position in orbital plane
    u = atan2(sqrt(1.0 - eph->e * eph->e) * sinE, cosE - eph->e) + eph->omg;
    r = Ak * (1.0 - eph->e * cosE);
    i = eph->i0 + eph->idot * tk;
    sin2u = sin(2.0 * u);
    cos2u = cos(2.0 * u);

    // Apply harmonic corrections
    u += eph->cus * sin2u + eph->cuc * cos2u;
    r += eph->crs * sin2u + eph->crc * cos2u;
    i += eph->cis * sin2u + eph->cic * cos2u;

    x = r * cos(u);
    y = r * sin(u);
    cosi = cos(i);

    // Special treatment for BeiDou GEO satellites
    if (sys == SYS_CMP && prn <= 5) {
        O = eph->OMG0 + eph->OMGd * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);

        xg = x * cosO - y * cosi * sinO;
        yg = x * sinO + y * cosi * cosO;
        zg = y * sin(i);

        sino = sin(omge * tk);
        coso = cos(omge * tk);

        rs[0] = xg * coso + yg * sino * COS_5 + zg * sino * SIN_5;
        rs[1] = -xg * sino + yg * coso * COS_5 + zg * coso * SIN_5;
        rs[2] = -yg * SIN_5 + zg * COS_5;
    }
    else {
        O = eph->OMG0 + (eph->OMGd - omge) * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        rs[0] = x * cosO - y * cosi * sinO;
        rs[1] = x * sinO + y * cosi * cosO;
        rs[2] = y * sin(i);
    }

    // Satellite clock bias correction
    tk = timediff(time, eph->toc);
    *dts = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;

    // Apply relativity correction
    *dts -= 2.0 * sqrt(mu * eph->A) * eph->e * sinE / SQR(CLIGHT);

    // Estimate position and clock error variance
    *var = var_uraeph(eph->sva);
}

// Select B2b ephemeris data closest to given time
extern eph_t* selB2beph(gtime_t time, int sat, int iodn, const nav_t* nav)
{
    double t, tmax, tmin;
    int i, j = -1;

    // Set maximum allowed time gap depending on satellite system
    switch (satsys(sat, NULL)) {
    case SYS_QZS: tmax = MAXDTOE_QZS + 1.0; break;
    case SYS_GAL: tmax = MAXDTOE_GAL + 1.0; break;
    case SYS_CMP: tmax = MAXDTOE_CMP + 1.0; break;
    default:      tmax = MAXDTOE + 1.0; break;
    }
    tmin = tmax + 1.0;

    // Search for the best matching ephemeris
    for (i = 0; i < nav->n; i++) {
        if (nav->eph[i].sat != sat) continue;
        if (iodn >= 0 && nav->eph[i].iodc != iodn) continue;
        if ((t = fabs(timediff(nav->eph[i].toe, time))) > tmax) continue;
        if (iodn >= 0) return nav->eph + i;
        if (t <= tmin) { j = i; tmin = t; }
    }
    if (iodn >= 0 || j < 0) return NULL;
    return nav->eph + j;
}

// Calculate satellite position, velocity and clock error for B2b service
extern int ephpos_B2b(gtime_t time, gtime_t teph, int sat, const nav_t* nav, int iode, double* rs, double* dts, double* var, int* svh)
{
    eph_t* eph;
    geph_t* geph;
    double rst[3], dtst[1], tt = 1E-3;
    int i, sys;

    sys = satsys(sat, NULL);
    *svh = -1;

    if (sys == SYS_CMP) {
        if (!(eph = selB2beph(teph, sat, iode, nav))) return 0;
        eph2pos_CNAV(time, eph, rs, dts, var);
        time = timeadd(time, tt);
        eph2pos_CNAV(time, eph, rst, dtst, var);
        *svh = eph->svh;
    }
    else if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS) {
        if (!(eph = selB2beph(teph, sat, iode, nav))) return 0;
        eph2pos(time, eph, rs, dts, var);
        time = timeadd(time, tt);
        eph2pos(time, eph, rst, dtst, var);
        *svh = eph->svh;
    }
    else if (sys == SYS_GLO) {
        // Future extension: GLONASS not handled here
    }
    else return 0;

    // Differentiate position and clock to compute velocity and clock drift
    for (i = 0; i < 3; i++) rs[i + 3] = (rst[i] - rs[i]) / tt;
    dts[1] = (dtst[0] - dts[0]) / tt;

    return 1;
}

extern int satpos_B2b(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    double* rs, double* dts, double* var, int* svh)
{
    const b2bsatp_t* b2bsatp;
    b2bsat_t b2bsat = nav->b2bsat;
    eph_t* eph;
    double t1, t2, t4, er[3], ea[3], ec[3], rc[3], dant[3] = { 0 }, tk;
    int i, sys, slot, index, num, count = 0, IodCorr4 = -1, ind = -1;
    int maskSlot[MaskNSAT];


    for (i = 0; i < nav->b2bsat.nsat; i++) {
        b2bsatp = nav->b2bsat.b2bsats + i;
        if (b2bsatp->sat == sat)break;
    }
    if (i >= nav->b2bsat.nsat) {
        // If the satellite is not found, compute satellite position using broadcast ephemeris
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
        *svh = -1;
        return 0;
    }
    if (b2bsatp->b2btype4.Iodp != b2bsatp->b2btype1.Iodp) { // If Iodp doesn't match, use older ephemeris data

        for (i = 0; i < nav->b2bsat_pre.nsat; i++) {
            b2bsatp = nav->b2bsat_pre.b2bsats + i;
            if (b2bsatp->sat == sat)break;
        }
        if (i >= nav->b2bsat_pre.nsat) {
            // If the satellite is not found, compute satellite position using broadcast ephemeris
            ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
            *svh = -1;
            return 0;
        }
        if (b2bsatp->b2btype4.Iodp != b2bsatp->b2btype1.Iodp)
        {
            ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
            *svh = -1;
            return 0;
        }
    }
    double OrbCorr[3] = { 0 };
    OrbCorr[0] = b2bsatp->b2btype2.OrbCorr[0];
    OrbCorr[1] = b2bsatp->b2btype2.OrbCorr[1];
    OrbCorr[2] = b2bsatp->b2btype2.OrbCorr[2];
    double C0 = b2bsatp->b2btype4.C0; // Clock bias correction term
    // IOD (Issue of Data) matching check
    if (b2bsatp->b2btype1.IodSsr != b2bsatp->b2btype2.IodSsr &&
        b2bsatp->b2btype1.IodSsr != b2bsatp->b2btype4.IodSsr) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);// If not matching, directly calculate satellite position
        *svh = -1;
        return 0;
    }
    if (b2bsatp->b2btype2.IodCorr != b2bsatp->b2btype4.IodCorr) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);// If not matching, directly calculate satellite position
        *svh = -1;
        return 0;
    }


    t2 = timediff(time, b2bsatp->b2btype2.t0);
    t4 = timediff(time, b2bsatp->b2btype4.t0);

    if (fabs(t2) > MAXAGEB2b || fabs(t4) > MAXAGEB2b_CLOCK) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh); // If the satellite data is too old, directly calculate satellite position
        *svh = -1;
        return 0;
    }

    //* satellite postion and clock by broadcast ephemeris */
    if (!ephpos_B2b(time, teph, sat, nav, b2bsatp->b2btype2.IODN, rs, dts, var, svh)) return 0; // If matching, compute satellite position


    /* satellite clock for gps, galileo and qzss */
    sys = satsys(sat, NULL);
    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP) {
        if (!(eph = selB2beph(teph, sat, b2bsatp->b2btype2.IODN, nav))) return 0;

        /* satellite clock by clock parameters */
        tk = timediff(time, eph->toc);
        dts[0] = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;
        dts[1] = eph->f1 + 2.0 * eph->f2 * tk;

        /* relativity correction */
        dts[0] -= 2.0 * dot(rs, rs + 3, 3) / CLIGHT / CLIGHT; // Apply relativity correction to clock
    }
    /* Calculate radial, along-track, and cross-track unit vectors in ECEF */
    if (!normv3(rs + 3, ea)) return 0;
    cross3(rs, rs + 3, rc);
    if (!normv3(rc, ec)) {
        *svh = -1;
        return 0;
    }
    cross3(ea, ec, er); // Calculate radial, cross, and along-track unit vectors for corrections

    /* Satellite antenna offset correction */
    satantoff(time, rs, sat, nav, dant);

    for (i = 0; i < 3; i++) {
        rs[i] += -(er[i] * OrbCorr[0] + ea[i] * OrbCorr[1] + ec[i] * OrbCorr[2]); // Apply satellite position corrections
    }
    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
    if (fabs(fabs(C0) - 26.2128) < 1E-6) {
        C0 = 0;
    }
    dts[0] -= C0 / CLIGHT;

    /* Variance by B2b URA (User Range Accuracy) */
    double Ura[2] = { b2bsatp->b2btype2.UraClass,b2bsatp->b2btype2.UraValue };
    *var = varUraB2b(Ura); // Calculate the variance using B2b URA class and value
    return 1;
}
extern double prange_dualfrequency(const obsd_t* obs, const nav_t* nav, double* var
    , const double* dantr, const double* dants, const prcopt_t* opt)
{
    double gamma, tgd1 = 0.0, tgd2 = 0.0, tgd = 0, t = 0, dcb1 = 0, dcb2 = 0, dcb = 0;
    double freq[NFREQ] = { 0 };
    int i;
    const b2bsatp_t* b2bsatp;
    // Calculate the frequencies for each signal from the observation data
    for (i = 0; i < NFREQ; i++)
    {
        freq[i] = sat2freq(obs->sat, obs->code[i], nav);
    }
    // Calculate the ionospheric-free combination coefficient: f1^2 / f2^2
    gamma = SQR(freq[0]) / SQR(freq[1]); /* f1^2/f2^2 */

    // Adjustments for double frequency (L1 and L2) pseudo-range values with 
    // satellite and receiver-specific delays.
    double P1 = obs->P[0] - dantr[0] - dants[0];
    double P2 = obs->P[1] - dantr[1] - dants[1];
    int sat = obs->sat, sys;

    // Find the satellite in the broadcast satellite data (b2bsat)
    for (i = 0; i < nav->b2bsat.nsat; i++) {
        b2bsatp = nav->b2bsat.b2bsats + i;
        if (b2bsatp->sat == sat)break;
    }
    // If the satellite is not found, calculate the ionospheric-free pseudo-range without 
    // any correction from DCB (Differential Code Bias).
    if (i >= nav->b2bsat.nsat) {
        return  ((P2 - gamma * P1) - (dcb2 - gamma * dcb1)) / (1.0 - gamma);
    }
    // Calculate the time difference between the observation time and the satellite time
    t = timediff(obs->time, b2bsatp->b2btype3.t0);
    // If the time difference is too large, use the simplified calculation for the ionospheric-free range
    if (fabs(t) > MAXAGEB2b_CBIAS) return  ((P2 - gamma * P1) - (dcb2 - gamma * dcb1)) / (1.0 - gamma);
    // Determine the satellite system type
    if (!(sys = satsys(sat, NULL))) return 0.0;
    *var = 0;
    // Check if the dual-frequency data is valid (i.e., non-zero values for both frequencies)
    if (P1 == 0 || P2 == 0)return 0.0; // Incomplete dual-frequency data, cannot process
    if (P1 == 0 || P2 == 0)return 0.0; // Incomplete dual-frequency data, cannot process
    // Check if the satellite¡¯s ionospheric correction (IodSsr) is valid
    if (b2bsatp->b2btype3.IodSsr != b2bsatp->b2btype2.IodSsr &&
        b2bsatp->b2btype3.IodSsr != b2bsatp->b2btype4.IodSsr &&
        b2bsatp->b2btype3.IodSsr != b2bsatp->b2btype1.IodSsr)return  ((P2 - gamma * P1) - (dcb2 - gamma * dcb1)) / (1.0 - gamma);

    // Specific handling for GPS system to account for differential code bias (DCB)
    if (sys == SYS_GPS) {
        // Loop through each code and apply corresponding DCB values
        for (int i = 0; i < 2; i++) {
            if (obs->code[i] == CODE_L1C) dcb = b2bsatp->b2btype3.SatDCB[0];//L1C/A
            if (obs->code[i] == CODE_L1P) dcb = b2bsatp->b2btype3.SatDCB[1];//L1P
            if (obs->code[i] == CODE_L1L) dcb = b2bsatp->b2btype3.SatDCB[4];//L1C(P)
            if (obs->code[i] == CODE_L1X) dcb = b2bsatp->b2btype3.SatDCB[5];//L1C(D+P)
            if (obs->code[i] == CODE_L2L) dcb = b2bsatp->b2btype3.SatDCB[7];//L2C(L)
            if (obs->code[i] == CODE_L2X) dcb = b2bsatp->b2btype3.SatDCB[8];//L2C(M+L)
            if (obs->code[i] == CODE_L5I) dcb = b2bsatp->b2btype3.SatDCB[11];//L5I
            if (obs->code[i] == CODE_L5Q) dcb = b2bsatp->b2btype3.SatDCB[12];//L5Q
            if (obs->code[i] == CODE_L5X) dcb = b2bsatp->b2btype3.SatDCB[12];//L5 I+Q
            // Store DCB values for both frequencies
            if (i == 0)  dcb1 = dcb;
            else dcb2 = dcb;
        }
        // Return the ionospheric-free pseudo-range corrected for DCB
        return ((P2 - gamma * P1) - (dcb2 - gamma * dcb1)) / (1.0 - gamma);
    }
    else if (sys == SYS_QZS) {
    }
    else if (sys == SYS_GLO) { /* G1-G2 */
        return (P2 - gamma * P1) / (1.0 - gamma);
    }
    // Specific handling for BeiDou (COMPASS) system
    else if (sys == SYS_CMP) {
        // Loop through each code and apply corresponding DCB values for BeiDou
        for (int i = 0; i < 2; i++) {
            if (obs->code[i] == CODE_L2I) dcb = b2bsatp->b2btype3.SatDCB[0];//B1I
            if (obs->code[i] == CODE_L1D) dcb = b2bsatp->b2btype3.SatDCB[1];//B1C(D)
            if (obs->code[i] == CODE_L1P) dcb = b2bsatp->b2btype3.SatDCB[2];//B1C(P)
            if (obs->code[i] == CODE_L5D) dcb = b2bsatp->b2btype3.SatDCB[4];//B2a(D)
            if (obs->code[i] == CODE_L5P) dcb = b2bsatp->b2btype3.SatDCB[5];//B2a(P)
            if (obs->code[i] == CODE_L7I) dcb = b2bsatp->b2btype3.SatDCB[7];//B2b-I
            if (obs->code[i] == CODE_L7Q) dcb = b2bsatp->b2btype3.SatDCB[8];//B2b-Q
            if (obs->code[i] == CODE_L6I) dcb = b2bsatp->b2btype3.SatDCB[12];//B3I
            // Store DCB values for both frequencies
            if (i == 0)  dcb1 = dcb;
            else dcb2 = dcb;
        }
        // Return the ionospheric-free pseudo-range corrected for DCB
        return ((P2 - gamma * P1) - (dcb2 - gamma * dcb1)) / (1.0 - gamma);
    }
}


