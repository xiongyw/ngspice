/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000 AlansFixes
**********/

#ifndef RES
#define RES

#include "ngspice/ifsim.h"
#include "ngspice/cktdefs.h"
#include "ngspice/gendefs.h"
#include "ngspice/complex.h"
#include "ngspice/noisedef.h"

/* definitions used to describe resistors */


/* information used to describe a single instance */

typedef struct sRESinstance {
    struct sRESmodel *RESmodPtr;            /* backpointer to model */
    struct sRESinstance *RESnextInstance;   /* pointer to next instance of
                                             * current model*/

    IFuid RESname;      /* pointer to character string naming this instance */
    int RESstate;       /* not used but needed for sructure consistency */
    const int RESposNode;     /* number of positive node of resistor */
    const int RESnegNode;     /* number of negative node of resistor */

    double REStemp;     /* temperature at which this resistor operates */
    double RESdtemp;    /* delta-temperature of a particular instance  */
    double RESconduct;  /* conductance at current analysis temperature */
    double RESresist;   /* resistance at temperature Tnom */
    double REScurrent;  /* The dc current in the resistor */
    /* serban */
    double RESacResist;             /* AC resistance, useful for fancy .ac analyses */
    double RESacConduct;            /* AC conductance */
    double RESwidth;                /* width of the resistor */
    double RESlength;               /* length of the resistor */
    double RESscale;                /* Scale factor */
    double RESm;                    /* Multiplicity factor for this instance */
    double REStc1;                  /* first temperature coefficient of resistors */
    double REStc2;                  /* second temperature coefficient of resistors */
    double REStce;                  /* exponential temperature coefficient of resistors */
    double RESbv_max;               /* Maximum resistor voltage */
    int    RESnoisy;                /* Set if the resistor generates noise */
    double RESeffNoiseArea;         /* effective resistor area for noise calculation */
    double *RESposPosPtr;           /* pointer to sparse matrix diagonal at
                                     * (positive,positive) */
    double *RESnegNegPtr;           /* pointer to sparse matrix diagonal at
                                     * (negative,negative) */
    double *RESposNegPtr;           /* pointer to sparse matrix offdiagonal at
                                     * (positive,negative) */
    double *RESnegPosPtr;           /* pointer to sparse matrix offdiagonal at
                                     * (negative,positive) */
    unsigned RESresGiven    : 1;    /* flag to indicate resistance was specified */
    unsigned RESwidthGiven  : 1;    /* flag to indicate width given */
    unsigned RESlengthGiven : 1;    /* flag to indicate length given */
    unsigned RESscaleGiven  : 1;    /* flag to indicate scale given */
    unsigned REStempGiven   : 1;    /* indicates temperature specified */
    unsigned RESdtempGiven  : 1;    /* indicates delta-temp specified  */
    /* serban */
    unsigned RESacresGiven  : 1;    /* indicates AC value specified */
    unsigned RESmGiven      : 1;    /* indicates M parameter specified */
    unsigned REStc1Given    : 1;    /* indicates tc1 parameter specified */
    unsigned REStc2Given    : 1;    /* indicates tc2 parameter specified */
    unsigned REStceGiven    : 1;    /* indicates tce parameter specified */
    unsigned RESnoisyGiven  : 1;    /* indicates if noisy is specified */
    unsigned RESbv_maxGiven : 1;    /* flags indicates maximum voltage is given */
    int    RESsenParmNo;            /* parameter # for sensitivity use;
                                     * set equal to  0 if not a design parameter*/

    /* indices to array of RES noise sources */

#define RESTHNOIZ  0     /* Thermal noise source */
#define RESFLNOIZ  1     /* Flicker noise source */
#define RESTOTNOIZ 2     /* Total noise          */

#define RESNSRCS   3     /* the number of RES noise sources */


#ifndef NONOISE
    double RESnVar[NSTATVARS][RESNSRCS];
#else /* NONOISE */
    double **RESnVar;
#endif /* NONOISE */

#ifdef KLU
    BindElement *RESposPosBinding ;
    BindElement *RESnegNegBinding ;
    BindElement *RESposNegBinding ;
    BindElement *RESnegPosBinding ;
#endif

} RESinstance ;

#ifdef USE_CUSPICE

/* 
 * D[] means double, I[] means integer;
 * each pointer is an array, and each instance has a slot in the array;
 * CPU/GPU version are identical, for copying back and forth.
 */


typedef struct sRESparamCPUstruct {
    double *REScpuPointersD [8] ;
    #define REStc1Array REScpuPointersD[0]
    #define REStc2Array REScpuPointersD[1]
    #define RESmArray REScpuPointersD[2]
    #define RESconductArray REScpuPointersD[3]
    #define REStempArray REScpuPointersD[4]
    #define RESdtempArray REScpuPointersD[5]
    #define REScurrentArray REScpuPointersD[6]
    #define RESgValueArray REScpuPointersD[7]

    int *REScpuPointersI [5] ;
    #define REStc1GivenArray REScpuPointersI[0]
    #define REStc2GivenArray REScpuPointersI[1]
    #define RESmGivenArray REScpuPointersI[2]
    #define RESposNodeArray REScpuPointersI[3]
    #define RESnegNodeArray REScpuPointersI[4]
} RESparamCPUstruct ;

typedef struct sRESparamGPUstruct {
    double *REScudaPointersD [8] ;
    #define d_REStc1Array REScudaPointersD[0]
    #define d_REStc2Array REScudaPointersD[1]
    #define d_RESmArray REScudaPointersD[2]
    #define d_RESconductArray REScudaPointersD[3]
    #define d_REStempArray REScudaPointersD[4]
    #define d_RESdtempArray REScudaPointersD[5]
    #define d_REScurrentArray REScudaPointersD[6]
    #define d_RESgValueArray REScudaPointersD[7]

    int *REScudaPointersI [5] ;
    #define d_REStc1GivenArray REScudaPointersI[0]
    #define d_REStc2GivenArray REScudaPointersI[1]
    #define d_RESmGivenArray REScudaPointersI[2]
    #define d_RESposNodeArray REScudaPointersI[3]
    #define d_RESnegNodeArray REScudaPointersI[4]
} RESparamGPUstruct ;
#endif

/* per model data */

typedef struct sRESmodel {       /* model structure for a resistor */
    int RESmodType; /* type index of this device type */
    struct sRESmodel *RESnextModel; /* pointer to next possible model in
                                     * linked list */
    RESinstance * RESinstances; /* pointer to list of instances that have this
                                 * model */
    IFuid RESmodName;       /* pointer to character string naming this model */

#ifdef USE_CUSPICE
    unsigned int has_cuda:1 ;
#endif

    /* --- end of generic struct GENmodel --- */

    double REStnom;         /* temperature at which resistance measured */
    double REStempCoeff1;   /* first temperature coefficient of resistors */
    double REStempCoeff2;   /* second temperature coefficient of resistors */
    double REStempCoeffe;   /* exponential temperature coefficient of resistors */
    double RESsheetRes;     /* sheet resistance of devices in ohms/square */
    double RESdefWidth;     /* default width of a resistor */
    double RESdefLength;    /* default length of a resistor */
    double RESnarrow;       /* amount by which device is narrower than drawn */
    double RESshort;        /* amount by which device is shorter than drawn */
    double RESfNcoef;       /* Flicker noise coefficient */
    double RESfNexp;        /* Flicker noise exponent */
    double RESres;          /* Default model resistance */
    double RESbv_max;       /* Maximum resistor voltage */
    double RESlf;           /* length exponent for noise calculation */
    double RESwf;           /* width exponent for noise calculation */
    double RESef;           /* frequncy exponent for noise calculation */
    unsigned REStnomGiven       :1; /* flag to indicate nominal temp. was given */
    unsigned REStc1Given        :1; /* flag to indicate tc1 was specified */
    unsigned REStc2Given        :1; /* flag to indicate tc2 was specified */
    unsigned REStceGiven        :1; /* flag to indicate tce was specified */
    unsigned RESsheetResGiven   :1; /* flag to indicate sheet resistance given*/
    unsigned RESdefWidthGiven   :1; /* flag to indicate default width given */
    unsigned RESdefLengthGiven  :1; /* flag to indicate default length given */
    unsigned RESnarrowGiven     :1; /* flag to indicate narrow effect given */
    unsigned RESshortGiven      :1; /* flag to indicate short effect given */
    unsigned RESfNcoefGiven     :1; /* flag to indicate kf given */
    unsigned RESfNexpGiven      :1; /* flag to indicate af given */
    unsigned RESresGiven        :1; /* flag to indicate model resistance given */
    unsigned RESbv_maxGiven     :1; /* flags indicates maximum voltage is given */
    unsigned RESlfGiven         :1; /* flags indicates lf is given */
    unsigned RESwfGiven         :1; /* flags indicates wf is given */
    unsigned RESefGiven         :1; /* flags indicates ef is given */

#ifdef USE_CUSPICE
    RESparamCPUstruct RESparamCPU ;
    RESparamGPUstruct RESparamGPU ;

    int offset ;
    int n_values ;
    int n_Ptr ;
    int *PositionVector ;
    int *d_PositionVector ;

    int n_instances ;
#endif

} RESmodel;

/* device parameters */
#define RES_RESIST 1
#define RES_WIDTH 2
#define RES_LENGTH 3
#define RES_CONDUCT 4
#define RES_RESIST_SENS 5
#define RES_CURRENT 6
#define RES_POWER 7
#define RES_TEMP 8
/* serban */
#define RES_ACRESIST 10
#define RES_ACCONDUCT 11
#define RES_M 12 /* pn */
#define RES_SCALE 13 /* pn */
#define RES_DTEMP 14 /* pn */
#define RES_NOISY 15 /* pn */
/* tanaka */
#define RES_TC1 16
#define RES_TC2 17
#define RES_BV_MAX 18
#define RES_TCE 19

/* model parameters */
#define RES_MOD_TC1 101
#define RES_MOD_TC2 102
#define RES_MOD_RSH 103
#define RES_MOD_DEFWIDTH 104
#define RES_MOD_DEFLENGTH 105
#define RES_MOD_NARROW 106
#define RES_MOD_R 107
#define RES_MOD_TNOM 108
#define RES_MOD_SHORT 109
#define RES_MOD_KF 110
#define RES_MOD_AF 111
#define RES_MOD_BV_MAX 112
#define RES_MOD_LF 113
#define RES_MOD_WF 114
#define RES_MOD_EF 115
#define RES_MOD_TCE 116

/* device questions */
#define RES_QUEST_SENS_REAL      201
#define RES_QUEST_SENS_IMAG      202
#define RES_QUEST_SENS_MAG       203
#define RES_QUEST_SENS_PH        204
#define RES_QUEST_SENS_CPLX      205
#define RES_QUEST_SENS_DC        206

/* model questions */

#include "resext.h"

#endif /*RES*/
