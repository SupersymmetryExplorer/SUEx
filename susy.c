#include <stdio.h>
#include <math.h>

#include "variables.h"
#include "readwrite.h"
#include "functions.h"
#include "electroweak.h"
#include "loop.h"
#include "numericx.h"
#include "radiative.h"
#include "initcond.h"
#include "higgs.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  This is the MAIN routine of the program, to be used as it is or to be called by any other routine,as discussed below.
//	 The routine has the following four important input control parameters:
//

//  IKNOWL: which sets the degree of control on the various parts of the algorithm. It has two possible values:
//
//  -- IKNOWL = 0: 	blind use of the program, no warning and other messages.
//  						default values are set for the control parameters and the program
//  						gives just the results from the physical input.
//
//  -- IKNOWL = 1:  	some warning/error messages are collected in the SuSpect2.out file
//  						(this is the recommended choice).
// ------------------------------------------------------------------------------------------------------------------------------------

//  INPUT: sets input (and output) control, it offers now 4 possibilities:
//
//  =0 : model and option input parameters ONLY read in file suspect2.in.
//			(output generated in both suspect2.out and SLHA format suspect2_lha.out)
//
//  =1 : define yourself IN suspect2_call.f all relevant input and parameters.
//			(i.e. NO reading of input files):
//  		see example of  input in accompanying file suspect2_call.f
//  		Maybe more convenient e.g. for scan over the model parameter space.
//  		(output generated in both suspect2.out and SLHA format suspect2_lha.out)
//
//  =2 : same as input =0 but read SLHA format input file: suspect2_lha.in
// 		(it writes also all output in the SLHA format file: suspect2_lha.out)
//
//  =11: same as input=1, but NO output file(s) suspect*.out generated
// ------------------------------------------------------------------------------------------------------------------------------------

//  ICHOICE: initializes the various options for the models to be considered, the degree of accuracy to be required, the features to be
//           included, etc. There are 10 possible choice at present and the options are described in more details in the input file:
//
//  NB: ICHOICE[..] superseded if using SLHA input file mode: in this case we follow SLHA standards and conventions (now adapted to ver. 2).
//
//  	-- ICHOICE[1]  : Choice of the model to be considered.
//  	-- ICHOICE[2]  : choice of perturbative order (1 or 2 loop) of the RGEs.
//  	-- ICHOICE[3]  : To impose or not the GUT scale.
//  	-- ICHOICE[4]  : Accuracy of the RGEs.
//  	-- ICHOICE[5]  : To impose or not the radiative EWSB.
//  	-- ICHOICE[6]  : Chose different input in the Higgs sector (MA,MU,Mhu,Mhd)
//  	-- ICHOICE[7]  : Choice of radiative corrections to (s)particles masses.
//  	-- ICHOICE[8]  : Choice of the EWSB scale
//  	-- ICHOICE[9]  : Accuracy of the physical spectrum calculation.
//  	-- ICHOICE[10] : Different options for calculation of Higgs boson masses.
//  	-- ICHOICE[11] : Running/pole H masses used in loops.
// ------------------------------------------------------------------------------------------------------------------------------------

//  ERRMESS: which provides a useful set of warning/error message flags, that are automatically written in the output file SUSPECT2.OUT:
//
//  	-- ERRMESS[i] = 0:  Everything is fine.
//  	-- ERRMESS[1] =-1:  Tachyon 3rd gen. sfermion from RGE
//  	-- ERRMESS[2] =-1:  Tachyon 1,2 gen. sfermion from RGE
//  	-- ERRMESS[3] =-1:  Tachyon A    (maybe temporary: see final mass)
//  	-- ERRMESS[4] =-1:  Tachyon 3rd gen. sfermion from mixing
//  	-- ERRMESS[5] =-1:  Mu inconsistent (or unstable) after many iterations
//  	-- ERRMESS[6] =-1:  Non-convergent mu from EWSB
//  	-- ERRMESS[7] =-1:  EWSB maybe inconsistent  (!but RG-improved only check)
//  	-- ERRMESS[8] =-1:  V_Higgs maybe UFB or CCB (!but RG-improved only check)
//  	-- ERRMESS[9] =-1:  Higgs boson masses are NaN
//  	-- ERRMESS[10]=-1:  RGE problems (non-pert and/or Landau poles)
// --------------------------------------------------------------------------------------------------------------------------------------

// ========================== //
//  The program starts here:  //
// ========================== //

void susy(int iknowl, int input, int ichoice[], double errmess[])
{
    int icount, iremember, inorc;
    int imod[3];

    unsigned int n = 31;

    double y[n + 1];
    double errnogo;
    double pizz_mz;
    double tachsqrc;
    double bup;
    double alfinv, alphas, mt, mc, qewsb, ehigh, mhu2, mhd2;
    double ma2;
    double dvdvu2, dvdvd2;
    double mapole2;

    double au1, ad1, al1, al0, ad0, au0, al10, ad10, au10;
    double mtaur0, msl0, mbr0, mtr0, msq0, mer0, mel0, mdr0, mur0, muq0, m10, m20, m30;

    int rmhalf;
    double rm0;

    FILE *fout, *fin;
    FILE *finlha;

    //
    //		Initializing various control parameters + other parameters:
    //  ----------------------------------------------------------------------------------------------

    for (int ierr = 1; ierr <= 10; ierr++)
        errmess[ierr] = 0.0; // Initialize all the error terms to 0.0

    errnogo = 0.0;
    irge = 0;
    iflop = 0;
    tachsqrc = 0.0;
    icount = 0;
    iremember = 0;

    for (int irg = 1; irg <= 31; irg++)
        y[irg] = 0.0;

    aaml = 0.0;
    ytauewsb = 0.0;
    ybewsb = 0.0;
    ytewsb = 0.0;
    pizz_mz = 0.0;

    pizzp = 0.0;
    inorc = 0;

    inonpert = 0; // added for non-pert pbs control

    bup = 0.0;
    sterr = 0.0;
    sberr = 0.0;
    stauerr = 0.0;
    stnuerr = 0.0;
    errma2z = 0.0;

    // further reinitializations added
    alsewsb = 0.0;
    g2ewsb = 0.0;
    g1ewsb = 0.0;
    vuewsb = 0.0;
    vdewsb = 0.0;

    aamh = 0.0;
    aamch = 0.0;

    marun = 0.0;
    mapole = 1.0; // initialization at 1st call (value later superseded)

    piaa = 0.0;
    tadba = 0.0;
    d2ma = 0.0;

    kmaflag = 0;
    imbmb = 0;

    dmc1 = 0.0;
    dmc2 = 0.0;
    dmn1 = 0.0;
    dmn2 = 0.0;
    dmn3 = 0.0;
    dmn4 = 0.0;
    mgluino = 0.0;

    dmst1 = 0.0;
    dmst2 = 0.0;
    dmsu1 = 0.0;
    dmsu2 = 0.0;
    dmsb1 = 0.0;
    dmsb2 = 0.0;
    dmsd1 = 0.0;
    dmsd2 = 0.0;
    dmsl1 = 0.0;
    dmsl2 = 0.0;
    dmse1 = 0.0;
    dmse2 = 0.0;
    dmsn1 = 0.0;
    dmsntau = 0.0;

    thetout = 0.0;
    thebout = 0.0;
    thelout = 0.0;

    dml = 0.0;
    dmh = 0.0;
    dmch = 0.0;
    alfa = 0.0;

    // open and read the input file
    //------------------------------------------------------------------------------------------------------

    //  read input:
    //  physical input parameters:

    if (input == 0)
    {
        fin = fopen("susy.in", "r");

        read_n_line(fin, 10);
        fscanf(fin, "%d", &ichoice[1]);
        read_n_line(fin, 4);
        fscanf(fin, "%d", &ichoice[2]);
        read_n_line(fin, 4);
        fscanf(fin, "%d", &ichoice[3]);
        read_n_line(fin, 4);
        fscanf(fin, "%d", &ichoice[4]);
        read_n_line(fin, 4);
        fscanf(fin, "%d", &ichoice[5]);
        read_n_line(fin, 4);
        fscanf(fin, "%d", &ichoice[6]);
        read_n_line(fin, 6);
        fscanf(fin, "%d", &ichoice[7]);
        read_n_line(fin, 4);
        fscanf(fin, "%d", &ichoice[8]);
        read_n_line(fin, 3);
        fscanf(fin, "%d", &ichoice[9]);
        read_n_line(fin, 6);
        fscanf(fin, "%d", &ichoice[10]);
        read_n_line(fin, 6);
        fscanf(fin, "%d", &ichoice[11]);
        read_n_line(fin, 4);
        fscanf(fin, "%lf %lf %lf %lf %lf", &alfinv, &alphas, &mt, &mbmb, &mtau);
        read_n_line(fin, 3);
        fscanf(fin, "%lf %lf", &ehigh, &qewsb);
        read_n_line(fin, 4);

        //
        //  Minimal Supergravity (mSUGRA) Input
        //  ----------------------------------------------------------------------------------------------------------------------------------
        if (ichoice[1] == 10)
        {
            fscanf(fin, "%lf %d %lf %lf %lf", &rm0, &rmhalf, &a0, &tgbeta, &sgnmu0);
            m0 = rm0;
            mhalf = rmhalf;
        }

        //
        // Gauge mediated supersymmetry breaking (GMSB)  input:
        //  ----------------------------------------------------------------------------------------------------------------------------------
        else if (ichoice[1] == 11)
        {
            read_n_line(fin, 4);
            fscanf(fin, "%lf %lf %lf %lf %d %d", &mgmmess, &mgmsusy, &tgbeta, &sgnmu0, &nl, &nq);
        }

        //
        // Anomaly mediated supersymmetry breaking (AMSB) input:
        //  ----------------------------------------------------------------------------------------------------------------------------------
        else if (ichoice[1] == 12)
        {
            read_n_line(fin, 8);
            fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &m32, &am0, &tgbeta, &sgnmu0, &cq, &cu, &cd, &cl, &ce, &chu, &chd);
        }

        //
        // i.e. non-universal arbitrary input case
        //  ----------------------------------------------------------------------------------------------------------------------------------
        else
        {
            read_n_line(fin, 12);
            fscanf(fin, "%lf %lf %lf %lf", &mhu2, &mhd2, &tgbeta, &sgnmu0);

            read_n_line(fin, 2);
            fscanf(fin, "%lf %lf %lf", &m1, &m2, &m3);

            read_n_line(fin, 2);
            fscanf(fin, "%lf %lf %lf %lf %lf", &msl, &mtaur, &msq, &mtr, &mbr);

            read_n_line(fin, 2);
            fscanf(fin, "%lf %lf %lf %lf %lf", &mel, &mer, &muq, &mur, &mdr);

            read_n_line(fin, 2);
            fscanf(fin, "%lf %lf %lf %lf %lf %lf", &al, &au, &ad, &al1, &au1, &ad1);

            read_n_line(fin, 3);
            fscanf(fin, "%lf %lf", &aama, &mu);

            mu = sgnmu0 * fabs(mu); // Add to avoid inconsistent user's input
        }

        fclose(fin);

        if (gf == 0.0)
            gf = 1.16639e-5; // only if not already defined
        if (mz == 0.0)
            mz = 91.187; // only if not already defined
    }

    int igut_in;

    if (ichoice[1] == 10 || ichoice[1] == 1)
        igut_in = 1;
    else
        igut_in = 0;

    double tgbet0, beta_z;
    int ihrcsave;
    int iaccsave;
    int inorge, kmaflag;
    int irgmax, irgsave;

    double cw2, sw, cw, rmtau, rmtau2, mb, alphas0, g32;
    //
    // Stage 1
    // -----------------------------------------------------------------------------------------------

    //
    // Case where MHu and MHd are not input. These values are initilized due to some precosionary
    // measure so that nothing get diverged in the middle of the program
    // ------------------------------------------------------------------------------------------------

    // Added reinitialization of mhu2,mhd2 for scans:

    if (ichoice[1] != 2 && ichoice[6] == 0)
    {
        mhu2 = 1.0e4;
        mhd2 = 1.0e4;
        // (nb ichoice[6]=0 -> ma,mu input thus these mhu2,mhd2 values are
        // irrelevant but only initialized for convergence of iteration control)
    }

    //
    // Stage 2
    // -----------------------------------------------------------------------------------------------

    //
    // An Introduction to Supersymmetry - Manuel Drees (arXiv:hep-ph/9611409v1)
    // tan(beta) = v_bar/v    Eq.(66)
    // -----------------------------------------------------------------------------------------------

    ihflag = ichoice[10]; // Higgs mass flag.
    ihrcsave = ihflag;
    ipolemz = ichoice[11]; // Higgs masses at loop-level at mZ
    tgbet0 = tgbeta;
    beta_z = atan(tgbet0);

    if (ichoice[1] == 12)
        rm0 = am0; // Choice of model : AMSB   (cMSSM)
                   //  blind use: assign protection default values to control parameters:
    if (ichoice[2] == 0)
        ichoice[2] = 11; // All the RGEs are at 1-loop order

    //  (i.e. 1-loop rge at first run by default, if 2-loop not chosen)
    //	 :essai        if(ichoice[3] == 0 && ichoice[1] != 2) ichoice[3]=1
    //  (i.e. gut scale always consitently recalculated as g1(gut)=g2(gut)										// Check what is trying to do

    if (ichoice[4] != 1 && ichoice[4] != 2) // RGE flag
        ichoice[4] = 1;
    iaccsave = ichoice[4];
    // (i.e. protections against wrong or undefined rge accuracy setup)

    if (ichoice[5] != 1) //	No radiative EWSB imposed (only in pMSSM) : 0
        ichoice[5] = 1;  // Consistent EWSB (automatic in cMSSMs)     : 1
                         // (i.e. always ewsb)

    if (ichoice[8] != 0) // Default EWSB scale=(mt_L*mt_R)^(1/2)       : 1
        ichoice[8] = 1;  // Arbitrary EWSB scale (to be given below)  : 0
                         // (i.e to be sure that not choosing default ewsb scale
                         //  =(m_t_l*m_t_r)^(1/2) is on purpose)
                         // choose frozen scale in rge parameters:

    kpole = ichoice[8];   // for some no rge purposes:
    inorge = ichoice[1];  // for special case where ma(pole) is really input:
    kmaflag = ichoice[6]; // choose susy r.c. options:

    if (ichoice[7] == 1) // SUSY radiative corrections to the (s)particles masses
        isfrc = 0;       // only mt,mb,mtau susy r.c.
    else if (ichoice[7] == 2)
        isfrc = 1; //  mt,mb,mtau  + (all) squarks + (all) gaugino susy rc:

    // optimize number of long (rg+ full spectrum) iterations
    irgmax = 50;
    irgsave = irgmax;

    //
    // Stage 3
    // -----------------------------------------------------------------------------------------------

    //
    // An Introduction to Supersymmetry - Manuel Drees (arXiv:hep-ph/9611409v1)
    // Signature of m1,m2 and m3    Eq.(63)
    // -----------------------------------------------------------------------------------------------

    if (ichoice[1] <= 2) // Signs of masses in different models
    {
        // one could have arbitrary m1,m2,m3 signs
        sgnm1 = m1 / fabs(m1);
        sgnm2 = m2 / fabs(m2);
        sgnm3 = m3 / fabs(m3);
    }
    else if (ichoice[1] == 10 || ichoice[1] == 11)
    {
        // msugra (or gmsb) case
        sgnm1 = 1.0;
        sgnm2 = 1.0;
        sgnm3 = 1.0;
    }

    else if (ichoice[1] == 12)
    {
        // amsb case
        sgnm1 = 1.0;
        sgnm2 = 1.0;
        sgnm3 = -1.0;
    }

    //
    // Stage 4
    // -----------------------------------------------------------------------------------------------

    //
    // A Supersymmetry Primer - Stephen P Martin (arXiv:hep-ph/9709356v6)
    //
    // -----------------------------------------------------------------------------------------------

    //  rename input parameters for slha  and other input choice purpose

    // tan(beta): the ratio of the vevs of the two.Higgs doublet fields.
    // mHu2 ,mHd2 : the Higgs mass parameters squared.
    // M1,M2,M3: the bino, wino and gluino mass parameters.
    // mq,muR,mdR,ml,meR: the first/second generation sfermion mass parameters.
    // mQ,mtR,mbR,mL ,mtauR: the third generation sfermion mass parameters.
    // Au,Ad,Ae: the first/second generation trilinear couplings.
    // At,Ab,Atau : the third generation trilinear couplings.

    double mhd20, mhu20, mu0;

    if (!(input != 2 && ichoice[1] == 10))
    {
        mhu20 = mhu2; // 6.3.1
        mhd20 = mhd2; // 6.3.1

        msl0 = msl;
        mtaur0 = mtaur;
        msq0 = msq;
        mtr0 = mtr;
        mbr0 = mbr;

        mel0 = mel;
        mer0 = mer;
        muq0 = muq;
        mur0 = mur;
        mdr0 = mdr;

        al0 = al; // 6.3.1
        au0 = au; // 6.3.1
        ad0 = ad; // 6.3.1
        al10 = al1;
        au10 = au1;
        ad10 = ad1;
        mu0 = mu;

        m10 = m1;
        m20 = m2;
        m30 = m3;
    }

    beta = atan(tgbeta);

    // some other basic parameter definitions

    cpi = 1.0 / (16 * pi * pi);

    //
    // Gf (in GeV^-1) : the Fermi constant determined from the muon decay
    // Mz (in GeV)    : the pole mass of the Z boson
    //--------------------------------------------------------------------------

    if (gf == 0.0)
        gf = 1.16639e-5; // only if not already defined
    sw2 = 0.22210;       // only starting value sw2_drbar calculated below

    if (mz == 0.0)
        mz = 91.187; // only if not already defined

    zm = mz;

    //
    // Guess starting point for susym , elow, ehigh scales:
    // Set E_Low to Mz
    //----------------------------------------------------------------------------
    delow = mz;
    if (ichoice[3] != 0)
        ehigh = 1.0e17;

    if (ichoice[1] == 10 || ichoice[1] == 12)
        susym = 0.5 * (rm0 + rmhalf) + mz;

    else if (ichoice[1] == 11)
    {
        susym = mz;
        rm0 = susym;
    }
    // nb this "rm0" is not m0, only used at 1rst iter to guess mu(0),b(0)

    else if (ichoice[1] == 1)
    {
        rm0 = (msl + mtaur + msq + mtr + mbr + mel + mer + muq + mur + mdr) / 10.0;
        susym = 0.5 * (rm0 + (m1 + m2 + m3) / 3.0) + mz;
    }

    double gut; // GUT time scale
    gut = ehigh;
    kunif = ichoice[3];
    wistep = 1.0e2;
    nf = 6.0;

    //
    //  mw, sw^2 msbar at z scale (values may be changed):
    //-----------------------------------------------------------------------------
    cw2 = 1.0 - sw2; // Sw2 : 0.22210 (initiallized before)
    sw = sqrt(sw2);
    cw = sqrt(cw2);
    mw = mz * cw;
    wm = mw;
    rmtau = mtau;
    rmtau2 = rmtau * rmtau;

    //  some saving
    mc = 1.42;
    mc0 = mc;
    mb0 = 4.9; // Value just used for very first initialization
    mb = mb0;
    mt0 = mt;
    mbpole = mb;
    mtpole = mt;
    mtaurun = mtau;
    mbrun = mb0;
    mtrun = mt0;

    // (initial values only! at 1st iter mrun = mpole)
    // passing from alpha_s(mz) msbar to alpha_s(mz) drbar:

    //
    //	Stage 5
    //--------------------------------------------------------------------------

    alphas0 = alphas;
    g32 = 4 * pi * alphas0 / (1.0 - (1.0 / 2) * alphas0 / (2 * pi));
    alphas = g32 / 4 / pi;

    // (nb value in fact used at first rg only, does not include susy etc r.c.)
    // passing from alpha(mz) msbar to alpha(mz) drbar:

    double e2, sw20, cw20, g120, g220, acc, rmbms;
    double rmb, rmb2, rmt2;
    int nloop;

    alpha = 1.0 / (alfinv - 1.0 / pi / 6.0);
    // (nb value used at first rg iteration only, does not include susy etc r.c.)

    //
    // Stage 6
    // Introduction to Elementary Particles : David Griffiths
    // Electroweak Unification : P346
    // --------------------------------------------------------------------------

    e2 = 4 * pi * alpha;
    sw20 = sw2;
    cw20 = 1.0 - sw20;
    g12 = e2 / cw20;
    g22 = e2 / sw20;
    g120 = g12;
    g220 = g22;

    acc = 1.e-8;
    nloop = 2;
    nnlo = 1;
    idrflag = 0;

    // Calaulate the Lambda
    // xlambda is the scale Lambda. Check in arXiv:hep-ph/06072009v2
    xlambda = xitla(nloop, alphas0, acc); // alphas0 is msbar

    // Consider the effect of masses : charm, bottom, top quark
    // Check in arXiv:hep-ph/06072009v2  Page : 17
    n0 = 5;
    alsini(acc);

    // Calcuttate the varying mass of the quarks
    // Ref: arXiv:hep-ph/9511344v1; Physical Review D, Vol 43, No 5; Physics LeRe~B 341 (1994) 73-83
    // --------------------------------------------------------------------------------------------------
    imbmb = 0;
    rmbms = runm(mz, 5);

    mb = mbpole;
    mb0 = mb;

    // rmbms is mb running(mz) in msbar scheme

    mc = runm(mz, 4); // Defining the charm quark mass

    // now defining running quark masses in drbar at z scale:

    rmb = rmbms * (1.0 - alphas / (3 * pi) - 23 * pow(alphas, 2) / (72 * pow(pi, 2)) + 3 * g22 / (128 * pow(pi, 2)) + 13 * g12 / (1152 * pow(pi, 2)));

    // rmb is mb running(mz) in drbar scheme (what is mostly used after)
    // xlambda=xitla(nloop,alphas0,acc)
    // call alsini(acc)

    rmb2 = rmb * rmb;
    rmtop = runm(mt, 6); // Top quark mass
    rmt2 = rmtop * rmtop;

    if (ichoice[1] == 2 && ichoice[6] == 1)
    {
        iremember = 1;
        ichoice[1] = 0; // a trick to simplify the bottom-up case
    }

    int nok, nbad, ifix;
    double eps;
    double mtlog, delmt;
    double temp_pizz[1] = {0.0}, temp_piww[1] = {0.0}, temp_piww0[1] = {0.0}, temp_m3z[1] = {0}, temp_alphadr[1] = {0};
    double pizz, piww, piww0;
    double m3z, alphadr;
    double pass_sw2[1], pass_alphadr[1], pass_alphas[1];
    double vd2, cbeta, sbeta, vu2, vd_mz, vu_mz;
    double su_deriv1, su_rkqc, su_deriv2;
    double ysave[32];
    double mtaugut, mbgut, mtgut;
    double db, rmu0, b0;
    double rmhu2old;
    double rmel, rmdr, rmur, rmuq, rmer;
    double rmu, bold, rmuold, b;
    double rmino1, rmino2, rmino3;
    double ewsb2;
    double c2beta, wm2, zm2;
    double rmst12, rmst22, rmsb12, rmsb22, rmstau12, rmstau22;
    double rmhd2old, sb2, cb2, mzdr2, madr2, rmhu2, rmhd2;
    double errhuold, errhdold, errstop;
    double madr2old;
    double r1, r2, r3, test1, test2, test3;
    double mhu2old;
    double alz, adz, auz, mtaurz, mslz, mbrz, mtrz, msqz, merz, mdrz, murz, muqz, melz, mglu;
    double delgino;
    double m1z, m2z, mtausave, mbsave, mtsave, b_mz, mu_mz;
    double msntau_mz, delmb;
    double delmtau, delmtop;
    double dal, dau, dad, dal1, dau1, dad1, dtgbeta;
    double dmhu2, dmhd2, dm1, dm2, dm3, dma, dmsl;
    double dmtaur, dmsq, dmtr, dmbr, dmel, dmer;
    double dmuq, dmur, dmdr, dmu;
    double errhu, errhd;

    //
    // Stage 7
    // long iteration (on rge + spectrum once defined) starts here:
    // ----------------------------------------------------------------------------------------------------------------------------

l44:
    irge = irge + 1;

    //  reinitialize at each rge iter higgs rc choice (1 or 2 loop):
    ihflag = ihrcsave;

    // reinitialize  error messages until last iteration:
    if (irge <= irgmax)
        for (int i = 1; i <= 10; i++)
            errmess[i] = 0.0; // No error till

    tbeta = tgbet0;

    // calculating s^2_w_drbar(mz), g1_drbar(mz), g2_drbar(mz) incl. susy r.c.:
    if (irge >= 2)
    {
        // (because at first call no susy physical masses etc are defined)
        // First need to compute q=pizz(mz), q=piww(mz):
        scale = mz;

        // Precision corrections in the minimal supersymmetric standard model
        //	By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang
        //	arXiv:0709.1075v1  [hep-ph]
        // Calculating the self energy of W and Z bosons. These will be required for calculating The running of the coupling constants

        su_pixx(sw2, sqrt(g22), sqrt(g12), tbeta, &pizz, &piww, &piww0, 0.0); // pixx with pole mt
        pizz_mz = pizz;

        // Precision corrections in the minimal supersymmetric standard model
        //		By D.M. Pierce, Jonathan A. Bagger, K.T. Matchev and Ren-jie Zhang
        //		arXiv:0709.1075v1  [hep-ph]
        // Now the more complete calculation of g1,g2,sw2 (mz) in drbar:
        su_runningcp(alphas0, mt, rmtop, m3z, tbeta, pizz, piww, piww0, &alphadr, &alphas, &sw2);

        e2 = 4 * pi * alphadr;
        cw2 = 1.0 - sw2;
        // following redef of sw etc added

        sw = sqrt(sw2);
        cw = sqrt(cw2);
        mw = sqrt((pow(mz, 2) + pizz) * cw2 - piww); // Cos^2(theta_w) = Mw^2/Mz^2
        wm = mw;
        g12 = e2 / cw2;
        g22 = e2 / sw2;
        g32 = 4 * pi * alphas;
    }

    //
    // Higgs vev at z scale: tbeta = vu/vd
    // (nb in our normalization mz = (g12+g22)/2*(vu2+vd2),
    // and there are no factors of sqrt(2.0) in the higgs doublet components
    // (cf ramond et al prd49(1994) 4882)
    // ----------------------------------------------------------------------------------------------------------------------------

    if (irge == 1)
        pizz = 0.0;
    else
    {
        su_pixx(sw2, sqrt(g22), sqrt(g12), tbeta, &pizz, &piww, &piww0, rmtop); // pizz with running mt
        pizz_mz = pizz;
    }

    if (isnan(pizz) || pow(mz, 2) + pizz <= 0.0)
    {
        //
        // protections added
        // non-pert or nan pb, uses tree-level values temporarily:
        // -------------------------------------------------------------------------------------------------------------------------

        pizz = 0.0;

        if (irge == irgmax)
            inonpert = -1;
    }

    //
    // Stage 8
    // ---------------------------------------------------------------------------------------------------------------------------

    //
    //  Renormalizing group study of the standard model and its extensions: The minimal supersymmetric standard model :
    //  D.J.Castano, E.J.Piard and P.Ramond
    //  Eq. (4.13)
    // ---------------------------------------------------------------------------------------------------------------------------

    // This is the first thing to be called for irge = 1

    vd2 = 2 * (mz * mz + pizz) / (g12 + g22) / (1.0 + tbeta * tbeta);
    cbeta = 1.0 / sqrt(1.0 + pow(tbeta, 2));
    sbeta = tbeta * cbeta;
    vu2 = vd2 * tbeta * tbeta;
    vd = sqrt(vd2);
    vu = sqrt(vu2);

    vd_mz = vd;
    vu_mz = vu;

    //
    // defining yukawa couplings at z scale:
    // ----------------------------------------------------------------------------------------------------------------------------

    if (irge == 1)
    {
        y[4] = mtau / vd; // qcd corrections to mt(mz) (yt(mz)=y[6]) in drbar including logs:
        mtlog = log(pow(mt / mz, 2));
        delmt = alphas / pi * (5.0 / 3 - mtlog) + pow(alphas, 2) * (0.8760 - 0.384 * mtlog + 0.038 * pow(mtlog, 2));

        y[6] = mt / vu * (1.0 - delmt);
        y[5] = rmb / vd;
    }

    //
    // Higgs vev at z scale: y[7]=ln vu, y[8]=ln vd
    // ----------------------------------------------------------------------------------------------------------------------------

    y[7] = 0.5 * log(vu2);
    y[8] = 0.5 * log(vd2);

    //
    // 1st stage: evolution of gauge + yukawa cpl from mz to gut:
    // for irge=1 (iter. 1) yukawa's determined from qcd corrections only
    // ----------------------------------------------------------------------------------------------------------------------------

    y[1] = 5.0 * g12 / 3.0; // (i.e usual su[5] normalisation of g1)
    y[2] = g22;
    y[3] = g32;

    //
    // set rge accuracy choices (3 different)
    // ----------------------------------------------------------------------------------------------------------------------------

    if (ichoice[4] == 0)
    {
        h1 = 0.20; // probably the stepsize of the integration process. Decreasing step size will make the process slower.
        eps = 1.e-3;
    }
    else if (ichoice[4] == 1)
    {
        h1 = 0.06;
        eps = 1.0e-3;
    }
    else if (ichoice[4] == 2)
    {
        h1 = 0.010;
        if (ichoice[1] == 0 || ichoice[1] == 2)
            h1 = 0.005; // more precise rge for pmssm
        eps = 2.0e-5;
    }

    if (ichoice[3] != 0 && irge == 1)
        ehigh = 1.0e17;

    //
    // Note ehigh = 1.e17 will be superseded by true unification scale (where y[1]=y[2] by def.):
    // ----------------------------------------------------------------------------------------------------------------------------

    double x1, x2;

    if (ichoice[1] == 0)
    {
        //
        // Case where only mass spectrum at EWSB scale is calculated:
        // It is then assumed that all mssm parameters are defined at EWSB scale, except tanbeta(mz). The EWSB scale is an input
        // arbitrarily chosen, and the only rge performed is to calculate the gauge+yukawa+vevs from their input values at mz
        // scale to their consistent values at ewsb scale.
        // -------------------------------------------------------------------------------------------------------------------------

        if (ichoice[8] == 0)
        {
            if (qewsb == 0.0)
                qewsb = 1.05 * zm;
        }

        //
        // (protections in case of badly chosen ewsb scale input in this case)
        // -------------------------------------------------------------------------------------------------------------------------

        else
        {
            if (irge == 1)
                qewsb = sqrt(msq * mtr);
            else
                qewsb = sqrt(msttr1 * msttr2);

            if (qewsb < mz)
                qewsb = mz + 1.0e-1; // Added protection
        }

        x1 = log(zm);
        x2 = log(qewsb);
    }

    else
    {
        //
        // means all other cases where rge is performed from mz to gut scales
        // ----------------------------------------------------------------------------------------------------------------------------

        if (ichoice[8] == 1)
        {
            if (irge == 1)
                qewsb = mz;
            else
            {
                qewsb = sqrt(msttr1 * msttr2);
                if (qewsb < mz)
                    qewsb = mz + 1.e-1; // added protection
            }
        }
        x1 = log(zm);
        x2 = log(1.0e20);

        // essai       if(ichoice[1] == 2 && ehigh != 0.0) x2=log(ehigh)

        if (ichoice[3] == 0 && ehigh != 0.0) // GUT scale imposed. Higher value has been provided.
            x2 = log(ehigh);
    }

    ifirst = 0;
    jfirst = 0;
    scale = qewsb;

    //
    // First step: run from mz to high scale with initial conditions
    // g_i(mz), yukawa_i, find gut scale etc.
    // ----------------------------------------------------------------------------------------------------------------------------

    if (ichoice[2] == 11) // All the RGEs are at 1-loop order (faster)
        su_odeint(y, n, x1, x2, eps, h1, 1.0e-8, 1);
    else if (ichoice[2] == 21) // 2-loop RGEs for gauge+Yukawas+gauginos
        su_odeint(y, n, x1, x2, eps, h1, 1.0e-8, 2);

    // Stage 9
    // ----------------------------------------------------------------------------------------------------------------------------

    // protection against rge num. pbs (landau poles, non-perturbativity):
    if (iflop == 1)
    {
        errmess[10] = -1.0;
        writing(); //	goto 801;	// this has to be clearly written
    }

    if (ichoice[1] == 0)
    {
        g1ewsb = sqrt(3 * y[1] / 5);
        g2ewsb = sqrt(y[2]);
        alsewsb = y[3] / 4 / pi;
        ytauewsb = y[4];
        ybewsb = y[5];
        ytewsb = y[6];
        vuewsb = exp(y[7]);
        vdewsb = exp(y[8]);
        tbeta = vuewsb / vdewsb;
        goto l880;
    }

// (exact) gauge (g1=g2) unif. if required:
l882:
    if (egut != 0.0 && ichoice[3] != 0)
    {
        ehigh = egut;
        for (int irg = 1; irg <= 31; irg++)
            y[irg] = ygut[irg];
        y[2] = y[1];
    }

    for (int i = 1; i <= 8; i++)
        ysave[i] = y[i];

    vu = exp(y[7]);
    vd = exp(y[8]);
    mtaugut = vd * y[4];
    mbgut = vd * y[5];
    mtgut = vu * y[6];

    if (ichoice[1] == 2 && irge == irgmax) // "bottom-up" RGE (arbitrary MSSM input at low energy)
    {
        mhu2 = y[12];
        mhd2 = y[13];
        mtaur = sqrt(y[14]);
        msl = sqrt(y[15]);
        mbr = sqrt(y[16]);
        mtr = sqrt(y[17]);
        msq = sqrt(y[18]);
        mer = sqrt(y[24]);
        mel = sqrt(y[25]);
        mdr = sqrt(y[26]);
        mur = sqrt(y[27]);
        muq = sqrt(y[28]);
        al = y[9];
        ad = y[10];
        au = y[11];
        al1 = y[29];
        ad1 = y[30];
        au1 = y[31];

        db = y[19];
        mu = sgnmu0 * exp(y[23]);
        m1 = sgnm1 * exp(y[20]);
        m2 = sgnm2 * exp(y[21]);
        m3 = sgnm3 * exp(y[22]);
        writing();
    }

    //
    // 2d stage: Evolution from high (gut) scale down to low energy
    // ----------------------------------------------------------------------------------------------------------------------------

    // now taking input rmu0,b0 values (!only guess initialization values)

    if (ichoice[6] == 0 && irge > 1)
    {
        mhu2 = y[12];
        mhd2 = y[13];
        // i.e. for ma_pole,mu(ewsb) input: in this case mhu2,mhd2(gut) should not
        // be reinitialized (except at first rge iteration where they are undefined)
    }

    //
    // guess mu(gut) value at first time run (later superseeded by ewsb mu)
    // this apply in particular in msugra or non-univ cases:
    // ----------------------------------------------------------------------------------------------------------------------------

    if (rm0 == 0.0)
        rm0 = 1.0e-4; // protection
    if (irge == 1)
        rmu0 = 1.1 * rm0;

l7:

    rmu0 = sgnmu0 * fabs(rmu0);

    // also guess value for b(gut):
    b0 = 2 * rm0;
    // set up boundary conditions at gut scale:
    // yukawa coupling (eventual) unification at gut scale:

    if (ichoice[3] >= 2)
    {
        y[5] = y[4];
        ysave[5] = y[5];
    }

    // (tau- b unification)
    if (ichoice[3] == 3)
    {
        y[6] = y[5];
        ysave[6] = y[6];
    }

    // (tau-b-top unification):
    // caution then tanbeta is constrained
    // ( not yet operative  )

    // Higgs initial vev at gut scale: fixed from evolution from z scale (see above)

l2:
    icount = icount + 1;

    // icount is a counter for things to be done only at first iteration
    errhu = 1.0e3;
    errhd = 1.0e3;
    ifix = 0;

    // errhu,hd,ifix are convergence control parameters for ichoice[6] = 0
    // if on different high-energy input (msugra, amsb,gmsb,..) starts here:
    //	ichoice[1] = 12;
l77:
    if (ichoice[1] == 1 || ichoice[1] == 10)
    {
        //  unconstrained mssm: general case; ! now includes sugra case with
        //  universality as special case in same algorithm

        y[9] = al0;
        y[10] = ad0;
        y[11] = au0;
        y[29] = al10;
        y[30] = ad10;
        y[31] = au10;
        y[12] = mhu2;
        y[13] = mhd2;
        y[14] = pow(mtaur0, 2);
        y[15] = pow(msl0, 2);
        y[16] = pow(mbr0, 2);
        y[17] = pow(mtr0, 2);
        y[18] = pow(msq0, 2);

        if (irge == 1)
            y[19] = b0;
        if (irge == 1)
            y[23] = log(fabs(rmu0));
        y[24] = pow(mer0, 2);
        y[25] = pow(mel0, 2);
        y[26] = pow(mdr0, 2);
        y[27] = pow(mur0, 2);
        y[28] = pow(muq0, 2);
        y[20] = log(fabs(m10));
        y[21] = log(fabs(m20));
        y[22] = log(fabs(m30));
    }

    //
    // Bottom-up rge case with soft(ewsb) input:
    // Initialize reasonable gut scale (guess) values at firt rge iteration to catch convergence
    // ----------------------------------------------------------------------------------------------------------------------------

    else if (ichoice[1] == 2 && irge == 1)
    {
        for (int j = 9; j <= 11; j++)
            y[j] = 0.0;

        for (int j = 29; j <= 31; j++)
            y[j] = 0.0;

        for (int kk = 12; kk <= 18; kk++)
            y[kk] = 100.0;

        for (int kk = 24; kk <= 28; kk++)
            y[kk] = 100.0;

        y[19] = b0;

        for (int l = 20; l <= 22; l++)
            y[l] = log(100.0);

        y[23] = 0.0;
    }

    //
    // AMSB case:  m3/2, m0,c_q, etc (coeffs of m0),sgn(mu0) input)
    // ----------------------------------------------------------------------------------------------------------------------------

    else if (ichoice[1] == 12)
    {
        double pass_m[3];

        pass_m[0] = m1;
        pass_m[1] = m2;
        pass_m[2] = m3;
        su_amsbsub(rm0, m32, cq, cu, cd, cl, ce, chu, chd, y, pass_m);

        m10 = pass_m[0];
        m20 = pass_m[1];
        m30 = pass_m[2];

        y[20] = log(fabs(m10));
        y[21] = log(fabs(m20));
        y[22] = log(fabs(m30));
        // remaining needed parameters:
        if (irge == 1)
            y[19] = b0;
        if (irge == 1)
            y[23] = log(fabs(rmu0));

        // forces radiative EW breaking (if was not chosen before:)
        ichoice[5] = 1;
    }

    double xewsb;
    x1 = log(ehigh);

    //
    // Generic end of running scale: determined "consistently" by default:
    //  - at mz scale for gauge+yukawas couplings, that are defined at mz,
    //  - at ewsb scale (!by default = sqrt(mst_l*mst_r)).
    //  for all others rg parameters
    // ----------------------------------------------------------------------------------------------------------------------------

    if (scale == 0.0)
        scale = mz + 1.0e-1; // protection against undefined
    xewsb = log(scale);
    x2 = log(mz);
    h1 = -h1;

l3:
    int issb, istab;
    issb = 0;
    istab = 0;
    ifirst = 0;
    jfirst = 0;

    if (ichoice[1] != 11)
    {
        //  rge is made in two steps from gut scale to ewsb; then mz
        if (ichoice[2] == 11)
            su_odeint(y, n, x1, xewsb, eps, h1, 1.e-8, 1);
        else if (ichoice[2] == 21)
            su_odeint(y, n, x1, xewsb, eps, h1, 1.e-8, 2);
    }
    //
    // This is GMSB case:
    // Input are mgmmess,mgmgsusy,nl,nq, sgn(mu) and tbeta) but intermediate (messenger) scale: mgmmess for rge of soft terms
    // ----------------------------------------------------------------------------------------------------------------------------

    else
    {
        double xint;
        if (irge == 1)
            y[19] = b0;
        if (irge == 1)
            y[23] = log(fabs(rmu0));

        xint = log(mgmmess);

        double pass_m[3];
        if (ichoice[2] == 11)
            su_odeint(y, n, x1, xint, eps, h1, 1.e-8, 1);
        else if (ichoice[2] == 21)
            su_odeint(y, n, x1, xint, eps, h1, 1.e-8, 2);

        pass_m[0] = m1;
        pass_m[1] = m2;
        pass_m[2] = m3;

        // - now input soft susy-breaking terms at messenger scale:
        su_gmsbsub(mgmmess, mgmsusy, nl, nq, y, pass_m);

        m10 = pass_m[0];
        m20 = pass_m[1];
        m30 = pass_m[2];

        y[20] = log(fabs(m10));
        y[21] = log(fabs(m20));
        y[22] = log(fabs(m30));

        for (int i = 29; i <= 31; i++)
            y[i] = 0.0;

        // next rge down to ewsb scale: forces as usual radiative ew breaking:

        ichoice[5] = 1;

        if (ichoice[2] == 11)
            su_odeint(y, n, xint, xewsb, eps, h1, 1.e-8, 1);
        else if (ichoice[2] == 21)
            su_odeint(y, n, xint, xewsb, eps, h1, 1.e-8, 2);
    }

    //
    // (last endif = end of non-univ/msugra/gmsb/amsb distinctions)
    // protection against big troubles in rge (e.g. landau poles):
    // ----------------------------------------------------------------------------------------------------------------------------

    if (iflop == 1)
    {
        errmess[10] = -1.0;
        writing();
    }

    //
    // check for problems (non-perturbativity or/and landau poles) in rge:
    // ----------------------------------------------------------------------------------------------------------------------------

    for (int i = 1; i <= 31; i++)
    {
        if (isnan(y[i]))
            errmess[10] = -1.0;
    }

    if (errmess[10] == -1.0)
        writing();

    if (ichoice[1] == 2)
    {
        //
        // New algorithm for EWSB soft terms input with bottom-up rge unconstrained mssm:
        // general case: note in this case al0 etc are supposed to be soft terms input values at low ewsb input scale
        // -------------------------------------------------------------------------------------------------------------------------
        y[9] = al0;
        y[10] = ad0;
        y[11] = au0;
        y[29] = al10;
        y[30] = ad10;
        y[31] = au10;

        if (ichoice[6] == 1 || irge == 1)
        {
            y[12] = mhu2;
            y[13] = mhd2;
        }

        // (nb otherwise it means that ma_pole and mu(ewsb) are input)
        y[14] = mtaur0 * mtaur0;
        y[15] = msl0 * mtaur0;
        y[16] = mbr0 * mtaur0;
        y[17] = mtr0 * mtaur0;
        y[18] = msq0 * mtaur0;

        if (irge == 1)
            y[19] = b0;
        if (irge == 1)
            y[23] = log(fabs(rmu0));

        y[24] = mer0 * mer0;
        y[25] = mel0 * mer0;
        y[26] = mdr0 * mer0;
        y[27] = mur0 * mer0;
        y[28] = muq0 * mer0;
        y[20] = log(fabs(m10));
        y[21] = log(fabs(m20));
        y[22] = log(fabs(m30));
    }

    vu = exp(y[7]);
    vd = exp(y[8]);

    // saving all rge parameters at ewsb scale:

    for (int ip = 1; ip <= 31; ip++)
        yewsb[ip] = y[ip];

    // saving also yukawas and others at ewsb scale:

l886:
    ytauewsb = y[4];
    ybewsb = y[5];
    ytewsb = y[6];
    alsewsb = y[3] / (4 * pi);
    g2ewsb = sqrt(y[2]);
    g1ewsb = sqrt(3 * y[1] / 5);
    vuewsb = exp(y[7]);
    vdewsb = exp(y[8]);
    tbeta = exp(y[7]) / exp(y[8]);
    atau = y[9];
    ab = y[10];
    atop = y[11];
    al1 = y[29];
    ad1 = y[30];
    au1 = y[31];
    rmhu2 = y[12];
    rmhd2 = y[13];

    // Stage 10
    // --------------------------------------------------------------------------------------------------------------------------

    //
    // Change (after 10 iter) of standard fixed point algorithm:
    // mhu_new = mhu_ewsb -> mhu_new = (1-c)*mhu_old + c*mhu_ewsb, c=0.3
    // This trick may help recovering convergence if on the boarder:
    // also increasing rge accuracy in this case:
    // ----------------------------------------------------------------------------------------------------------------------------

    if (irge >= 10)
    {
        rmhu2 = 0.7 * rmhu2old + 0.3 * rmhu2;
        ichoice[4] = 2; //  (i.e. also increasing rge accuracy in this case)
    }

    for (int kk = 14; kk <= 18; kk++)
    {
        if (y[kk] < 0.0)
        {
            if (irge == irgmax)
                errmess[1] = -1.0;
            if (iknowl == 2)
            {
                printf("\nbad input: one m^2(3rd gen. sf) <0 from rge");
                printf("\nmaybe temporary due to iteration. wait and see");
            }
        }
    }

    for (int kk = 24; kk <= 28; kk++)
    {
        if (y[kk] < 0.0)
        {
            if (irge == irgmax)
                errmess[2] = -1.0;
            if (iknowl == 2)
            {
                printf("bad input: one m^2(1,2 gen. sf) <0 from rge");
                printf("maybe temporary due to iteration. wait and see");
            }
        }
    }

    if (errmess[1] == -1.0 || errmess[2] == -1.0)
        writing();

    rmtaur = sqrt(y[14]);
    rml = sqrt(y[15]);
    rmbr = sqrt(y[16]);
    rmtr = sqrt(y[17]);
    rmq = sqrt(y[18]);
    rmer = sqrt(y[24]);
    rmel = sqrt(y[25]);
    rmdr = sqrt(y[26]);
    rmur = sqrt(y[27]);
    rmuq = sqrt(y[28]);

    // Modif (temporary, until final conv) rescue in case tachyon rge sf

    if (irge < irgmax) // protections against nan
    {
        if (y[14] < 0.0)
            rmtaur = 1.0;
        if (y[15] < 0.0)
            rml = 1.0;
        if (y[16] < 0.0)
            rmbr = 1.0;
        if (y[17] < 0.0)
            rmtr = 1.0;
        if (y[18] < 0.0)
            rmq = 1.0;

        if (y[24] < 0.0)
            rmer = 1.0;
        if (y[25] < 0.0)
            rmel = 1.0;
        if (y[26] < 0.0)
            rmdr = 1.0;
        if (y[27] < 0.0)
            rmur = 1.0;
        if (y[28] < 0.0)
            rmuq = 1.0;
    }
    else
    {
        if (y[14] < 0.0)
            errmess[1] = -1.0;
        if (y[15] < 0.0)
            errmess[1] = -1.0;
        if (y[16] < 0.0)
            errmess[1] = -1.0;
        if (y[17] < 0.0)
            errmess[1] = -1.0;
        if (y[18] < 0.0)
            errmess[1] = -1.0;

        if (y[24] < 0.0)
            errmess[2] = -1.0;
        if (y[25] < 0.0)
            errmess[2] = -1.0;
        if (y[26] < 0.0)
            errmess[2] = -1.0;
        if (y[27] < 0.0)
            errmess[2] = -1.0;
        if (y[28] < 0.0)
            errmess[2] = -1.0;

        if (errmess[1] == -1.0 || errmess[2] == -1.0)
            writing();
    }

    b = y[19];

    if (ichoice[1] == 2 && ichoice[6] == 0) // if mu(EWSB) input
    {
        rmu0 = mu;
        y[23] = log(fabs(mu));
    }

    rmu = sgn(rmu0) * exp(y[23]);

    if (ichoice[6] == 0)
        rmu = mu; // if mu(EWSB) input

    if (ichoice[5] == 1)
    {
        bold = b;
        rmuold = 1.0;
    }
    else             // means no radiative ew required
        printf(" "); // check what should be the command ; //continue;

    rmino1 = sgnm1 * exp(y[20]);
    rmino2 = sgnm2 * exp(y[21]);
    rmino3 = sgnm3 * exp(y[22]);

    //  interface with higgs mass spectrum calculations:

    ihflag = ichoice[10];

    msl = rml;
    mtaur = rmtaur;
    msq = rmq;
    mtr = rmtr;
    mbr = rmbr;

    mel = rmel;
    mer = rmer;
    muq = rmuq;
    mur = rmur;
    mdr = rmdr;

    al = atau;
    au = atop;
    ad = ab;
    mu = rmu;

    m1 = rmino1;
    m2 = rmino2;
    m3 = rmino3;

    // Tan beta at the relevant (EWSB) scale: tbeta
    // (extracted from common/su_tbewsb/vuewsb,vdewsb)

    tbeta = vuewsb / vdewsb;
    beta = atan(tbeta);

    // set EWSB scale if not defaut value (used in rge, v_eff, and susy r.c):

    if (ichoice[8] == 0)
    {
        scale = qewsb;
        ewsb2 = qewsb * qewsb;
    }
    else if (ichoice[8] == 1)
    {
        // default ewsb scale:
        if (msttr1 != 0.0 && irge > 1)
        {
            scale = sqrt(msttr1 * msttr2);
        }
        else
        {
            scale = sqrt(msq * mtr);
        }

        if (scale < mz) // Added protections
        {
            scale = mz + 1.0e-1;

            if (scale < sqrt(msq * mtr))
            {
                scale = sqrt(msq * mtr);
                qewsb = scale;
            }
        }

        ewsb2 = pow(scale, 2);
    }

    // Stage 11
    // -------------------------------------------------------------------------------------------------------------------------

    // Gaugino masses:
    su_gaugino(mu, m1, m2, m3, beta, a, gmc, gmn, xmn);

    if (inonpert == -1 && irge == irgmax)
    {
        errmess[10] = -1.0;
        writing();
    }

    dmc1 = gmc[1];
    dmc2 = gmc[2];
    dmn1 = xmn[1];
    dmn2 = xmn[2];
    dmn3 = xmn[3];
    dmn4 = xmn[4];

    //
    // Set up the conditions for radiative sym. break. and stability:
    //--------------------------------------------------------------------------------------------------------------------------------

    cbeta = 1.0 / sqrt(1.0 + tbeta * tbeta);
    sbeta = tbeta * cbeta;
    c2beta = cbeta * cbeta - sbeta * sbeta;
    wm2 = wm * wm;
    zm2 = zm * zm;

    if (ichoice[6] == 0)
    {
        //
        // Input is ma_pole!, mu(ewsb). consistent m^2_hu, m^2_d from EWSB with iteration.
        // -------------------------------------------------------------------------------------------------------------------------

    l66:

        ifix = ifix + 1;
        inonpert = 0;

        if (ichoice[1] == 0)
        {
            ewsb2 = pow(qewsb, 2);
            ytewsb = rmtop / vu;
        }
        // Gaugino masses

        if (ichoice[7] == 2 && irge == irgmax)
            inorc = 1;

        su_gaugino(mu, m1, m2, m3, beta, a, gmc, gmn, xmn);

        if (inonpert == -1 && irge == irgmax)
        {
            errmess[10] = -1.0;
            writing();
        }

        // Sfermion masses
        su_sfermion(msq, mtr, mbr, msl, mtaur, al, au, ad, mu, gmst, msb, gmsl, gmsu, gmsd, gmse, gmsn);

        if (sterr == -1.0 || sberr == -1.0 || stauerr == -1.0 || stnuerr == -1.0)
        {
            // means there is really the tachyonic sfermion mass problem:
            // can't even calculate higgs spectrum etc, so has to stop really.
            errmess[4] = -1.0;
            writing(); // goto Ll801;
        }

        // Higgs masses
        // The user's ma input value is used in higgs spectrum calc.:

        ma2 = aama * aama;
        mapole2 = ma2; // In that case the input is really ma_pole

        if (ewsb2 < mz * mz)
            ewsb2 = pow(qewsb, 2);

        su_susycp(tgbeta);

        if (inonpert == -1 && irge == irgmax)
        {
            errmess[10] = -1.0;
            writing();
        }

        alfa = a;

        //
        // Call one-loop effective potential corrections to mh^2_1,2:
        // dvdvd2, dvdvu2 are d(v_eff)/d(vd^2) and d(v_eff)/d(vu^2) which add corrections to m^2_phid (rmhd2) and m^2_phiu (rmhu2) resp.

        rmtaur = mtaur;
        rml = msl;
        rmbr = mbr;
        rmtr = mtr;
        rmq = msq;
        atau = al;
        ab = ad;
        atop = au;
        rmst12 = pow(msttr1, 2);
        rmst22 = pow(msttr2, 2);
        rmsb12 = pow(msbtr1, 2);
        rmsb22 = pow(msbtr2, 2);
        rmstau12 = pow(gmsl[1], 2);
        rmstau22 = pow(gmsl[2], 2);
        dmsu1 = gmsu[1];
        dmsu2 = gmsu[2];
        dmsd1 = gmsd[1];
        dmsd2 = gmsd[2];
        dmse1 = gmse[1];
        dmse2 = gmse[2];
        dmsn1 = gmsn[1];
        dmsntau = gmsn[3];
        dmc1 = gmc[1];
        dmc2 = gmc[2];
        dmn1 = xmn[1];
        dmn2 = xmn[2];
        dmn3 = xmn[3];
        dmn4 = xmn[4];
        rmu = mu;

        ewsb2 = pow(scale, 2);

        if (ytewsb == 0.0)
            ytewsb = rmtop / vu;
        if (ytauewsb == 0.0)
            ytauewsb = rmtau / vu;
        if (ybewsb == 0.0)
            ybewsb = rmb / vu;

        su_vloop2(ewsb2, mu, au, ad, al, &dvdvd2, &dvdvu2, &pizz);

        if (inonpert == -1 && irge == irgmax)
        {
            errmess[10] = -1.0;
            writing();
        }

        if (ifix == 1)
        {
            dvdvd2 = 0.0;
            dvdvu2 = 0.0;
            rmhu2old = 0.0;
            rmhd2old = 0.0;
        }

        sb2 = pow(sin(beta), 2);
        cb2 = pow(cos(beta), 2);
        mzdr2 = pow(mz, 2) + pizz;

        madr2 = mapole2 + piaa - tadba - d2ma;

        rmhu2 = (cb2 - sb2) * mzdr2 / 2 + cb2 * madr2 - mu * mu - dvdvu2;
        rmhd2 = (sb2 - cb2) * mzdr2 / 2 + sb2 * madr2 - mu * mu - dvdvd2;

        // (note: -dvdvi2: to get "tree-level" values of m^2_hu, m^2_hd,
        // thus without v_eff loop corrections)

        mhu2 = rmhu2; // may be dmhu2,dmhd2
        mhd2 = rmhd2;
        b = (rmhd2 + rmhu2 + dvdvd2 + dvdvu2 + 2 * pow(rmu, 2)) * sin(2 * beta) / (2 * rmu);

        // - to be compared to previous m^2_hu,hd values:

        errhuold = errhu;
        errhdold = errhd;
        errhu = (rmhu2 - rmhu2old) / rmhu2;
        errhd = (rmhd2 - rmhd2old) / rmhd2;
        errstop = 1.0e-2;

        if (ifix > 1)
            errstop = 1.0e-4;

        if (fabs(errhu) > errstop && fabs(errhd) > errstop)
        {
            rmhu2old = rmhu2;
            rmhd2old = rmhd2;
            goto l66;
        }

        y[12] = rmhu2;
        y[13] = rmhd2;

        if (ichoice[6] == 0)
        {
            // stop long (rge) iterations on spectrum when xx % accuracy reached:
            // (usually needs ~ 3-4 iterations). nb conv. test is made on ma_run(ewsb)
            if (irge == 1)
                madr2old = 0.0;

            if (ichoice[9] <= 1)
                if (fabs(1.0 - madr2old / madr2) < 2.0e-2)
                    irgmax = irge;
                else if (fabs(1.0 - madr2old / madr2) < 2.0e-4)
                    irgmax = irge;
            madr2old = madr2;
        }
    }

    // Now comes ichoice[6] != 0 i.e:
    // input parameters m_hu, m_hd
    // consistent mu, b from ewsb conditions

    else
    {
        //
        // stop long iterations on spectrum when xx % accuracy reached:
        // (usually needs ~ 3-4 iterations)

        if (ichoice[1] != 2)
        {
            if (irge == 1)
                rmhu2old = 0.0;
            if (ichoice[9] <= 1)
            { // ! 1% accuracy
                if (fabs(1.0 - rmhu2old / rmhu2) < 2.0e-2)
                    irgmax = irge;
            }
            else
            { // ! 0.01% accuracy
                if (fabs(1.0 - rmhu2old / rmhu2) < 2.0e-4)
                    irgmax = irge;
            }
            rmhu2old = rmhu2;
            if (irge == irgsave)
                errmess[5] = -1.0;
        }

        // --- algorithm to find a consistent mu with v_eff corrections:
        double errmu, errb;
        errb = 1.0e5;
        errmu = 1.0e5;

        if (ichoice[5] == 1)
        {
            // i.e. we want ewsb to determine mu and b
            dvdvd2 = 0.0;
            dvdvu2 = 0.0;
            ifix = 0;
        l80:
            ifix = ifix + 1;
            inonpert = 0;
            mu = rmu;

            su_gaugino(mu, m1, m2, m3, beta, a, gmc, gmn, xmn);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            dmc1 = gmc[1];
            dmc2 = gmc[2];
            dmn1 = xmn[1];
            dmn2 = xmn[2];
            dmn3 = xmn[3];
            dmn4 = xmn[4];

            // equation for ma_run:

            if (irge == 1)
                pizz = 0.0;

            ma2 = (rmhu2 + dvdvu2 - rmhd2 - dvdvd2) / cos(2 * beta) - zm * zm - pizz;
            double dmhu2 = rmhu2;
            double dmhd2 = rmhd2;
            double masave;

            if (ma2 >= 0.0)
            {
                aama = sqrt(ma2);
                masave = aama;
                errmess[3] = 0.0;
            }
            else
            {
                // allows for temporary ma^2 < 0 (before ewsb converges)
                // and attempt to retrieve a correct ma via a correct mu etc.
                // gives approximate ma_run(ewsb) values just so that calculation
                // (ewsb iteration) can proceed for a while:
                aama = 1.1;

                if (!(isnan(mapole)) && mapole != 0.0)
                    aama = mapole;

                masave = aama;
                su_gaugino(mu, m1, m2, m3, beta, a, gmc, gmn, xmn);

                if (inonpert == -1 && irge == irgmax)
                {
                    errmess[10] = -1.0;
                    writing();
                }
            }

            //
            //  now calculate sfermion masses and mixing angle:
            //

            su_sfermion(msq, mtr, mbr, msl, mtaur, al, au, ad, mu, gmst, msb, gmsl, gmsu, gmsd, gmse, gmsn);

            if (sterr == -1.0 || sberr == -1.0 || stauerr == -1.0 || stnuerr == -1.0)
            {
                // means there is really the tachyonic sfermion mass problem:
                // can't even calculate higgs spectrum etc, so has to stop really.
                errmess[4] = -1.0;

                if (iknowl == 2)
                    printf(" caution: m^2_sf < 0! . has been changed to 0 ");

                writing();
            }

            if (tachsqrc == -1.0)
            {
                errmess[4] = -1.0;
                writing();
            }
            // otherwise (= no tachyonic sfermions) calculate higgs mass spectrum:

            su_susycp(tbeta);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            // protection against nan higgs that could occurs despite previous protec.
            if (isnan(aaml) || isnan(aamh) || isnan(aamch))
            {
                errmess[9] = -1.0;
                writing();
            }

            if (aaml == 0.0 || aaml > 1.0e10 || aamh > 1.0e10)
            {
                if (irge == irgmax)
                {
                    errmess[9] = -1.0;
                    writing();
                }
            }

            double rmst12 = pow(msttr1, 2);
            double rmst22 = pow(msttr2, 2);
            double rmsb12 = pow(msbtr1, 2);
            double rmsb22 = pow(msbtr2, 2);
            double rmstau12 = pow(gmsl[1], 2);
            double rmstau22 = pow(gmsl[2], 2);
            dmsu1 = gmsu[1];
            dmsu2 = gmsu[2];
            dmsd1 = gmsd[1];
            dmsd2 = gmsd[2];
            dmse1 = gmse[1];
            dmse2 = gmse[2];
            dmsn1 = gmsn[1];
            dmsntau = gmsn[3];
            alfa = a;

            //  call one-loop effective potential corrections to mh^2_1,2:
            //
            //  dvdvd2, dvdvu2 are d(v_eff)/d(vd^2) and d(v_eff)/d(vu^2) which
            //  add corrections to m^2_hd (rmhd2) and m^2_hu (rmhu2) resp.
            // -----------------------------------------------------------------------------

            if (ytewsb == 0.0)
                ytewsb = rmtop / vu;
            if (ytauewsb == 0.0)
                ytauewsb = rmtau / vu;
            if (ybewsb == 0.0)
                ybewsb = rmb / vu;

            su_vloop2(ewsb2, mu, au, ad, al, &dvdvd2, &dvdvu2, &pizz);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            if (isnan(dvdvd2) || isnan(dvdvu2))
            {
                if (irge == irgmax && ifix != 1)
                {
                    errmess[3] = -1.0;
                    writing();
                }
                else
                {
                    // maybe due to uncorrect spectrum at 1rst iter., give it a chance
                    if (isnan(dvdvd2))
                        dvdvd2 = 0.0;
                    if (isnan(dvdvu2))
                        dvdvu2 = 0.0;
                }
            }

            // now the radiative breaking conditions define true mu(mz):
            //
            // tree-level ewsb conditions as (first time!) mu guess:

            double rmu2;
            if (ifix == 1)
                rmu2 = (rmhd2 - (rmhu2)*pow(tbeta, 2)) / (pow(tbeta, 2) - 1.0) - (zm * zm + pizz) / 2.0;
            else
                rmu2 = (rmhd2 + dvdvd2 - (rmhu2 + dvdvu2) * pow(tbeta, 2)) / (pow(tbeta, 2) - 1.0) - (zm * zm + pizz) / 2.0;

            if (rmu2 <= 0.0)
            {
                if (iknowl == 2)
                    printf("warning: mu^2(ewsb) <0 (may be temporary)");

                if (irge == irgmax && ifix >= 5)
                {
                    // consider the mu^2 < 0 problem irremediable:
                    errmess[6] = -1.0;
                    writing();
                }
                else
                {
                    // take approximate mu "rg" =f(ma,m_hu,mhd) to attempt to retrieve
                    // ewsb convergence:
                    if (pow(aama, 2) - rmhu2 - rmhd2 > 0.0)
                        rmu = sgn(rmu0) * sqrt((aama * aama - rmhu2 - rmhd2) / 2);
                    else
                        // take arbitrary small mu to attempt to retrieve ewsb convergence:
                        rmu = sgn(rmu0) * 10.0;
                }
                rmu = rmu / 2.0;
            }

            else
            {
                rmu = sgn(rmu0) * sqrt(rmu2);
                //  !! added: change (after 10 iter) of standard fixed point algorithm:
                //   mu_new = mu_ewsb -> mu_new = (1-c)* mu_old + c*mu_ewsb, c=0.3
                //   to try recovering convergence if on the boarder:
                if (ifix >= 10)
                    rmu = 0.70 * rmuold + 0.30 * rmu;
                mu = rmu;
            }

            //
            //  and true b(ewsb):
            //  tree-level ewsb conditions as first time mu guess:

            if (ifix == 1)
                b = (rmhd2 + rmhu2 + 2 * pow(rmu, 2)) * sbeta * cbeta / rmu;
            else
                b = (rmhu2 + dvdvu2 + rmhd2 + dvdvd2 + 2 * rmu * rmu) * sbeta * cbeta / rmu;

            //
            // - to be compared to evolved mu values:

            double errmuold;
            errmuold = errmu;
            errmu = (rmu - rmuold) / rmuold;

            //  i.e. considers as unconvergent mu from ewsb either:
            //    -inaccurate (> 1e-4) convergence;
            //    - more than 5 tolerated iterations if ma^2 was in fact <0,
            //  so that convergence is around fake ma,mu
            //  since ma was articifially forced = mz temporarily in that case
            // ====================================================================================

            if (!(fabs(errmu) < 5.e-5 && ma2 > 0.0 && rmu2 > 0.0 || ifix == 20))
            {
                if (ma2 <= 0.0 && ifix == 5)
                    goto l81;

                // !!added to get out if really unconvergent ewsb:
                if (fabs(errmu) > fabs(errmuold) && ifix > 15)
                {
                    if (irge == irgmax)
                    {
                        errmess[6] = -1.0;
                        writing();
                    }
                    else
                        goto l81;
                }

                rmuold = rmu;

                goto l80;
            }

            //  ( end of the iterative loop on consistent mu,b )
        }

    l81:

        // (previous endif = end of the choice m_hu,mhd or ma,mu input)

        if (ma2 <= 0.0 && ifix == 5 && irge == irgmax)
        {
            errmess[6] = -1.0;
            errmess[3] = -1.0;

            if (iknowl == 2)
                printf("Consistent ewsb unconvergent below 1d-4");
        }
        if (ifix == 20 && irge == irgmax)
        {
            errmess[6] = -1.0;
            errmess[3] = -1.0;

            if (iknowl == 2)
                printf("Consistent ewsb unconvergent below 1d-4");
        }

        if (ichoice[1] == 1 && ifix == 20)
        {
            errmess[6] = -1.0;
            if (iknowl == 2)
                printf("Consistent ewsb unconvergent below 1e-4");
        }
    }

    //    control ssb v stability (naive rg improved checks of ufb/ccb):

    r1 = rmhd2 + dvdvd2 + mu * mu;
    r2 = rmhu2 + dvdvu2 + mu * mu;
    r3 = b * mu;
    test1 = r1 * r2 - r3 * r3;
    test2 = ma2 + 2 * r3;
    test3 = ma2 - 2 * r3;

    if (ichoice[5] == 1)
    {
        if (test1 >= 0.0 && irge == irgmax)
        {
            errmess[7] = -1.0;
            if (iknowl == 2)
                printf("warning!: ew sym. break may be not realized");
        }

        if (test2 < 0.0 || test3 < 0.0 && irge == irgmax)
        {
            errmess[8] = -1.0;

            if (iknowl == 2)
                printf("warning: potential maybe unbounded from below");
        }

        // ccb (simplest!) constraints, checked at ewsb scale:
        //
        // https://arxiv.org/pdf/hep-ph/0211331
        // Check Eq 40 for CCB condition
        // -------------------------------------------------------------------------

        double ccbt, ccbb, ccbtau, ccbu, ccbd, ccbl;

        if (irge == irgmax)
        {
            ccbt = pow(atop, 2) - 3 * (pow(msq, 2) + pow(mtr, 2) + rmhu2 + pow(rmu, 2));
            ccbb = pow(ab, 2) - 3 * (pow(msq, 2) + pow(mbr, 2) + rmhd2 + pow(rmu, 2));
            ccbtau = pow(atau, 2) - 3 * (pow(msl, 2) + pow(mtaur, 2) + rmhd2 + pow(rmu, 2));
            ccbu = pow(au1, 2) - 3 * (pow(muq, 2) + pow(mur, 2) + rmhu2 + pow(rmu, 2));
            ccbd = pow(ad1, 2) - 3 * (pow(muq, 2) + pow(mdr, 2) + rmhd2 + pow(rmu, 2));
            ccbl = pow(al1, 2) - 3 * (pow(mel, 2) + pow(mer, 2) + rmhd2 + pow(rmu, 2));

            if (ccbt > 0.0 || ccbb > 0.0 || ccbtau > 0.0)
                // ! these are points which do not pass those naive ccb constraints
                errmess[8] = -1.0;

            if (ccbu > 0.0 || ccbd > 0.0 || ccbl > 0.0)
                errmess[8] = -1.0;
        }
    }

    else
    {
        double rmu2;
        // means no radiative ew required
        // now b = y[19] and mu =exp(y(23)) are determined from ew breaking
        // (however not radiative breaking in this case)
        rmu2 = (rmhd2 + dvdvd2 - (rmhu2 + dvdvu2) * pow(tgbeta, 2)) / (pow(tgbeta, 2) - 1.0) - pow(zm, 2) / 2.0;

        if (rmu2 <= 0.0)
        {
            if (iknowl == 2)
            {
                printf(" warning: initial rmu0(high) inconsistent.");
                printf(" has been changed");
            }

            rmu0 = rmu0 / 2;

            for (int i = 1; i <= 8; i++)
                y[i] = ysave[i];

            x2 = x1;
            h1 = -h1;
            goto l7;
        }

        // check how to skip the following part
        rmu = sgn(rmu0) * sqrt(rmu2);

        // - .. b(mz):
        b = (rmhd2 + dvdvd2 + rmhu2 + dvdvu2 + 2 * rmu2) * sbeta * cbeta / rmu;
        //    control of ssb and v stability scales:
        r1 = rmhd2 + dvdvd2 + rmu2;
        r2 = rmhu2 + dvdvu2 + rmu2;
        r3 = b * rmu;

        test1 = r1 * r2 - r3 * r3;
        test2 = r1 + r2 + 2 * r3;
        test3 = r1 + r2 - 2 * r3;

        if (test1 > 0.0)
        {
            errmess[7] = -1.0;

            if (iknowl == 2)
            {
                printf("warning: m^2(hu),m^2(hd) inconsistent with ewsb");
                printf("%d %d", rmhu2, rmhd2);
                printf("have been changed ");
            }

            mhu2 = 1.5 * mhu2;
            mhd2 = mhd2;

            for (int i = 1; i <= 8; i++)
                y[i] = ysave[i];

            x2 = x1;
            h1 = -h1;
            goto l7;
        }

        if (test2 < 0.0 || test3 < 0.0)
        {
            errmess[8] = -1.0;

            if (iknowl == 2)
            {
                printf(" warning: potential unbounded from below! ");
                printf(" m^2(hu),m^2(hd) values been changed ");
            }

            mhu2 = 1.5 * mhu2;
            mhd2 = mhd2;

            for (int i = 1; i <= 8; i++)
                y[i] = ysave[i];

            x2 = x1;
            h1 = -h1;
            goto l7;
        }
    }

    if (ichoice[5] != 1)
    {
        aama = sqrt(rmhu2 + rmhd2 + 2 * rmu * rmu);
        //  calculate higgs mass spectrum
        if (ewsb2 < mz * mz)
            ewsb2 = qewsb * qewsb;

        if (ytewsb == 0.0)
            ytewsb = rmtop / vu;

        su_susycp(tgbeta);

        if (inonpert == -1 && irge == irgmax)
        {
            errmess[10] = -1.0;
            writing();
        }

        //  calculate sfermion masses and mixing angle:
        su_sfermion(msq, mtr, mbr, msl, mtaur, al, au, ad, mu, gmst, msb, gmsl, gmsu, gmsd, gmse, gmsn);
    }

    //
    // special case of unconstrained mssm with low-en input:
    // =====================================================

l880:
    if (ichoice[1] == 0)
    {
        //  case of the pmssm (unconstrained mssm, low-en input)

        if (ichoice[6] == 0)
        {
            //   input parameter of the pmssm  is ma_pole , mu(ewsb)
            // stop long iterations on spectrum when xx % accuracy reached:
            // (usually needs ~ 3-4 iterations)

            if (irge == 1)
                mhu2old = 0.0;
            if (irge < irgmax)
            {
                if (ichoice[9] <= 1)
                {
                    if (fabs(1.0 - mhu2old / mhu2) < 2.0e-2)
                        irgmax = irge; //! 1% accuracy
                }
                else
                {
                    if (fabs(1.0 - mhu2old / mhu2) < 2.0e-4)
                        irgmax = irge; //! .01%
                }
            }

            mhu2old = mhu2;
            if (irge == irgsave)
                errmess[5] = -1.0;
            //   gaugino masses
        l881:
            beta = atan(tbeta);
            if (ichoice[7] == 2 && irge == irgmax)
                inorc = 1;

            su_gaugino(mu, m1, m2, m3, beta, a, gmc, gmn, xmn);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            //
            //  sfermion masses

            su_sfermion(msq, mtr, mbr, msl, mtaur, al, au, ad, mu, gmst, msb, gmsl, gmsu, gmsd, gmse, gmsn);

            if (sterr == -1.0 || sberr == -1.0 || stauerr == -1.0 || stnuerr == -1.0)
            {
                // means there is really the tachyonic sfermion mass problem:
                // can't even calculate higgs spectrum etc, so has to stop really.
                errmess[4] = -1.0;
                writing();
            }

            //
            //   higgs masses

            ma2 = aama * aama;
            mapole2 = ma2; // in that case the input is really ma_pole

            if (ewsb2 < mz * mz)
                ewsb2 = qewsb * qewsb;

            su_susycp(tgbeta);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            alfa = a;

            // Check of ewsb in this parametrization:
            // Note we include in the ewsb consistency relations all the
            // v_eff contributions +loop: indeed, it is consistent with the fact
            // that all higgs masses are calculated with 1- +2-loop contributions:
            // -------------------------------------------------------------------------

            rmtaur = mtaur;
            rml = msl;
            rmbr = mbr;
            rmtr = mtr;
            rmq = msq;
            atau = al;
            ab = ad;
            atop = au;

            double rmst12 = pow(msttr1, 2);
            double rmst22 = pow(msttr2, 2);
            double rmsb12 = pow(msbtr1, 2);
            double rmsb22 = pow(msbtr2, 2);
            double rmstau12 = pow(gmsl[1], 2);
            double rmstau22 = pow(gmsl[2], 2);

            rmt2 = pow(mtrun, 2);
            rmtop = mtrun;
            rmb2 = pow(mbrun, 2);
            rmtau2 = pow(mtaurun, 2);
            dmsu1 = gmsu[1];
            dmsu2 = gmsu[2];
            dmsd1 = gmsd[1];
            dmsd2 = gmsd[2];
            dmse1 = gmse[1];
            dmse2 = gmse[2];
            dmsn1 = gmsn[1];
            dmsntau = gmsn[3];
            alfa = a;
            dmc1 = gmc[1];
            dmc2 = gmc[2];
            dmn1 = xmn[1];
            dmn2 = xmn[2];
            dmn3 = xmn[3];
            dmn4 = xmn[4];
            rmu = mu;
            ewsb2 = pow(scale, 2);

            su_vloop2(ewsb2, mu, au, ad, al, &dvdvd2, &dvdvu2, &pizz);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            double sb2 = pow(sin(beta), 2);
            double cb2 = pow(cos(beta), 2);
            double mzdr2 = mz * mz + pizz;
            double madr2 = mapole2 + piaa - tadba - d2ma;
            rmhu2 = (cb2 - sb2) * mzdr2 / 2 + cb2 * madr2 - pow(mu, 2) - dvdvu2;
            rmhd2 = (sb2 - cb2) * mzdr2 / 2 + sb2 * madr2 - pow(mu, 2) - dvdvd2;

            // (note: -dvdvi2: to get "tree-level" values of m^2_hu, m^2_hd,
            // thus without v_eff loop corrections)

            double dmhu2 = rmhu2;
            double dmhd2 = rmhd2;
            b = (rmhd2 + rmhu2 + dvdvd2 + dvdvu2 + 2 * pow(rmu, 2)) * sin(2 * beta) / (2 * rmu);

            // control of ssb and v stability scales:
            r1 = rmhd2 + dvdvd2 + rmu * rmu;
            r2 = rmhu2 + dvdvu2 + rmu * rmu;
            r3 = b * rmu;
            test1 = r1 * r2 - r3 * r3;
            test2 = r1 + r2 + 2 * r3;
            test3 = r1 + r2 - 2 * r3;
            mhu2 = rmhu2;
            mhd2 = rmhd2;

            if (ichoice[1] == 0 || ichoice[1] == 2)
            {
                rmino1 = m1;
                rmino2 = m2;
                rmino3 = m3;
                sgnm1 = m1 / fabs(m1);
                sgnm2 = m2 / fabs(m2);
                sgnm3 = m3 / fabs(m3);
            }
            else
            {
                rmino1 = sgnm1 * exp(y[20]);
                rmino2 = sgnm2 * exp(y[21]);
                rmino3 = sgnm3 * exp(y[22]);
            }

            if (test1 >= 0.0)
                errmess[7] = -1.0;

            if (test2 < 0.0 || test3 < 0.0)
                errmess[8] = -1.0;

            if (ichoice[1] == 2 && irge == irgmax)
                writing();
        }

        //
        // Endif of ichoice[1] pmssm
        // --------------------------------------------------
        else
        {
            //==========================
            //   input parameter of pmssm  is mhd2,mhu2:
            ichoice[5] = 1;
            rmhd2 = mhd2;
            rmhu2 = mhu2;

            if (irge == 1)
            {
                pizz = 0.0;
                dvdvd2 = 0.0;
                dvdvu2 = 0.0;
            }

            ewsb2 = pow(qewsb, 2);

            if (irge >= 2)
            {
                su_vloop2(ewsb2, mu, au, ad, al, &dvdvd2, &dvdvu2, &pizz);
                if (inonpert == -1 && irge == irgmax)
                {
                    errmess[10] = -1.0;
                    writing();
                }
            }
            double rmu2;
            ma2 = (rmhu2 + dvdvu2 - rmhd2 - dvdvd2) / cos(2 * beta) - zm * zm - pizz;

            //
            // --- algorithm to find a consistent mu with v_eff corrections:
            // --- the radiative breaking conditions define true mu(mz):
            //
            tgbeta = tbeta;
            rmu2 = (rmhd2 + dvdvd2 - (rmhu2 + dvdvu2) * pow(tbeta, 2)) / (pow(tbeta, 2) - 1.0) - (pow(zm, 2) + pizz) / 2.0;

            if (rmu2 <= 0.0)
            {
                if (iknowl == 2)
                {
                    printf(" caution: initial m^2_hu,hd inconsistent");
                    printf(" their values were changed so that mu^2 >=0! ");
                }
                //  find the minimal values of m^2_hu,hd to guarantee mu^2 >0,ma>0:
                rmhu2 = (1.e-6 + mz * mz / 2) * (1.0 - pow(tgbeta, 2)) / (1.0 + pow(tgbeta, 2)) + (aama * aama - 2 * 1.0e-6) / (1.0 + pow(tgbeta, 2));
                rmhd2 = -rmhu2;
                rmu = sgnmu0 * 1.0e-6;
                rmu2 = rmu * rmu;
            }
            else
            {
                // rmu^2 >0 from the input
                rmu = sgnmu0 * sqrt(rmu2);
                rmu2 = rmu * rmu;
            }

            double rmu2old;
            // stop long iterations on spectrum when xx % accuracy reached:
            // (usually needs ~ 3-4 iterations)
            if (irge == 1)
                rmu2old = 0.0;

            if (ichoice[9] <= 1)
            {
                if (fabs(1.0 - rmu2old / rmu2) < 0.02)
                    irgmax = irge;
            }
            else
            {
                if (fabs(1.0 - rmu2old / rmu2) < 0.0002)
                    irgmax = irge; //!.002->.0002
            }

            rmu2old = rmu2;

            if (irge == irgsave)
                errmess[5] = -1.0;
            //
            // - .. and true b(mz):

            b = (rmhd2 + dvdvd2 + rmhu2 + dvdvu2 + 2 * rmu2) * sbeta * cbeta / rmu;
            mu = rmu;

            su_gaugino(mu, m1, m2, m3, beta, a, gmc, gmn, xmn);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            //  calculate sfermion masses and mixing angle:
            su_sfermion(msq, mtr, mbr, msl, mtaur, al, au, ad, mu, gmst, msb, gmsl, gmsu, gmsd, gmse, gmsn);

            if (sterr == -1.0 || sberr == -1.0 || stauerr == -1.0 || stnuerr == -1.0)
            {
                // means there is really the tachyonic sfermion mass problem:
                // can't even calculate higgs spectrum etc, so has to stop really.
                errmess[4] = -1.0;
                writing();
            }

            //
            if (ma2 >= 0.0)
                aama = sqrt(ma2);
            else
            {
                aama = 1.0e-6;
                errmess[3] = -1.0;
            }

            if (ewsb2 < mz * mz)
                ewsb2 = pow(qewsb, 2);

            su_susycp(tgbeta);

            if (inonpert == -1 && irge == irgmax)
            {
                errmess[10] = -1.0;
                writing();
            }

            alfa = a;

            //    control of ssb and v stability scales:

            r1 = rmhd2 + dvdvd2 + rmu2;
            r2 = rmhu2 + dvdvu2 + rmu2;
            r3 = b * rmu;
            test1 = r1 * r2 - r3 * r3;
            test2 = r1 + r2 + 2 * r3;
            test3 = r1 + r2 - 2 * r3;

            if (test1 >= 0.0)
                errmess[7] = -1.0;

            if (test2 < 0.0 || test3 < 0.0)
                errmess[8] = -1.0;

            if (ichoice[1] == 0 || ichoice[1] == 2)
            {
                rmino1 = m1;
                rmino2 = m2;
                rmino3 = m3;
                sgnm1 = m1 / fabs(m1);
                sgnm2 = m2 / fabs(m2);
                sgnm3 = m3 / fabs(m3);
            }
            else
            {
                rmino1 = sgnm1 * exp(y[20]);
                rmino2 = sgnm2 * exp(y[21]);
                rmino3 = sgnm3 * exp(y[22]);
            }
        }
    }

    if (irge == irgmax && iremember == 1)
        ichoice[1] = 2;

    // (a trick to simplify the case of bottom-up evol with mhu,mhd input)
    if (ichoice[1] == 2 && irge == irgmax)
    {
        //
        // unconstrained mssm runned up to high scale
        //  now the final run from q_ewsb to q_final:
        x1 = log(qewsb);
        x2 = log(ehigh);
        h1 = h1 * (x2 - x1) / fabs(x2 - x1);
        //
        if (ichoice[2] == 11)
            su_odeint(y, n, x1, x2, eps, h1, 1.e-8, 1);
        else if (ichoice[2] == 21)
            su_odeint(y, n, x1, x2, eps, h1, 1.e-8, 2);

        goto l882;
    }

    //     ****************************************************************
    //      susy radiative corrections to tau,b,t and sparticle masses:
    //     ****************************************************************
    // 		recovering all rge parameter values at mz scale:

l884:
    if (ichoice[1] != 0 && irge >= 2)
    {
        y[19] = b;
        y[23] = log(fabs(mu));
        xewsb = log(qewsb); //!! added to be consistent with new
                            //!! protections for tachyonic sfermions

        if (ichoice[2] == 11)
            su_odeint(y, n, xewsb, x2, eps, h1, 1.e-8, 1);
        else if (ichoice[2] == 21)
            su_odeint(y, n, xewsb, x2, eps, h1, 1.e-8, 2);

        vu = vu_mz;
        vd = vd_mz;

        rmtau = y[4] * vd;
        rmb = y[5] * vd;
        rmtop = y[6] * vu;
        mtaurun = rmtau;
        mbrun = rmb;
        mtrun = rmtop;
    }
    else if (ichoice[1] == 0)
    {
        y[9] = al;
        y[10] = ad;
        y[11] = au;
        y[29] = al1;
        y[30] = ad1;
        y[31] = au1;
        y[12] = mhu2;
        y[13] = mhd2;
        y[14] = mtaur * mtaur;
        y[15] = msl * msl;
        y[16] = mbr * mbr;
        y[17] = mtr * mtr;
        y[18] = msq * msq;
        y[19] = b;
        y[23] = log(fabs(mu));
        y[24] = mer * mer;
        y[25] = mel * mel;
        y[26] = mdr * mdr;
        y[27] = mur * mur;
        y[28] = muq * muq;
        y[20] = log(fabs(m1));
        y[21] = log(fabs(m2));
        y[22] = log(fabs(m3));
        x1 = log(qewsb);

        x2 = log(mz);
        h1 = -h1;

        if (ichoice[2] == 11)
            su_odeint(y, n, x1, x2, eps, h1, 1.e-8, 1);
        else if (ichoice[2] == 21)
            su_odeint(y, n, x1, x2, eps, h1, 1.e-8, 2);

        vu = vu_mz;
        vd = vd_mz;

        rmtau = y[4] * vd;
        rmb = y[5] * vd;
        rmtop = y[6] * vu;
        mtaurun = rmtau;
        mbrun = rmb;
        mtrun = rmtop;
    }

    //
    // Incorporating leading susy rc to gluino mass:
    // -----------------------------------------------------------
    if (ichoice[7] == 2)
    {
        su_ginocr(alsewsb, m3, mb0, mt0, &delgino);
        mgluino = fabs(m3) / (1.0 - delgino / fabs(m3));
    }
    else
        mgluino = fabs(m3);

    mglu = mgluino;

    if (ichoice[7] >= 1)
    {
        //======  incorporating mb,mt,mtau corrections:
        // first redefining all needed soft etc parameters now at mz scale:
        alz = y[9];
        adz = y[10];
        auz = y[11];
        mtaurz = sqrt(y[14]);
        mslz = sqrt(y[15]);
        mbrz = sqrt(y[16]);
        mtrz = sqrt(y[17]);
        msqz = sqrt(y[18]);
        merz = sqrt(y[24]);
        melz = sqrt(y[25]);
        mdrz = sqrt(y[26]);
        murz = sqrt(y[27]);
        muqz = sqrt(y[28]);
    }

    // Modif (temporary, until final conv) rescue in case tachyon rge sf
    if (irge < irgmax) //! protections against nan
    {
        if (y[14] < 0.0)
            mtaurz = 1.0;
        if (y[15] < 0.0)
            mslz = 1.0;
        if (y[16] < 0.0)
            mbrz = 1.0;
        if (y[17] < 0.0)
            mtrz = 1.0;
        if (y[18] < 0.0)
            msqz = 1.0;

        if (y[24] < 0.0)
            merz = 1.0;
        if (y[25] < 0.0)
            melz = 1.0;
        if (y[26] < 0.0)
            mdrz = 1.0;
        if (y[27] < 0.0)
            murz = 1.0;
        if (y[28] < 0.0)
            muqz = 1.0;
    }
    else
    {
        if (y[14] < 0.0)
            errmess[1] = -1.0;
        if (y[15] < 0.0)
            errmess[1] = -1.0;
        if (y[16] < 0.0)
            errmess[1] = -1.0;
        if (y[17] < 0.0)
            errmess[1] = -1.0;
        if (y[18] < 0.0)
            errmess[1] = -1.0;

        if (y[24] < 0.0)
            errmess[2] = -1.0;
        if (y[25] < 0.0)
            errmess[2] = -1.0;
        if (y[26] < 0.0)
            errmess[2] = -1.0;
        if (y[27] < 0.0)
            errmess[2] = -1.0;
        if (y[28] < 0.0)
            errmess[2] = -1.0;

        if (errmess[1] == -1.0 || errmess[2] == -1.0)
            writing();
    }

    mu_mz = sgnmu0 * exp(y[23]);
    b_mz = y[19];
    m1z = sgnm1 * exp(y[20]);
    m2z = sgnm2 * exp(y[21]);
    m3z = sgnm3 * exp(y[22]);

    if (irge == 1)
    {
        mtausave = rmtau;
        mbsave = rmb;
        mtsave = rmtop;
    }

    // calculating all sfermion parameters at mz scale:
    su_sfbpmz(pizz_mz, msqz, mtrz, mbrz, mslz, mtaurz, muqz, murz, mdrz, melz, merz, alz, auz, adz, mu_mz, b_mz, tgbet0, rmtau, rmb, rmtop);

    if (sterr == -1.0 || sberr == -1.0 || stauerr == -1.0 || stnuerr == -1.0)
    {
        // means there is really the tachyonic sfermion mass problem at q=mz
        errmess[4] = -1.0;
        if (errma2z == -1.0)
            // put error flag: ma^2(mz)<0 at last iter, considered irremediable
            errmess[3] = errma2z;
        writing();
    }

    if (errma2z == -1.0)
    {
        // stop/ put error flag: ma^2(mz)<0 at last iter, considered irremediable
        errmess[3] = errma2z;
        writing();
    }

    su_bmsusycr(alphas, mb, rmtop, rmb, y[6], tgbet0, m2z, m3z, msqz, mtrz, mbrz, auz, adz, mu_mz, &delmb);
    // now susy rc to tau and top  masses:

    msntau_mz = sqrt(pow(mslz, 2) + 0.50 * (pow(mz, 2) + pizz_mz) * cos(2 * beta_z));

    if (isnan(msntau_mz))
        msntau_mz = 1.0; //! protection

    delmtau = su_taumscr(tgbet0, mu_mz, m2z, msntau_mz); // ! changed

    su_topmscr(alphas, mt, mb0, rmtop, rmb, y[6], y[5], tgbet0, m3z, msqz, mtrz, mbrz, auz, adz, mu_mz, &delmtop);

    //  nb: susy rc to quark masses redefines their respective yukawas
    // (we assume the top, b, tau pole masses do not change, within exp.acc.)
    if (irge < irgmax)
    {
        //  redefining running mtau,mb,mtop masses and yuk. cplgs at z scale:
        //  modif in mb resummations (since 2.11 version):
        //  for t,b we have generically: m(pole) = m(run,q) * (1 +cr_qcd(q)+cr_susy(q) )
        //  from which we want to extract e.g. mb(run,mz).
        // 1) no resummation for mtop: (mt = mt_pole,delmtop = cr_qcd(mt)+cr_susy(mt)
        // i.e. delmtop contains all corrections):

        rmtop = mtpole * (1.0 + delmtop);

        // similarly for mtau:
        rmtau = mtau * (1.0 + delmtau);

        // 2) now for mb: note that in eqs. below: rmb is mb(run,mz)(qcd+susy);
        //  delmb =  cr_susy(mz)only, as cr_qcd(mz) is already taken into account before
        //  also resummation is made for mb which may be relevant for large tb
        rmb = mbsave / (1.0 + delmb);

        y[4] = rmtau / vd;
        y[5] = rmb / vd;
        y[6] = rmtop / vu;
    }

    //
    //  Now this will redefine yukawas at high scale as well:
    // ---------------------------------------------------------------

    mtaurun = rmtau;
    mbrun = rmb;
    mtrun = rmtop;

    if (irge < irgmax)
    {
        // saving some parameters:
        dmhu2 = rmhu2;
        dmhd2 = rmhd2;
        dm1 = rmino1;
        dm2 = rmino2;
        dm3 = rmino3;
        dtgbeta = tgbeta;
        dma = aama;
        dml = ml;
        dmh = aamh;
        dmch = aamch;

        dmc1 = gmc[1];
        dmc2 = gmc[2];
        dmn1 = xmn[1];
        dmn2 = xmn[2];
        dmn3 = xmn[3];
        dmn4 = xmn[4];

        dmst1 = gmst[1];
        dmst2 = gmst[2];
        dmsu1 = gmsu[1];
        dmsu2 = gmsu[2];
        dmsb1 = msb[1];
        dmsb2 = msb[2];
        dmsd1 = gmsd[1];
        dmsd2 = gmsd[2];
        dmsl1 = gmsl[1];
        dmsl2 = gmsl[2];
        dmse1 = gmse[1];
        dmse2 = gmse[2];
        dmsn1 = gmsn[1];
        dmsntau = gmsn[3];

        dmsl = msl;
        dmtaur = mtaur;
        dmsq = msq;
        dmtr = mtr;
        dmbr = mbr;
        dmel = mel;
        dmer = mer;
        dmuq = muq;
        dmur = mur;
        dmdr = mdr;
        dal = al;
        dau = au;
        dad = ad;
        dal1 = al1;
        dau1 = au1;
        dad1 = ad1;
        dma = aama;
        dmu = mu;

        goto l44;
    }
    else
    {
        double mtcr, mbcr, mtaucr;
        //==========
        // means that no rc are required
        mtcr = mt;
        mbcr = mb;
        mtaucr = mtau;
    }

    //
    // last thing: calculating now the r.c to chargino, neutralino masses:
    // ---------------------------------------------------------------------
    if (ichoice[7] == 2)
        inorc = 1;

    su_gaugino(mu, m1, m2, m3, beta, a, gmc, gmn, xmn);

    dmc1 = gmc[1];
    dmc2 = gmc[2];
    dmn1 = xmn[1];
    dmn2 = xmn[2];
    dmn3 = xmn[3];
    dmn4 = xmn[4];

//
// Now comes the writing in the outputs part.
//------------------------------------------------------------------
l801:
    //
    // additional theoretical and experimental limits checks (g-2 etc)
    //------------------------------------------------------------------
    errnogo = errmess[4] + errmess[9] + errmess[10];

    if (errnogo == 0.0)
    {
        //
        // 1) the rho parameter (su[2]_custodial breaking at loop-level):
        //---------------------------------------------------------------
        crho = 0.0;
        su_delrho(mt, gmst, msb, gmsl, gmsn[3], thetout, thebout, thelout, &crho);

        //  2) g_mu -2 sm + susy contributions:
        double vv;
        su_gminus2(mel, mer, al1, mu, tgbeta, gmc[1], gmc[2], dxmn); //,u[3],vv[3],z[5], gmuon);
        //  3) what follow is for interface with b-> s gamma calculation:

        int imod_bs, io_bs;
        double bsdeltp, bsvkm, bsl, bsthet, bstheb, xsuh;
        imod_bs = 2;
        io_bs = 1;
        bsdeltp = 0.9;
        bsvkm = 0.95;
        bsl = 0.105;
        // (re)define st,sb mixing to match b->s gamma routine conventions:
        // = flip angles def so that m_sf_1 > m_sf_2 (bsg conventions)

        bsthet = (thetout - pi / 2) / pi;
        bstheb = (thebout - pi / 2) / pi;

        double xsvl, xsul;
        xsuh = min1(gmst[2], mgluino, gmsu[1], gmsd[1]);
        xsul = max(gmst[1], gmc[1]);
        xsvl = min1(gmst[1], msb[1], gmsu[1], gmsd[1]);
        xsvl = (xsvl < mgluino) ? xsvl : mgluino;

        int inlosusy, ihv;
        if (xsvl >= 400.0)
        {
            inlosusy = 1;
            ihv = 1;
        }
        else if (fabs(bsthet) < 0.10 && xsuh > 2 * xsul)
        {
            inlosusy = 1;
            ihv = 0;
        }
        else
        {
            inlosusy = 0;
            ihv = 0;
        }
        double bsgchm[3], bsgflag, mmm2, ubsg[3][3], vbsg[3][3], ierr, c70, c71, c80, c81, ee, rbox;

        bsgchm[1] = gmc[2];
        bsgchm[2] = gmc[1];
        bsgflag = 0.0;

        chargino(tgbeta, gmc[1], mu, mmm2, bsgchm, ierr);

        matching(imod_bs, io_bs, inlosusy, ihv, mw, alphas0, mt, aamch, tgbeta, gmst[1], gmst[2], bsthet, msb[1], msb[2], bstheb, gmsd[1], mgluino, au, ad, rmu, bsgchm, &c70, &c80, &c71, &c81, &ee, &rbox, &ierr);
        su_bsg(alphas0, mt, mbpole - mc0, mc0 / mbpole, alfinv, mw, rmb, rmb, bsvkm, bsl, bsdeltp, io_bs, c70, c71, c80, c81, ee, rbox, &brsg);

        //
        // 4) calculating some fine-tuning parameters for info
        su_finetune(mu, tgbeta, rmhd2, rmhu2, &czmu, &czbmu, &ctmu, &ctbmu);
    }

    //
    // saving final soft etc parameters and output masses:
    // special case:
    //------------------------------------------------------------------------------------

    if (ichoice[1] == 2)
    {
        rmhu2 = y[12];
        rmhd2 = y[13];
    }

    if (ichoice[1] != 2)
    {
        dmhu2 = rmhu2;
        dmhd2 = rmhd2;
        dm1 = rmino1;
        dm2 = rmino2;
        dm3 = rmino3;
        dtgbeta = tgbeta;
        dmsl = msl;
        dmtaur = mtaur;
        dmsq = msq;
        dmtr = mtr;
        dmbr = mbr;
        dmel = mel;
        dmer = mer;
        dmuq = muq;
        dmur = mur;
        dmdr = mdr;
        dal = al;
        dau = au;
        dad = ad;
        dal1 = al1;
        dau1 = au1;
        dad1 = ad1;
        dmu = mu;
    }

    dma = ma;
    dml = ml;
    dmh = mh;
    dmch = mch;
    dmc1 = gmc[1];
    dmc2 = gmc[2];
    dmn1 = xmn[1];
    dmn2 = xmn[2];
    dmn3 = xmn[3];
    dmn4 = xmn[4];
    dmst1 = gmst[1];
    dmst2 = gmst[2];
    dmsu1 = gmsu[1];
    dmsu2 = gmsu[2];
    dmsb1 = msb[1];
    dmsb2 = msb[2];
    dmsd1 = gmsd[1];
    dmsd2 = gmsd[2];
    dmsl1 = gmsl[1];
    dmsl2 = gmsl[2];
    dmse1 = gmse[1];
    dmse2 = gmse[2];
    dmsn1 = gmsn[1];
    dmsntau = gmsn[3];

    // SUSPECT OUTPUT WRITING (in SUSPECT2.out)

    if (errmess[1] == -1.0 || errmess[2] == -1.0 || errmess[4] == -1.0 || errmess[6] == -1.0 || errmess[9] == -1.0 || errmess[10] == -1.0)
        printf("\nCAUTION UNRELIABLE OUTPUT! check errmess below");

    FILE *fpout;
    fpout = fopen("suspect2.out", "w");

    if (ichoice[1] == 10)
    {
        fprintf(fpout, "\n        SUSPECT2.4 OUTPUT: MSUGRA CASE");
        fprintf(fpout, "\n        ------------------------------\n");
    }
    else if (ichoice[1] == 11)
    {
        fprintf(fpout, "\n		 SUSPECT2.4 OUTPUT: GMSB CASE");
        fprintf(fpout, "\n		 ----------------------------\n");
    }
    else if (ichoice[1] == 12)
    {
        fprintf(fpout, "\n		 SUSPECT2.4 OUTPUT: AMSB CASE");
        fprintf(fpout, "\n		 ----------------------------");
    }
    else
    {
        fprintf(fpout, "\n		 SUSPECT2.4 OUTPUT: pMSSM CASE");
        fprintf(fpout, "\n		 -----------------------------\n");
    }

    if (ichoice[1] == 0)
    {
        fprintf(fpout, "\n		 Spectrum calculation only at low (EWSB) energy scale");
        fprintf(fpout, "\n		 ----------------------------------------------------\n");
    }
    if (ichoice[1] == 2)
    {
        fprintf(fpout, "\n		 Bottom-up: RGE from low (EWSB) to GUT energy scale");
        fprintf(fpout, "\n		 --------------------------------------------------\n");
    }

    fprintf(fpout, "\n		 Input values:  ");
    fprintf(fpout, "\n		 ---------------\n");

    if (ichoice[1] == 10)
    {
        fprintf(fpout, "\n		 m_0 \t m_1/2 \t A_0 \t tan(beta) \t sign(mu)");
        fprintf(fpout, "\n      %e \t   %e  \t  %e \t  %e     \t  %e\n", rm0, rmhalf, a0, tgbet0, sgnmu0);
    }
    else if (ichoice[1] == 11)
    {
        fprintf(fpout, "\n		 M_mess \t M_susy \t nl \t nq \t tan(beta) \t sign(mu)");
        fprintf(fpout, "\n       %e   \t   %e   \t %e \t %e \t  %e       \t %e \n", mgmmess, mgmsusy, nl, nq, tgbet0, sgnmu0);
    }

    else if (ichoice[1] == 12)
    {
        fprintf(fpout, "\n		 M_3/2 \t m_0 \t tan(beta) \t sign(mu)");
        fprintf(fpout, "\n     %e   \t  %e  \t  %e  \t  %e ", m32, am0, tgbet0, sgnmu0);

        fprintf(fpout, "\n		 cQ \t cuR \t cdR \t cL \t ceR \t cHu \t cHd");
        fprintf(fpout, "\n     %e \t %e  \t  %e \t %e \t %e  \t %e  \t %e ", cq, cu, cd, cl, ce, chu, chd);
    }

    fprintf(fpout, "\n     M_top \t mb_mb \t M_tau \t 1/alpha \t sw**2(M_Z) \t alpha_S");
    fprintf(fpout, "\n			%e  \t  %e   \t  %e   \t   %e    \t   %e       \t    %e  ", mt, mbmb, mtau, alfinv, sw20, alphas0);

    if (ichoice[1] != 0)
    {
        fprintf(fpout, "\n     M_GUT \t M_EWSB \t E_LOW \t (input or ouput scales)");
        if (ichoice[3] == 0)
            fprintf(fpout, "\n     %e \t %e \t %e \n", ehigh, sqrt(ewsb2), delow);

        else if (ichoice[3] == 1)
            fprintf(fpout, "\n     %e \t %e \t %e \n", egut, sqrt(ewsb2), delow);
    }

    if (ichoice[1] == 1)
    {
        fprintf(fpout, "\n     Input non-universal soft terms at M_GUT");
        fprintf(fpout, "\n     ---------------------------------------");
    }
    if (ichoice[1] == 0 || ichoice[1] == 2)
    {
        fprintf(fpout, "\n     Input non-universal soft terms at M_EWSB");
        fprintf(fpout, "\n     ----------------------------------------");
    }

    if (ichoice[1] == 0 || ichoice[1] == 1)
    {
        if (ichoice[6] == 0)
        {
            fprintf(fpout, "\n     Q_EWSB \t mu \t M_A \t tan(beta) \t sign(mu)");
            fprintf(fpout, "\n		    %e  \t %e \t %e  \t  %e     \t  %e", qewsb, mu0, ma, tbeta, sgnmu0);
        }
        else if (ichoice[6] == 1)
        {
            fprintf(fpout, "\n     Q_EWSB \t M^2_Hu \t M^2_Hd \t tan(beta) \t sign(mu)");
            fprintf(fpout, "\n		    %e  \t %e \t %e  \t  %e     \t  %e", qewsb, mhu20, mhd20, tbeta, sgnmu0);
        }

        fprintf(fpout, "\n     M_1 \t M_2 \t M_3");
        fprintf(fpout, "\n      %e \t %e \t %e ", m10, m20, m30);

        fprintf(fpout, "\n     m_eR \t m_eL \t m_dR \t m_uR \t m_qL");
        fprintf(fpout, "\n       %e \t  %e  \t  %e  \t  %e  \t  %e", mer0, mel0, mdr0, mur0, muq0);

        fprintf(fpout, "\n     m_tauR \t m_tauL \t m_bR \t m_tR \t m_QL");
        fprintf(fpout, "\n        %e  \t  %e    \t   %e \t  %e  \t  %e ", mtaur0, msl0, mbr0, mtr0, msq0);

        fprintf(fpout, "\n     Atau \t Abottom \t Atop \t Al \t Ad \t Au");
        fprintf(fpout, "\n       %e \t   %e    \t  %e  \t %e \t %e \t %e", al0, ad0, au0, al10, ad10, au10);
    }

    if (ichoice[1] == 1 || ichoice[1] == 10)
    {
        fprintf(fpout, "\n     Fermion masses and gauge couplings: Q=HIGH/EWSB");
        fprintf(fpout, "\n     -----------------------------------------------\n");
        fprintf(fpout, "\n     M_top \t M_bot \t M_tau \t g1 \t g2 \t g3");
        fprintf(fpout, "\n       %e  \t  %e   \t  %e   \t 5e \t %e \t %e", mtgut, mbgut, mtaugut, sqrt(ysave[1]), sqrt(ysave[2]), sqrt(ysave[3]));
        fprintf(fpout, "\n       %e %e %e %e %e %e", ytewsb * vuewsb, ybewsb * vdewsb, ytauewsb * vdewsb, sqrt(5. / 3.) * g1ewsb, g2ewsb, sqrt(4 * pi * alsewsb));
    }

    else
    {
        fprintf(fpout, "\n     Fermion masses and gauge couplings: Q=EWSB");
        fprintf(fpout, "\n     ------------------------------------------\n");
        fprintf(fpout, "\n   M_top','M_bot \t M_tau \t g1 \t g2 \t g3");
        fprintf(fpout, "\n   %e %e %e %e %e %e", ytewsb * vuewsb, ybewsb * vdewsb, ytauewsb * vdewsb, sqrt(5. / 3.) * g1ewsb, g2ewsb, sqrt(4 * pi * alsewsb));
    }

    if (ichoice[1] != 0)
    {
        fprintf(fpout, "mu parameter and soft terms at M_EWSB:");
        fprintf(fpout, "--------------------------------------");
        fprintf(fpout, "\n   mu \t B \t M^2_Hu \t M^2_Hd");
        fprintf(fpout, "\n   %e %e %e %e\n", rmu, b, rmhu2, rmhd2);

        fprintf(fpout, "   M_1 \t M_2 \t M_3");
        fprintf(fpout, "\n  %e %e %e ", m1, m2, m3);

        fprintf(fpout, "\n   m_tauR \t m_tauL \t m_bR \t m_tR \t m_QL");
        fprintf(fpout, "\n %e %e %e %e %e", rmtaur, rml, rmbr, rmtr, rmq);

        fprintf(fpout, "\n   m_eR \t m_eL \t m_dR \t m_uR \t m_qL");
        fprintf(fpout, "\n    %e %e %e %e %e", rmer, rmel, rmdr, rmur, rmuq);

        fprintf(fpout, "\n   Atau \t Abottom \t Atop \t Al \t Ad \t Au");
        fprintf(fpout, "\n    %e %e %e %e %e %e", al, ad, au, al1, ad1, au1);
    }

    if (ichoice[1] == 2)
    {
        fprintf(fpout, "mu parameter and soft terms at M_GUT:");
        fprintf(fpout, "--------------------------------------");
        fprintf(fpout, "\n   mu \t B \t M^2_Hu \t M^2_Hd");
        fprintf(fpout, "\n %e %e %e %e\n", mugut, 90.0, mhu2gut, mhd2gut);

        fprintf(fpout, "   M_1 \t M_2 \t M_3");
        fprintf(fpout, "\n %e %e %e\n", m1gut, m2gut, m3gut);

        fprintf(fpout, "\n   m_tauR \t m_tauL \t m_bR \t m_tR \t m_QL");
        fprintf(fpout, "\n %e %e %e %e %e\n", mtaurgut, mslgut, mbrgut, mtrgut, msqgut);

        fprintf(fpout, "\n   m_eR \t m_eL \t m_dR \t m_uR \t m_qL");
        fprintf(fpout, "\n %e %e %e %e %e\n", mergut, melgut, mdrgut, murgut, muqgut);

        fprintf(fpout, "\n   Atau \t Abottom \t Atop \t Al \t Ad \t Au");
        fprintf(fpout, "\n %e %e %e %e %e %e\n", algut, adgut, augut, al1gut, ad1gut, au1gut);
    }

    fprintf(fpout, "\nMass matrices and mixing angles:");
    fprintf(fpout, "\n--------------------------------");
    fprintf(fpout, "\n   tan(beta) \t alpha_(h,H)");
    fprintf(fpout, "\n   %e %e\n", tbeta, alfa);

    fprintf(fpout, "\n   thet_tau \t thet_b \t thet_t");
    fprintf(fpout, "\n   %e %e %e\n", thel, theb, thet);

    fprintf(fpout, "\n\n   Z[i][j]");
    fprintf(fpout, "\n  %e %e %e %e\n", z[1][1], z[1][2], z[1][3], z[1][4]);
    fprintf(fpout, "\n  %e %e %e %e\n", z[2][1], z[2][2], z[2][3], z[2][4]);
    fprintf(fpout, "\n  %e %e %e %e\n", z[3][1], z[3][2], z[3][3], z[3][4]);
    fprintf(fpout, "\n  %e %e %e %e\n", z[4][1], z[4][2], z[4][3], z[4][4]);

    fprintf(fpout, "\n   U(i,j) \t V(i,j)");
    fprintf(fpout, "\n %e %e %e %e", u[1][1], u[1][2], v[1][1], v[1][2]);
    fprintf(fpout, "\n %e %e %e %e", u[2][1], u[2][2], v[2][1], v[2][2]);

    fprintf(fpout, "Final Higgs and SUSY particle masses: ");
    fprintf(fpout, "------------------------------------- ");

    if (ma2 > 0.0)
    {
        fprintf(fpout, "\n   h   \t H \t A \t H+");
        fprintf(fpout, "\n %e %e %e %e", ml, mh, ma, mch);
    }
    else
        fprintf(fpout, "\nMA**2 <0! NO further Higgs masses calculated\n");

    fprintf(fpout, "\n   chi+_1 \t chi+_2 \t chi0_1 \t chi0_2 \t chi0_3 \t chi0_4");
    fprintf(fpout, "\n   %e %e %e %e %e %e", gmc[1], gmc[2], xmn[1], xmn[2], xmn[3], xmn[4]);

    fprintf(fpout, "\n   gluino");
    fprintf(fpout, "\n   %e", mgluino);

    fprintf(fpout, "\n   stop_1 \t stop_2 \t sup_1 \t sup_2");
    fprintf(fpout, "\n   %e %e %e %e", gmst[1], gmst[2], gmsu[1], gmsu[2]);

    fprintf(fpout, "\n   sbot_1 \t sbot_2 \t sdown_1 \t sdown_2");
    fprintf(fpout, "\n   %e %e %e %e", msb[1], msb[2], gmsd[1], gmsd[2]);

    fprintf(fpout, "\n   stau_1 \t stau_2 \t snutau \t selec_1 \t selec_2 \t snuelec");
    fprintf(fpout, "\n   %e %e %e %e %e %e", gmsl[1], gmsl[2], gmsn[3], gmse[1], gmse[2], gmsn[1]);

    fprintf(fpout, "Low-energy/LEP precision parameter values:");
    fprintf(fpout, "\n   Delta_rho \t g_mu -2 \t Br(b->s gamma)");
    fprintf(fpout, "\n %e %e %e", crho, gmuon, brsg);

    fprintf(fpout, "Fine-tuning values for info: fine-tuned if >>1");
    fprintf(fpout, "dmZ^2/mZ^2(mu^2) dmZ^2/mZ^2(B.mu) dmt/mt(mu) \t dmt/mt(B.mu)");
    fprintf(fpout, " %e %e %e %e", czmu, czbmu, ctmu, ctbmu);

l1000:
    if (iknowl != 0)
    {
        fprintf(fpout, "\nWarning/Error Flags: errmess(1)-(10):");
        fprintf(fpout, "\n-------------------------------------");
        for (int i = 1; i <= 10; i++)
            fprintf(fpout, "\nError(%d) : %e", i, errmess[i]);
        fprintf(fpout, "\n---------------------------------");
        fprintf(fpout, "\nerrmess(i)= 0: Everything is fine.");
        fprintf(fpout, "\nerrmess(1)=-1: tachyon 3rd gen. sfermion from RGE");
        fprintf(fpout, "\nerrmess(2)=-1: tachyon 1,2 gen. sfermion from RGE");
        fprintf(fpout, "\nerrmess(3)=-1: tachyon A    (maybe temporary: see final mass) ");
        fprintf(fpout, "\nerrmess(4)=-1: tachyon 3rd gen. sfermion from mixing");
        fprintf(fpout, "\nerrmess(5)=-1: mu unstable after many iter");
        fprintf(fpout, "\nerrmess(6)=-1: non-convergent mu from EWSB ");
        fprintf(fpout, "\nerrmess(7)=-1: EWSB maybe inconsistent (!but RG-improved only check) ");
        fprintf(fpout, "\nerrmess(8)=-1: V_eff maybe UFB or CCB (!but RG-improved only check) ");
        fprintf(fpout, "\nerrmess(9)=-1: Higgs boson masses are NaN ");
        fprintf(fpout, "\nerrmess(10)=-1: RGE problems (non-pert and/or Landau poles)");

        if (errmess[1] == -1.0)
        {
            fprintf(fpout, "\n   Bad input: one m^2(3rd gen. sf) <0 from RGE ");
            fprintf(fpout, "\n   maybe artefact of algorithm, see final result");
        }

        if (errmess[2] == -1.0)
        {
            fprintf(fpout, "\n   Bad input: one m^2(1,2 gen. sf) <0 from RGE ");
            fprintf(fpout, "\n   maybe artefact of algorithm, see final result");
        }

        if (errmess[1] == -1.0 || errmess[2] == -1.0)
        {
            fprintf(fpout, " Tachyonic RGE: UNRELIABLE OUTPUT! ");
        }

        if (errmess[3] == -1.0)
        {
            fprintf(fpout, "\n   Warning:  MA^2(Q) <0 at a scale MZ<Q<EWSB ! ");
            fprintf(fpout, "\n   check final results ");
        }

        if (errmess[4] == -1.0)
        {
            fprintf(fpout, "\n   STOP: one tachyonic m^2(3rd gen. sf) <0 ");
            fprintf(fpout, "\n   UNRELIABLE OUTPUT! ");
        }

        if (errmess[5] == -1.0)
            fprintf(fpout, "\n Warning: MU unstable after many iter");

        if (errmess[6] == -1.0)
            fprintf(fpout, "\nWARNING: EWSB unconvergent after 20 iter.");

        if (errmess[7] == -1.0)
        {
            fprintf(fpout, "\nEW Sym. Break. may be not realized ");
            fprintf(fpout, "\n(however from naive tree-level analysis) ");
        }

        if (errmess[8] == -1.0)
        {
            fprintf(fpout, "\nPotential may be unbounded from below ");
            fprintf(fpout, "\n(however from naive tree-level analysis) ");
        }

        if (errmess[9] == -1.0)
            fprintf(fpout, "\nPROBLEM: some Higgs masses are NaN! ");

        if (errmess[10] == -1.0)
            fprintf(fpout, "\nSTOP: non-pert. R.C., or Landau pole in RGE!");
    }

    fclose(fpout);
    printf("\n\n RUN TERMINATED : OUTPUT in suspect2.out");
}

int main()
{
    int ichoice[12];
    double errmess[31];

    susy(1, 0, ichoice, errmess);
}