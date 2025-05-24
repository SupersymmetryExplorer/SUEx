#ifndef _VARIABLES_H_
#define _VARIABLES_H_

/*************************************************************************/
/*                   Defining Constants									 */
/*************************************************************************/
static double pi = 3.14159265359;

extern int kpole;
extern int irge, irgmax, ifix, isfrc, inorc; 
extern int iflop;		 
extern double tachsqrc; 

extern double mst1bp, mst2bp, msn1bp, msta1bp, msta2bp;
extern double msu1, msu2, msd1, msd2, mse1, mse2, msn1, msntau, msta1, msta2, msb1, msb2, mst1, mst2; 
extern double msb1sfbp, msb2sfbp, mst1sfbp, mst2sfbp;
extern double thett, thetb, thetl; 

extern double m1, mtaur, mu, m2, m3;							  
extern double ytauewsb, ybewsb, ytewsb, alsewsb, g2ewsb, g1ewsb; 
extern double pizzp;											  
extern double sterr, sberr, stauerr, stnuerr;					  
extern double vuewsb, vdewsb;									  
extern double kmaflag, piaa, tadba, d2ma;						  

extern double ma, ml, mh, mch, marun; 
extern double errma2z;    //  ! added for ma^2(mz) <0 control

// the u(1), su(2), su(3) soft susy-breaking gaugino masses:

//   charginos 1,2 masses, neutralinos 1-4 masses, gluino mass
extern double dmc1, dmc2, dmn1, dmn2, dmn3, dmn4, mgluino;

//  stop 1,2 and sup 1,2 = scharm 1,2 masses
extern double dmst1, dmst2, dmsu1, dmsu2;

//  sbottom 1,2 and sdown 1,2 = sstrange 1,2 masses
extern double dmsb1, dmsb2, dmsd1, dmsd2;

//  stau 1,2 ; selectron (=smuon) 1,2; sneut_e,mu, sneut_tau masses
extern double dmsl1, dmsl2, dmse1, dmse2, dmsn1, dmsntau;

//  light, heavy, charged higgs masses, neutral (h,h) mix angle alpha
extern double dml, dmh, dmch, alfa;

// stop, sbottom, stau mixing angles


extern double mbmb, imbmb;

extern double inonpert;
extern double mapole, mchpole;
extern double ipolemz;
extern double scale;
extern double czmu, czbmu, ctmu, ctbmu;

extern double rmtop, susym, egut;
extern int kunif;
extern double gf, alpha, mz, mw;
extern double nf, cpi, zm, wm, tbeta;
extern double xlambda, mc0, mb0, mt0;
extern double mtau, mbpole, mtpole;
extern double mtaurun, mbrun, mtrun;
extern double ytau, yb, yt, ytop;			   
extern double yewsb[32];
extern double msbtr1, msbtr2, msttr1, msttr2; 
extern double istflip, isbflip;			   


extern int n0;
//********************************************************************* //
//  "standard model" input parameters (couplings and fermion masses):   //
//********************************************************************* //

//  dmt,dmtau are pole masses but dmb is mb(mb)_msbar
extern double alfinv, alphas, mt, mb, mc;

//  rg evolution scale parameters ewsb scale, high and low rge ends:
extern double delow;

//   mssm parameters of the scalar sector:
extern double mhu2, mhd2;

//   the soft-susy breaking slepton mass terms (3d and then 1/2 gen.):
extern double msl;

//   the soft-susy breaking squark mass terms (3d and then 1/2 gen.):
extern double msq, mtr, mbr;

//   the soft-susy breaking trilinear couplings (3d and then 1/2 gen.):
extern double al, au, ad;
extern double al1, au1, ad1;
extern double aama, aaml, aamh, aamch;

// low-energy contrained parameter values: rho-1, g_mu-2, br(b->s gamma):
extern double crho, gmuon, brsg;

// gut scale mssm parameters output:
extern double mhu2gut, mhd2gut, magut, mugut;
extern double m1gut, m2gut, m3gut;
extern double mslgut, mtaurgut, melgut, mergut;
extern double msqgut, mtrgut, mbrgut, muqgut, murgut, mdrgut;
extern double algut, augut, adgut;
extern double al1gut, au1gut, ad1gut;
extern double sgnmu0, tgbeta;

//  msugra case input parameters:
extern double m0, a0;
extern int mhalf;

//  gmsb case input parameters:
extern double mgmmess, mgmsusy;
extern int nl, nq;
//  amsb case input parameters:
extern double m32, am0, cq, cu, cd, cl, ce, chu, chd;

extern double mel, mer, muq, mur, mdr;
extern double gmn[5], xmn[5], gmc[3], gmst[3], msb[4], gmsl[3], gmsu[3], gmsd[3], gmse[3], gmsn[5];
extern double bcoup, a, gat, gab, glt, glb, ght, ghb, ghvv, glvv;
extern double beta, adum;
extern double gcen[3][3], gctb[3][3], glee[3][3], gltt[3][3], glbb[3][3], ghee[3][3], ghtt[3][3], ghbb[3][3], gatt, gabb, gaee;
extern double ac1[3][3], ac2[3][3], ac3[3][3], an1[5][5], an2[5][5], an3[5][5], acnl[3][5], acnr[3][5];
extern double vu, vd, atop, ab, atau, rmllt, rmllb, rmlltau, rmrrt, rmrrb, rmrrtau;
extern double rmtaur, rml, rmbr, rmtr, rmq;
extern double g12, g22, sw2;
extern double sgnm1, sgnm2, sgnm3;
extern double u[3][3], v[3][3], z[5][5], dxmn[5];

// ***************************************************************************
//  			            commons internal to the routine
//
//	  ("internal" means that normally the user does not have to care about
//    any parameters defined by the commons etc below: in particular none
// of these commons below should be necessary for interface with other codes)
// ***************************************************************************

extern double wistep, h1;
extern double jfirst, ygut[32];
extern int ifirst;
extern double nnlo, idrflag;
extern int ihflag;

// ***************************************************************************
//  			  Other commons to the routine
// ***************************************************************************

extern double smin_warn, extpar_warn, muma_warn, algo_warn;
extern double b0, b1, b[13];
extern double f[21];
extern double meldum, merdum, muqdum, murdum, mdrdum;
extern double madr2save;
extern double la1, la2, la3, la4, la5, la6, la7;

// The values of Lambda's are stored
//
// On the running coupling constant in QCD
// M. Prosperi, M. Raciti and C. Simolo
// arXiv:hep-ph/0607209v2
// ------------------------------------------------------------------------------------------------------------
//
// For 1-loop and 2-loop calculations respectively
// ------------------------------------------------------------------------------------------------------------
extern double xlb1[7], xlb2[7];

// ***************************************************************************
//    These variables are common to the functions : SU_PIXX(), SU_TOPMSCR()
// ***************************************************************************
extern double hml, hmh, hma, hmch, halfa; 
extern double mhrunz, mlrunz, marunz, mchrunz, alpharunz;				 
extern double thet, theb, thel;
extern double thetbp, thebbp, thelbp;
extern double ama0;
extern double thetout, thebout, thelout;
extern double marunpwsb, mlrunpwsb, mhrunpwsb, mchrunpwsb, alfarunpwsb; 
extern double stoplc, stophc, glc, xstopc, amuc; 
extern double asc, tc, zm1, ioc;				 

#endif