#ifndef _VARIABLES_H_
#define _VARIABLES_H_

#ifdef __cplusplus
extern "C" {
#endif

/*************************************************************************/
/*                   Defining Constants											 */
/*************************************************************************/
//static double pi = 3.14159265359;

int kpole;
int irge, irgmax, ifix, isfrc, inorc; // su_strc

// int irge;						// su_strc
int iflop;		 // su_good/
double tachsqrc; // su_tachyrc

double mst1bp, mst2bp, msn1bp, msta1bp, msta2bp;											   // su_bpew
double msu1, msu2, msd1, msd2, mse1, mse2, msn1, msntau, msta1, msta2, msb1, msb2, mst1, mst2; // su_bpew

double msb1sfbp, msb2sfbp, mst1sfbp, mst2sfbp;

double thett, thetb, thetl; // su_bpew

double m1, mtaur, mu, m2, m3;							  //   common/su_break/   //,msl,msq,mtr,mbr
double ytauewsb, ybewsb, ytewsb, alsewsb, g2ewsb, g1ewsb; // su_yukaewsb
double pizzp;											  // run_p
double sterr, sberr, stauerr, stnuerr;					  // su_errsf/
double vuewsb, vdewsb;									  //  su_tbewsb/
double kmaflag, piaa, tadba, d2ma;						  // common/su_mainput/,kmaflag

double ma, ml, mh, mch, marun; // su_hmass

double errma2z; // su_errma   //  ! added for ma^2(mz) <0 control

//   the u(1), su(2), su(3) soft susy-breaking gaugino masses:
// double dm1,dm2,dm3;   					//su_mssmgpar/dm1,dm2,dm3

//   charginos 1,2 masses, neutralinos 1-4 masses, gluino mass
double dmc1, dmc2, dmn1, dmn2, dmn3, dmn4, mgluino; // su_outginos

//  stop 1,2 and sup 1,2 = scharm 1,2 masses
double dmst1, dmst2, dmsu1, dmsu2; //  su_outsqu

//  sbottom 1,2 and sdown 1,2 = sstrange 1,2 masses
double dmsb1, dmsb2, dmsd1, dmsd2; // su_outsqd

//  stau 1,2 ; selectron (=smuon) 1,2; sneut_e,mu, sneut_tau masses
double dmsl1, dmsl2, dmse1, dmse2, dmsn1, dmsntau; // su_outslep

//  light, heavy, charged higgs masses, neutral (h,h) mix angle alpha
double dml, dmh, dmch, alfa; // su_outhiggs

//  stop, sbottom, stau mixing angles
// double thet,theb,thel;						//su_outmix

double mbmb, imbmb; // su_mbmb //  added for mb(mb) input

double inonpert;				 //   common/su_nonpert///  ! added for non-pert problems
double mapole, mchpole;			 //   common/pietro/
double ipolemz;					 //   common/su_polemz/
double scale;					 //   common/su_renscale/
double czmu, czbmu, ctmu, ctbmu; //   common/su_ftune/

double rmtop, susym, egut;			   // su_sthresh
int kunif;							   // su_gunif
double gf, alpha, mz, mw;			   // su_param
double nf, cpi, zm, wm, tbeta;		   // su_cte
double xlambda, mc0, mb0, mt0;		   // su_als
double mtau, mbpole, mtpole;		   // su_fmasses  //
double mtaurun, mbrun, mtrun;		   // su_runmasses
double ytau, yb, yt, ytop;			   // su_yuka
double yewsb[32];					   // su_allewsb   // 10 set by hand. check the correct on
double msbtr1, msbtr2, msttr1, msttr2; // su_treesfer
double istflip, isbflip;			   // common/su_mixflip/

// double msq,mtr,mbr;

int n0;
// ********************************************************************* //
//   "standard model" input parameters (couplings and fermion masses):   //
//	********************************************************************* //

//  dmt,dmtau are pole masses but dmb is mb(mb)_msbar
double alfinv, alphas, mt, mb, mc; // su_smpar	//,mtau,sw2

//  rg evolution scale parameters ewsb scale, high and low rge ends:
double delow; // su_rgscal        //dqewsb,dehigh,

//   mssm parameters of the scalar sector:
double mhu2, mhd2; // su_mssmhpar  // ,ma,mu;

//   the soft-susy breaking slepton mass terms (3d and then 1/2 gen.):
double msl; // su_mssmslep //,mtaur,mel,mer

//   the soft-susy breaking squark mass terms (3d and then 1/2 gen.):
double msq, mtr, mbr; // su_mssmsqua  //,mur,mdr,muq

//   the soft-susy breaking trilinear couplings (3d and then 1/2 gen.):
double al, au, ad;	  // su_atri3
double al1, au1, ad1; // su_atri12
double aama, aaml, aamh, aamch;

// low-energy contrained parameter values: rho-1, g_mu-2, br(b->s gamma):
double crho, gmuon, brsg; // su_lowen

//   gut scale mssm parameters output:
double mhu2gut, mhd2gut, magut, mugut;				   // su_mssmhgut
double m1gut, m2gut, m3gut;							   // su_mssmggut
double mslgut, mtaurgut, melgut, mergut;			   // su_mssmslgut
double msqgut, mtrgut, mbrgut, muqgut, murgut, mdrgut; // su_mssmsqgut
double algut, augut, adgut;							   // su_a3gut
double al1gut, au1gut, ad1gut;						   // su_a12gut

//  tan(beta) and sign(mu)
double sgnmu0, tgbeta; // su_radewsb

//  msugra case input parameters:
double m0, a0; // su_msugra
int mhalf;

//  gmsb case input parameters:
double mgmmess, mgmsusy; // su_gmsb
int nl, nq;
//  amsb case input parameters:
double m32, am0, cq, cu, cd, cl, ce, chu, chd; // su_amsb

double mel, mer, muq, mur, mdr;																							 // su_break
double gmn[5], xmn[5], gmc[3], gmst[3], msb[4], gmsl[3], gmsu[3], gmsd[3], gmse[3], gmsn[5];							 // su_smass    //check the dimensions correctly
double bcoup, a, gat, gab, glt, glb, ght, ghb, ghvv, glvv;																 // su_hcoup
double beta, adum;																										 //	su_hmix
double gcen[3][3], gctb[3][3], glee[3][3], gltt[3][3], glbb[3][3], ghee[3][3], ghtt[3][3], ghbb[3][3], gatt, gabb, gaee; // su_cplhsf
double ac1[3][3], ac2[3][3], ac3[3][3], an1[5][5], an2[5][5], an3[5][5], acnl[3][5], acnr[3][5];						 // su_cplhino
double vu, vd, atop, ab, atau, rmllt, rmllb, rmlltau, rmrrt, rmrrb, rmrrtau;											 // su_cteloop
double rmtaur, rml, rmbr, rmtr, rmq;																					 // su_soft
double g12, g22, sw2;																									 // su_cpl
double sgnm1, sgnm2, sgnm3;																								 // su_sgnm123
double u[3][3], v[3][3], z[5][5], dxmn[5];																				 // su_matino

// ***************************************************************************
//  			            commons internal to the routine
//
//	  ("internal" means that normally the user does not have to care about
//    any parameters defined by the commons etc below: in particular none
// of these commons below should be necessary for interface with other codes)
// ***************************************************************************

double wistep, h1;		 // su_stepwi
double jfirst, ygut[32]; // su_stegut
int ifirst;				 // su_stegut
double nnlo, idrflag;	 // su_qcdflag
int ihflag;				 // su_hflag

// ***************************************************************************
//  			            other commons to the routine
// ***************************************************************************

double smin_warn, extpar_warn, muma_warn, algo_warn; // su_slha_warn

double b0, b1, b[13];						   // su_bernou
double f[21];								   // su_facto
double meldum, merdum, muqdum, murdum, mdrdum; // su_break2
double madr2save;							   // su_savemar
double la1, la2, la3, la4, la5, la6, la7;	   // common/hself_hdec/

//   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  //

// ode_path

// The values of Lambda's are stored
//
//	On the running coupling constant in QCD
// M. Prosperi, M. Raciti and C. Simolo
// arXiv:hep-ph/0607209v2
// ------------------------------------------------------------------------------------------------------------
//
// For 1-loop and 2-loop calculations respectively
// ------------------------------------------------------------------------------------------------------------
double xlb1[7], xlb2[7]; // su_alslam

// double mudim2; //, divergence, lambda2;					//su_cutoff

// ***************************************************************************
//  	These variables are common to the functions : SU_PIXX(), SU_TOPMSCR()
// ***************************************************************************
double hml = 0.0, hmh = 0.0, hma = 0.0, hmch = 0.0, halfa = 0.0; // SU_higgsrunz
double mhrunz, mlrunz, marunz, mchrunz, alpharunz;				 // SU_higgsrunz

////////////////////////////////////

double thet, theb, thel;
double thetbp, thebbp, thelbp;
double ama0 = 0.0;

double thetout, thebout, thelout; // SU_outmix

double marunpwsb = 0.0, mlrunpwsb = 0.0, mhrunpwsb = 0.0, mchrunpwsb = 0.0, alfarunpwsb = 0.0; // SU_runhiggsewsb

double stoplc, stophc, glc, xstopc, amuc; // st
double asc, tc, zm1, ioc;				  // alp
// double nf,cpi,zm,wm,tbeta;

#ifdef __cplusplus
}
#endif

#endif