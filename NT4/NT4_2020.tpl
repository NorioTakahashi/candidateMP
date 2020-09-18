// Management Procedure Simple "NT4" modified from "NT3" and "NT2" (based on "NT1") code
//
// by norio takahashi
// start  : Apr  04, 2018 (as "NT1")
// end    : Jul  26, 2018 - modified to "NT4", just added tuning year's condition to "NT3"
// modified: Sep 09, 2020 - to compare with Cape Town Procedure based on base19
//
// This MP utilizes the code developed by Rich Hillary (CSIRO) for the empirical POP index based on CKMR data.
// The original source code is "C1GT1CK1a.tpl" ptovided by Rich on 180522.
// See the following references about the empirical POP index:
// Hillary R, Preece A, Davies C (2016) Methods for data generation in projections.
// CCSBT-OMMP/1609/07 (CCSBT-ESC/1609/BGD06, CCSBT-OMMP/1806/BGD01)
//
//
DATA_SECTION
//  !!ad_comm::change_datafile_name("NT4_ctrl.dat"); // set parameter values for MP
  // POP target part parameters
  init_number    ssb_pop_target  // target value of SSB comparable to POP index
  init_int            t_pop      // number of years over which POP indices (SSB) are averaged

  init_int    tune_year  // tuning year

  // average POP <= target POP case
  init_number    k1_cpue  // gain parameter when cpue slope < 0 (decrease)
  init_number    k2_cpue  // gain parameter when cpue slope >= 0 (increase)
  init_int       t1_cpue  // number of years over which slope of CPUEs(age4+) is estimated

  // average POP > target POP case
  init_number    k3_cpue  // gain parameter when cpue slope < 0 (decrease)
  init_number    k4_cpue  // gain parameter when cpue slope >= 0 (increase)
  init_int       t2_cpue  // number of years over which slope of CPUEs(age4+) is estimated

  // GT(age2) part parameters
  // used as limit
  init_number    k1_gt_limit   // smoothing parameter when delta_gt_limit < 1
  init_int       t_gt_limit    // number of years over which GT estimates(age2) are averaged
  init_number    n_age2_limit  // limit value of number of age 2 fish

  // parameters for continuity parts
  // for CPUE-based TAC
  init_number   alpha1
  init_number   alpha2

  // for final TAC
  init_number   beta1
  init_number   beta2

  // switch for which TAC from parts of HCR is used
  init_int    swit_cpue       // switch for TAC from cpue slope part of HCR is used
  init_int    swit_gt_limit   // switch for TAC from GT limit part of HCR is used

  init_number    max_change_up    // maximum TAC increase
  init_number    max_change_down  // maximum TAC decrease
  init_number    min_change       // minimum TAC increase
  init_number    maxTAC           // maximum TAC (for capping TAC)

  init_int    debug_write  // debug write? yes = 1, no = 0


  int    first_cpue_yr   // first year of historical cpue data
  int    last_cpue_yr    // last year of historical cpue data

  !!first_cpue_yr = 1969;
  !!last_cpue_yr = 2019;

  init_vector    hist_cpue(first_cpue_yr,last_cpue_yr)  // historical age4+ cpue

  // we now change files to get simulated data from sbtOMdata file(cpoied from Rich's test code)
  !!ad_comm::change_datafile_name("sbtOMdata");
  init_int first_yr
  init_int last_yr
  init_int implementation_yr // actual quota implementation year.
  init_int current_yr //"current" is implementation year minus nlag (extra lag from decision to implementation)
  init_vector quota(first_yr-1,implementation_yr-1) // historicaly catches
  init_vector proj_cpuen(first_yr,current_yr-2) // projected CPUE (in numbers)
  init_vector proj_cpueb(first_yr,current_yr-2) // projected CPUE (in biomass)
  init_vector proj_asurv(first_yr,current_yr-1) // projected AS (in relative biomass)
  init_vector proj_gtn(first_yr,last_yr)    // GT index
  init_vector proj_cvgtn(first_yr,last_yr)  // GT index CV
  init_vector proj_gtrec(first_yr,last_yr)  // GT index #matches

  // get the number of POP and HSP data (f/ Rich's code)
  !!ad_comm::change_datafile_name("POP_nsim.dat");
  init_int nPOPdat
  !!ad_comm::change_datafile_name("POP_proj.dat");
  init_imatrix POPproj(1,nPOPdat,1,5)

  // get the historical CKMR data (f/ Rich's code)
  // POPs first
  !!ad_comm::change_datafile_name("POP_hist.dat");
  init_int nPOPdat2
  init_imatrix POPhist(1,nPOPdat2,1,5)
  //!!cout << nPOPdat << " " << current_yr << " " << POPproj(nPOPdat) << endl;


  vector obs_cpue(first_cpue_yr,current_yr-2)
  vector obs_gtn(2016,current_yr-2)
  vector obs_cvgtn(2016,current_yr-2)
  vector obs_gtrec(2016,current_yr-2)
  vector yrs(first_cpue_yr,current_yr-2)

  // required vectors/matrices for MP later on (f/ Rich's code)
//  vector w(2016,current_yr-2);

  vector POPind(2002,current_yr-5)
  matrix Mja(2002,last_yr-5,2006,last_yr)
  matrix Rja(2002,last_yr-5,2006,last_yr)
  matrix wja(2002,last_yr-5,2006,last_yr)


PARAMETER_SECTION
  objective_function_value f

 LOC_CALCS

  // primary variables
  double  mu_pop;      // average POP index over most recent years
  double  t_cpue;      // work: number of years over which slope of CPUEs(age4+) is estimated
  double  slope_cpue;  // work: slope of cpue for age4+
  double  k_cpue;      // work: gain parameter for cpue slope
  double  mu_gt;       // average GT estimate over most recent years
//  double  slope_gt;   // slope of GT estimates
//  double  mean_pop;   // recent mean of POP/HSP
//  double  delta_wrk;  // work variable

  int  flg_gt_limit;  // flag of whether GT limit is used

  double muCK; // f/ Rich's code

  double  TAC_cpue;     // TAC by cpue slope
  double  TAC1_cpue, TAC2_cpue;
  double  wt_cpue, wt;
  double  TAC_gt_limit; // TAC by GT limit
//  double  TAC_gt;       // TAC by GT slope
//  double  TAC_pop;      // TAC by POP/HSP
  double  tmp_tac;      // temporary TAC before constraints
  double  tac_change;   // amount of tac change
  double  TAC;          // final specified TAC for current year


//  // for Rich's code
//  int y;
//  double lamGT;
//  double muGT;
//  double mut;

  // variables (f/ Rich's code)
  int a,c,cc,i,y;
  double wsum;

  // initialize vector and matrix variables for POP index calculation
  // (f/ Rich's code)
  POPind.initialize();
  Mja.initialize();
  Rja.initialize();
  wja.initialize();
//  w.initialize();

  // concatenate the CPUE and GT index data
  for (int ii = first_cpue_yr; ii <= current_yr - 2; ii++) {

    yrs(ii) = ii;

    if(ii <= last_cpue_yr) {
      obs_cpue(ii) = hist_cpue(ii);  //age4+
    }

    if(ii > last_cpue_yr) {
      obs_cpue(ii) = proj_cpuen(ii);  //age4+
    }

    // copied Rich's example MP code
    if(ii == 2016){  // hard-wired (norio: from GeneTagging_Data_GT2016_update_2019-04-30.xlsx)
      obs_gtn(ii) = 2271416.;
      obs_cvgtn(ii) = 0.224;
      obs_gtrec(ii) = 20.;
    }
    if(ii == 2017){  // hard-wired (norio: from GeneTagging2017_Data_2019-04-30.xlsx)
      obs_gtn(ii) = 1154020.;
      obs_cvgtn(ii) =  0.122;
      obs_gtrec(ii) = 67.;
    }
    if(ii == 2018){  // hard-wired from GT estimate in 2020
      obs_gtn(ii) = 1140000.;
      obs_cvgtn(ii) =  0.12;
      obs_gtrec(ii) = 66.;
    }
    if(ii == 2019){  // assumed for the OMMP10, June 2019
      obs_gtn(ii) = -11.;
      obs_cvgtn(ii) =  -11.;
      obs_gtrec(ii) = 0.;
    }

    if(ii > 2019){
      obs_gtn(ii) = proj_gtn(ii);
      obs_cvgtn(ii) = proj_cvgtn(ii);
      obs_gtrec(ii) = proj_gtrec(ii);
    }

  }


  // ******** CK POP-index (f/ Rich's code) ********

  int cmin = 2002;          // min cohort (norio comment)
  int cmax = current_yr-5;  // max cohort
  int ymin = 2006;          // min adult capture year
  int ymax = current_yr;    // max adult capture year

  ///////////////
  // POP index //
  ///////////////

  // aggregate across adult capture ages

  for(int n=1;n<=nPOPdat2;n++) {

     c = POPhist(n,1);
     y = POPhist(n,2);
     a = POPhist(n,3);
     Mja(c,y) += POPhist(n,5);
     Rja(c,y) += POPhist(n,4);

  }

  for(int n=1;n<=nPOPdat;n++) {

     c = POPproj(n,1);
     y = POPproj(n,2);
     a = POPproj(n,3);
     Mja(c,y) += POPproj(n,5);
     Rja(c,y) += POPproj(n,4);

  }

  // create the weighting for the index and sum the index

  for(c=cmin;c<=cmax;c++) {

     wsum = sum(Rja(c));
     for(y=ymin;y<=ymax;y++) wja(c,y) = Rja(c,y)/wsum;
     for(y=ymin;y<=ymax;y++) {

       if(Rja(c,y) > 0) POPind(c) += (Mja(c,y)/Rja(c,y)) * wja(c,y);

     }

  }

  //cout << POPind << endl;

  //Iref = sum(POPind(c1,c2))/double(c2-c1+1);
  muCK = 0.;
  //for(y=cmax-tau3+1;y<=cmax;y++) muCK += POPind(y)/double(tau3);
  for(y=cmax-t_pop+1 ;y<=cmax;y++) muCK += POPind(y)/double(t_pop);

  //cout << current_yr << " " << Iref << " " << muCK << " " << delCK << " " << delCPUE << " " << delGT << endl;
  //cout << current_yr << " " << muCK << " " << endl;

  mu_pop = muCK;


  // ******** CPUE(age4+) slope part of HCR ********

  int yr1, yr2;

  if(current_yr <= tune_year){

    yr1 = current_yr-2 - t1_cpue + 1;
    yr2 = current_yr-2;

    slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );

    k_cpue = slope_cpue < 0 ? k1_cpue : k2_cpue;

    TAC_cpue = quota(implementation_yr-1) * (1. + k_cpue * slope_cpue);

  }else if(current_yr > tune_year){

    if(mu_pop <= (alpha1*ssb_pop_target)){

      yr1 = current_yr-2 - t1_cpue + 1;
      yr2 = current_yr-2;

      slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );

      k_cpue = slope_cpue < 0 ? k1_cpue : k2_cpue;

      TAC_cpue = quota(implementation_yr-1) * (1. + k_cpue * slope_cpue);

    }else if(mu_pop >= (alpha2*ssb_pop_target)){

      yr1 = current_yr-2 - t2_cpue + 1;
      yr2 = current_yr-2;

      slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );

      k_cpue = slope_cpue < 0 ? k3_cpue : k4_cpue;

      TAC_cpue = quota(implementation_yr-1) * (1. + k_cpue * slope_cpue);

    }else if((alpha1*ssb_pop_target) < mu_pop && mu_pop < (alpha2*ssb_pop_target)){

      // calc TAC1_cpue
      yr1 = current_yr-2 - t1_cpue + 1;
      yr2 = current_yr-2;

      slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );

      k_cpue = slope_cpue < 0 ? k1_cpue : k2_cpue;

      TAC1_cpue = quota(implementation_yr-1) * (1. + k_cpue * slope_cpue);

      // calc TAC2_cpue
      yr1 = current_yr-2 - t2_cpue + 1;
      yr2 = current_yr-2;

      slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );

      k_cpue = slope_cpue < 0 ? k3_cpue : k4_cpue;

      TAC2_cpue = quota(implementation_yr-1) * (1. + k_cpue * slope_cpue);

      // calc weight
      wt_cpue = (mu_pop/ssb_pop_target - alpha1)/(alpha2 - alpha1);

      // calc weighted average TAC
      TAC_cpue = wt_cpue*TAC2_cpue + (1. - wt_cpue)*TAC1_cpue;

    } // end of ssb_pop_target if

  } // end of tune_year if


  // ******** GT(age2) part of HCR ********

  // (1) used as limit

  //----calculate ratio of (recent mean GT estimate)/(limit value of N_age2)

  // first, need to remove data of when no esitimate is available ("-11" is embedded)
  // the following code was based on the code (C1GT2.tpl) provided by Rich

  yr1 = current_yr - t_gt_limit - 1 < first_yr-2 ? first_yr-2 : current_yr - t_gt_limit - 1;
  yr2 = current_yr-2;

  yr1 = yr1 < 2016 ? 2016 : yr1;

  int cnt;
  cnt = 0;

  mu_gt = 0.;

  for(int ii=yr1; ii<=yr2; ii++) {

    mu_gt = obs_gtn(ii) != -11. ? mu_gt+obs_gtn(ii) : mu_gt;
    cnt = obs_gtn(ii) != -11. ? cnt+1 : cnt;

  }

  mu_gt /= cnt;

  mu_gt = mu_gt < 1.e-6? n_age2_limit : mu_gt;  // if mu_gt = 0 then set n_age2_limit to avoid error and set TAC by CPUE only

  TAC_gt_limit = quota(implementation_yr-1) * k1_gt_limit * (mu_gt/n_age2_limit)*(mu_gt/n_age2_limit);

  //flg_gt_limit = mu_gt < n_age2_limit ? 1 : 0;


  // ******** set TAC ********
  
  // use TAC from CPUE slope part of HCR only
  if(swit_cpue==1 && swit_gt_limit==0) tmp_tac = TAC_cpue;

  // use TAC from GT limit part of HCR only
  if(swit_gt_limit==1 && swit_cpue==0){

    tmp_tac = flg_gt_limit==1 ? TAC_gt_limit : quota(implementation_yr-1);

  }

  // use TAC from CPUE slope and GT limit parts of HCR
  if(swit_cpue==1 && swit_gt_limit==1){

    //tmp_tac = flg_gt_limit==1 ? min(TAC_gt_limit, TAC_cpue) : TAC_cpue;

    if(mu_gt <= (beta1*n_age2_limit)){

      tmp_tac = min(TAC_gt_limit, TAC_cpue);

    }else if(mu_gt >= (beta2*n_age2_limit)){

      tmp_tac = TAC_cpue;

    }else if((beta1*n_age2_limit) < mu_gt && mu_gt < (beta2*n_age2_limit)){

      wt = (mu_gt/n_age2_limit - beta1)/(beta2 - beta1);

      tmp_tac = wt*TAC_cpue + (1. - wt)*min(TAC_gt_limit, TAC_cpue);
      
    }

  }


  tac_change = tmp_tac - quota(implementation_yr-1);

  if( fabs(tac_change) <= min_change ){

    TAC = quota(implementation_yr-1);

  }
  else if( tac_change > max_change_up ) {

    TAC = quota(implementation_yr-1) + max_change_up;

  }
  else if( tac_change < -1 * max_change_down) {

    TAC = quota(implementation_yr-1) - max_change_down;

  }
  else{

    TAC = tmp_tac;

  }

  if(TAC < 0.) TAC = 0.;
  if(TAC > maxTAC) TAC = maxTAC;

  //----output for sbtproj.exe (copied from Rich's test code)
  ofstream ofs("tacfile");
  ofs << implementation_yr << " " << TAC << endl;


  //----debug write?
  if(debug_write){
    DebugWrite1(slope_cpue, mu_gt, mu_pop);
  }

  exit(1);
 END_CALCS

PROCEDURE_SECTION
  f=0.;


FUNCTION void DebugWrite1(double slope_cpue, double mu_gt, double mu_pop)

  ofstream fdebug;
  int ii;


  fdebug.open("debug_write_ctrl.txt", ios::trunc);

  fdebug << "---- check NT*_ctrl.dat inputs ----" << endl;

  fdebug << "ssb_pop_target = " << ssb_pop_target << endl;
  fdebug << "t_pop          = " << t_pop          << endl;

  fdebug << "tune_year      = " << tune_year      << endl;

  fdebug << "k1_cpue        = " << k1_cpue        << endl;
  fdebug << "k2_cpue        = " << k2_cpue        << endl;
  fdebug << "t1_cpue        = " << t1_cpue        << endl;

  fdebug << "k3_cpue        = " << k3_cpue        << endl;
  fdebug << "k4_cpue        = " << k4_cpue        << endl;
  fdebug << "t2_cpue        = " << t2_cpue        << endl;

  fdebug << "k1_gt_limit    = " << k1_gt_limit    << endl;
  fdebug << "t_gt_limit     = " << t_gt_limit     << endl;
  fdebug << "n_age2_limit   = " << n_age2_limit   << endl;

  fdebug << "alpha1         = " << alpha1         << endl;
  fdebug << "alpha2         = " << alpha2         << endl;
  fdebug << "beta1          = " << beta1          << endl;
  fdebug << "beta2          = " << beta2          << endl;

  fdebug << "swit_cpue     = " << swit_cpue     << endl;
  fdebug << "swit_gt_limit = " << swit_gt_limit << endl;

  fdebug << "max_change_up   = " << max_change_up   << endl;
  fdebug << "max_change_down = " << max_change_down << endl;
  fdebug << "min_change      = " << min_change      << endl;
  fdebug << "maxTAC      = " << maxTAC      << endl;

  fdebug << "debug_write  = " << debug_write  << endl;

  for(ii = first_cpue_yr; ii <= last_cpue_yr; ii++){
    fdebug << "hist_cpue(" << ii << ") = " << hist_cpue(ii) << endl;
  }

  fdebug.close();


  fdebug.open("debug_write_cpue_gt.txt", ios::trunc);

  fdebug << "---- check values of cpue, GT, GT CV, POP set to variables ----" << endl;
  for(ii = first_cpue_yr; ii <= current_yr-2; ii++){
    fdebug << "obs_cpue(" << ii << ") = " << obs_cpue(ii) << endl;
  }

  fdebug << endl;
  for(ii = 2016; ii <= current_yr-2; ii++){
    fdebug << "obs_gtn(" << ii << ") = " << obs_gtn(ii) << endl;
  }

  fdebug << endl;
  for(ii = 2016; ii <= current_yr-2; ii++){
    fdebug << "obs_cvgtn(" << ii << ") = " << obs_cvgtn(ii) << endl;
  }

  fdebug << endl;
  for(ii = 2002; ii <= current_yr-5; ii++){
    fdebug << "POPind(" << ii << ") = " << POPind(ii) << endl;
  }

  fdebug.close();


  fdebug.open("debug_write_primary1.txt", ios::app);

  fdebug << "---- check values set to primary variables - 1 ----" << endl;
  fdebug << "current_yr = " << current_yr << endl;
  fdebug << "  slope_cpue = " << slope_cpue << endl;
  fdebug << "  mu_gt      = " << mu_gt      << endl;
  fdebug << "  mu_pop     = " << mu_pop     << endl;
  fdebug << endl;

  fdebug.close();

  return;
