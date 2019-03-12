// Management Procedure Simple "NT4" modified from "NT3" and "NT2" (based on "NT1") code
//
// by norio takahashi
// start  : Apr  04, 2018 (as "NT1")
// end    : Jul  26, 2018 - modified to "NT4", just added tuning year's condition to "NT3"
//
//
DATA_SECTION
  !!ad_comm::change_datafile_name("NT4_ctrl.dat"); // set parameter values for MP
  // POP target part parameters
  init_number    ssb_pop_target  // target value of SSB comparable to POP index
  init_int            t_pop                 // number of years over which POP indices (SSB) are averaged

  init_int    tune_year  // tuning year

  // average POP <= target POP case
  init_number    k1_cpue  // gain parameter when cpue slope < 0 (decrease)
  init_number    k2_cpue  // gain parameter when cpue slope >= 0 (increase)
  init_int            t1_cpue  // number of years over which slope of CPUEs(age4+) is estimated

  // average POP > target POP case
  init_number    k3_cpue  // gain parameter when cpue slope < 0 (decrease)
  init_number    k4_cpue  // gain parameter when cpue slope >= 0 (increase)
  init_int            t2_cpue   // number of years over which slope of CPUEs(age4+) is estimated

  // GT(age2) part parameters
  // used as limit
  init_number    k1_gt_limit     // smoothing parameter when delta_gt_limit < 1
  init_int            t_gt_limit        // number of years over which GT estimates(age2) are averaged
  init_number    n_age2_limit  // limit value of number of age 2 fish

  // switch for which TAC from parts of HCR is used
  init_int    swit_cpue        // switch for TAC from cpue slope part of HCR is used
  init_int    swit_gt_limit   // switch for TAC from GT limit part of HCR is used

  init_number    max_change_up       // maximum TAC increase
  init_number    max_change_down  // maximum TAC decrease
  init_number    min_change             // minimum TAC increase

  init_int    debug_write  // debug write? yes = 1, no = 0


  int    first_cpue_yr   // first year of historical cpue data
  int    last_cpue_yr    // last year of historical cpue data

  !!first_cpue_yr = 1969;
  !!last_cpue_yr = 2016;

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
  double  mu_pop;     // average POP index over most recent years
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
  double  TAC_gt_limit; // TAC by GT limit
//  double  TAC_gt;       // TAC by GT slope
//  double  TAC_pop;      // TAC by POP/HSP
  double  tmp_tac;      // temporary TAC before constraints
  double  tac_change;   // amount of tac change
  double  TAC;         // final specified TAC for current year


//  // for Rich's code
//  int y;
//  double lamGT;
//  double muGT;
//  double mut;

  // variables (f/ Rich's code)
  int a,c,cc,i,y;
  double wsum;

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
    if(ii == 2016){
      obs_gtn(ii) = 2417786.; // hard-wired (norio corrected according to Ann's e-mail 180531)
      obs_cvgtn(ii) = 0.2132; // hard-wired (norio added f/ estimates Excel file)
      obs_gtrec(ii) = 22.;     // hard-wired
    }

    if(ii > 2016){
      obs_gtn(ii) = proj_gtn(ii);
      obs_cvgtn(ii) = proj_cvgtn(ii);
      obs_gtrec(ii) = proj_gtrec(ii);
    }

  }


  // ******** CPUE(age4+) slope part of HCR ********

  // CK POP-index (f/ Rich's code)

  int cmin = 2002;
  int cmax = current_yr-5;
  int ymin = 2006;
  int ymax = current_yr;

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

//  if(mean_pop >  ssb_pop_target) TAC_pop = quota(implementation_yr-1) * (1. + k1_pop_target*delta_wrk);
//  if(mean_pop <= ssb_pop_target) TAC_pop = quota(implementation_yr-1) * (1. + k2_pop_target*delta_wrk);


  // ---- calculate average CPUE over most recent years

  int yr1, yr2;

  if( mu_pop <= ssb_pop_target ){

    yr1 = current_yr-2 - t1_cpue + 1;
    yr2 = current_yr-2;

    slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /
                 sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );

    k_cpue = slope_cpue < 0 ? k1_cpue : k2_cpue;

  }
  else if( mu_pop > ssb_pop_target ){

    if(current_yr > tune_year){
      yr1 = current_yr-2 - t2_cpue + 1;
      yr2 = current_yr-2;
    }
    else{
      yr1 = current_yr-2 - t1_cpue + 1;
      yr2 = current_yr-2;
    }

    slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /
                 sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );

    if(current_yr > tune_year){
      k_cpue = slope_cpue < 0 ? k3_cpue : k4_cpue;
    }
    else{
      k_cpue = slope_cpue < 0 ? k1_cpue : k2_cpue;
    }

  }

  TAC_cpue = quota(implementation_yr-1) * (1. + k_cpue * slope_cpue);


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

  TAC_gt_limit = quota(implementation_yr-1) * k1_gt_limit * (mu_gt/n_age2_limit)*(mu_gt/n_age2_limit);

  flg_gt_limit = mu_gt < n_age2_limit ? 1 : 0;


  // (2) used as slope

  // ---- calculate GT slope

  // added based on the code (C1GT2.tpl) provided by Rich
//  int n = current_yr-2;
//  int m = current_yr - t_gt - 1 < first_yr-2 ? first_yr-2 : current_yr - t_gt - 1;
//  m = m < 2016 ? 2016 : m;
//  double cnt,sum1,sum2;
//  cnt = sum1 = sum2 = 0.;
//  // minimum number of years to calculate trend over is 3?
//  if(n-m > 3) {
//  		sum1 = sum2 = muGT = mut = 0.;
//  		for(y=m;y<=n;y++) {
//
//  			muGT = obs_gtn(y) != -11. ? muGT+log(obs_gtn(y)) : muGT;
//			mut = obs_gtn(y) != -11. ? mut+double(y) : mut;
//			cnt = obs_gtn(y) != -11. ? cnt+1. : cnt;
//
//		}
//
//		mut /= cnt;
//		muGT /= cnt;
//		for(y=m;y<=n;y++) sum1 = obs_gtn(y) != -11. ? sum1+(log(obs_gtn(y))-muGT)*(double(y)-mut) : sum1;
//		for(y=m;y<=n;y++) sum2 = obs_gtn(y) != -11. ? sum2+(double(y)-mut)*(double(y)-mut) : sum2;
//		lamGT = sum1/sum2;
//
//  }
//  else lamGT = 0.;
//
//  slope_gt = lamGT;
//
//  double k_gt = slope_gt <0 ? k1_gt : k2_gt;
//
//  TAC_gt = quota(implementation_yr-1) * (1. + k_gt * slope_gt);


//  // ---- POP target part of HCR
//
//  // CK POP-index (f/ Rich's code)
//
//  int cmin = 2002;
//  int cmax = current_yr-5;
//  int ymin = 2006;
//  int ymax = current_yr;
//
//  ///////////////
//  // POP index //
//  ///////////////
//
//  // aggregate across adult capture ages
//
//  for(int n=1;n<=nPOPdat2;n++) {
//
//     c = POPhist(n,1);
//     y = POPhist(n,2);
//     a = POPhist(n,3);
//     Mja(c,y) += POPhist(n,5);
//     Rja(c,y) += POPhist(n,4);
//
//  }
//
//  for(int n=1;n<=nPOPdat;n++) {
//
//     c = POPproj(n,1);
//     y = POPproj(n,2);
//     a = POPproj(n,3);
//     Mja(c,y) += POPproj(n,5);
//     Rja(c,y) += POPproj(n,4);
//
//  }
//
//  // create the weighting for the index and sum the index
//
//  for(c=cmin;c<=cmax;c++) {
//
//     wsum = sum(Rja(c));
//     for(y=ymin;y<=ymax;y++) wja(c,y) = Rja(c,y)/wsum;
//     for(y=ymin;y<=ymax;y++) {
//
//       if(Rja(c,y) > 0) POPind(c) += (Mja(c,y)/Rja(c,y)) * wja(c,y);
//
//     }
//
//  }
//
//  //cout << POPind << endl;
//
//  //Iref = sum(POPind(c1,c2))/double(c2-c1+1);
//  muCK = 0.;
//  //for(y=cmax-tau3+1;y<=cmax;y++) muCK += POPind(y)/double(tau3);
//  for(y=cmax-t_pop+1 ;y<=cmax;y++) muCK += POPind(y)/double(t_pop);
//
//  //cout << current_yr << " " << Iref << " " << muCK << " " << delCK << " " << delCPUE << " " << delGT << endl;
//  //cout << current_yr << " " << muCK << " " << endl;
//
//
//  //----calculate ratio of (recent mean POP estimate)/(limit value of POP)
//
//  mean_pop = muCK;
//
//  delta_wrk = (ssb_pop_target - mean_pop)/mean_pop;
//
//  if(mean_pop >  ssb_pop_target) TAC_pop = quota(implementation_yr-1) * (1. + k1_pop_target*delta_wrk);
//  if(mean_pop <= ssb_pop_target) TAC_pop = quota(implementation_yr-1) * (1. + k2_pop_target*delta_wrk);


  // ----set TAC

  // use TAC from CPUE slope part of HCR only
  if(swit_cpue==1 && swit_gt_limit==0) tmp_tac = TAC_cpue;

  // use TAC from GT limit part of HCR only
  if(swit_gt_limit==1 && swit_cpue==0){

    tmp_tac = flg_gt_limit==1 ? TAC_gt_limit : quota(implementation_yr-1);

  }

  // use TAC from CPUE slope and GT limit parts of HCR
  if(swit_cpue==1 && swit_gt_limit==1){

    tmp_tac = flg_gt_limit==1 ? min(TAC_gt_limit, TAC_cpue) : TAC_cpue;

  }

//  // use TAC from GT slope part of HCR only
//  if(swit_gt==1 && (swit_cpue+swit_gt_limit+swit_pop)==0) tmp_tac = TAC_gt;

//  // use TAC from GT limit and GT slope parts of HCR only
//  if(swit_gt_limit==1 && swit_gt==1 && (swit_cpue+swit_pop)==0){
//
//    tmp_tac = flg_gt_limit==1 ? min(TAC_gt_limit, TAC_gt) : TAC_gt;
//
//  }

//  // use TAC from CPUE slope and GT slope parts of HCR only
//  if(swit_cpue==1 && swit_gt==1 && (swit_gt_limit+swit_pop)==0) (TAC_cpue + TAC_gt)/2.;

//  // use TAC from POP target part of HCR only
//  if(swit_pop==1 && (swit_cpue+swit_gt_limit+swit_gt)==0) tmp_tac = TAC_pop;

//  // use TAC from cpue slope, GT limit, and GT slope parts of HCR only
//  if(swit_cpue==1 && swit_gt_limit==1 && swit_gt==1 && swit_pop==0){
//
//    tmp_tac = (TAC_cpue + TAC_gt)/2.;
//    tmp_tac = flg_gt_limit==1 ? min(TAC_gt_limit, tmp_tac) : tmp_tac;
//
//  }

//  // use TAC from cpue slope, GT limit, GT slope, and POP target parts of HCR
//  if(swit_cpue==1 && swit_gt_limit==1 && swit_gt==1 && swit_pop==1){
//
//    tmp_tac = (TAC_cpue + TAC_gt + TAC_pop)/3.;
//    tmp_tac = flg_gt_limit==1 ? min(TAC_gt_limit, tmp_tac) : tmp_tac;
//
//  }


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


  //----output for sbtproj.exe (copied from Rich's test code
  ofstream ofs("tacfile");
  ofs << implementation_yr << " " << TAC << endl;


  //----debug write?
  if(debug_write){
    DebugWrite1(slope_cpue, mu_gt, mu_pop);
//    DebugWrite1(slope_cpue, mean_gt/n_age2_limit, slope_gt, mean_pop);
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

//  fdebug << "k1_gt        = " << k1_gt        << endl;
//  fdebug << "k2_gt        = " << k2_gt        << endl;
//  fdebug << "t_gt         = " << t_gt         << endl;

//  fdebug << "k1_pop_target  = " << k1_pop_target  << endl;
//  fdebug << "k2_pop_target  = " << k2_pop_target  << endl;
//  fdebug << "ssb_pop_target = " << ssb_pop_target << endl;
//  fdebug << "t_pop          = " << t_pop          << endl;

  fdebug << "swit_cpue     = " << swit_cpue     << endl;
  fdebug << "swit_gt_limit = " << swit_gt_limit << endl;
//  fdebug << "swit_gt       = " << swit_gt       << endl;
//  fdebug << "swit_pop      = " << swit_pop      << endl;

  fdebug << "max_change_up   = " << max_change_up   << endl;
  fdebug << "max_change_down = " << max_change_down << endl;
  fdebug << "min_change      = " << min_change      << endl;

  fdebug << "debug_write  = " << debug_write  << endl;

  for(ii = first_cpue_yr; ii <= last_cpue_yr; ii++){
    fdebug << "hist_cpue(" << ii << ") = " << hist_cpue(ii) << endl;
  }

  fdebug.close();


  fdebug.open("debug_write_cpue_gt.txt", ios::trunc);
//  fdebug.open("debug_write_cpue_gt_pop.txt", ios::trunc);

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
//  fdebug << "mean_gt/n_age2_limit = " << ratio_gt_limit << endl;
//  fdebug << "slope_gt             = " << slope_gt       << endl;
//  fdebug << "mean_pop             = " << mean_pop       << endl;
  fdebug << endl;

  fdebug.close();

  return;
